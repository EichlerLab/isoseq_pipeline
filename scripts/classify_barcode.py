__author__ = 'etseng@pacb.com'


import os, sys
import os.path as op
import math
import re
import logging
import multiprocessing
from collections import defaultdict, namedtuple
from pbcore.util.Process import backticks
from pbcore.io.FastaIO import FastaReader, FastaWriter
from pbtools.pbtranscript.PBTranscriptException import PBTranscriptException
from pbtools.pbtranscript.io.DOMIO import DOMReader
from pbtools.pbtranscript.io.ReadAnnotation import ReadAnnotation
from pbtools.pbtranscript.io.Summary import ClassifySummary
from pbtools.pbtranscript.Utils import revcmp, realpath, \
    generateChunkedFN, cat_files, real_upath, ln


PBMATRIXFN = "PBMATRIX.txt"
PRIMERFN = "primers.fa"
PRIMERFRONTENDFN = "primers.front_end.fa"
PRIMERCHIMERAFN = "primers.chimera.fa"
PRIMERREPORTFN = "primer_info.csv"
FRONTENDDOMFN = "hmmer.front_end.dom"
FLCHIMERADOMFN = "hmmer.fl.chimera.dom"
NFLCHIMERADOMFN = "hmmer.nfl.chimera.dom"
CLASSIFYSUMMARY = "classify_summary.txt"


# ChimeraDetectionOptions:
# Minimum length to output a (trimmed) sequence.
# Minimum phmmer score for primer hit.
# Minimum distance a primer has to be from end of sequence.
# Maximum distance between adjacent primer hits to consider as chimera.
# Search primers within windows of length primer_search_window.
# Apply chimera detection on non-full-length reads.
ChimeraDetectionOptions = namedtuple("ChimeraDetectionOptions",
                                     ("min_seq_len min_score min_dist_from_end max_adjacent_hit_dist " +
                                      "primer_search_window detect_chimera_nfl"))


class PBRead(object):

    """Class for PacBio read."""

    def __init__(self, read):
        self.sequence = read.sequence
        self.name = read.name
        self.isCCS = False
        self.start, self.end = None, None
        self.movie, self.zmw = None, None

        # pattern: m....../1123/ccs
        m = re.search(r"(.+)/(\d+)/ccs", self.name)
        if m is None:
            # pattern: m...../1123 (alternative ccs)
            m = re.search(r"(.+)/(\d+)$", self.name)

        if m is not None:
            self.isCCS = True
            self.movie, self.zmw = m.groups()[0], int(m.groups()[1])
            self.start, self.end = 0, len(self.sequence)
        else:  # pattern: m...../1123/23_450
            m = re.search(r"(.+)/(\d+)/(\d+)_(\d+)", self.name)
            if m is not None:
                self.movie = m.groups()[0]
                self.zmw, self.start, self.end = \
                    [int(x) for x in m.groups()[1:]]
            else:
                raise ValueError("Unsupported PacBio read {r}".
                                 format(r=self.name))


class ClassifierException(PBTranscriptException):

    """
    Exception class for Classifier.
    """

    def __init__(self, msg):
        print msg # Necessary because PBTranscriptException isn't printing the error message
        PBTranscriptException.__init__(self, "classify", msg)


class Classifier(object):

    """
    Class for classifying reads based on whether they are full length and
    have their 5' primer, 3' primer and poly A tail seen, trim primers and
    PolyA tails from reads, and finally determine whether the trimmed reads
    are chimeras.
    """

    def __init__(self, reads_fn="ccs.fasta", out_dir="classifyOut/",
                 out_reads_fn="isoseq_draft.fasta", primer_fn_forward=None, primer_fn_reverse=None,
                 primer_report_fn=None, summary_fn=None,
                 cpus=1, change_read_id=True,
                 opts=ChimeraDetectionOptions(50, 10, 100, 50, 150, False),
                 out_nfl_fn=None, out_flnc_fn=None,
                 ignore_polyA=False, keep_primer=False, reuse_dom=False):
        self.reads_fn = realpath(reads_fn)
        self.out_dir = realpath(out_dir)
        self.cpus = cpus
        self.change_read_id = change_read_id
        self.chimera_detection_opts = opts
        self.ignore_polyA = ignore_polyA
        self.keep_primer = keep_primer # if True, primers are not removed (useful for targeted)
        self.reuse_dom = reuse_dom

        # for now, the barcoded primer files must be given!
        assert primer_fn_forward is not None
        assert primer_fn_reverse is not None
        self.primer_fn_forward = primer_fn_forward
        self.primer_fn_reverse = primer_fn_reverse
        # The output fasta file.
        self.out_all_reads_fn = realpath(out_reads_fn)

        # Intermediate output fasta file before chimera detection.
        #     trimmed full-length reads: fl.trimmed.fasta
        # and
        #     trimmed non-full-length reads: nfl.trimmed.fasta
        self._trimmed_fl_reads_fn = op.join(self.out_dir, "fl.trimmed.fasta")
        self._trimmed_nfl_reads_fn = op.join(self.out_dir, "nfl.trimmed.fasta")

        self.primer_front_back_fn = op.join(self.out_dir, PRIMERFRONTENDFN)
        self.primer_chimera_fn = op.join(self.out_dir, PRIMERCHIMERAFN)

        # The output primer file: primer_info.csv
        self.primer_report_fn = primer_report_fn \
            if primer_report_fn is not None else \
            ".".join(out_reads_fn.split('.')[:-1]) + "." + PRIMERREPORTFN
        # primer reports for nfl reads before chimera detection. Note that
        # chimera detection is not necessary for nfl reads.
        self._primer_report_nfl_fn = op.join(self.out_dir,
                                             "primer_report.nfl.csv")
        # primer reports for fl reads after chimera detection. Note that
        # chimera detection is required for fl reads.
        self._primer_report_fl_fn = op.join(self.out_dir,
                                            "primer_report.fl.csv")

        # The matrix file: PBMATRIX.txt
        self.pbmatrix_fn = op.join(self.data_dir, PBMATRIXFN)

        # The output phmmer Dom file for trimming primers: hmmer.front_end.dom
        self.out_front_back_dom_fn = op.join(self.out_dir, FRONTENDDOMFN)
        # The output phmmer Dom file for chimera detection:
        #     hmmer.fl.chimera.dom and hmmer.nfl.chimera.dom
        self.out_trimmed_fl_dom_fn = op.join(self.out_dir, FLCHIMERADOMFN)
        self.out_trimmed_nfl_dom_fn = op.join(self.out_dir, NFLCHIMERADOMFN)

        self.chunked_front_back_reads_fns = None
        self.chunked_front_back_dom_fns = None

        #self.chunked_trimmed_reads_fns = None
        #self.chunked_trimmed_reads_dom_fns = None

        # The summary file: *.classify_summary.txt
        self.summary = ClassifySummary()
        self.summary_fn = summary_fn if summary_fn is not None else \
            ".".join(out_reads_fn.split('.')[:-1]) + \
            "." + CLASSIFYSUMMARY

        self.out_nfl_fn = realpath(out_nfl_fn) if out_nfl_fn is not None \
            else op.join(self.out_dir, "nfl.fasta")
        self.out_nflnc_fn = op.join(self.out_dir, "nflnc.fasta")
        self.out_nflc_fn = op.join(self.out_dir, "nflc.fasta")

        self.out_flnc_fn = realpath(out_flnc_fn) if out_flnc_fn is not None \
            else op.join(self.out_dir, "flnc.fasta")
        self.out_flc_fn = op.join(self.out_dir, "flc.fasta")

    def __str__(self):
        return ("reads_fn={0}\n".format(self.reads_fn) +
                "primer_fn={0}\n".format(self.primer_fn) +
                "out_all_reads_fn={0}\n".format(self.out_all_reads_fn) +
                "pbmatrix_fn={0}\n".format(self.pbmatrix_fn) +
                "out_front_back_dom_fn={0}\n".
                format(self.out_front_back_dom_fn))

    @property
    def data_dir(self):
        """Return the data dir which has primers.fa and PBMATRIX.txt."""
        return op.join(op.dirname(op.realpath(__file__)), "data")

    def _validate_inputs(self, reads_fn, primer_fn_forward, primer_fn_reverse, pbmatrix_fn):
        """Validate whether input files and required data files all exist."""
        logging.info("Checking input files.")
        if not op.exists(reads_fn):
            raise ClassifierException(
                "Unable to find reads file: {fn}".format(fn=reads_fn))
        if not op.exists(primer_fn_forward):
            raise ClassifierException(
                "Unable to find forward primer file: {fn}".format(fn=primer_fn_forward))
        if not op.exists(primer_fn_reverse):
            raise ClassifierException(
                "Unable to find reverse primer file: {fn}".format(fn=primer_fn_reverse))
        if not op.exists(pbmatrix_fn):
            raise ClassifierException(
                "Unable to find matrix file for PacBio reads: {fn}".
                format(fn=pbmatrix_fn))

    def _checkPhmmer(self):
        """Check phmmer can be called successfully."""
        logging.info("checking for phmmer existence.")
        _output, errCode, errMsg = backticks("phmmer -h > /dev/null")
        if errCode != 0:
            raise ClassifierException("Unable to invoke phmmer.\n{e}".
                                      format(e=errMsg))

    def _processPrimers(self, primer_fn_forward, primer_fn_reverse, window_size, primer_out_fn,
                        revcmp_primers=False):
        """
        Do basic sanity checks that:
        (1) all primers in forward start with f_xxx and are unique
        (2) all primers in reverse start with r_xxx and are unique
        (3) check that no forward primers appear in reverse primers (no symmetry)
        (4) write the primers (f_xxx, f_xxx_revcmp, r_xxx, r_xxx_revcmp) all to one primer file
        """
        def sanity_check_primers(reader, prefix):
            """
            Go through the primers, check that the prefix exists and all seqs are unique
            """
            primers = {} # primer -> sequence, but can also contain the revcmp version with _revcmp suffix
            for r in reader:
                if not r.name.startswith(prefix):
                    errMsg = "Primer should start with %s, but saw: %s" % (prefix, r.name)
                    raise ClassifierException(errMsg)
                if len(r.sequence) > window_size:
                    errMsg = "Primer {n} has length {l} which is longer than {k}.".\
                     format(n=r.name, l=len(r.sequence), k=window_size)
                    logging.error(errMsg)
                    raise ClassifierException(errMsg)
                ss = r.sequence.upper()
                if ss in primers.itervalues():
                    errMsg = "Duplicate sequences found for", ss
                    raise ClassifierException(errMsg)
                primers[r.name.strip()] = r.sequence
                # revcmp not needed becuz phmmer does both strands apparently...
                #primers[r.name.strip() + "_revcmp"] = revcmp(r.sequence)
            return primers


        logging.info("Process primers for {case}.".
                     format(case=("finding primers" if not revcmp_primers
                                  else "detecting chimeras")))
        reader_f = FastaReader(primer_fn_forward)
        reader_r = FastaReader(primer_fn_reverse)

        primers_f = sanity_check_primers(reader_f, prefix="f_")
        primers_r = sanity_check_primers(reader_r, prefix="r_")

        reader_f.close()
        reader_r.close()

        same_seqs = set(primers_f.values()).intersection(primers_r.values())
        if len(same_seqs) > 0:
            errMsg = "Identical sequences found in both Forward/Reverse!\n"
            errMsg += "\n".join(same_seqs)
            raise ClassifierException(errMsg)

        # Write Fi and reverse-complemented Ri to primer_out_fn
        with open(primer_out_fn, 'w') as f:
            for (name, seq) in primers_f.iteritems():
                f.write(">{n}\n{s}\n".format(n=name, s=seq))
            for (name, seq) in primers_r.iteritems():
                f.write(">{n}\n{s}\n".format(n=name, s=revcmp(seq)))
        return primers_f.keys() + primers_r.keys()

    @property
    def numReads(self):
        """Return the number of reads in reads_fn."""
        cmd = "grep -c '>' {r}".format(r=real_upath(self.reads_fn))
        output, errCode, errMsg = backticks(cmd)
        if errCode != 0:
            raise ClassifierException(
                "Error reading file {r}:{e}".
                format(r=self.reads_fn, e=str(errMsg)))
        return int(output[0])

    def _chunkReads(self, reads_fn, reads_per_chunk, chunked_reads_fns,
                    extract_front_back_only=True, window_size=100):
        """Split reads within reads_fn into multiple chunks each containing
        at most 'reads_per_chunk' reads, save to files in 'chunked_reads_fns'.
        If extract_front_back_only is true, extract the first and the last
        'window_size' bases and save them as readname_front and readname_back.
        Otherwise, copy read names and sequences entirely.
        """
        logging.debug("Split {f} into ".format(f=reads_fn) +
                      "{n} chunks, ".format(n=len(chunked_reads_fns)) +
                      "each containing at most {n} reads.".
                      format(n=reads_per_chunk))
        if extract_front_back_only:
            logging.debug("Extract exactly {k} bases from front" +
                          " and end of each read.".format(k=window_size))

        freader = FastaReader(reads_fn)
        chunkIndex = -1
        fwriter = None
        for i, read in enumerate(freader):
            if i % reads_per_chunk == 0:
                chunkIndex += 1
                if fwriter is not None:
                    fwriter.close()
                    fwriter = None
                fwriter = open(chunked_reads_fns[chunkIndex], 'w')
            rcseq = revcmp(read.sequence)
            if extract_front_back_only:
                fwriter.write(">{n}_front\n{s}\n>{n}_back\n{rcs}\n".format(
                              n=read.name, s=read.sequence[:window_size],
                              rcs=rcseq[:window_size]))
            else:
                fwriter.write(">{n}\n{s}\n".format(n=read.name,
                                                   s=read.sequence))

        if fwriter is not None:
            fwriter.close()

    def _startPhmmers(self, chunked_reads_fns, chunked_dom_fns,
                      out_dom_fn, primer_fn, pbmatrix_fn):
        """Run phmmers on chunked reads files in 'chunked_reads_fns' and
        generate chunked dom files as listed in 'chunked_dom_fns', finally
        concatenate dom files to 'out_dom_fn'."""
        logging.info("Start to launch phmmer on chunked reads.")
        jobs = []
        for reads_fn, domFN in zip(chunked_reads_fns, chunked_dom_fns):
            p = multiprocessing.Process(
                target=self._phmmer,
                args=(reads_fn, domFN, primer_fn, pbmatrix_fn))
            jobs.append((p, domFN))
            p.start()

        for p, domFN in jobs:
            p.join()
            cmd = "cat {0} >> {1}".format(real_upath(domFN),
                                          real_upath(out_dom_fn))
            _output, errCode, errMsg = backticks(cmd)
            if errCode != 0:
                raise ClassifierException(
                    "Error concatenating dom files: {e}".
                    format(e=str(errMsg)))

        self._cleanup(chunked_reads_fns)
        self._cleanup(chunked_dom_fns)

    def _phmmer(self, reads_fn, domFN, primer_fn, pbmaxtrixFN):
        """Invoke phmmer once."""
        cmd = "phmmer --cpu 1 --domtblout {d} --noali --domE 1 ".\
              format(d=real_upath(domFN)) + \
              "--mxfile {m} ".format(m=real_upath(pbmaxtrixFN)) + \
              "--popen 0.07 --pextend 0.07 {r} {p} > /dev/null".\
              format(r=real_upath(reads_fn), p=real_upath(primer_fn))
        logging.debug("Calling phmmer: {cmd}".format(cmd=cmd))
        _output, errCode, errMsg = backticks(cmd)
        if (errCode != 0):
            raise ClassifierException(
                "Error calling phmmer: {e}.".format(e=str(errMsg)))

    def _getBestFrontBackRecord(self, domFN, min_score):
        """Parses DOM output from phmmer and fill in best_of_front, best_of_back
           bestOf: sequence id ---> DOMRecord
        """
        logging.info("Get the best front & back primer hits.")
        # bestOf_ = {} # key: sid --> (score, primer
        best_of_front = defaultdict(lambda: None)
        best_of_back = defaultdict(lambda: None)

        reader = DOMReader(domFN)
        for r in reader:
            # allow missing adapter
            if r.sStart > 48 or r.pStart > 48:
                continue

            if r.score < min_score:
                continue

            # ex: sid m160213_091647_42134_c100957952550000001823213806221633_s1_p0/54497/ccs_front
            # ex: pid f_G11
            if r.sid.endswith('_front'):  # _front
                bestOf = best_of_front
                r.sid = r.sid[:-6]
            elif r.sid.endswith('_back'):  # _back
                bestOf = best_of_back
                r.sid = r.sid[:-5]
            else:
                raise ClassifierException(
                    "Unable to parse a read {r} in phmmer dom file {f}.".
                    format(r=r.sid, f=domFN))
            if r.sid not in bestOf:
                bestOf[r.sid] = {}
            if (r.pid not in bestOf[r.sid]) or \
                (bestOf[r.sid][r.pid].score < r.score):
                bestOf[r.sid][r.pid] = r
        return (best_of_front, best_of_back)

    def _getChimeraRecord(self, domFN, opts):
        """Parses phmmer DOM output from trimmed reads for chimera
           detection, return DOMRecord of suspicious chimeras, which
           have primer hits in the MIDDLE of the sequence.
        """
        logging.info("Identify chimera records from {f}.".
                     format(f=domFN))
        # sid --> list of DOMRecord with primer hits in the middle
        # of sequence.
        suspicous_hits = defaultdict(lambda: [])
        reader = DOMReader(domFN)
        for r in reader:
            # A hit has to be in the middle of sequence, and with
            # decent score.
            if r.sStart > opts.min_dist_from_end and \
               r.sEnd < r.sLen - opts.min_dist_from_end and \
               r.score > opts.min_score:
                suspicous_hits[r.sid].append(r)
        return suspicous_hits

    def _updateChimeraInfo(self, suspicous_hits, in_read_fn, out_nc_fn,
                           out_c_fn, primer_report_fn,
                           write_report_header=True):
        """
        in_read_fn --- a fasta of full-length reads or a fasta of
                       non-full-length reads.
        For each full-length read in in_read_fn FASTA file, detect whether
        it is chimeric or not, and write its annotation to
        primer_report_fn.
        Return:
            (num_nc, num_c, num_nc_bases, num_c_bases)
        """
        logging.debug("Update chimera info for reads in {f} ".
                      format(f=in_read_fn))
        logging.debug("Write primer report to {rpt}".
                      format(rpt=primer_report_fn))

        num_nc, num_c, num_nc_bases, num_c_bases = 0, 0, 0, 0
        with FastaReader(in_read_fn) as reader, \
                FastaWriter(out_nc_fn) as writer, \
                FastaWriter(out_c_fn) as writer_chimera, \
                open(primer_report_fn, 'w') as reporter:
            if write_report_header:
                reporter.write(ReadAnnotation.header(delimiter=",") + "\n")
            for r in reader:
                # e.g. r.name="movie/zmw/0_100_CCS fiveend=1;threeend=100;"
                readid = r.name.split()[0]
                annotation = ReadAnnotation.fromString(r.name,
                                                       ignore_polyA=self.ignore_polyA)
                if readid not in suspicous_hits:  # Non-chimeric reads
                    # Primer of a primer-trimmed read can not be None.
                    # assert(annotation.primer is not None)
                    annotation.chimera = 0
                    num_nc += 1
                    num_nc_bases += len(r.sequence)
                    writer.writeRecord(annotation.toAnnotation(),
                                       r.sequence)
                else:  # chimeric reads
                    annotation.chimera = 1
                    num_c += 1
                    num_c_bases += len(r.sequence)
                    writer_chimera.writeRecord(annotation.toAnnotation(),
                                               r.sequence)

                reporter.write(annotation.toReportRecord(delimitor=",") + "\n")
            return (num_nc, num_c, num_nc_bases, num_c_bases)

    def _findPolyA(self, seq, min_a_num=8, three_start=None):
        """
        Find poly A tail, which has at least 'min_a_num' A bases and at most
        two non-A bases in 3' of sequence. Return index of the very first base,
        if a polyA tail is found; otherwise, return -1.
        """
        polyA = 'A' * min_a_num
        offset = 50
        startEnd = three_start - offset if three_start is not None \
            else len(seq) - offset
        # search within the last <offset> bp
        i = seq.rfind(polyA, startEnd)
        if i > 0:
            nonA = 0
            # backtrace to the front of polyA, allowing only 2 max non-A
            while i >= 0:
                nonA += (seq[i] != 'A')
                if nonA > 2:
                    break
                i -= 1
            return i + 1
        else:
            return -1

    def _pickBestPrimerCombo(self, dFront, dBack, primer_names, min_score):
        """Pick up best primer combo.

        best_of_front/Back: {read_id: {primer_name:DOMRecord}}
        If the read is '+' strand: then front -> f_A1, back -> r_B2
        else: front -> r_B2, back -> f_A1
        Returns: (primer_combo, strand, DOM rec for f, DOM rec for r)
        """
        logging.debug("dFront={0}".format(dFront))
        logging.debug("dBack={0}".format(dBack))

        # simply pick the best one in front
        best_front = None
        if dFront is not None:
            for rec in dFront.itervalues():
                if rec.score >= min_score and (best_front is None or rec.score > best_front.score):
                    best_front = rec

        best_back = None
        if dBack is not None:
            for rec in dBack.itervalues():
                if rec.score >= min_score and (best_back is None or rec.score > best_back.score):
                    best_back = rec

        if best_front is None:
            if best_back is None:
                return ("NA", '?', None, None)
            elif best_back.pid.startswith('f_'):
                return (best_back.pid[2:]+"+NA", '-', best_back, None)
            else:
                return ("NA+"+best_back.pid[2:], '+', None, best_back)
        if best_back is None:
            if best_front.pid.startswith('f_'):
                return (best_front.pid[2:]+"+NA", '+', best_front, None)
            else:
                return ("NA+"+best_front.pid[2:], '-', None, best_front)

        if best_front.pid.startswith('f_'):
            if best_back.pid.startswith('r_'):
                return (best_front.pid[2:]+"+"+best_back.pid[2:], '+', best_front, best_back)
            else: # conflict! strand unresolved
                return (best_front.pid[2:]+"+"+best_back.pid[2:], '?', best_front, best_back)
        else:
            if best_back.pid.startswith('f_'):
                return (best_back.pid[2:]+"+"+best_front.pid[2:], '-', best_back, best_front)
            else:
                return (best_front.pid[2:]+"+"+best_back.pid[2:], '?', best_front, best_back)


    def _trimBarCode(self, reads_fn, out_fl_reads_fn, out_nfl_reads_fn,
                     primer_report_nfl_fn,
                     best_of_front, best_of_back, primer_names,
                     min_seq_len, min_score, change_read_id,
                     ignore_polyA, keep_primer):
        """Trim bar code from reads in 'reads_fn', annotate each read,
        indicating:
            whether its 5' primer, 3' primer and polyA tail are seen,
            start positions of its 5' primer, 3' primer and polyA tail,
            and primer info.
        , save non-full-length reads to 'out_nfl_reads_fn',
        , save full-length reads to 'out_fl_reads_fn', which can later be
        used in chimera detection
        , write primer info of nfl reads to _primer_report_nfl_fn.

        Note that chimera detection is not necessary for nfl reads, but
        is required for fl reads. So we only write primer info for nfl here
        and will write primer info for fl reads when chimera detection
        is done.

        best_of_front/Back: {read_id: {primer_name:DOMRecord}}
        min_seq_len: minimum length to output a read.
        min_score: minimum score to output a read.
        change_read_id: if True, change read ids to 'movie/zmw/start_end'.
        """
        logging.info("Trim bar code away from reads.")
        logging.debug("Writing full-length trimmed reads to {f}".
                      format(f=out_fl_reads_fn))
        logging.debug("Writing non-full-length trimmed reads to {f}".
                      format(f=out_nfl_reads_fn))
        logging.debug("Writing primer reports before chimera detection to {f}".
                      format(f=primer_report_nfl_fn))

        with FastaReader(reads_fn) as fareader, \
                FastaWriter(out_nfl_reads_fn) as nfl_fawriter, \
                FastaWriter(out_fl_reads_fn) as fl_fawriter, \
                open(primer_report_nfl_fn, 'w') as reporter:
            for read in fareader:
                self.summary.num_reads += 1  # number of ROI reads
                pbread = PBRead(read)
                logging.debug("Pick up best primer combo for {r}".
                              format(r=read.name))
                primerName, strand, fw, rc = self._pickBestPrimerCombo(
                    best_of_front[read.name], best_of_back[read.name],
                    primer_names, min_score)
                logging.debug("read={0}\n".format(read.name) +
                              "strand={0} fw={1} rc={2}".
                              format(strand, fw, rc))

                if (strand == '?') or (fw is None and rc is None):
                    # No primer seen in this sequence, classified
                    # as non-full-length
                    newName = pbread.name
                    if change_read_id:
                        newName = "{m}/{z}/{s1}_{e1}{isccs}".format(
                                  m=pbread.movie, z=pbread.zmw,
                                  s1=pbread.start, e1=pbread.end,
                                  isccs=("_CCS" if pbread.isCCS else ""))
                    annotation = ReadAnnotation(ID=newName, primer=primerName)
                    # Write reports of nfl reads
                    reporter.write(annotation.toReportRecord(delimitor=",") + "\n")
                    if len(read.sequence) >= min_seq_len:
                        # output non-full-length reads to nfl.trimmed.fasta
                        nfl_fawriter.writeRecord(annotation.toAnnotation(),
                                                 read.sequence)
                        self.summary.num_nfl += 1
                    else:
                        self.summary.num_filtered_short_reads += 1
                    continue
                seq = read.sequence if strand == "+" else revcmp(read.sequence)
                five_end, three_start = None, None
                if fw is not None:
                    five_end = fw.sEnd
                    self.summary.num_5_seen += 1
                if rc is not None:
                    three_start = len(seq) - rc.sEnd
                    self.summary.num_3_seen += 1

                s, e = pbread.start, pbread.end
                # Try to find polyA tail in read
                polyAPos = self._findPolyA(seq, three_start=three_start)
                if polyAPos >= 0 and not ignore_polyA:  # polyA found and not to ignore it
                    if not keep_primer:
                        seq = seq[:polyAPos]
                        e1 = s + polyAPos if strand == "+" else e - polyAPos
                    else:
                        e1 = e if strand == '+' else s
                    self.summary.num_polyA_seen += 1
                elif three_start is not None:  # polyA not found but 3' found
                    if not keep_primer:
                        seq = seq[:three_start]
                        e1 = s + three_start if strand == "+" else e - three_start
                    else:
                        e1 = e if strand == '+' else s
                else: # polyA not found and 3' not found
                    e1 =  e if strand == "+" else s

                if five_end is not None:
                    if not keep_primer:
                        seq = seq[five_end:]
                        s1 = s + five_end if strand == "+" else e - five_end
                    else:
                        s1 = s if strand == '+' else e
                else:
                    s1 = s if strand == "+" else e

                newName = pbread.name
                if change_read_id:
                    newName = "{m}/{z}/{s1}_{e1}{isccs}".format(
                        m=pbread.movie, z=pbread.zmw, s1=s1, e1=e1,
                        isccs=("_CCS" if pbread.isCCS else ""))
                # Create an annotation
                annotation = ReadAnnotation(ID=newName, strand=strand,
                                            fiveend=five_end, polyAend=polyAPos,
                                            threeend=three_start, primer=primerName,
                                            ignore_polyA=ignore_polyA)

                # Write reports for nfl reads
                if annotation.isFullLength is not True:
                    reporter.write(annotation.toReportRecord(delimitor=",") + "\n")

                if len(seq) >= min_seq_len:
                    if annotation.isFullLength is True:
                        # Write long full-length reads
                        fl_fawriter.writeRecord(annotation.toAnnotation(), seq)
                        self.summary.num_fl += 1
                    else:
                        # Write long non-full-length reads.
                        nfl_fawriter.writeRecord(annotation.toAnnotation(), seq)
                        self.summary.num_nfl += 1
                else:
                    self.summary.num_filtered_short_reads += 1

    def _validate_outputs(self, out_dir, out_all_reads_fn):
        """Validate and create output directory."""
        logging.info("Creating output directory {d}.".format(d=out_dir))
        if op.exists(out_dir):
            logging.warn("Output directory {d} already exists.".
                         format(d=out_dir))
        else:
            os.mkdir(out_dir)
        if op.exists(out_all_reads_fn):
            logging.warn("Existing output file {f} will be overwritten.".
                         format(f=out_all_reads_fn))

    def _cleanup(self, fileList):
        """Remove files in the list if they exist."""
        if fileList is None:
            return

        logging.debug("Clean up intermediate files: {fs}".
                      format(fs=",".join(fileList)))
        for f in fileList:
            if op.exists(f):
                os.remove(f)

    def runPrimerTrimmer(self):
        """Run PHMMER to identify barcodes and trim them away.
        (1) create forward/reverse primers
        (2) copy input with just the first/last k bases
        (3) run phmmer
        (4) parse phmmer DOM output, trim barcodes and output summary
        """
        logging.info("Start to find and trim 3'/5' primers and polyAs.")
        # Sanity check input primers and create forward/reverse primers
        # for primer detection.
        primer_names = self._processPrimers(
            primer_fn_forward=self.primer_fn_forward,
            primer_fn_reverse=self.primer_fn_reverse,
            window_size=self.chimera_detection_opts.primer_search_window,
            primer_out_fn=self.primer_front_back_fn,
            revcmp_primers=False)

        logging.info("reuse_dom = {0}".format(self.reuse_dom))
        if op.exists(self.out_front_back_dom_fn) and self.reuse_dom:
            logging.warn("Primer detection output already exists. Parsing {0}".
                         format(self.out_front_back_dom_fn))
        else:
            # Split reads in reads_fn into smaller chunks.
            num_chunks = max(min(self.cpus, self.numReads), 1)
            reads_per_chunk = int(math.ceil(self.numReads / (float(num_chunks))))
            num_chunks = int(math.ceil(self.numReads / float(reads_per_chunk)))

            logging.debug("Split reads into {n} chunks".format(n=num_chunks))
            # Divide input reads into smaller chunks and extract only
            # the front and the end segment from each read.
            self.chunked_front_back_reads_fns = generateChunkedFN(self.out_dir,
                                                                  "in.front_end.fa_split", num_chunks)

            # Dom output of phmmer for the above front/end sequences.
            self.chunked_front_back_dom_fns = generateChunkedFN(self.out_dir,
                                                                "out.front_end.hmmer_split", num_chunks)

            # Split reads within 'reads_fn' into 'num_chunks' chunks, and only
            # extract the front and end segment from each read.
            window_size = self.chimera_detection_opts.primer_search_window
            self._chunkReads(reads_fn=self.reads_fn,
                             reads_per_chunk=reads_per_chunk,
                             chunked_reads_fns=self.chunked_front_back_reads_fns,
                             extract_front_back_only=True,
                             window_size=window_size)

            # Start n='num_chunks' phmmer.
            self._startPhmmers(
                chunked_reads_fns=self.chunked_front_back_reads_fns,
                chunked_dom_fns=self.chunked_front_back_dom_fns,
                out_dom_fn=self.out_front_back_dom_fn,
                primer_fn=self.primer_front_back_fn,
                pbmatrix_fn=self.pbmatrix_fn)

        # Parse dome file, and return dictionary of front & back.
        best_of_front, best_of_back = self._getBestFrontBackRecord(
            self.out_front_back_dom_fn, self.chimera_detection_opts.min_score)


        # Trim bar code away
        self._trimBarCode(reads_fn=self.reads_fn,
                          out_fl_reads_fn=self._trimmed_fl_reads_fn,
                          out_nfl_reads_fn=self._trimmed_nfl_reads_fn,
                          primer_report_nfl_fn=self._primer_report_nfl_fn,
                          best_of_front=best_of_front,
                          best_of_back=best_of_back,
                          primer_names=primer_names,
                          min_seq_len=self.chimera_detection_opts.min_seq_len,
                          min_score=self.chimera_detection_opts.min_score,
                          change_read_id=self.change_read_id,
                          ignore_polyA=self.ignore_polyA,
                          keep_primer=self.keep_primer)

        # Clean intemediate files: chunked reads files and chunked dom files.
        self._cleanup(self.chunked_front_back_reads_fns)
        self._cleanup(self.chunked_front_back_dom_fns)
        logging.info("Done with finding and trimming primers and polyAs.")

    def _detect_chimera(self, in_fasta, out_nc_fasta, out_c_fasta,
                        primer_report_fn, out_dom, num_reads, job_name):
        """Detect chimeric reads from in_fasta, call phmmer to generate a
        dom file (out_dom), save non-chimeric reads to out_nc_fasta and
        chimeric reads to out_c_fasta.
            in_fasta --- either a fasta of trimmed fl reads, or a fasta of
                         trimmed nfl reads.
            out_nc_fasta --- an output fasta of non-chimeric reads
            out_c_fasta --- an output fasta of chimeric reads
            primer_report_fn --- an output primer report
            out_dom --- phmmer output
            num_reads --- number of reads in in_fasta
            job_name --- either 'fl' or 'nfl'
        Return:
            (num_nc, num_c, num_nc_bases, num_c_bases)
        """
        if op.exists(out_dom) and self.reuse_dom:
            logging.warn("Chimera detection output already exists. Parse {o}.".
                         format(o=out_dom))
        else:
            num_chunks = max(min(num_reads, self.cpus), 1)
            reads_per_chunk = int(math.ceil(num_reads / float(num_chunks)))
            num_chunks = int(math.ceil(num_reads / float(reads_per_chunk)))

            chunked_reads_fns = generateChunkedFN(self.out_dir,
                                                  "in.{n}.trimmed.fa_split".format(n=job_name), num_chunks)

            chunked_dom_fns = generateChunkedFN(self.out_dir,
                                                "out.{n}.trimmed.hmmer_split".format(n=job_name), num_chunks)

            self._chunkReads(reads_fn=in_fasta,
                             reads_per_chunk=reads_per_chunk,
                             chunked_reads_fns=chunked_reads_fns,
                             extract_front_back_only=False)

            self._startPhmmers(chunked_reads_fns=chunked_reads_fns,
                               chunked_dom_fns=chunked_dom_fns,
                               out_dom_fn=out_dom,
                               primer_fn=self.primer_chimera_fn,
                               pbmatrix_fn=self.pbmatrix_fn)

        suspicous_hits = self._getChimeraRecord(out_dom,
                                                self.chimera_detection_opts)

        # Update chimera information
        (num_nc, num_c, num_nc_bases, num_c_bases) = \
            self._updateChimeraInfo(suspicous_hits=suspicous_hits,
                                    in_read_fn=in_fasta,
                                    out_nc_fn=out_nc_fasta,
                                    out_c_fn=out_c_fasta,
                                    primer_report_fn=primer_report_fn,
                                    write_report_header=True if job_name == "fl" else False)

        return (num_nc, num_c, num_nc_bases, num_c_bases)

    def runChimeraDetector(self):
        """Call chimera detection on full-length reads, and non-full-length
        reads if required."""
        # Create forward/reverse primers for chimera detection.
        self._processPrimers(
            primer_fn_forward=self.primer_fn_forward,
            primer_fn_reverse=self.primer_fn_reverse,
            window_size=self.chimera_detection_opts.primer_search_window,
            primer_out_fn=self.primer_chimera_fn,
            revcmp_primers=True)

        # Detect chimeras among full-length reads, separate flnc reads and
        # flc reads.
        logging.info("Detect chimeric reads from trimmed full-length reads.")
        (self.summary.num_flnc, self.summary.num_flc,
         self.summary.num_flnc_bases, _x) = \
            self._detect_chimera(in_fasta=self._trimmed_fl_reads_fn,
                                 out_nc_fasta=self.out_flnc_fn,
                                 out_c_fasta=self.out_flc_fn,
                                 primer_report_fn=self._primer_report_fl_fn,
                                 out_dom=self.out_trimmed_fl_dom_fn,
                                 num_reads=self.summary.num_fl,
                                 job_name="fl")
        assert(self.summary.num_fl == self.summary.num_flnc +
               self.summary.num_flc)
        logging.info("Done with chimera detection on trimmed full-length " +
                     "reads.")

        # Detect chimeras among non-full-length reads if required, separate
        # nflnc reads and nflc reads, rewrite self.primer_report_nfl_fn.
        if self.chimera_detection_opts.detect_chimera_nfl is True:
            logging.info("Detect chimeric reads from trimmed non-full-length " +
                         "reads.")
            (self.summary.num_nflnc, self.summary.num_nflc, _x, _y) = \
                self._detect_chimera(in_fasta=self._trimmed_nfl_reads_fn,
                                     out_nc_fasta=self.out_nflnc_fn,
                                     out_c_fasta=self.out_nflc_fn,
                                     primer_report_fn=self._primer_report_nfl_fn,
                                     out_dom=self.out_trimmed_nfl_dom_fn,
                                     num_reads=self.summary.num_nfl,
                                     job_name="nfl")
            assert(self.summary.num_nfl == self.summary.num_nflnc +
                   self.summary.num_nflc)
            logging.info("Done with chimera detection on trimmed " +
                         "non-full-length reads.")

            # Concatenate out_nflnc_fn and out_nflc_fn as out_nfl_fn
            cat_files(src=[self.out_nflnc_fn, self.out_nflc_fn],
                      dst=self.out_nfl_fn)
            # Concatenate out_flnc and out_nflnc to make out_all_reads_fn
            cat_files(src=[self.out_flnc_fn, self.out_nflnc_fn],
                      dst=self.out_all_reads_fn)

        else:
            # Soft link _trimmed_nfl_reads_fn as out_nfl_fn
            ln(self._trimmed_nfl_reads_fn, self.out_nfl_fn)
            # Concatenate out_flnc and out_nfl to make out_all_reads_fn
            cat_files(src=[self.out_flnc_fn, self.out_nfl_fn],
                      dst=self.out_all_reads_fn)

        # primer info of fl/nfl reads reported to _primer_report_fl_fn
        # and _primer_report_nfl_fn, concatenate them in order to make
        # a full report: primer_report_fn.
        cat_files(src=[self._primer_report_fl_fn, self._primer_report_nfl_fn],
                  dst=self.primer_report_fn)

        # Delete intermediate files.
        self._cleanup([self._primer_report_nfl_fn,
                       self._primer_report_fl_fn])

    def run(self):
        """Classify/annotate reads according to 5' primer seen,
        3' primer seen, polyA seen, chimera (concatenation of two
        or multiple transcripts with primers seen in the middle of
        a read)
        (1) Create and validate input/output
        (2) Check phmmer is runnable
        (3) Find primers using phmmer and trim away primers and polyAs
        (4) Detect chimeras from trimmed reads
        """
        # Validate input files and required data files.
        self._validate_inputs(self.reads_fn, self.primer_fn_forward, self.primer_fn_reverse, self.pbmatrix_fn)

        # Validate and create output dir.
        self._validate_outputs(self.out_dir, self.out_all_reads_fn)

        # Sanity check phmmer can be called successfully.
        self._checkPhmmer()

        # Find and trim primers and polyAs.
        self.runPrimerTrimmer()

        # Check whether no fl reads detected.
        no_flnc_errMsg = "No full-length non-chimeric reads detected."
        if self.summary.num_fl == 0:
            logging.error(no_flnc_errMsg)
            raise ClassifierException(no_flnc_errMsg)

        # Detect chimeras and generate primer reports.
        #print >> sys.stderr, "TURNING OFF CHIMERA DETECTOR FOR NOW"
        self.runChimeraDetector()

        try:
            # Write summary.
            logging.info("Writing report to {f}".format(f=self.summary_fn))
            self.summary.write(self.summary_fn)
        except ZeroDivisionError:
            logging.error(no_flnc_errMsg)
            raise ClassifierException(no_flnc_errMsg)

        return 0


if __name__ == "__main__":
    obj = Classifier()
    obj.run()
