from optparse import OptionParser
from Bio import SeqIO
import re
import sys

'''
 Takes as input a 
1) a gff created by the collapse pipeline
 For example:
chr1	PacBio	transcript	206186183	206204423	.	-	.	gene_id "PB.1"; transcript_id "PB.1.1";
chr1	PacBio	exon	206186183	206187373	.	-	.	gene_id "PB.1"; transcript_id "PB.1.1";
chr1	PacBio	exon	206195752	206195876	.	-	.	gene_id "PB.1"; transcript_id "PB.1.1";
chr1	PacBio	exon	206199807	206199884	.	-	.	gene_id "PB.1"; transcript_id "PB.1.1";
chr1	PacBio	exon	206203844	206204061	.	-	.	gene_id "PB.1"; transcript_id "PB.1.1";
chr1	PacBio	exon	206204252	206204423	.	-	.	gene_id "PB.1"; transcript_id "PB.1.1";
chr1	PacBio	transcript	206186183	206204584	.	-	.	gene_id "PB.1"; transcript_id "PB.1.2";
chr1	PacBio	exon	206186183	206187373	.	-	.	gene_id "PB.1"; transcript_id "PB.1.2";
chr1	PacBio	exon	206195752	206195876	.	-	.	gene_id "PB.1"; transcript_id "PB.1.2";
chr1	PacBio	exon	206199807	206199884	.	-	.	gene_id "PB.1"; transcript_id "PB.1.2";
chr1	PacBio	exon	206204252	206204584	.	-	.	gene_id "PB.1"; transcript_id "PB.1.2";

2) a fasta file created by applying bedtools getfasta on that same file


Output: a. a fasta file with correct sequences

'''


regionRe = re.compile('(.*):(\d+)\-(\d+)')

def ParseRegion(regionStr):
    regionM = regionRe.match(regionStr)
    if (regionM is None):
        print "Malformatted region " + regionStr
        sys.exit(0)
    m = regionM.groups()
    start = int(m[1])
    end   = int(m[2])
    chrom = m[0]
    return (chrom,start,end)


def ParseGffLine(gff_line):
    lineS = line.strip().split("\t")
    line_type = lineS[2]
    chrom = lineS[0]
    start = int(lineS[3])-1
    end = int(lineS[4])
    strand = lineS[6]
    name = lineS[8]
    return (line_type,chrom,start,end,strand,name)

if __name__=="__main__" :
    opts=OptionParser()
    opts.add_option('', '--gff_infile', dest='gff_infile', default=None) 
    opts.add_option('', '--fasta_infile', dest='fasta_infile', default=None) 
    opts.add_option('', '--outfile', dest='outfile', default=None) 
    opts.add_option('', '--ignore_strand', action="store_true", dest='ignore_strand', default=False)
    (o,args)=opts.parse_args()    

    gff_inF = open(o.gff_infile, 'r')
    fasta_generator = SeqIO.parse(o.fasta_infile, "fasta")

    output_handle = open(o.outfile, "w")
    chrom = ""
    start = -1
    end = -1
    transcript_strand = ""
    first_flag = True
    for line in gff_inF :
        line_type,chrom,start,end,strand,name = ParseGffLine(line)
        fasta_seq = fasta_generator.next()
        fasta_chrom, fasta_start, fasta_end = ParseRegion(fasta_seq.id)
        
        # making sure the fasta seq matches the gff line
        if (fasta_chrom != chrom) or (fasta_start != start) or (fasta_end != end) :
            print "ERROR", fasta_seq.id, line
            sys.exit(0)

        if line_type == 'transcript' :
            if first_flag :
                first_flag = False
            else:
                if (not o.ignore_strand) and (transcript_strand == "-") :
                    transcript_seq.seq = transcript_seq.seq.reverse_complement() #[::-1]
                SeqIO.write([transcript_seq], output_handle, "fasta")
            transcript_seq = fasta_seq
            transcript_seq.seq = ""
            transcript_strand = strand
#            transcript_seq.id += " " + name
#            print name
            transcript_seq.id = name.replace(" ","_") + transcript_seq.id
        elif line_type == 'exon' :
            transcript_seq.seq += fasta_seq.seq
        else :
            print "ERROR unrecognized line_type", line_type
            sys.exit(0)

    if len(transcript_seq.seq) > 0 :
        SeqIO.write([transcript_seq], output_handle, "fasta")
    gff_inF.close()
    output_handle.close()
            


            


            
            



