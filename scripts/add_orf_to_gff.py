from __future__ import print_function

from optparse import OptionParser
from Bio import SeqIO
import re
import sys

'''
input : 
1) gff file, e.g. :
/net/eichler/vol21/projects/human_primate_brain_ips_transcriptomes/nobackups/pacbio_data/cluster_by_region/test3_24Oct2016/merged_trimmed/adult_brain/mapping_cluster_output/collapse/collapsed/FAM72A.all_quivered_hq.100_30_0.99.sorted.bam.yes5merge.collapsed.gff

2) fasta file with angel results, e.g. :
/net/eichler/vol21/projects/human_primate_brain_ips_transcriptomes/nobackups/pacbio_data/cluster_by_region/test3_24Oct2016/merged_trimmed/adult_brain/mapping_cluster_output/collapse/collapsed/fasta_of_gff/FAM72A.all_quivered_hq.100_30_0.99.sorted.bam.yes5merge.collapsed.gff.fasta.final.cds

3) The original fasta file on which angel was run on. This file is needed in order to get the length of each such original sequence so that "-" strand ORFs can be analysed.


output:
gff with orf

method :
1. read fasta file ids to a dictionary name2cds .  for example: PB.1.3 -> (870,1319)
2. go over the gff input file, for every exon that is within the range of the cds add a CDS line
'''

debug =True

# returns pb_id, strand, start, end
def getAngelNameAndCdsRegion(angel_fasta_id) :
    regionM = re.compile('.*(PB\.\d+\.\d+).* pos:(\d+)-(\d+)').match(angel_fasta_id)
    if regionM is None :
        print("Malformatted region " + angel_fasta_id)
        sys.exit(0)
    m = regionM.groups()
    return m[0], int(m[1]), int(m[2])


def getAngelName(fasta_id):
    regionM = re.compile('.*(PB\.\d+\.\d+).*').match(fasta_id)
    if regionM is None :
        print("Malformatted region " + fasta_id)
        sys.exit(0)
    m = regionM.groups()
    return m[0]


def ParseAngelGffLine(angel_gff_line):
    lineS = line.strip().split("\t")
    line_type = lineS[2]
    chrom = lineS[0]
    start = int(lineS[3])
    end = int(lineS[4])
    strand = lineS[6]
    name = lineS[8]
    angel_nameM = re.compile('.*(PB\.\d+\.\d+).*').match(name) 
    if angel_nameM is None :
        print("Malformatted gff name", name)
        sys.exit(0)
    angel_name = angel_nameM.groups()[0]
    return (line_type,chrom,start,end,strand,name,angel_name,lineS)


def parseAngelFastaIds(orf_fasta_infile) :
    name2cds = {}
    input_handle = open(orf_fasta_infile, 'r') 
    for seq in SeqIO.parse(input_handle, "fasta") :
#        print seq.id
#        print seq.description
        name,cds_start,cds_end = getAngelNameAndCdsRegion(seq.description)
        if name in name2cds :
            print("ERROR, name was already seen before:", name, ". Using only the longest ORF.")
            if cds_end - cds_start < name2cds[name][1] - name2cds[name][0] :
                continue
        name2cds[name] = (cds_start, cds_end)

    input_handle.close()
    return name2cds
        
def parseFastaLengths(transcripts_fasta_infile) :
    name2len = {}
    input_handle = open(transcripts_fasta_infile, 'r') 
    for seq in SeqIO.parse(input_handle, "fasta") :
        name = getAngelName(seq.description)
        if name in name2len :
            print("ERROR, name was already seen before", name)
        name2len[name] = len(seq.seq)
    
#    print name2len
    input_handle.close()
    return name2len

if __name__=="__main__" :
    opts=OptionParser()
    opts.add_option('', '--gff_infile', dest='gff_infile', default=None) 
    opts.add_option('', '--orf_fasta_infile', dest='orf_fasta_infile', default=None) 
    opts.add_option('', '--transcripts_fasta_infile', dest='transcripts_fasta_infile', default=None) 

    opts.add_option('', '--outfile', dest='outfile', default=None) 
    (o,args)=opts.parse_args()    


    print("*** reading orf_fasta_infile ***")
    name2cds = parseAngelFastaIds(o.orf_fasta_infile)
#    print name2cds
    
    print("*** reading transcripts_fasta_infile ***")
    name2len = parseFastaLengths(o.transcripts_fasta_infile)
    
    cds_start = -1
    cds_end = -1
    cds_len = -1

    print("*** reading gff ***")

    gff_inF = open(o.gff_infile, 'r')
    outF = open(o.outfile, 'w')
    for line in gff_inF :
        outF.write(line)
        line_type,chrom,start,end,strand,name,angel_name,lineS = ParseAngelGffLine(line)
        if line_type == 'transcript' :
            if angel_name not in name2cds :
                print("skipping transcript, name in gff does not exist in input fasta", angel_name)
                cds_start = -1
                cds_end = -1
                continue
            cds_start, cds_end = name2cds[angel_name]
            cds_len = cds_end - cds_start
            if strand == "-" :
                # convert from seq position to real coordinates based on transcript start
                genomic_cds_start = start + name2len[angel_name] - cds_end
            else :
                genomic_cds_start = -1
            cds_end = -1
                
        if line_type == 'exon' :
            exon_length = end - start + 1
            if cds_start <= exon_length:
                if strand == "+":
                    # convert from seq position to real coordinates based on transcript start
                    genomic_cds_start = cds_start + start - 1
                exon_cds_start = max(start, genomic_cds_start)
                cds_end = exon_cds_start + cds_len
                exon_cds_end = min(end, cds_end)
                if exon_cds_start < exon_cds_end:
                    lineS[3] = str(exon_cds_start)
                    lineS[4] = str(exon_cds_end)
                    lineS[2] = "CDS"
                    outF.write("\t".join(lineS) + "\n")
                    exon_cds_len = exon_cds_end - exon_cds_start + 1
                    cds_len -= exon_cds_len
            cds_start = max(0, cds_start - exon_length)

    
    outF.close()
    gff_inF.close()
            
    

