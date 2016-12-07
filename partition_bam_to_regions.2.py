import sys # argv
import os
import pysam # samtools for python
from optparse import OptionParser


def parse_read(read) :
    return read.reference_id, read.reference_start, read.reference_end


if __name__=="__main__" :
    opts=OptionParser()
    opts.add_option('', '--inFile', dest='inFileName', default=None)  # sam\bam file should be sorted 
    opts.add_option('', '--outDir', dest='outDir', default=None)
    opts.add_option('', '--outBedFile', dest='outBedFile', default=None)
    opts.add_option('', '--min_mapq', dest='min_mapq', type ='int', default=0)  
    

    (o,args)=opts.parse_args()

    bamfile = pysam.Samfile(o.inFileName, "rb") # r=read, b=bamfile
    if o.outDir :
        if not o.outDir.endswith("/"):
            o.outDir += "/"
        if not os.path.exists(o.outDir) :
            os.makedirs(o.outDir)
        tmpOutFileName = o.outDir + "cur_region.bam"
        nonMappedOutFileName = o.outDir + "non_mapped.bam"
        outF = pysam.Samfile(tmpOutFileName , "wb", template=bamfile)
        non_mapped_outF = pysam.Samfile(nonMappedOutFileName , "wb", header=bamfile.header, template=bamfile)

    if o.outBedFile :
        outBed = open(o.outBedFile, "w")

    firstLineFlag = True
    region_start = None
    region_end = None
    region_chr = None

    for read in bamfile :
        # read not mapped
        if read.reference_id == -1 :
            if o.outDir :
                non_mapped_outF.write(read)
            continue

        # read mapq < cutoff
        if read.mapq < o.min_mapq :
            continue

        # first read in the file
        if firstLineFlag :
            region_chrom_id, region_start, region_end = parse_read(read)
            firstLineFlag = False
            region_count = 1
            if o.outDir :
                outF.write(read)
            continue
        
        if (read.reference_id == region_chrom_id) and (read.reference_start <= region_end) :
            region_end = read.reference_end
            region_count += 1
            if o.outDir :
                outF.write(read)
        else :
            regionName = [bamfile.getrname(region_chrom_id), str(region_start), str(region_end)]
            if o.outDir :
                regionFileName = "_".join(regionName) + ".bam"
                outF.close()
                os.rename(tmpOutFileName, o.outDir + regionFileName)
                outF = pysam.Samfile(tmpOutFileName , "wb", header=bamfile.header, template=bamfile)
                outF.write(read)

            if o.outBedFile :
                outBed.write("\t".join(regionName) + "\t" + str(region_count) + "\n")

            region_count = 1
            
            region_chrom_id, region_start, region_end = parse_read(read)


    regionName = [bamfile.getrname(region_chrom_id), str(region_start), str(region_end)]
    if o.outDir :
        outF.close()
        regionFileName = "_".join(regionName) + ".bam"
        os.rename(tmpOutFileName, o.outDir + regionFileName)
        non_mapped_outF.close()
               
    if o.outBedFile :
        outBed.write("\t".join(regionName) + "\t" + str(region_count) + "\n")
        outBed.close()
            
        
