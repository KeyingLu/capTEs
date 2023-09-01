import os
import re
import sys
import argparse
# srcpath=sys.argv[0]
# srcdir=os.path.dirname(os.path.dirname(os.path.abspath(srcpath)))
# sys.path.append('%s/pyPackages' % srcdir)
import pandas as pd
import numpy as np



    
def parseArgs():
    # dir of source code
    srcpath=sys.argv[0]
    srcdir=os.path.dirname(os.path.abspath(srcpath))
    softwaredir = os.path.dirname(srcdir)
    #
    parser = argparse.ArgumentParser(description=' preprocessing pipeline with capTEs data')
    parser.add_argument('-id', '--ID', type=str, required=True, help="Sample name")
    parser.add_argument('-i','--input', type=str, required=True, help="input fastq file")
    parser.add_argument('-o','--outdir',type=os.path.abspath, required=True, help="outdir path")
    parser.add_argument('-r','--reference',type=str, required=True, 
                        help="The path to the reference file. File must be indexed by samtools faidx")
    parser.add_argument('-g','--gtf',type=str, required=True, 
                        help="The path to the gtf file")
    parser.add_argument('-b','--TEbed',type=str, required=True, 
                        help="The path to the TE bed file, format must be: chromosome\tstart\tend\tsubfamily\tClass/family\n"
                             "example: chr1 107 155 AluYk3 SINE/Alu")
    #
    args = parser.parse_args()
    args.softwaredir = softwaredir
    return args


def trimming_and_filtering(ID, fastq, softwaredir):
    out_fq = "%s.trimming.Q7.L300.fq.gz" % ID
    cmd = "porechop -i {fastq} -t 40 | \
        {softwaredir}/citeTools/NanoFilt -q 7 -l 300 --headcrop 0 --tailcrop 0 | " \
          "gzip > {out_fq}".format(fastq=fastq, out_fq=out_fq)
    os.system(cmd)
    return(out_fq)

    
def minimap2_run(ID, fq, ref, softwaredir):
    bam = "%s.sort.bam" % ID
    sam = bam[:-3]+ "sam"
    cmd = "minimap2 \
          -ax splice --MD {ref} {fq} > {sam}".format(ref=ref, fq=fq, sam=sam)
    os.system(cmd)
    cmd = "samtools sort -@ 100 {sam} -o {bam}".format(
          sam=sam, bam=bam)
    os.system(cmd)
    cmd = "samtools index {bam}".format(bam=bam)
    os.system(cmd)
    os.system("rm -rf %s" % sam)
    return(bam)
    
    
def minimap2_transcriptom_run(ID, fq, ref, gtf, softwaredir, transcript_fa='NULL'):
    bam = "%s.transcriptom.sort.bam" % ID
    # prepare transcript fa without intron
    if transcript_fa == 'NULL':
        transcript_fa = "transcript.fa"
        cmd = "gffread -w {outFa} -g {refFa} {gtf}".format(outFa=transcript_fa, refFa=ref, gtf=gtf)
        os.system(cmd)
    # mapping
    sam = bam[:-3] + "sam"
    cmd = "minimap2 \
        -ax map-ont --MD {ref} {fq} > {sam}".format(ref=transcript_fa, fq=fq, sam=sam)
    os.system(cmd)
    cmd = "samtools sort -@ 100 {sam} -o {bam}".format(
        sam=sam, bam=bam)
    os.system(cmd)
    cmd = "samtools index {bam}".format(bam=bam)
    os.system(cmd)
    os.system("rm -rf %s" % sam)
    # unique bam
    uniqueBam = bam[:-3] + 'unique.bam'
    cmd = "samtools view -h  %s | \
        awk '$1~/@/ || $5==60 {print $0}' | samtools view -b  > %s" % (bam, uniqueBam)
    os.system(cmd)
    bed = bam[:-3] + 'unique.bed'
    cmd = "bedtools bamtobed -i {bam} > {bed}".format(
        bam=uniqueBam, bed=bed)
    os.system(cmd)
    return(uniqueBam)


def RM_run(fq, softwaredir):
    fa = fq[:-5] + 'fa'
    os.system("%s/citeTools/seqkit fq2fa %s -j 50 -o %s" % (softwaredir, fq, fa))
    cmd = "{softwaredir}/citeTools/RepeatMasker/RepeatMasker \
        -engine rmblast \
        -trf_prgm {softwaredir}/citeTools/TRF-4.09.1/build/src \
        -libdir {softwaredir}/citeTools/RepeatMasker/Libraries \
        -rmblast_dir {softwaredir}/citeTools/rmblast-2.11.0/bin \
        -pa 50 -a -species human {fa} -dir ./".format(
        softwaredir=softwaredir, fa=fa)
    print(cmd)
    os.system(cmd) 


    
def produce_TE_gtf(TE_bed):
    df = pd.read_csv(TE_bed, sep="\t")
    df.columns = ["chr", "start", "end", "repeat_name", "class_id", "strand"]
    repeat_names = list(df["repeat_name"].value_counts().index)
    repeat_dup = pd.DataFrame(np.repeat(1, len(repeat_names)))
    repeat_dup.index = repeat_names
    repeat_dup.columns = ["Times"]
    out_gtf = "TE.gtf"
    with open(out_gtf, "w") as fout:
        for index, row in df.iterrows():
            chrom = row["chr"]
            start = str(row["start"])
            end = str(row["end"])
            strand = row["strand"]
            if strand == "C":
                strand = "-" 
            gene_id = row["repeat_name"]
            transcript_id = gene_id + "_dup" + str(repeat_dup.loc[gene_id]["Times"])
            repeat_dup.loc[gene_id]["Times"]  += 1
            family_id = row["class_id"].split("/")[1]
            class_id = row["class_id"].split("/")[0]
            anno = 'gene_id "%s"; transcript_id "%s"; family_id "%s"; class_id "%s";' % (
                transcript_id, transcript_id, family_id, class_id)
            string = [chrom, "hg38_rmsk", "exon", start, end, ".", strand, ".", anno]
            fout.write("\t".join(string) + "\n")
    return(out_gtf)

 

def featureCounts_run(ID, bam, gtf, softwaredir):
    out = "%s.TE.featureCounts" % ID
    cmd = "featureCounts -a {gtf} \
            -o {out} {bam} \
            -F GTF -t exon -g transcript_id \
            -f -O --minOverlap 20 -M --fraction \
            -s 0 -T 50 -L --maxMOp 10 \
            --tmpDir ./".format(gtf=gtf, out=out, bam=bam)
    os.system(cmd)

    
def main():
    args=parseArgs()
    ID = args.ID
    fastq = os.path.abspath(args.input)
    outdir = os.path.abspath(args.outdir)
    ref = os.path.abspath(args.reference)
    gtf = os.path.abspath(args.gtf)
    TE_bed = os.path.abspath(args.TEbed) # eg: chr1     107     155     AluYk3  SINE/Alu  
    softwaredir = args.softwaredir
    os.chdir(outdir) 
    clean_fq = trimming_and_filtering(ID, fastq, softwaredir)
    bam1 = minimap2_run(ID, clean_fq, ref, softwaredir)
    bam2 = minimap2_transcriptom_run(ID, clean_fq, ref, gtf, softwaredir)
    RM_run(clean_fq, softwaredir)
    TE_gtf = produce_TE_gtf(TE_bed)
    featureCounts_run(ID, bam1, TE_gtf, softwaredir)
    

    
if __name__=="__main__":
    main()