import os
import sys
import argparse




    
def parseArgs():
    # dir of source code
    srcpath=sys.argv[0]
    srcdir=os.path.dirname(os.path.abspath(srcpath))
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
                        help="The path to the TE bed file, format must be chromosome\tstart\tend\tsubfamily\tClass/family."
                             "example: chr1 107 155 AluYk3 SINE/Alu")
    #
    args = parser.parse_args()
    args.srcdir=srcdir
    return args




def main():
    args = parseArgs()
    ID = args.ID
    fastq = os.path.abspath(args.input)
    outdir = os.path.abspath(args.outdir)
    ref = os.path.abspath(args.reference)
    gtf = os.path.abspath(args.gtf)
    TE_bed = os.path.abspath(args.TEbed)
    srcdir = args.srcdir
    outdir = '%s/%s' % (outdir, ID)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    os.chdir(outdir) 
    cmd = 'python {srcdir}/scripts/preparing_files.py \
        -id {ID} -i {fastq} -o {outdir} -r {ref} -g {gtf} -b {TEbed}'.format(
        srcdir=srcdir, ID=ID, fastq=fastq, outdir=outdir, ref=ref, gtf=gtf, TEbed=TE_bed)
    os.system(cmd)
    cmd = 'python {srcdir}/scripts/parseAdapter.py -i {fastq} -o {outdir}'.format(srcdir=srcdir, fastq=fastq, outdir=outdir)
    os.system(cmd)
    cmd = 'Rscript {srcdir}/scripts/Target_efficiency.r {ID} {gtf}'.format(srcdir=srcdir, ID=ID, gtf=gtf)
    os.system(cmd)
    cmd = 'python {srcdir}/scripts/cutSites.py -id {ID}'.format(srcdir=srcdir, ID=ID)
    os.system(cmd)

    

    
if __name__=="__main__":
    main()