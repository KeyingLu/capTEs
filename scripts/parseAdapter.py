# parseAdapter.py
import os
import re
import sys
import argparse

def parseArgs():
    # dir of source code
    srcpath=sys.argv[0]
    srcdir=os.path.dirname(os.path.abspath(srcpath))
    softwaredir = os.path.dirname(srcdir)
    # parser
    parser = argparse.ArgumentParser(description='the content of adapters in origin fastqs')
    parser.add_argument('-i', '--input', type=str, required=False,
            help="input fastq file")
    parser.add_argument('-o','--outdir',type=os.path.abspath, required=False,
            help="outdir path")
    args = parser.parse_args()
    args.softwaredir = softwaredir
    return args
    
    
def get_trimmed_log(fastq, outDir, softwaredir):
    os.chdir(outDir)
    cmd = "porechop \
            -i {fastq} \
            --check_reads 10000 -t 30 \
            --verbosity 3 \
            --end_threshold 75 \
            -o tmp.fastq > score75.log".format(fastq=fastq)
    os.system(cmd)
    cmd = "porechop \
            -i {fastq} \
            --check_reads 10000 -t 30 \
            --verbosity 3 \
            --end_threshold 95 \
            -o tmp.fastq > score95.log".format(fastq=fastq)
    os.system(cmd)
    os.system('rm -rf tmp.fastq')
    return(['%s/score75.log' % outDir, '%s/score95.log' % outDir])

    
def porechop_trimmed_info(log):
    outfile = log + ".txt"
    fin = open(log, "r")
    lines = fin.readlines()
    fin.close()
    ###
    readName_index = []
    for line in lines:
        if re.search("had adapters trimmed from their start", line):
            break
        if re.search(r'runid', line):
            index = lines.index(line)
            readName_index.append(index)
    ###
    read_trimmed_info = {} # readname, start_Smart-seq3, start_SQK-NSK007, end_Smart-seq3, end_SQK-NSK007
    for i in range(0, len(readName_index)):
        index = readName_index[i]  
        if i != (len(readName_index)-1):
            index_latter = readName_index[i+1]
            trimmed_info = lines[index:index_latter]
        else:
            trimmed_info=lines[index:(index+15)]
        ####
        # initialization
        read_name = trimmed_info[0].split()[0]
        read_trimmed_info[read_name] = {}
        start_index = 0
        end_index = 0
        for n in ["start_Smart-seq3", "start_SQK-NSK007", "end_Smart-seq3", "end_SQK-NSK007"]:
            read_trimmed_info[read_name][n] = 0
        ##
        # start alignments
        for info in trimmed_info:
            if re.search(r'start alignments', info):
                start_index = trimmed_info.index(info)
        # end alignments
        for info in trimmed_info:
            if re.search(r'end alignments', info):
                end_index = trimmed_info.index(info)
        # for start 
        if start_index != 0:
            if end_index != 0:
                sub = trimmed_info[start_index:end_index]  
            else:
                sub = trimmed_info[start_index:]
            for info in sub:
                if re.search(r'Smart-seq3', info):
                    read_trimmed_info[read_name]["start_Smart-seq3"] = 1
                if re.search(r'SQK-NSK007', info):
                    read_trimmed_info[read_name]["start_SQK-NSK007"] = 1 
        # for end
        if end_index != 0 :
            sub = trimmed_info[end_index:]
            for info in sub:
                if re.search(r'Smart-seq3', info):
                    read_trimmed_info[read_name]["end_Smart-seq3"] = 1
                if re.search(r'SQK-NSK007', info):
                    read_trimmed_info[read_name]["end_SQK-NSK007"] = 1 
    with open(outfile, "w") as fout:
        fout.write("\t".join(["read_name", "start_Smart-seq3", "start_SQK-NSK007", "end_Smart-seq3", "end_SQK-NSK007"]) + "\n")
        for read_name in read_trimmed_info.keys():
            fout.write("%s\t%s\t%s\t%s\t%s\n" % (
                read_name, 
                str(read_trimmed_info[read_name]["start_Smart-seq3"]),
                str(read_trimmed_info[read_name]["start_SQK-NSK007"]),
                str(read_trimmed_info[read_name]["end_Smart-seq3"]),
                str(read_trimmed_info[read_name]["end_SQK-NSK007"]),  
            ))


def main():
    args = parseArgs()
    fastq = args.input
    outDir = args.outdir
    softwaredir = args.softwaredir
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    logs = get_trimmed_log(fastq, outDir, softwaredir)
    for log in logs:
         porechop_trimmed_info(log)
    

if __name__=="__main__":
    main()