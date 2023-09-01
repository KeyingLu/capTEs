import os
import re
import sys
from multiprocessing import Pool
import argparse
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Align


    
def parseArgs():
    # dir of source code
    srcpath=sys.argv[0]
    srcdir=os.path.dirname(os.path.abspath(srcpath))
    #
    parser = argparse.ArgumentParser(description='To  cut site preferences')
    parser.add_argument('-id', '--ID', type=str, required=True, help="Sample name")
    #
    args = parser.parse_args()
    args.srcdir = srcdir
    return(args)



# function
def cut_sites(ID, mobile_element_name, mobile_element_fasta="NULL", sequence="NULL", start=0, winwidth = 80, min_score = 100):
    fastq = "%s.trimming.Q7.L300.fq.gz" % ID
    out1 = "read_start.fq"
    out2 = "read_end.fq"
    if not os.path.exists(out1):
        cmd = "bash %s/SeqExtract.sh %s %s %s %s" % (srcdir, fastq, out1, out2, 80)
        os.system(cmd)
    ### Store mobile element sequence and reverse complement
    mobile_element = ""
    if sequence != "NULL":
        mobile_element = SeqRecord(Seq(sequence), id=mobile_element_name)
    if mobile_element_fasta != "NULL":
        mobile_element = SeqIO.read(mobile_element_fasta, "fasta")
    mobile_element_seq = mobile_element.seq
    rev_comp_mobile_element_seq = mobile_element_seq.reverse_complement()  # TCCGCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCGGCC
    mobile_element_len = len(mobile_element_seq)
    ### Create Aligner
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 5
    aligner.mismatch_score = -5
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -5
    ### Align reads to mobiole element sequence
    read_fastq = out1
    outfile =  "read_start%s_cutsite_%s.txt" % (start, mobile_element_name)
    with open(outfile, "w") as results_out:
        i = 1  # Keep track of how many reads have been processed
        #
        for seq in SeqIO.parse(read_fastq, "fastq"):
            # if i % 50 == 0:
                # print(str(i) + " reads processed")
            results_out.write(seq.id + "\n")
            #
            # Use Bio.Align package to do pairwise alignments
            #
            forward_alignments = aligner.align(seq.seq[start:(winwidth+start)], mobile_element_seq)
            if not any(a.score >= min_score for a in forward_alignments):
                results_out.write("No alignments with score >= 100 on the forward strand" + "\n")
            else:
                if ((forward_alignments[0].aligned)[0][0][0]) == 0:
                    results_out.write("Forward alignment score: " + "\t" + str(forward_alignments[0].score) + "\n")
                    results_out.write(str(forward_alignments[0].aligned) + "\n")
                    results_out.write(seq.id + ":\t" + "Forward alignment starting mobile element base:" + "\t" + str(
                        (forward_alignments[0].aligned)[1][0][0] + 1) + "\t" + mobile_element_name + "\n")
                else:
                    results_out.write("Forward strand alignment does not begin at start of read" + "\n")
            #
            rev_comp_alignments = aligner.align(seq.seq[start:(winwidth+start)], rev_comp_mobile_element_seq)
            if not any(a.score >= min_score for a in rev_comp_alignments):
                results_out.write("No alignments with score >= %s on the reverse strand" % min_score + "\n")
                results_out.write("\n")
            else:
                if ((rev_comp_alignments[0].aligned)[0][0][0]) == 0:
                    results_out.write(
                        "Reverse strand alignment score: " + "\t" + str(rev_comp_alignments[0].score) + "\n")
                    results_out.write(str(rev_comp_alignments[0].aligned) + "\n")
                    results_out.write(seq.id + ":\t" + 
                        "Reverse alignment starting mobile element base:" + "\t" + str(mobile_element_len - (
                            (rev_comp_alignments[0].aligned)[1][0][0])) + "\t" + mobile_element_name + "\n")
                    results_out.write(
                        "\n")  # subtracting the first alignment base from the size of the element to get the starting base with respect to the forward mobile element sequence
                else:
                    results_out.write("Reverse strand alignment does not begin at start of read" + "\n")
                    results_out.write("\n")
            i = i + 1



            
if __name__=="__main__":
    args = parseArgs()
    ID = args.ID
    srcdir = args.srcdir
    for i in range(0,16): # for barcode samples    
        cut_sites(ID, "Alu", "NULL", 
                    "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGT",
                    i, 30, 80)
        # # print(ID + ": Alu finished\n")
        # mobile_element_fasta = "/NAS/wg_looking/gRNA_ONT/database/hg38reps.L1HS.fa"
        # print(ID + ": L1 finished\n")
        # pool.apply_async(cut_sites, 
        #                  args=(ID, "L1", mobile_element_fasta, 
        #                        "NULL",
        #                        i, 80, 100))
    #
    cmd = 'Rscript {srcdir}/cutSites_plot.r {ID}'.format(srcdir=srcdir, ID=ID)
    os.system(cmd)
    os.system('rm -rf read_start*_cutsite_Alu*.txt')
