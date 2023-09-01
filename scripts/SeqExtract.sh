#!/bin/bash

fastq=$1
start_out=$2
end_out=$3
READ_CHECK=$4


gunzip -c $fastq | gawk -v TL="$READ_CHECK" \
    -v OUT1="$start_out"\
    -v OUT2="$end_out"\
   'NR%2==1{
     print $0 > OUT1;
     print $0 > OUT2;
   }
   NR%2==0{
     print substr($0, 1, TL) > OUT1;
     print substr($0, length($0) - TL + 1, TL)  > OUT2;
   }'
