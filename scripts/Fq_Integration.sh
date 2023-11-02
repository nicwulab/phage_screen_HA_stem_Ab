#!/usr/bin/bash

OUT=$(grep "^OUTPUT" Files_loc |awk '$1="";{print $NF}')

# integrate all fq into one file after add a prefix of reads name
for fq in $(grep "Fastq" Files_loc |awk '$1="";{print}'); do
    Pref=$(echo $fq| awk -F"/" '{print $NF}'| awk -F"_" '{print $1"_"$2}')
    sed  "s/^@m/@$Pref\-m/" $fq >> $OUT/All.fq
done

# Trim
## here is some problem with this program. So, I could only trimming them one by one
cutadapt -j 60 -g TTCAGttcaggaggaatttaaaatgaaaaagac  -o trimmed1.fastq --untrimmed-output untri1.fq $OUT/All.fq
cutadapt -j 60 -g AGCTATGACCCACTCTTTCAACAGTCTTATCGTCATCG  -o trimmed2.fastq --untrimmed-output untri2.fq  untri1.fq
cat trimmed1.fastq trimmed2.fastq > head_out.fq
cutadapt -j 60 -a CGATGACGATAAGACTGTTGAAAGAGTGGGTCATAGCT  -o trimmed1.fastq --untrimmed-output untri3.fq head_out.fq
cutadapt -j 60 -a gtctttttcattttaaattcctcctgaaCTGAA  -o trimmed2.fastq --untrimmed-output untri4.fq  untri3.fq
cat trimmed1.fastq trimmed2.fastq > HeadTail_out.fq
cat untri2.fq untri4.fq > $OUT/HeadTail_failed.fq

mv HeadTail_out.fq $OUT
rm untri* trimmed* head_out.fq

# Mark the duplicated reads
seqkit rmdup $OUT/HeadTail_out.fq  -s -i -o $OUT/clean.fastqc.gz -D Result/duplicates.txt

# Counts Total
grep "^@P"  $OUT/HeadTail_out.fq | sed 's/@//'| awk -F"-" '{print $1}'|sort |uniq -c | awk '{print "Total",$2,$1}' > Result/Counts_total.csv

# Count the uniq
zgrep "^@P"  $OUT/clean.fastqc.gz | sed 's/@//'| awk -F"-" '{print $1}'|sort |uniq -c | awk '{print "Unique",$2,$1}' > Result/Counts_uniq.csv


# Count the duplicated part
python3 scripts/duplciated_counts.py
head Result/Dupli_count.json | sed 's/}//;s/{//;s///;s/"//g;s/, /,/g;s/://g'| tr "," "\n" > Result/Counts_duplciate.csv
awk '{print $2}'  Result/Counts_duplciate.csv |sort|uniq -c| awk '{print "Repeat",$2,$1}' > Result/Counts_repeat.csv


# calculate the ratio and exponatial regression
python3 scripts/counts_cal.py

# grep the reads from top list (101 reads)
for i in $(cat Result/ExpRe_100_list.txt  Result/Ratio_100_list.txt | sort|uniq); do
    grep -A1 $(awk -v NUM=$(echo $i+1|bc) 'NR==NUM' Result/duplicates.txt | awk '{print $2}'| sed 's/,//') $OUT/HeadTail_out.fq | sed "s/^@/>$i:/" >> Top_list.fa &
done

for i in $(cat Result/Wt_100_list.txt | sort|uniq| sort -n); do
    grep -A1 $(awk -v NUM=$(echo $i+1|bc) 'NR==NUM' Result/duplicates.txt | awk '{print $2}'| sed 's/,//') $OUT/HeadTail_out.fq | sed "s/^@/>$i:/" >> Wt_Top100_ratio.fa &
done
for i in $(cat Result/Cm_100_list.txt | sort|uniq| sort -n); do
    grep -A1 $(awk -v NUM=$(echo $i+1|bc) 'NR==NUM' Result/duplicates.txt | awk '{print $2}'| sed 's/,//') $OUT/HeadTail_out.fq | sed "s/^@/>$i:/" >> Cm_Top100_ratio.fa &
done

# blast the link seq
seqkit fq2fa $OUT/clean.fastqc.gz -o $OUT/clean.fa
blastn -query linker.fa -out blast.out -db $OUT/blastdb/All -outfmt "6 qacc sacc evalue sstart send" -evalue 1e-5 -max_target_seqs 1156158 -num_threads 60 -max_hsps 1
python scripts/split_fq.py 


# blast
pyir -m 60  $OUT/clean_split.fa --outfmt tsv -o $OUT/clean


# do the same to the duplicated reads:
## extract the duplciated reads 
seqkit fq2fa $OUT/HeadTail_out.fq -o $OUT/HeadTail_out.fa
awk '$1="";$2="";{print $0}' Result/duplicates.txt | sed 's/,//g;s/^ //'| tr ' ' '\n'| sed '/^$/d' > Result/Duplicated_list.list

seqkit grep -f Result/Duplicated_list.list $OUT/HeadTail_out.fa -o $OUT/Duplicated.fa
# makeblastdb -parse_seqids -dbtype nucl -in ../Duplicated.fa -out Dupli
blastn -query linker.fa -out blast_du.out -db $OUT/blastdb/Dupli -outfmt "6 qacc sacc evalue sstart send" -evalue 1e-5 -max_target_seqs 1156158 -num_threads 60 -max_hsps 1
python scri

pyir -m 60  $OUT/Duplicated_split.fa --outfmt tsv -o $OUT/Duplicated


# family counts
zcat $OUT/clean.tsv.gz $OUT/Duplicated.tsv.gz | awk -F'\t' '{print $10,$1}'|sed 's/IG//g'|awk -F"-m" '{print $1}'| awk '{print $2,$1}'| awk -F"*" '{print $1}'| sort| uniq -c|sort -n| awk '{print $2,$3,$1}' > Result/V_class_family_count.tsv

# Plot the family counts
Rscript C_class_plot.R