#!/usr/bin/bash

OUT=$(grep "^OUTPUT" Files_loc |awk '$1="";{print $NF}')

# integrate all fq into one file after add a prefix of reads name
for fq in $(grep "Fastq" Files_loc |awk '$1="";{print}'); do
    Pref=$(echo $fq| awk -F"/" '{print $NF}'| awk -F"_" '{print $1"_"$2}')
    sed  "s/^@m/@$Pref\-m/" $fq >> $OUT/All.fq
done

# Mark the duplicated reads
seqkit rmdup $OUT/All.fq  -s -i -o $OUT/clean.fastqc.gz -D Result/duplicates.txt

# Counts Total
grep "^@P"  $OUT/All.fq | sed 's/@//'| awk -F"-" '{print $1}'|sort |uniq -c | awk '{print "Total",$2,$1}' > Result/Counts_total.csv

# Count the uniq
zgrep "^@P"  $OUT/clean.fastq.gz | sed 's/@//'| awk -F"-" '{print $1}'|sort |uniq -c | awk '{print "Unique",$2,$1}' > Result/Counts_uniq.csv


# Count the duplicated part
python script/duplicated_counts.py
head Result/Dupli_count.json | sed 's/}//;s/{//;s///;s/"//g;s/, /,/g;s/://g'| tr "," "\n" > Result/Counts_duplciate.csv
awk '{print $2}'  Result/Counts_duplciate.csv |sort|uniq -c| awk '{print "Repeat",$2,$1}' > Result/Counts_repeat.csv

