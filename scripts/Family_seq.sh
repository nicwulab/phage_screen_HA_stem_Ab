
FAMILY=HV1-69
OUT=Result/$FAMILY
Num=200
mkdir $OUT

# extract the sequence based on their families
#zcat PacBio/clean.tsv.gz | awk -v FAMILY=IG$FAMILY -F"\t" '$95==FAMILY{print ">"$1"\n"$23}' > $OUT/$FAMILY\_domain.fa # the domain from the 

#zcat PacBio/clean.tsv.gz | awk -v FAMILY=IG$FAMILY -F"\t" '$95==FAMILY{print ">"$1"\n"$2}' > $OUT/$FAMILY\_split.fa # the main part for this seq


# for complit seuquence
#zcat PacBio/clean.tsv.gz | awk -v FAMILY=IG$FAMILY -F"\t" '$95==FAMILY{print $1}'| sed 's/_dn//;s/_up//'  > $OUT/IG$FAMILY.list
#seqkit grep -f $OUT/IG$FAMILY.list  PacBio/HeadTail_out.fa -o $OUT/IG$FAMILY\_uniq_complete.fa


# Top 100 $FAMILY
zcat PacBio/Duplicated.tsv.gz | awk -v FAMILY=IG$FAMILY -F"\t" '$95==FAMILY{print $1}'| sed 's/_dn//;s/_up//'  > $OUT/IG$FAMILY\_duplic.list
cat $OUT/IG$FAMILY\_duplic.list $OUT/IG$FAMILY.list | awk '{print "N "$1}'| sort| uniq > $OUT/IG$FAMILY\_all.list


cat $OUT/IG$FAMILY\_all.list Result/ID_reads.txt  | awk 'seen[$2]++'| awk -F '-' '{print $1}'| sort| uniq -c| awk '{print $1,$2,$3}'| sort -n -r > $OUT/IG$FAMILY\_ReadsID_counts.csv


Rscript scripts/Family_cal.R $FAMILY $Num $OUT

for i in $(cat $OUT/IG$FAMILY\_Alltop$Num.list ); do grep "^$i " Result/ID_reads.txt | head -n 1 ;done | awk '{print $2}'> tmp

seqkit grep -f tmp PacBio/HeadTail_out.fa -o $OUT/IG$FAMILY\_Alltop$Num.fasta
# add the ID at the head of each seq
for i in $(grep ">" $OUT/IG$FAMILY\_Alltop$Num.fasta | sed 's/>//'); do TMP=$(grep $i Result/ID_reads.txt|head -n 1); ID=$(echo $TMP| awk '{print $1}'); echo $ID $i; sed -i "s=$i=$ID:$i=" $OUT/IG$FAMILY\_Alltop$Num.fasta;done 

zcat PacBio/clean.tsv.gz | head -n 1 > $OUT/IG$FAMILY.tsv
zgrep -E $(cat tmp | tr "\n" "|"| sed 's/|$/\n/') PacBio/clean.tsv.gz >> $OUT/IG$FAMILY.tsv
