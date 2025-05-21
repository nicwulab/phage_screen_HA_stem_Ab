# miniHA_Abs_Phage_Screen

# Method

Analysis of PacBio sequencing data Circular consensus sequences (CCSs) were generated from the raw subreads using SMRTLink v13.0, setting the parameters to require 99.9% accuracy and a minimum of 3 passes. CSSs in FASTQ format were parsed using the SeqIO module in BioPython and filtered based on the base calling quality score, where any read with more than five nucleotides of phred quality score <40 were removed. (I am not sure about this part, the reads we got are well processed normal fastq reads) The adapter sequences were then identified on each read and trimmed from the Fab sequences. Reads that did not have the complete adapter sequences were also removed.
The filtered reads were then use seqkit 2.7.0 [PMC5051824] to remove convert it into fasta and remove the duplicates. Unique sequences aligned to the reference Fab sequences by blast+ 2.12.0+ [PMC2803857]. Frequency (F) of a Fab i of a given sample s was computed for each replicate as follows:
F_(i,s)=(〖readcount〗_(i,s)+1)/(∑_s▒( 〖readcount〗_(i,s)+1))   (1)
A pseudocount of 1 was added to each mutant to avoid division by zero in subsequent steps. We then calculate the enrichment (E) of a Fab i of a given replicate k after the phage display panning as follows:
E_(i,k)=F_(post-selection,i,k)/F_(pre-selection,i,k)    (2)


# Integration

After integration, the total number of reads are:

| Sample     | Counts  |
|------------|--------|
| P0_Ab| 377573 |
| P2_cm| 330712 |
| P2_wt| 310041 |
| P3_cm| 397649 |
| P3_wt| 357229 |

377573 + 330712 + 310041 + 397649 + 357229 = 1778204



# Flank sequence trim


Ps-r: 5’-3’ TTCAGttcaggaggaatttaaaatgaaaaagac
ps-f: 5’-3’ AGCTATGACCCACTCTTTCAACAGTCTTATCGTCATCG

Reads without Flank sequence:

TTCAGttcaggaggaatttaaaatgaaaaagac
TTCAGTTCAGGAGGAATTTAAAATGAACCTATTGCCTACGGCAGCCGCTGGATTGT

AGCTATGACCCACTCTTTCAACAGTCTTATCGTCATCG
AGCTATGACCCACTCAAGCCCCTTCCCTGGAGCCTGGCGGACCCAG


## Trimmed result

Untrimmed: 8494
Trimmed:  1,764,710
~ 0.47%

# Mark duplication

`[INFO] 608562 duplicated records removed`

# Statistic result

## total
|Sample|Counts|
|:-|:-|
 P0_Ab| 375990
 P2_cm| 329090
 P2_wt| 307941
 P3_cm| 395761
 P3_wt| 355928

## Duplicated reads

|Reads ID| Sample | Counts|
|:-|:-|:-|
0 |P0_Ab| 3
0| P2_cm| 324
0| P2_wt| 32
0| P3_cm| 36503
0| P3_wt| 75
1| P2_cm| 43
1| P2_wt| 90
1| P3_cm| 1989
1| P3_wt| 34175
2| P2_cm| 579

# ratio matrix

![](Picture/First_100_bar.png)

# Find the link sequence

`ttctagataattaattaggaggaatttaaaatgaaatacctattgcctacggcagccgctggattgttattactcgctgcccaaccagccatggcc`
Found the link sequence in 1154584 reads from 1156148 reads.

# split the reads 

split the reads into up & down independent sequences if the link sequences if the link is found. If not, keep the full sequence.

Result: 2310732 sequence in library

# blast the reads to the IG library

```
Splitting input fasta file ../PacBio/clean_split.fa
2,245,021 sequences successfully split into 2310 pieces
Starting process pool using 62 processors
2310723seq [44:49, 859.14seq/s]                              
2,245,021 sequences processed in 2,695.56 seconds, 832 sequences / s
Zipping up final output
Analysis complete, result file: Result/clean.tsv.gz
```
Result: 6.2G tsv file

# Separate the HV-chain out and Removing the domain which contains the stop codon

Because we mainly focus on the Heavy chain in this stage, separated file could dramatically reducing the calculation time.

critirias for filtering:
- less than 10 aa on the head
- aa seq longer than 60
- no stop codon on the v-domain

||Before|After|Clean|
|:-|:-|:-|:-|
Unique|2310724|1155627|942051|
Duplicated| 1216025|607813|538878|

# Families

![](Picture/Family_all_bar.png)

# Seqlogo
![](Picture/Kabat_all_duplicated.png)

```
options:
  -h, --help            show this help message and exit
  -i INPUT [INPUT ...], -I INPUT [INPUT ...], --input INPUT [INPUT ...]
                        input the annotated blast results. It is a list. You could input the Cleared results only. But if you want
                        to cont the ratio, please add the duplicated annotation result, too.
  -o OUTPUT, -O OUTPUT, --output OUTPUT
                        Please remember that the output is a directory
  -fv VFAMILY [VFAMILY ...], -FV VFAMILY [VFAMILY ...], --vfamily VFAMILY [VFAMILY ...]
  -t [TOP], -T [TOP], --top [TOP]
  -s SAMPLE [SAMPLE ...], -S SAMPLE [SAMPLE ...], --sample SAMPLE [SAMPLE ...]
  -p PSTAGE [PSTAGE ...], -P PSTAGE [PSTAGE ...], --pstage PSTAGE [PSTAGE ...]

Usage:
python scripts/I_Can_Do_all.py -i {Input annotated results}  -O '{Output directory} -s {samples} -p {stages} -t {number for top}

Example:
python scripts/I_Can_Do_all.py -i PacBio/clean.tsv.gz PacBio/Duplicated.tsv.gz  -O Result/SeqLogo -s wt cm -p P2 P3 -t 100 -f IGHV1-69 IGHV3-21
```

```bash
python scripts/I_Can_Do_all.py -i PacBio/clean.tsv.gz PacBio/Duplicated.tsv.gz  -O Result/SeqLogo -s 'Ab'
```



## Grep the overlap results for Tomas


```bash

for i in $(cat PacBio/20231110_Tomas.list ); do
    zgrep  $i PacBio/clean.tsv.gz > tmp
    IDS=$(awk '{print $1}' tmp | sed 's/_up//;s/_dn//'| tr '\n' ',')
    NUM=$(wc -l tmp| awk '{print $1}')
    COUNTS=$(awk -F"-" '{print $1}' tmp | sort|uniq -c| awk '{print $2":"$1}'| tr '\n'  ',')
    echo -e $i"\t"$NUM"\t"$COUNTS"\t"$IDS
done > Result/20231110_Tomas_Overlap.list
```

```bash

zcat PacBio/clean_IGH.tsv.gz| head -n 1 > ../HA_Abs/Ab_epitope/data/Top_500.tsv
zgrep -E $(head -n 500  Result/duplicates_IGH_AA.txt | awk '{print $2}'| sed 's/,//'| tr '\n' '|'| sed 's/|$/\n/') PacBio/clean_IGH.tsv.gz >> ../HA_Abs/Ab_epitope/data/Top_500.tsv
```





## Kabat Number and SeqLogo

![](Picture/Kabat_all_duplicated.png)
