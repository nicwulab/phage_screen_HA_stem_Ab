# miniHA_Abs_Phage_Screen

## General Pipeline
- [x] Integration of Fastq files
- [x] Flank sequence trim
- [x] Removal of redundancy and marking of duplications
- [x] Top list for the AntiBody 
- [ ] Blast analysis of L/H chain
- [ ] Removal of low-quality bases in L/H region
- [x] Counting of reads
- [ ] Analysis of variations


# Integration

After integration, the total number of reads are:

| Sample     | Counts  |
|------------|--------|
| P0_Ab| 377573 |
| P2_cm| 330712 |
| P2_wt| 310041 |
| P3_cm| 397649 |
| P3_wt| 357229 |

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


