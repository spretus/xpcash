
# u2os cells are female, therefore dont use chrY!
# LS paper doesnt incl chrY on the genome graph.

#YJ: "About which reference to use, I used “GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz”  based on Heng Li’s blog (author of Samtools): https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use. Here you can also find his explanations about some advantage/ disadvantage of different references."

#"If you map reads to GRCh38 or hg38, use the following:
#ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

mkdir hg38_bwa_fgcz
 mv hg38* hg38_bwa_fgcz/
 

# manually dl'd ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz to local then scp'd to /data/myanco/accessory_files/hg38/
myanco@login0:/data/myanco/accessory_files$ mkdir hg38

scp GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz myanco@cluster.s3it.uzh.ch:/data/myanco/accessory_files/hg38/

# /cluster/home/yanjiang/genome_index/bwa_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna 

gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
# fna looks like a fasta file.

# 
grep 'chr' GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
# >chr1  AC:CM000663.2  gi:568336023  LN:248956422  rl:Chromosome  M5:6aef897c3d6ff0c78aff06ac189178dd  AS:GRCh38
# ...
# >chrEBV  AC:AJ507799.2  gi:86261677  LN:171823  rl:decoy  M5:6743bd63b3ff2b5b8985d8933c53290a  SP:Human_herpesvirus_4  tp:circular


myanco@login1:/data/myanco/accessory_files/hg38_bwa_fgcz/hg38_BWAIndex$ grep 'chrEBV' genome.fa
# not found..is the herpes virus.

grep 'chrUn_KI270582v1' genome.fa
# >chrUn_KI270582v1


grep 'chr' /data/myanco/accessory_files/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna | sort -k1,1 | uniq | wc -l
# 195
nano HengLi.hg38ref.chrs.bed
# cp all header lines, containing chr id, from heng li's genome ref into a 1-col list.

grep 'chr' /data/myanco/accessory_files/hg38_bwa_fgcz/hg38_BWAIndex/genome.fa | sort -k1,1 | uniq | wc -l
# 455

#grep's -w option will find only lines that contain your target word as a complete word.
grep -w -f /data/myanco/accessory_files/hg38/HengLi.hg38ref.chrs.bed /data/myanco/accessory_files/hg38_bwa_fgcz/hg38_BWAIndex/genome.fa | sort -k 1,1V | uniq | wc -l
# 194   # only the herpes virus line is missing

# count bp:
grep -v ">" /data/myanco/accessory_files/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna | wc | awk '{print $3-$1}'
# 3,099,922,541

grep -v ">" /data/myanco/accessory_files/hg38_bwa_fgcz/hg38_BWAIndex/genome.fa | wc | awk '{print $3-$1}'
# 3,209,286,105


grep 'chr' /data/myanco/accessory_files/hg38_bwa_fgcz/hg38_BWAIndex/genome.fa | sort -k1,1 | uniq >chrs.bwaidx.hg38.tmp

grep -w -f /data/myanco/accessory_files/hg38/chrs.bwaidx.hg38.tmp /data/myanco/accessory_files/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna | sort -k 1,1V | uniq | wc -l
# 194


# NOT the same genome at all! hg38 has more chr's and more bp. Use Heng Li's.

# https://www.biostars.org/p/468917/
# GenCode version basically the same as Heng Li version.
# the one exception is the herpes virus



# https://www.gencodegenes.org/human/
Content	Regions	Description	Download
#Comprehensive gene annotation	CHR	
#It contains the comprehensive gene annotation on the reference chromosomes only
#This is the main annotation file for most users

#GTF GFF3
#Comprehensive gene annotation	ALL	
#It contains the comprehensive gene annotation on the reference chromosomes, scaffolds, assembly patches and alternate loci (haplotypes)
#This is a superset of the main annotation file

# 
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz

 gunzip gencode.v39.annotation.gtf.gz

head gencode.v39.annotation.gtf
##description: evidence-based annotation of the human genome (GRCh38), version 39 (Ensembl 105)
##provider: GENCODE
##contact: gencode-help@ebi.ac.uk
##format: gtf
##date: 2021-09-02
chr1	HAVANA	gene	11869	14409	.	+	.	gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; level 2; hgnc_id "HGNC:37102"; havana_gene "OTTHUMG00000000961.2";
chr1	HAVANA	transcript	11869	14409	.	+	.	gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";



# What about chromosome sizes and bwa alignment files:

# for chr sizes can derive these:
# for line in file_of_fna_headers; do grep find chr in fna file, grep -v exclude header line | wc | awk print $3-$1 ;done

#!/bin/bash

while IFS= read -r line; do
  echo "tester: $line"
done < "$1"

# ==========
#!/bin/bash
while IFS='' read -r LINE || [ -n "${LINE}" ]; do
    echo "processing line: ${LINE}" >>output.tmp
done < ~/Documents/NaegeliLab/file.txt


# ======
sh test.sh
cat output.tmp
processing line: a
processing line: v
processing line: x
processing line: w

# but how to grab the lines until the next header line?
# https://unix.stackexchange.com/questions/21076/how-to-show-lines-after-each-grep-match-until-other-specific-match

grep 'chr' GCA_000001405.15_GRCh38_no_alt_analysis_set.fna | sed -e '1d;$d'


# ==========
#!/bin/bash
while IFS='' read -r LINE || [ -n "${LINE}" ]; do
#    echo "processing line: ${LINE}" >>output.tmp
# will skip last line = herpes virus bc it doesnt end with a new line char.
	grep ${LINE} /data/myanco/accessory_files/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna 
	
done < /data/myanco/accessory_files/hg38/HengLi.hg38ref.chrs.bed




grep -v ">" /data/myanco/accessory_files/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna | wc | awk '{print $3-$1}'



GCA_000001405.15_GRCh38_no_alt_analysis_set.fna


# calculate chr sizes from fasta ( or fna): https://www.biostars.org/p/173963/#174150

samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna


An fai index file is a text file consisting of lines each with five TAB-delimited columns for a FASTA file and six for FASTQ:

# http://www.htslib.org/doc/faidx.html

#NAME	Name of this reference sequence
#LENGTH	Total length of this reference sequence, in bases
#OFFSET	Offset in the FASTA/FASTQ file of this sequence's first base
#LINEBASES	The number of bases on each line
#LINEWIDTH	The number of bytes in each line, including the newline
#QUALOFFSET	Offset of sequence's first quality within the FASTQ file


cut -f 1,2 GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai >GCA_hg38.chrom_chr.sizes

# 195 GCA_hg38.chrom_chr.sizes

# generate bwa index files from this genome?

awk '{total +=$2} END {print total/1e9}' GCA_hg38.chrom_chr.sizes
# 3.09992  Gb

awk '{total +=$2} END {print total/1e9}' ../hg38_bwa_fgcz/hg38.chrom_chr.sizes
# 3.08827   # this doesnt even match the total nt counted from the fasta file?!




grep -v ">" /data/myanco/accessory_files/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna | wc | awk '{print $3-$1}'
# 3,099,922,541

grep -v ">" /data/myanco/accessory_files/hg38_bwa_fgcz/hg38_BWAIndex/genome.fa | wc | awk '{print $3-$1}'
# 3,209,286,105

# bwa aligner used for xpc chip? h3k4me3 chip? cpd-seq?


# redo tads, intertads:
bedtools complement -i 360min_rep12cat.TADs.mrg.hg38.bed -g /data/myanco/accessory_files/hg38/GCA_hg38.chrom_chr.sizes | sort -k1,1V -k2,2n >iTADcatrep12.bedtcompl.hg38.bed

# calculate WITHOUT chrY since u2os are female derived:
#grep '>' my.fa | grep -v 'chrU' | sed 's|[>,]||g' | xargs samtools faidx my.fa > my.filtered.fa

#grep '>' GCA_000001405.15_GRCh38_no_alt_analysis_set.fna | grep -v 'chrY' | sed 's|[>,]||g' | xargs samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna >GCA_000001405.15_GRCh38_no_alt_analysis_set.nochrY.fna


#grep '>' GCA_000001405.15_GRCh38_no_alt_analysis_set.fna | grep -v 'chrY' | sed 's|[>,]||g' | xargs samtools faidx my.fa > my.filtered.fa

# seqkit install:

conda install -n myenv -c bioconda seqkit


seqkit grep -vrp "^chrY" /data/myanco/accessory_files/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna > /data/myanco/accessory_files/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.nochrY.fna

grep 'chrY'  GCA_000001405.15_GRCh38_no_alt_analysis_set.nochrY.fna
# nothin
grep 'chrY'  GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
#>chrY  AC:CM000686.2  gi:568336000  LN:57227415  rl:Chromosome  M5:ce3e31103314a704255f3cd90369ecce  AS:GRCh38  hm:10001-2781479,56887903-57217415
#>chrY_KI270740v1_random  AC:KI270740.1  gi:568335373  LN:37240  rg:chrY  rl:unlocalized  M5:69e42252aead509bf56f1ea6fda91405  AS:GRCh38



grep -v ">" /data/myanco/accessory_files/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.nochrY.fna | wc | awk '{print $3-$1}'
# 3,042,657,886 bp

samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.nochrY.fna

cut -f 1,2 GCA_000001405.15_GRCh38_no_alt_analysis_set.nochrY.fna.fai >GCA_hg38.chrom_chr.nochrY.sizes

# 193 GCA_hg38.chrom_chr.nochrY.sizes

# generate bwa index files from this genome?

awk '{total +=$2} END {print total/1e9}' GCA_hg38.chrom_chr.nochrY.sizes
# 3.04266 Gb

mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GCA_000001405.15_GRCh38_no_alt_analysis_set.withchrY.fna

mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai GCA_000001405.15_GRCh38_no_alt_analysis_set.withchrY.fna.fai

mv GCA_hg38.chrom_chr.sizes GCA_hg38.chrom_chr.withchrY.sizes


mv GCA_000001405.15_GRCh38_no_alt_analysis_set.withchrY.fna GCA_000001405.15_GRCh38_no_alt_analysis_set.withchrY.noncanon.fna

mv GCA_000001405.15_GRCh38_no_alt_analysis_set.withchrY.fna.fai GCA_000001405.15_GRCh38_no_alt_analysis_set.withchrY.noncanon.fna.fai

mv GCA_hg38.chrom_chr.withchrY.sizes GCA_hg38.chrom_chr.withchrY.noncanon.sizes



mv GCA_000001405.15_GRCh38_no_alt_analysis_set.nochrY.fna GCA_000001405.15_GRCh38_no_alt_analysis_set.nochrY.noncanon.fna

mv GCA_000001405.15_GRCh38_no_alt_analysis_set.nochrY.fna.fai GCA_000001405.15_GRCh38_no_alt_analysis_set.nochrY.noncanon.fna.fai

mv GCA_hg38.chrom_chr.nochrY.sizes GCA_hg38.chrom_chr.nochrY.noncanon.sizes

# canonical chr's only:
chrM, GL, KI, EBV



grep -v 'chrM' GCA_hg38.chrom_chr.withchrY.sizes | grep -v 'GL' | grep -v 'KI' | grep -v 'EBV'
# this captured everyone.
grep -v '^chrM\|GL\|KI\|EBV' GCA_hg38.chrom_chr.withchrY.sizes


#seqkit grep -vrp "^chrM\|^chrUn\|^chr" /data/myanco/accessory_files/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna > /data/myanco/accessory_files/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.nochrY.fna

#seqkit grep -vrp '^chrM\|GL\|KI\|EBV' GCA_000001405.15_GRCh38_no_alt_analysis_set.withchrY.noncanon.fna >GCA_000001405.15_GRCh38_no_alt_analysis_set.withchrY.canon.fna
# not excluding anything also when grep search is in double quotes.

seqkit grep -vrp "^chrM" GCA_000001405.15_GRCh38_no_alt_analysis_set.withchrY.noncanon.fna | seqkit grep -vrp "GL" | seqkit grep -vrp "KI" | seqkit grep -vrp "EBV" >GCA_000001405.15_GRCh38_no_alt_analysis_set.withchrY.canon.fna

samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.withchrY.canon.fna

cut -f 1,2 GCA_000001405.15_GRCh38_no_alt_analysis_set.withchrY.canon.fna.fai >GCA_hg38.chrom_chr.withchrY.canon.sizes

grep -v ">" /data/myanco/accessory_files/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.withchrY.canon.fna | wc | awk '{print $3-$1}'
# 3,088,269,832 bp
awk '{total +=$2} END {print total/1e9}' GCA_hg38.chrom_chr.withchrY.canon.sizes
# 3.08827 Gb



seqkit grep -vrp "^chrM" GCA_000001405.15_GRCh38_no_alt_analysis_set.nochrY.noncanon.fna | seqkit grep -vrp "GL" | seqkit grep -vrp "KI" | seqkit grep -vrp "EBV" >GCA_000001405.15_GRCh38_no_alt_analysis_set.nochrY.canon.fna

samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.nochrY.canon.fna

cut -f 1,2 GCA_000001405.15_GRCh38_no_alt_analysis_set.nochrY.canon.fna.fai >GCA_hg38.chrom_chr.nochrY.canon.sizes


grep -v ">" /data/myanco/accessory_files/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.nochrY.canon.fna | wc | awk '{print $3-$1}'
# 3,031,042,417 bp
awk '{total +=$2} END {print total/1e9}' GCA_hg38.chrom_chr.nochrY.canon.sizes
# 3.03104 Gb



# ========================================================================================

# now define regions based on the most appropriate reference genome...
# use canonical chr's only and no chrY to define expression- and Hi-C based feats: bc thy are from u2os cells and because the tad coords are only to canonical chr's




360min_rep1_40kb_TAD.bed
360min_rep2_40kb_TAD.bed

interTADrep1.bedtcompl.bed
interTADrep2.bedtcompl.bed



awk '{print $3-$2}' 360min_rep1_40kb_TAD.bed | awk '{total +=$1} END {print total/1e6}'
# 2576.6 Mb
awk '{print $3-$2}' 360min_rep2_40kb_TAD.bed | awk '{total +=$1} END {print total/1e6}'
#2567.68

 awk '{print $3-$2}' interTADrep1.bedtcompl.bed  | awk '{total +=$1} END {print total/1e6}'
#509.061
awk '{print $3-$2}' interTADrep2.bedtcompl.bed  | awk '{total +=$1} END {print total/1e6}'
#520.701

# merge the TAD repls and THEN derive the interTADs? = conservative estimate of which regions are TADs.

bedtools intersect -a  360min_rep1_40kb_TAD.bed -b 360min_rep2_40kb_TAD.bed -v | wc -l
# 3
bedtools intersect -b  360min_rep1_40kb_TAD.bed -a 360min_rep2_40kb_TAD.bed -v | wc -l
# 10

# 3-10 tads dont intersect.

wc -l 360min_rep*_40kb_TAD.bed
# 1939 360min_rep1_40kb_TAD.bed
# 2103 360min_rep2_40kb_TAD.bed
# 4042 total

bedtools intersect -a  360min_rep1_40kb_TAD.bed -b 360min_rep2_40kb_TAD.bed | sort -k1,1V -k2,2n | bedtools merge -i stdin | wc -l
#1512

cat 360min_rep1_40kb_TAD.bed 360min_rep2_40kb_TAD.bed | sort -k1,1V -k2,2n | bedtools merge -i stdin | wc -l
# 1000

cat 360min_rep1_40kb_TAD.bed 360min_rep2_40kb_TAD.bed | sort -k1,1V -k2,2n | bedtools merge -i stdin | awk '{print $3-$2}' | awk '{total +=$1} END {print total/1e6}'
# 2612.52

# currenty each sep one covers 2576.6 Mb and 2567.68.
# merging reduces nr bins while keeping total distance covered roughly the same.

cat 360min_rep1_40kb_TAD.bed 360min_rep2_40kb_TAD.bed | sort -k1,1V -k2,2n | bedtools merge -i stdin >360min_rep12cat.TADs.mrg.hg38.bed




#bedtools complement -i <(sort -k1,1V -k2,2n 360min_rep1_40kb_TAD.bed) -g ../accessory_files/hg19.chrom_chr.sizes >interTADrep1.bedtcompl.bed

bedtools complement -i /data/myanco/accessory_files/HiC_U2OS_hg38/360min_rep12cat.TADs.mrg.hg38.bed -g /data/myanco/accessory_files/hg38/GCA_hg38.chrom_chr.nochrY.canon.sizes | wc -l
# 1021
bedtools complement -i /data/myanco/accessory_files/HiC_U2OS_hg38/360min_rep12cat.TADs.mrg.hg38.bed -g /data/myanco/accessory_files/hg38/GCA_hg38.chrom_chr.nochrY.canon.sizes | bedtools merge -i stdin | wc -l
#1021

bedtools complement -i /data/myanco/accessory_files/HiC_U2OS_hg38/360min_rep12cat.TADs.mrg.hg38.bed -g /data/myanco/accessory_files/hg38/GCA_hg38.chrom_chr.nochrY.canon.sizes | awk '{print $3-$2}' | awk '{total +=$1} END{print total/1e6}'
# 418.522 Mb


# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02108-x
# Proximity ligation has revealed the existence of chromosomal topologically associating domains (TADs) of 10e5 to 10e6 bp in size [9]. TADs are flanked by boundaries that are discovered as loci where there is a sharp break from preferential left-ward interactions to preferential right-ward interactions. In 2012, Dixon et al. thus modelled Hi-C data to detect 1723 human TAD boundaries [10, 11]. Similarly, the Blueprint consortium reported between 2800 and 3741 TADs in primary human blood cell types [12]. 


bedtools complement -i /data/myanco/accessory_files/HiC_U2OS_hg38/360min_rep12cat.TADs.mrg.hg38.bed -g /data/myanco/accessory_files/hg38/GCA_hg38.chrom_chr.nochrY.canon.sizes | sort -k1,1V -k2,2n >/data/myanco/accessory_files/HiC_U2OS_hg38/iTADcatrep12.bedtcompl.hg38-nochrY-canon.bed

awk '{print $3-$2}'  iTADcatrep12.bedtcompl.hg38-nochrY-canon.bed | awk '{total +=$1} END {print total/1e6}'
#418.522 Mb

bedtools merge -i iTADcatrep12.bedtcompl.hg38-nochrY-canon.bed | awk '{print $3-$2}' | awk '{total +=$1} END {print total/1e6}'
#418.522 Mb


wc -l iTADcatrep12.bedtcompl.hg38-nochrY-canon.bed
bedtools merge -i iTADcatrep12.bedtcompl.hg38-nochrY-canon.bed | wc -l
# both 1021

 awk '{total+=$2}END{print total/1e6}' /data/myanco/accessory_files/hg38/GCA_hg38.chrom_chr.nochrY.canon.sizes
#3031.04

=2612.52+418.522 =  3031.042 Mb or the total genome, canonial only and no chrY.
=2612.52/3031.042
=418.522/3031.042

# 0.8619214118445077
# 0.1380785881554924

# 86% of the genome is covered by tads and 14% by intertads.
# TADs and interTADs location:
/data/myanco/accessory_files/HiC_U2OS_hg38/360min_rep12cat.TADs.mrg.hg38.bed
/data/myanco/accessory_files/HiC_U2OS_hg38/iTADcatrep12.bedtcompl.hg38-nochrY-canon.bed



# ======================================================================================
# RNA-seq dataset used? Kang et al 2021 (same place where got the TAD coords) has RNA-seq U2OS too, but in what format? 
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141137
Supplementary file	Size	Download	File type/resource
GSE141137_NormalizedFilteredData.txt.gz	1.8 Mb	(ftp)(http)	TXT
# GSE141137_WaveExpressionClusters.txt.gz	520.2 Kb	(ftp)(http)	TXT  # cant be dl'd

 scp GSE141137_NormalizedFilteredData.txt.gz myanco@cluster.s3it.uzh.ch:/scratch/myanco/cpdsq/
 scp GSE141137_WaveExpressionClusters.txt.gz myanco@cluster.s3it.uzh.ch:/scratch/myanco/cpdsq/
 
 wc -l *txt
#   9835 GSE141137_NormalizedFilteredData.txt
#  12195 GSE141137_WaveExpressionClusters.txt
 
 head -n 2 GSE141137_NormalizedFilteredData.txt
	tags.HK193_S1...11705407.0.total.	tags.HK194_S2...14286360.0.total.	tags.HK195_S3...12338979.0.total.	tags.HK196_S4...14244682.0.total.	tags.HK197_S5...13086284.0.total.	tags.HK198_S6...13875615.0.total.	tags.HK199_S7...13026745.0.total.	tags.HK200_S8...20645841.0.total.	tags.HK201_S9...18684454.0.total.	tags.HK202_S10...16889948.0.total.	tags.HK203_S11...18886637.0.total.	tags.HK204_S12...12043294.0.total.	tags.HK205_S13...13328696.0.total.	tags.HK206_S14...7163236.0.total.	tags.HK207_S15...14752734.0.total.	tags.HK208_S16...13868016.0.total.	tags.HK209_S17...16299561.0.total.	tags.HK210_S18...17269124.0.total.	tags.HK211_S19...17814615.0.total.	tags.HK212_S20...20044347.0.total.	tags.HK213_S21...20961382.0.total.	tags.HK214_S22...17027984.0.total.
AFF4|AF5Q31|CHOPS|MCEF|-|5q31.1|protein-coding	2.18605736847299	5.03014075932146	1.99951116317672	2.98616862249085	1.90193019439316	3.40995377904997	4.34029628004216	3.81515812700679	5.2793983439447	3.93041993762174	4.72051801858514	2.19444039336689	4.70683692865605	1.49358827357924	3.21022270347137	2.12989292679379	3.59481561568513	4.72003404458232	3.61546340769388	5.58790269154836	4.16039251611991	4.41252597708562

head GSE141137_WaveExpressionClusters.txt
Annotation/Divergence	0 min	35 min	60 min	90 min	120 min	180 min	Asynch	Class
FOSB	4.196452695	3.009984053	8.073862357	7.891916689	7.496067371	4.857452648	3.247340967	1
CKS2	4.072638095	4.568441791	6.779246978	6.842474276	7.229590005	5.323462581	6.197954034	1

# doesnt look relevant, also u2os 360min not included.

# https://www.biostars.org/p/468917/
# GenCode version basically the same as Heng Li version.


# https://www.biostars.org/p/328824/
#"And all Heng Li says is: "If you map reads to GRCh37 or hg19, use hs37-1kg"
hs37 is a special genome reference prepared for 1000 genomes project by this method. You can find that data here.

Ultimately GENCODE is the organization project responsible for managing human/mouse genome data. They provide the authoritative genome data that is used by everyone including NCBI/UCSC/Ensembl.


#The b37 is hs37-1kg and does not include only the "25 longest sequences from GRCh37 (1-22,X,Y,MT)" but it is a 1000 Genome convention that includes: -The 24 "relatively complete" chromosomal sequences (named "1" to "22", "X" and "Y") downloaded individually from ENSEMBL. -The GRCh37.p2 (rCRS) mitochondrial sequence (named "MT") downloaded from MITOMAP or NCBI. -The unlocalized sequences, which were named after their accession numbers, such as "GL000191.1", "GL000194.1", etc. -The unplaced sequences, which were named after their accession numbers, such as "GL000211.1", "GL000241.1", etc. Only the alternate loci were not included in the b37 dataset.

#hs37d5 (known also as b37 + decoy) was released by The 1000 Genomes Project (Phase II), which introduced additional sequence (BAC/fosmid clones, HuRef contigs, Epstein-Barr Virus genome) to the b37 reference to help reduce false positives for mapping. Note that this one uses the primary assembly of GRCh37.p4 (not the one of GRCh37 w/o patches).

#As for hs37 (without -1kg) I think it is generated only by bwakit in BWA and according to their manual it corresponds to b37+EBV (Epstein-Barr Virus genome). EBV genome is also found in hs37d5 and GRCh38 and it is included because it is used in molecular biology for transformations and because it naturally infects B cells in ~90% of the world population.


# same coords for GENCODE39 as for hs37-1kg???
# 




 

# current rna seq data set is u2os cells but they are Fucci-flavored
# grep 'Fucci' *sh
#211115_cpdsq.sh:Here, we provide a resource consisting of mapped transcriptomes in unsynchronized HeLa and U2OS cancer cells sorted for cell cycle phase by Fucci reporter expression.

# Akan et al 2012:
#We sequenced six Illumina paired-end lanes for the osteosarcoma (U2OS) cell line, and five for each of the other two cell lines, glioblastoma (U251) and epidermoid carcinoma (A431). In total, there were 16 lanes, amounting to 1.23 billion paired-end reads. The data are publicly available [ERP001947] [16].  --> only fq's
# 

#Genes
# 692 genes with high or low expression in U2OS relative to other cell lines from the CCLE Cell Line Gene Expression Profiles dataset.
# https://maayanlab.cloud/Harmonizome/gene_set/U2OS/CCLE+Cell+Line+Gene+Expression+Profiles

# 692 total listed, 310 active 382 inactive = 45% and 55% active and inactive.

# by contrast from single-cell rna-seq bostrom et al, there were 18k uniq genes i counted as positive:
  gunzip -cd RNAsq.U2OS_G1_1.coords.gz | wc -l
#       61643   incl sep lines for each alterante transcript?
# 18,191 geneNames.PosReadCount.u2osG1_1.uniq.bed

# Chang et al Cell 2014:
https://www.cell.com/action/showPdf?pii=S2211-1247%2814%2900499-9
#We first examined p53 binding in U2OS cells 6 hr after UV treatment, which we confirmed to induce signatures of apoptosis and responses to DNA damage
# Our findings indicate that p53 may directly activate at least 151 annotated genes in U2OS cells in response to UV damage. This number substantially exceeds known direct gene targets of p53.

https://www.biostars.org/p/238/
# I want to know which of these genes are "active" (or in other words: are likely to produce enough protein products to have an effect). I'm not interested in them being differentially expressed or X-fold over- or under-expressed. All I want is the classification of them being likely "on" or "off".
# So far I log-transformed (basis 10) the RMA score and centered them (subtracted the median). I called all genes which had a transformed score <0 as being inactive and scores >0 as being active.


# You're right in thinking that your methodology isn't a very good representation of the system. mRNAs (and their protein products) have a huge dynamic range. Some are going to be expressed constantly at extremely low levels, and at the other extremes, you'll have genes that are highly expressed, but only for a short period of time. Taking the median level as the dividing line between on and off is going to give you huge numbers of false negatives (genes that are actually being transcribed and translated, but that you'll classify as "off")
#I'd look at what the background noise level is, then run some stats to determine which probes give you signal significantly above that level. Any gene meeting that criteria should probably be considered "on". I suspect that may not divide the set as nicely as you'd hope, though.
# I agree, using the background level as noise and using that as not-expressed at all sounds like a good approach


# determine if Bostrom et al G1_1 and G1_2 are repls ?
#The raw read data files, Read Counts and RPKM values are available as a GEO submission (https://www.ncbi.nlm.nih.gov/geo/, #GSE104736).
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104736

#GSM2806897	HeLa_Fucci_G1_1
#GSM2806898	HeLa_Fucci_G1_2
#GSM2806899	HeLa_Fucci_G1_3

#GSM2806900	HeLa_Fucci_S_1
#GSM2806901	HeLa_Fucci_S_2
#GSM2806902	HeLa_Fucci_S_3

#GSM2806903	HeLa_Fucci_G2_1
#GSM2806904	HeLa_Fucci_G2_2
#GSM2806905	HeLa_Fucci_G2_3

#GSM2806906	U2OS_Fucci_G1_1
#GSM2806907	U2OS_Fucci_G1_2

#GSM2806908	U2OS_Fucci_S_1
#GSM2806909	U2OS_Fucci_S_2

#GSM2806910	U2OS_Fucci_G2_1
#GSM2806911	U2OS_Fucci_G2_2

# HeLa cells were sorted at three timepoints, while U2OS cells were sorted at two timepoints. Each time into three groups, categorized as "G1", "S", and "G2".

# thus the _1, _2 (and _3 for HeLa) seem to be not replicates but rather 2 distnct timepoints!


# 211214_u2osSTRANDEDgenicregion
# RNA-seq pipe: 211115_cpdsq.sh

# RNA-seq in U2OS, single-cell:, G1_1 phase: 

# U2OS, single-cell RNA-seq:, G1_1 phase: Böström et al. 2017, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5724894/
	# dl'd the list of genes: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104736

		# 		GSE104736_ReadCounts.csv

scp ~/Downloads/GSE104736_ReadCounts.csv myanco@cluster.s3it.uzh.ch:/data/myanco/accessory_files/hg38/genesU2OS/

# 28,378 GSE104736_ReadCounts.csv

awk -F ',' '{print NF;exit}' GSE104736_ReadCounts.csv
#16

#U2OS G1_1 genes are in col11.
sed 's/,/\t/g' GSE104736_ReadCounts.csv | cut -f 11 | head
#U2OS_Fucci_G1_1

# filter out weird lines like those referring to months...


# doesnt matter which genome build, coords not given, just the gene names.
# grabbed the active and inactive gene bodies and TSS from this list by searching ucsc table browser for the gene name abbrev's. --> hg38 on genome browser...! UCSC not same genome build as Heng Li?? Should it matter for this part as the genes are likely to be found only on the canonical chr's?
# 

# filter out weird entries like dates instead of gene names:
awk 'NR>1' GSE104736_ReadCounts.csv | sed 's/,/\t/g' | cut -f 1,11 | grep -v 'mar' | grep -v 'oct' | grep -v 'sep' | grep -v 'dec' | awk '{if ($2>0) print $1}' | sed 's/_/\t/g' | cut -f 1 | sort -k1,1 | uniq | sed 's/-/\t/g' | cut -f 1 | sort -k1,1 | uniq >geneNames.PosReadCount.u2osG1_1.uniq.bed

# histogram of read counts: divide into further categories? define a threshhold?
# on external HD, cpdsq folder:
	#U2OS-G1_1.RNAsq.counts.bed
	#RNAsq.U2OS_G1_1.coords.gz
	#geneNames.PosReadCount.u2osG1_1.bed
	#geneNames.PosReadCount.u2osG1_1.uniq.bed


df <- read.table("/Volumes/My Book/Mish/cpdsq/U2OS-G1_1.RNAsq.counts.bed",header=TRUE)
hist(df$Count, 1000)
# huge peak at 0.
# exclude the 0s:
df1 <- subset(df,Count>0)
colnames(df1) <- c("ReadCount")
# still huge peak near 0. x-axis goes up to 200k, limit it:
 hist(df1$ReadCount,500000,xlim=c(0,10))
# highest at 1, then 2

hist(log(df1$ReadCount),5000,xlim=c(0,10))


# try taking the standard deviation and using that to call a thrshhold of active vs inactive genes?
# Erin does DGE across samples, but doesnt look at individual datasets.
# before getting the DGE she filters out genes with low read counts: ie she asks for them to have at least 10 reads across multiple samples!

# check if this dataset includes any biological replicates? G1_2 = ??

#To categorize genes with similar expression patterns, the spectrum of relationships was divided into 6 categories based on θ-value (Fig 1D, Panel G in S2 Fig). Each θ-category corresponds to a binary classification of gene expression over the three phases relative to each other. Categories 1, 3 and 5 represent patterns of gene expression that peak in a single cell cycle phase: G1, S or G2/M. Categories 2, 4 and 6 represent a composite pattern where gene expression is high in two consecutive cell cycle phases and low in the third phase, i.e., high in G1+S, S+G2/M and G2/M+G1.



# # connect the gene id from the u2os single-cell rna-seq paper to something in table browser's list of transcription coordinates that i've been using to split the genome into genes, TSS, and intergenic regions.


	# For the intersection via UCSC, it is based on gene id, thus strand info not used. For everything else, pay attention to strand!!!! Force strandedness bedtools subtract, bedtools merge, etc.

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz

 gunzip gencode.v39.annotation.gtf.gz

## not able to determine chr sizes..
# but people do seem to mix coords from hs37 and GENCODE genomes:
# https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-76578/protocols/
# llumina Casava1.7 software used for basecalling. The reads were mapped to human hs37 reference genome (including human rRNA sequence) using Tophat2 with default parameters. The RPKM values were calculated for all genes with annotation in Gencode V17 including mRNA and snoRNAs. 

# check hg38 other refs:
paste -d'\t' <(sort -k1,1V -k2,2n  /data/myanco/accessory_files/hg38_bwa_fgcz/hg38.chrom_chr.sizes | cut -f 2 ) <(sort -k1,1V -k2,2n /data/myanco/accessory_files/hg38/GCA_hg38.chrom_chr.withchrY.canon.sizes | cut -f 2) | awk '{print $2-$1}' | sort -k1,1g | uniq
 # same canonical chr sizes.



# http://genome.ucsc.edu/cgi-bin/hgTables\
#Connected Tables and Joining Fields

#Paste In Identifiers for UCSC Genes: input is : \
	#  /Volumes/My\ Book/Mish/cpdsq/geneNames.PosReadCount.u2osG1_1.uniq.bed
	
# output: RNAsq.U2OS_G1_1.coords.GENCODDEv39.gz

#group:  Genes and Gene Predictions   track:     GENCODE V39
#group:  knownGene

# GENCODE V39 same as Heng Liu's genome?
# GENCODE V39:

# output filename: 
	# RNAsq.U2OS_G1_1.coords.gz
	
#	output format: alll fields from selected table
	


#Note: 1130 of the 18191 given identifiers have no match in table knownGene, field name or in alias table kgAlias, field alias. Try the "describe table schema" button for more information about the table and field.
#10 example missing identifier(s):
#AD026
#AHSA2
#ALMS1P

# 1130/18191
# 0.06211863009180364    # 6% of the gene id's supplied by the paper can't be found...

# with hg19 and applying some filters by trial and error, i got it down to :
# ucsc upload:
#Note: 1117 of the 18191 given identifiers have no match in table knownGene, field name or in alias table kgAlias, field alias. Try the "describe table schema" button for more information about the table and field.
#10 example missing identifier(s):
#AADACP1  --> does not show up in error list for hg38.


# Bostrom et al is 2017 published. no note about which genome build they used.
# hg38 came out 2013, they may have used it ..
# 1117/18191
# 	 =0.0614039909845528
# also 6% not found for hg19.
	
mv RNAsq.U2OS_G1_1.coords.GENCODDEv39 RNAsq.U2OS_G1_1.coords.GENCODDEv39.bed
gzip RNAsq.U2OS_G1_1.coords.GENCODDEv39.bed
mv ~/Downloads/RNAsq.U2OS_G1_1.coords.GENCODDEv39.bed.gz /Volumes/My\ Book/Mish/cpdsq/
	

gunzip -cd RNAsq.U2OS_G1_1.coords.GENCODDEv39.bed.gz | wc -l
# 165,772. for hg19 is 61,643 ... why so much higher in hg39?!


# gunzip -cd RNAsq.U2OS_G1_1.coords.gz | wc -l
#   61643
# alternate transcripts are present...multiple per gene.


#scp RNAsq.U2OS_G1_1.coords.gz myanco@cluster.s3it.uzh.ch:/scratch/myanco/cpdsq/
#gunzip -cd RNAsq.U2OS_G1_1.coords.gz | awk 'NR>1' | cut -f 1  | sort -k1,1 | uniq | wc -l
#61642
gunzip -cd RNAsq.U2OS_G1_1.coords.gz| awk 'NR>1' | cut -f 2 | sort -k1,1V | uniq -c
# far fewer noncanon chrs than hg38.
# 45 chrs total of which 23 are canonical and 22 are noncano / chrY.

#gunzip -cd RNAsq.U2OS_G1_1.coords.gz| awk 'NR>1' | cut -f 2 | sort -k1,1V | uniq -c | sed 's/ /\t/g' | grep 'random\|hap\|chrY\|Un'



scp RNAsq.U2OS_G1_1.coords.GENCODDEv39.bed.gz myanco@cluster.s3it.uzh.ch:/scratch/myanco/cpdsq/
# any genes on noncanon chr's?
gunzip -cd RNAsq.U2OS_G1_1.coords.GENCODDEv39.bed.gz | awk 'NR>1' | cut -f 1 | sort -k1,1V | uniq -c
# yup a bunch.
gunzip -cd RNAsq.U2OS_G1_1.coords.GENCODDEv39.bed.gz | awk 'NR>1' | cut -f 1 | sort -k1,1V | uniq -c | grep 'alt\|fix\|Y\|random'
# | tr '[ ]' '[\t]' | sed 's/ /\t/g'
# total 277 lines without filtering, 254 lines of which are nncan and chY.


# 11781 hits total to these noncanonical chr's plus chrY.


gunzip -cd RNAsq.U2OS_G1_1.coords.GENCODDEv39.bed.gz | awk 'NR>1' | cut -f 1 | sort -k1,1V | uniq -c | grep -v 'alt\|fix\|Y\|random'
# 153,990  hits, far outweigh 12k on noncan/chrY.


#gunzip -cd RNAsq.U2OS_G1_1.coords.gz | awk 'NR>1' | cut -f 1  | sort -k1,1 | sed 's/\./\t/g' | cut -f 1 | sort -k1,1 | uniq | wc -l
#61641 ...

## align ID 12th col seems to be same as 1st col..
#gunzip -cd RNAsq.U2OS_G1_1.coords.gz | cut -f 12 | sed 's/\./\t/g' | cut -f 1| sort -k1,1n | uniq | wc -l
# 61,642


/data/myanco/accessory_files/regions/genic.hg19.bed
awk 'NR>1' /data/myanco/accessory_files/regions/genic.hg19.bed | wc -l
# 82960

# filter out the lines to noncanonical chromosomes??!

# https://gatk.broadinstitute.org/hc/en-us/articles/360035890951-Human-genome-reference-builds-GRCh38-or-hg38-b37-hg19

# The latest build of the human reference genome, officially named GRCh38 (for Genome Research Consortium human build 38) but commonly nicknamed Hg38 (for Human genome build 38), greatly expanded the repertoire of ALT contigs. These represent alternate haplotypes and have a significant impact on our power to detect and analyze genomic variation that is specific to populations that carry alternate haplotypes. Read on further below for more details.
# 
# We strongly recommend switching to GRCh38/hg38 if you are working with human sequence data. In addition to adding many alternate contigs, GRCh38 corrects thousands of small sequencing artifacts that cause false SNPs and indels to be called when using the GRCh37 assembly (b37/Hg19). It also includes synthetic centromeric sequence and updates non-nuclear genomic sequence.

# Unlocalized sequences (known to belong on a specific chromosome but with unknown order or orientation) are identified by the _random suffix.

# Unplaced sequences (chromosome of origin unknown) are identified by the chrU_ prefix.

# # The GRCh38 ALT contigs are recognizable by their _alt suffix; they amount to a total of 109Mb in length and span 60Mb of the primary assembly. Alternate contig sequences can be novel to highly diverged or nearly identical to corresponding primary assembly sequence. Sequences that are highly diverged from the primary assembly only contribute a few million bases. Most subsequences of ALT contigs are fairly similar to the primary assembly. This means that if we align sequence reads to GRCh38+ALT blindly, then we obtain many multi-mapping reads with zero mapping quality. Since many GATK tools have a ZeroMappingQuality filter, we will then miss variants corresponding to such loci.



# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4347522/
# Typically, Hi-C data is mapped to a known high-quality genome sequence and is used to answer questions regarding the 3D organization of genomes. However, it has recently been shown in a number of studies that Hi-C data can be useful to learn about the 1D arrangement of the genome sequence and thus solve a number of outstanding problems in the field of genome assembly
#Interactions of loci located in different nuclei are less frequent than those in the same nucleus. This principle seems obvious, but has important implications. In microbiome studies, which analyze large mixed populations of different species, high-throughput sequencing typically yields a large set of contigs, yet it is difficult to establish which contigs belong to the same genome. Using Hi-C data, we can determine that if two contigs interact frequently in 3D they are likely to belong to the same genome with high probability (Burton et al. 2014; Beitel et al. 2014).

#Interactions of loci located on different chromosomes are less frequent than those in the same chromosome. As discussed above, this pattern is both strong and ubiquitous. When performing de novo genome scaffolding, we can thus use Hi-C data to determine that contigs that interact frequently are likely to belong to the same chromosome 

#Interactions of loci located far from each other along a chromosome are less frequent than loci that are near each other. Using Hi-C data, we can arrange contigs which belong to the same chromosome such that strongly interaction contigs are positioned next to each other (Burton et al. 2013; Kaplan and Dekker 2013).



# filter out the non-canonical chr hits; grab the TSS ±1.5kb around the start TSS coordinate.

# 

#Can the coding strand of one gene be the template strand for a different gene?
#DNA is double-stranded, but only one strand serves as a template for transcription at any given time. This template strand is called the noncoding strand. ... In most organisms, the strand of DNA that serves as the template for one gene may be the nontemplate strand for other genes within the same chromosome.

# bedtools complement to find intergenic regions: retain strand info?
# https://www.biostars.org/p/117925/

# Taking the complement of both strands separately worked for me.
grep -v '+$' sorted_gff_file | bedtools complement -i stdin -g genome_file | sed 's/$/\t-/' > complement_file

grep '+$' sorted_gff_file | bedtools complement -i stdin -g genome_file | sed 's/$/\t+/' >> complement_file



 gunzip -cd RNAsq.U2OS_G1_1.coords.GENCODDEv39.bed.gz | head -n 2
#chrom	chromStart	chromEnd	name	score	strand	thickStart	thickEnd	reserved	blockCount	blockSizes	chromStarts	name2	cdsStartStat	cdsEndStat	exonFrames	type	geneName	geneName2	geneType	transcriptClass	source	transcriptType	tag	level	tier

#chr22	31468281	31496108	ENST00000397520.1	0	-	31468281	31488718	789624	4	21,74,123,156,	0,3562,20341,27671,	uc062djl.1	none	none	2,0,0,-1,	none	EIF4ENIF1	B1AKL6	none	coding	havana_homo_sapiens	protein_coding	cds_end_NF,mRNA_end_NF,ncRNA_host	2	all

# col1-3 = chr, start, stop of each gene; col 6 = strand, col7-8 = thickStart thickEnd 
# thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). When there is no thick part, thickStart and thickEnd are usually set to the chromStart position.
# thickEnd - The ending position at which the feature is drawn thickly (for example the stop codon in gene displays).

# Schema for GENCODE V39 - GENCODE V39
#thickStart	166022214	Start of where display should be thick (start codon)
#thickEnd	166022214	End of where display should be thick (stop codon)


gunzip -cd /scratch/myanco/cpdsq/RNAsq.U2OS_G1_1.coords.GENCODDEv39.bed.gz | awk 'NR>1' | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"}'



gunzip -cd /scratch/myanco/cpdsq/RNAsq.U2OS_G1_1.coords.GENCODDEv39.bed.gz | cut -f 1-3,6-8 | awk '{print $3-$2"\t"$6-$5}' | awk '{print $1-$2}' | awk 'NR>1' | wc -l
#165771
gunzip -cd /scratch/myanco/cpdsq/RNAsq.U2OS_G1_1.coords.GENCODDEv39.bed.gz | cut -f 1-3,6-8 | awk '{print $3-$2"\t"$6-$5}' | awk '{print $1-$2}' | awk 'NR>1' | awk '{if ($1>0) print $0}' | wc -l
#162123

# for about 3k entries, the length of the coding region is longer than the length of the transcript..

# get GENCODE v39 genome-wide gene coords from UCSC. # genic.hg19.bed = hg19 version

#sun 1/17/21 --> #for the intergenic regions, i dl'd ucsc table browser track group: genes and gene predictions, table: knownGene --> gwide.
#
assembly: hg38
group:   Genes and Gene Predictions  track:     GENCODE V39
table:   knownGene

# dont upload any list just get everyone
# output format: all fields from selected table
# output file name: knownGenes.UCSC.hg38.bed

 scp knownGenes.UCSC.hg38.bed myanco@cluster.s3it.uzh.ch:/scratch/myanco/cpdsq/
 
 cp /scratch/myanco/cpdsq/knownGenes.UCSC.hg38.bed /data/myanco/accessory_files/hg38/
 
 # STOPPED HERE 3/30/22
 
 
gunzip -cd /scratch/myanco/cpdsq/RNAsq.U2OS_G1_1.coords.GENCODDEv39.bed.gz | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6}' | 

# all the gene id's start with ENST --> also found in the UCSC file?


# active U2OS genes including their full coordinates: do a bedtools intersect with the known genes list, retain the lines from the known genes list that (1) overlap and (2) dont overlap, then extract the tss coordinates and subtract them from the coding sequence coords.


# knownGenes file , knownGenes.UCSC.hg38.bed,  has thickStart but not tx start.
# re-dl ...
# Table Bwoser: hg38, Genes and Gene Predictions, trck=GENCODE v39, table = knownGene

rm knownGenes.UCSC.hg38.bed

scp knowngene.hg38.GENCODEv39.gz myanco@cluster.s3it.uzh.ch:/data/myanco/accessory_files/hg38/

gunzip -cd knowngene.hg38.GENCODEv39.gz  | head
 # #name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	proteinID	alignID
 # ENST00000456328.2	chr1	+	11868	14409	11868	11868	3	11868,12612,13220,	12227,12721,14409,		uc286dmu.1




# u2os-activgenes.hg38
bedtools intersected list of u2os genes, //cpdsq/geneNames.PosReadCount.u2osG1_1.uniq.bed, with ucsc gencode v39 / hg38 knownGene table --> intersect: pull down lines with any overlap.
# output: u2os-activgenes.hg38.bed.gz

# next intersect: lines with no overlap.
# create one BED record per: Whole Gene

# output: u2os-INactivgenes.hg38.bed.gz --> doesnt work, empty..

gunzip -cd knowngene.hg38.GENCODEv39.gz | wc -l
# 266065

gunzip -cd knowngene.hg38.GENCODEv39.gz | cut -f 2,4-5 | sort -k1,1V -k2,2n | bedtools merge -i stdin | wc -l
# 37,630


gunzip -cd u2os-activgenes.hg38.bed.gz | wc -l
#  165772
gunzip -cd u2os-activgenes.hg38.bed.gz | cut -f 2,4-5 | sort -k1,1V -k2,2n | bedtools merge -i stdin | wc -l
#   15,350


# bedt int should yield about 266-166 or 100k lines of genes w no overlap, assuming each line is a nonoverlapping region of the genome.

bedtools intersect -a <(gunzip -cd knowngene.hg38.GENCODEv39.gz | awk 'NR>1' | awk '{print $2"\t"$4"\t"$5"\t"$2":"$6"-"$7}' | sort -k1,1V -k2,2n) -b <(gunzip -cd u2os-activgenes.hg38.bed.gz | awk 'NR>1' | awk '{print $2"\t"$4"\t"$5"\t"$2":"$6"-"$7}' | sort -k1,1V -k2,2n) -v >knowngene.hg38.notoverlapwU2OSactivGeneTxStartStop.bed
# 68k lines not merged, 26k when merged.
# u2os active genes covers 10,127 Mb when not merged, 1839Mb when merged.

# OR do bedt subtract.
#bedtools subtract -a <(gunzip -cd knowngene.hg38.GENCODEv39.gz | awk 'NR>1' | awk '{print $2"\t"$4"\t"$5"\t"$2":"$6"-"$7}' | sort -k1,1V -k2,2n ) -b <(gunzip -cd u2os-activgenes.hg38.bed.gz | awk 'NR>1' | awk '{print $2"\t"$4"\t"$5"\t"$2":"$6"-"$7}' | sort -k1,1V -k2,2n) -v >knowngene.hg38.notoverlapwU2OSactivGeneTxStartStop.bed

rm knowngene.hg38.notoverlapwU2OSactivGeneTxStartStop.bed



# TSS: ±1`.5kb from the tx start site
# how to know whether to start from the txstart or txend coord? --> need strand info!!
# put the u2os coords of u2os activ genes back into bedt overlap w knowngenes and retain strand info?

#	-s	Require same strandedness.  That is, only report hits in B
#		that overlap A on the _same_ strand.
#		- By default, overlaps are reported without respect to strand.

#	-S	Require different strandedness.  That is, only report hits in B
#		that overlap A on the _opposite_ strand.
#		- By default, overlaps are reported without respect to strand.




gunzip -cd u2os-activgenes.hg38.bed.gz | head


scp -r ~/Downloads/u2os-activgenes.hg38.bed.gz myanco@cluster.s3it.uzh.ch:/data/myanco/accessory_files/hg38/



# find regions that do not overlap with hg38:

gunzip -cd knowngene.hg38.GENCODEv39.gz | cut -f 2,4-7 | head
#chrom	txStart	txEnd	cdsStart	cdsEnd
#chr1	11868	14409	11868	11868

gunzip -cd knowngene.hg38.GENCODEv39.gz | awk '{print $2"\t"$4"\t"$5"\t"$2":"$6"-"$7}' | awk 'NR>1' | head
# chr1	11868	14409	chr1:11868-11868
# chr1	12009	13670	chr1:12009-12009

bedtools intersect -a <(gunzip -cd knowngene.hg38.GENCODEv39.gz | awk '{print $2"\t"$4"\t"$5"\t"$2":"$6"-"$7}' | awk 'NR>1' | sort -k1,1V -k2,2n) -b <(gunzip -cd u2os-activgenes.hg38.bed.gz | awk '{print $2"\t"$4"\t"$5"\t"$2":"$6"-"$7}' | awk 'NR>1' | sort -k1,1V -k2,2n) -v
# 68k lines
# if merge the ouptput it becomes 26k lines

gunzip -cd knowngene.hg38.GENCODEv39.gz | wc -l
# 266k lines, when merged by txStart-End it is 37k lines.
gunzip -cd u2os-activgenes.hg38.bed.gz | wc -l
# 166k lines, when merged by txStart-End it is 15k lines


# should the final bedt int -v output be bedt merged? diff lines correspond to diff possibl overlapping transcripts...?

#with bedt merging: total 37k , 15k u2os active genes, 26k regions that dont overlap.


bedtools intersect -a <(gunzip -cd knowngene.hg38.GENCODEv39.gz | awk '{print $2"\t"$4"\t"$5"\t"$3"\t"$2":"$6"-"$7}' | awk 'NR>1' | sort -k1,1V -k2,2n) -b <(gunzip -cd u2os-activgenes.hg38.bed.gz | awk '{print $2"\t"$4"\t"$5"\t"$3"\t"$2":"$6"-"$7}' | awk 'NR>1' | sort -k1,1V -k2,2n) -v
# retain strand info for the sake of getting the TSS?
# but dnt factor it into the bedt overlap? if a u2os gene doesnt overlap the known genes regardless of strand..

bedtools intersect -a <(gunzip -cd knowngene.hg38.GENCODEv39.gz | awk '{print $2"\t"$4"\t"$5"\t"$3"\t"$2":"$6"-"$7}' | awk 'NR>1' | sort -k1,1V -k2,2n) -b <(gunzip -cd u2os-activgenes.hg38.bed.gz | awk '{print $2"\t"$4"\t"$5"\t"$3"\t"$2":"$6"-"$7}' | awk 'NR>1' | sort -k1,1V -k2,2n) -v >knowngene.hg38.notoverlapwU2OSactivGeneTxStartStop.bed

# are u2os activ genes and u2os nonoverlap w u2os cative gen contain no ovlps?
bedtools intersect -a knowngene.hg38.notoverlapwU2OSactivGeneTxStartStop.bed -b <(gunzip -cd u2os-activgenes.hg38.bed.gz | awk '{print $2"\t"$4"\t"$5"\t"$3"\t"$2":"$6"-"$7}')

bedtools intersect -b knowngene.hg38.notoverlapwU2OSactivGeneTxStartStop.bed -a <(gunzip -cd u2os-activgenes.hg38.bed.gz | awk '{print $2"\t"$4"\t"$5"\t"$3"\t"$2":"$6"-"$7}')
# no overlaps in either direction.

awk '{if ($4!="-") print $0}' knowngene.hg38.notoverlapwU2OSactivGeneTxStartStop.bed | wc -l 
# 35807
awk '{if ($4=="-") print $0}' knowngene.hg38.notoverlapwU2OSactivGeneTxStartStop.bed | wc -l
# 32417



gunzip -cd knowngene.hg38.GENCODEv39.gz | awk '{print $2"\t"$4"\t"$5"\t"$3"\t"$2":"$6"-"$7}' | awk '{if ($4!="-") print $0}' | wc -l
# 136183

gunzip -cd knowngene.hg38.GENCODEv39.gz | awk '{print $2"\t"$4"\t"$5"\t"$3"\t"$2":"$6"-"$7}' | awk '{if ($4=="-") print $0}' | wc -l
# 129882
gunzip -cd knowngene.hg38.GENCODEv39.gz | awk '{print $2"\t"$4"\t"$5"\t"$3"\t"$2":"$6"-"$7}' | awk '{if ($4!="+") print $0}' | wc -l
# 129882

gunzip -cd knowngene.hg38.GENCODEv39.gz | wc -l
# 266065 total line count



awk '{if ($4=="+") print $1"\t"$2-1500"\t"$2+1500"\t"$4"\t"$5}' knowngene.hg38.notoverlapwU2OSactivGeneTxStartStop.bed | awk '{if ($2>0) print $0}' | sort -k1,1V -k2,2n >inactvTSS.plus-strand.tmp
# about 35 lines have neg start coords.

awk '{if ($4=="-") print $1"\t"$3-1500"\t"$3+1500"\t"$4"\t"$5}' knowngene.hg38.notoverlapwU2OSactivGeneTxStartStop.bed | awk '{if ($2>0) print $0}' | sort -k1,1V -k2,2n >inactvTSS.neg-strand.tmp


gunzip -cd u2os-activgenes.hg38.bed.gz | awk 'NR>1' | awk '{print $2"\t"$4"\t"$5"\t"$3"\t"$2":"$6"-"$7}' | awk '{if ($4=="+") print $1"\t"$2-1500"\t"$2+1500"\t"$4"\t"$5}' | cut -f 2 | sort -k1,1g | head
# about 80 start TSS coords are neg...
# 2 start TSS coords from neg strand are neg.

gunzip -cd u2os-activgenes.hg38.bed.gz | awk 'NR>1' | awk '{print $2"\t"$4"\t"$5"\t"$3"\t"$2":"$6"-"$7}' | awk '{if ($4=="+") print $1"\t"$2-1500"\t"$2+1500"\t"$4"\t"$5}' | awk '{if ($2>0) print $0}' | sort -k1,1V -k2,2n >actvTSS.plus-strand.tmp


gunzip -cd u2os-activgenes.hg38.bed.gz | awk 'NR>1' | awk '{print $2"\t"$4"\t"$5"\t"$3"\t"$2":"$6"-"$7}' | awk '{if ($4=="-") print $1"\t"$3-1500"\t"$3+1500"\t"$4"\t"$5}' | awk '{if ($2>0) print $0}' | sort -k1,1V -k2,2n >actvTSS.neg-strand.tmp


wc -l inactvTSS.plus-strand.tmp
#35776
awk '{if ($2>0) print $0}' inactvTSS.plus-strand.tmp | bedtools merge -i stdin | wc -l
# 18246


wc -l inactvTSS.neg-strand.tmp
# 32410
bedtools merge -i inactvTSS.neg-strand.tmp | wc -l
# 17864


# 4/8/22
# subtract the 3kb TSS regions from the transcripton start thru stop coords to get the gene bodies..

# dont include strand info?
gunzip -cd u2os-activgenes.hg38.bed.gz | awk 'NR>1' | awk '{print $2"\t"$4"\t"$5"\t"$2":"$6"-"$7}' | sort -k1,1V -k2,2n >u2os-actvGenes.hg38.bed4

# 
gene bodies: subtract out all tss not just active tss??
# or: strand info retain, subtract out tss on same strand as gene?
# to make the regions mutually exclusive, would need to subtract out all the TSS.

cat inactvTSS.plus-strand.tmp inactvTSS.neg-strand.tmp actvTSS.plus-strand.tmp actvTSS.neg-strand.tmp | sort -k1,1V -k2,2n >allTSS.3kb.tmp

cut -f 1-3 allTSS.3kb.tmp | bedtools merge -i stdin | awk '{print $3-$2}' | awk '{total+=$1}END{print total/1e6}'
# 277 Mb

cut -f 1-3 allTSS.3kb.tmp | awk '{print $3-$2}' | awk '{total+=$1}END{print total/1e6}'
# 701 Mb without merging or combining across strands.

cut -f 1-3 allTSS.3kb.tmp | uniq | awk '{print $3-$2}' | awk '{total+=$1}END{print total/1e6}'
# 610 Mb if tkae unique bins


# gne bodies, subtract out tss: dont use any options w bedt int, no -wa or -wb or -wao.-> this should fragment up the a bins by overlapped regions.
# test:
nano a.tmp
chr1	1	40
chr1	100	120


nano b.tmp
chr1	2	3
chr1	39	110

expect:
chr1	1	2
chr1	4	38
chr1	111	120


bedtools intersect -a a.tmp -b b.tmp
#chr1	2	3
#chr1	39	40
#chr1	100	110

bedtools complement -i b.tmp -g a.tmp
#***** WARNING: chr1:39-110 exceeds the length of chromosome (chr1)
#chr1	0	2
#chr1	3	39
#chr1	0	100



bedtools subtract -a a.tmp -b b.tmp
chr1	1	2
chr1	3	39
chr1	110	120




bedtools subtract -a u2os-actvGenes.hg38.bed4 -b <( cut -f 1-3 allTSS.3kb.tmp | uniq ) | head
#chr1	13509	14409	chr1:11868-11868
#chr1	13509	13670	chr1:12009-12009
#chr1	14403	28070	chr1:14403-14403
#chr1	818870	819837	chr1:817370-817370
#chr1	829320	830077	chr1:825137-825137

# 332,543 lines; some clearly merge-able

bedtools subtract -a u2os-actvGenes.hg38.bed4 -b <( cut -f 1-3 allTSS.3kb.tmp | uniq ) | sort -k1,1V -k2,2n | bedtools merge -i stdin | wc -l
# 45,018

bedtools subtract -a u2os-actvGenes.hg38.bed4 -b <( cut -f 1-3 allTSS.3kb.tmp | uniq ) | sort -k1,1V -k2,2n | bedtools merge -i stdin | awk '{print $3-$2}' | awk '{total+=$1}END{print total/1e6}'
# 1106.5
# u2os active gene bodies cover 1106 Mb


bedtools subtract -a u2os-actvGenes.hg38.bed4 -b <( cut -f 1-3 allTSS.3kb.tmp | uniq ) | sort -k1,1V -k2,2n | bedtools merge -i stdin >u2os-activGeneBods.hg38.bed



bedtools subtract -a knowngene.hg38.notoverlapwU2OSactivGeneTxStartStop.bed -b <(cut -f 1-3 allTSS.3kb.tmp | uniq) | sort -k1,1V -k2,2n | head
# chr1	34553	34573	-	chr1:34553-34553
#chr1	59097	61448	+	chr1:57597-57597


bedtools subtract -a knowngene.hg38.notoverlapwU2OSactivGeneTxStartStop.bed -b <(cut -f 1-3 allTSS.3kb.tmp | uniq) | sort -k1,1V -k2,2n | wc -l
# 92,255

bedtools subtract -a knowngene.hg38.notoverlapwU2OSactivGeneTxStartStop.bed -b <(cut -f 1-3 allTSS.3kb.tmp | uniq) | sort -k1,1V -k2,2n | bedtools merge -i stdin | wc -l
# 22,,238 --> about half nr frm active gne bodies.
# hg19: also with this derivation, active gene bodies dominated inactive gene bodies in U2OS.
bedtools subtract -a knowngene.hg38.notoverlapwU2OSactivGeneTxStartStop.bed -b <(cut -f 1-3 allTSS.3kb.tmp | uniq) | sort -k1,1V -k2,2n | bedtools merge -i stdin | awk '{print $3-$2}' | awk '{total+=$1}END{print total/1e6}'
# 439.355 Mb

bedtools subtract -a knowngene.hg38.notoverlapwU2OSactivGeneTxStartStop.bed -b <(cut -f 1-3 allTSS.3kb.tmp | uniq) | sort -k1,1V -k2,2n | bedtools merge -i stdin >u2os-inactvGeneBods.hg38.bed

cat actvTSS.plus-strand.tmp actvTSS.neg-strand.tmp | sort -k1,1V -k2,2n | bedtools merge -i stdin >u2os-actvTSS.3kb.hg38.bed
# 165,698 lines; bedt merge 45,560
# 170.08 Mb covered by active TSS:
awk '{print $3-$2}' u2os-actvTSS.3kb.hg38.bed | awk '{total+=$1}END{print total/1e6}'
#170.08

cat inactvTSS.plus-strand.tmp inactvTSS.neg-strand.tmp | sort -k1,1V -k2,2n | bedtools merge -i stdin >u2os-inactvTSS.3kb.hg38.bed
# 68,186 lines; 34,045 lines bedt merged.
# 112 Mb covered by inactive TSS
awk '{print $3-$2}' u2os-inactvTSS.3kb.hg38.bed | awk '{total+=$1}END{print total/1e6}'
#112.784


# check that all categories are mututally exclsiv:
bedtools intersect -a u2os-activGeneBods.hg38.bed -b u2os-inactvGeneBods.hg38.bed | head #no overlps, same reciprocall

bedtools intersect -a u2os-activGeneBods.hg38.bed -b u2os-actvTSS.3kb.hg38.bed | head 
# no ovlps, same reciprocal.

bedtools intersect -b u2os-activGeneBods.hg38.bed -a u2os-inactvTSS.3kb.hg38.bed | head 

bedtools intersect -a u2os-inactvGeneBods.hg38.bed -b u2os-inactvTSS.3kb.hg38.bed | head 

bedtools intersect -b u2os-inactvGeneBods.hg38.bed -a u2os-actvTSS.3kb.hg38.bed | head 


# no overlaps.
# bedt mergeable? --> nobody is mergeable.
mv u2os-activGeneBods.hg38.bed u2os-actvGenBods.hg38.mrg.bed
mv u2os-inactvGeneBods.hg38.bed u2os-inactvGenBods.hg38.mrg.bed
mv u2os-actvTSS.3kb.hg38.bed u2os-actvTSS.3kb.hg38.mrg.bed
mv u2os-inactvTSS.3kb.hg38.bed u2os-inactvTSS.3kb.hg38.mrg.bed


awk '{print $3-$2}' u2os-actvTSS.3kb.hg38.mrg.bed | sort -k1,1g | tail
# longest merged actv tss region is 33kb, gene-dense region?
awk '{print $3-$2}' u2os-inactvTSS.3kb.hg38.mrg.bed | sort -k1,1g | tail
#75kb

	# next: concat every genic region and then bedt subtract from rest of hg38?
# or do sequential several subtraction steps

cat u2os-actvGenBods.hg38.mrg.bed u2os-inactvGenBods.hg38.mrg.bed u2os-actvTSS.3kb.hg38.mrg.bed u2os-inactvTSS.3kb.hg38.mrg.bed | wc -l
# 146,861
cat u2os-actvGenBods.hg38.mrg.bed u2os-inactvGenBods.hg38.mrg.bed u2os-actvTSS.3kb.hg38.mrg.bed u2os-inactvTSS.3kb.hg38.mrg.bed | sort -k1,1V -k2,2n | bedtools merge -i stdin | wc -l
# 34399
# increase distance to 1 for bedt merge :
cat u2os-actvGenBods.hg38.mrg.bed u2os-inactvGenBods.hg38.mrg.bed u2os-actvTSS.3kb.hg38.mrg.bed u2os-inactvTSS.3kb.hg38.mrg.bed | sort -k1,1V -k2,2n | bedtools merge -i stdin -d -1 | wc -l
# 144,466

# filter out noncanon chromosoems from genomic features!!!
for i in *mrg.bed; do mv -- "$i" "${i/mrg.bed/noncanon.mrg.bed}"; done

cat u2os-actvGenBods.hg38.noncanon.mrg.bed u2os-inactvGenBods.hg38.noncanon.mrg.bed u2os-actvTSS.3kb.hg38.noncanon.mrg.bed  u2os-inactvTSS.3kb.hg38.noncanon.mrg.bed  | cut -f 1 | sort -k1,1V | uniq

for f in *noncanon.mrg.bed; do grep -v 'fix\|chrUn\|chrM\|alt\|random' $f | cut -f 1 | sort -k1,1V | uniq ; done
# only canon left.

for f in *noncanon.mrg.bed; do grep -v 'fix\|chrUn\|chrM\|alt\|random\|chrY' $f >$f.canon; done

for i in *noncanon.mrg.bed.canon; do mv -- "$i" "${i/noncanon.mrg.bed.canon/canon.mrg.bed}"; done

for f in *.canon.mrg.bed; do cut -f 1 $f | uniq | wc -l; done
#23
#23
#23
#23


# or bedt complement 
bedtools complement -i <(cat u2os-actvGenBods.hg38.canon.mrg.bed u2os-inactvGenBods.hg38.canon.mrg.bed u2os-actvTSS.3kb.hg38.canon.mrg.bed u2os-inactvTSS.3kb.hg38.canon.mrg.bed | sort -k1,1V -k2,2n) -g /data/myanco/accessory_files/hg38/GCA_hg38.chrom_chr.nochrY.canon.sizes | wc -l
# 29,948

bedtools complement -i <(cat u2os-actvGenBods.hg38.canon.mrg.bed u2os-inactvGenBods.hg38.canon.mrg.bed u2os-actvTSS.3kb.hg38.canon.mrg.bed u2os-inactvTSS.3kb.hg38.canon.mrg.bed | sort -k1,1V -k2,2n) -g /data/myanco/accessory_files/hg38/GCA_hg38.chrom_chr.nochrY.canon.sizes | bedtools merge -i stdin | wc -l
 #29,948

bedtools complement -i <(cat u2os-actvGenBods.hg38.canon.mrg.bed u2os-inactvGenBods.hg38.canon.mrg.bed u2os-actvTSS.3kb.hg38.canon.mrg.bed u2os-inactvTSS.3kb.hg38.canon.mrg.bed | sort -k1,1V -k2,2n) -g /data/myanco/accessory_files/hg38/GCA_hg38.chrom_chr.nochrY.canon.sizes | sort -k1,1V -k2,2n | awk '{print $3-$2}' | sort -k1,1g | awk '{total +=$1}END{print total/1e6}'
# 1293.97 Mb intergenic, some are just 1bp...up to 21765977 bp or 21Mb.


bedtools complement -i <(cat u2os-actvGenBods.hg38.canon.mrg.bed u2os-inactvGenBods.hg38.canon.mrg.bed u2os-actvTSS.3kb.hg38.canon.mrg.bed u2os-inactvTSS.3kb.hg38.canon.mrg.bed | sort -k1,1V -k2,2n) -g /data/myanco/accessory_files/hg38/GCA_hg38.chrom_chr.nochrY.canon.sizes | sort -k1,1V -k2,2n >u2os-igenic.hg38.canon.mrg.bed

# cat up all the regions and bedt merge: shld have just 23 lins
cat u2os-actvGenBods.hg38.canon.mrg.bed u2os-inactvGenBods.hg38.canon.mrg.bed u2os-actvTSS.3kb.hg38.canon.mrg.bed u2os-inactvTSS.3kb.hg38.canon.mrg.bed u2os-igenic.hg38.canon.mrg.bed | sort -k1,1V -k2,2n | bedtools merge -i stdin | wc -l
# 23

cat u2os-actvGenBods.hg38.canon.mrg.bed u2os-inactvGenBods.hg38.canon.mrg.bed u2os-actvTSS.3kb.hg38.canon.mrg.bed u2os-inactvTSS.3kb.hg38.canon.mrg.bed u2os-igenic.hg38.canon.mrg.bed | sort -k1,1V -k2,2n | bedtools merge -i stdin | awk '{print $3-$2}' | awk '{total +=$1}END{print total/1e6}'
# 3031.04

awk '{print $3-$2}' GCA_hg38.chrom_chr.nochrY.canon.sizes | awk '{total+=$1}END{print total/1e6}'
# -3031.04

# check that the igenic doesnt overlap w the tss or gene bods:
bedtools intersect -b u2os-igenic.hg38.canon.mrg.bed -a u2os-actvGenBods.hg38.canon.mrg.bed 

bedtools intersect -b u2os-igenic.hg38.canon.mrg.bed -a u2os-inactvGenBods.hg38.canon.mrg.bed 

bedtools intersect -a u2os-igenic.hg38.canon.mrg.bed -b u2os-inactvTSS.3kb.hg38.canon.mrg.bed 

bedtools intersect -a u2os-igenic.hg38.canon.mrg.bed -b u2os-actvTSS.3kb.hg38.canon.mrg.bed 
# everyoe is mututally exclvs.


# divide igenic into tads vs itads:
 cat 360min_rep12cat.TADs.mrg.hg38.bed iTADcatrep12.bedtcompl.hg38-nochrY-canon.bed  | sort -k1,1V -k2,2n | bedtools merge -i stdin | wc -l
#23

bedtools subtract -a u2os-igenic.hg38.canon.mrg.bed -b /data/myanco/accessory_files/HiC_U2OS_hg38/iTADcatrep12.bedtcompl.hg38-nochrY-canon.bed | sort -k1,1V -k2,2n | wc -l
# 25,926 lines ;cant be merged.

bedtools subtract -a u2os-igenic.hg38.canon.mrg.bed -b /data/myanco/accessory_files/HiC_U2OS_hg38/360min_rep12cat.TADs.mrg.hg38.bed | sort -k1,1V -k2,2n | wc -l
# 4,615 lines; cant be further merged.


bedtools subtract -a u2os-igenic.hg38.canon.mrg.bed -b /data/myanco/accessory_files/HiC_U2OS_hg38/iTADcatrep12.bedtcompl.hg38-nochrY-canon.bed | sort -k1,1V -k2,2n >u2os-igenic.euchr.hg38.canon.mrg.bed

bedtools subtract -a u2os-igenic.hg38.canon.mrg.bed -b /data/myanco/accessory_files/HiC_U2OS_hg38/360min_rep12cat.TADs.mrg.hg38.bed | sort -k1,1V -k2,2n >u2os-igenic.heterochr.hg38.canon.mrg.bed

awk '{print $3-$2}' u2os-igenic.hg38.canon.mrg.bed | awk '{total +=$1}END{print total/1e6}'
#1293.97

cat u2os-igenic.euchr.hg38.canon.mrg.bed u2os-igenic.heterochr.hg38.canon.mrg.bed | awk '{print $3-$2}' | awk '{total+=$1}END{print total/1e6}'
#  1293.97

bedtools intersect -b u2os-igenic.euchr.hg38.canon.mrg.bed -a u2os-igenic.heterochr.hg38.canon.mrg.bed 


bedtools intersect -b u2os-igenic.euchr.hg38.canon.mrg.bed -a u2os-actvGenBods.hg38.canon.mrg.bed 

bedtools intersect -b u2os-igenic.heterochr.hg38.canon.mrg.bed -a u2os-actvGenBods.hg38.canon.mrg.bed 


cat u2os-igenic.euchr.hg38.canon.mrg.bed u2os-igenic.heterochr.hg38.canon.mrg.bed u2os-inactvTSS.3kb.hg38.canon.mrg.bed u2os-inactvGenBods.hg38.canon.mrg.bed u2os-actvTSS.3kb.hg38.canon.mrg.bed u2os-actvGenBods.hg38.canon.mrg.bed | sort -k1,1V -k2,2n | bedtools merge -i stdin | wc -l
# 23
cat u2os-igenic.euchr.hg38.canon.mrg.bed u2os-igenic.heterochr.hg38.canon.mrg.bed u2os-inactvTSS.3kb.hg38.canon.mrg.bed u2os-inactvGenBods.hg38.canon.mrg.bed u2os-actvTSS.3kb.hg38.canon.mrg.bed u2os-actvGenBods.hg38.canon.mrg.bed | sort -k1,1V -k2,2n | bedtools merge -i stdin | awk '{print $3-$2}' | awk '{total+=$1}END{print total/1e6}'
# 3031.04

awk '{print $3-$2}' GCA_hg38.chrom_chr.nochrY.canon.sizes | awk '{total+=$1}END{print total/1e6}'
# -3031.04

awk '{print $3-$2}' u2os-igenic.hg38.canon.mrg.bed | awk '{total+=$1}END{print total/1e6}'
1293.97

cat u2os-actvGenBods.hg38.canon.mrg.bed u2os-actvTSS.3kb.hg38.canon.mrg.bed u2os-inactvGenBods.hg38.canon.mrg.bed u2os-inactvTSS.3kb.hg38.canon.mrg.bed |  sort -k1,1V -k2,2n | bedtools merge -i stdin | wc -l
#29925; without merging 134505
cat u2os-actvGenBods.hg38.canon.mrg.bed u2os-actvTSS.3kb.hg38.canon.mrg.bed u2os-inactvGenBods.hg38.canon.mrg.bed u2os-inactvTSS.3kb.hg38.canon.mrg.bed |  sort -k1,1V -k2,2n | bedtools merge -i stdin | awk '{print $3-$2}' | awk '{total +=$1}END{print total/1e6}'  # 1737.07   # is exactly 100 pct minus pctg taken up by igenic.
cat u2os-actvGenBods.hg38.canon.mrg.bed u2os-actvTSS.3kb.hg38.canon.mrg.bed u2os-inactvGenBods.hg38.canon.mrg.bed u2os-inactvTSS.3kb.hg38.canon.mrg.bed |  sort -k1,1V -k2,2n | awk '{print $3-$2}' | awk '{total +=$1}END{print total/1e6}'   # slghtly higher witout merging, 1742.07 ...


awk '{print $3-$2}' u2os-igenic.euchr.hg38.canon.mrg.bed

# u2os-genic.fts.lengths.txt
feature	Mb	PctgOfGenm
wholegenome	3031.04	100
genic	1737.07	57.309372360641895
ignc	1293.97	42.690627639358105

igncEuchr	
igncHetchr	
actvTSS	
actvGenBods	
inactvTSS	
inactvGenBods	



# ========================
# active + inactive genes
# ========================
	# U2OS active gene coordinates:
gunzip -cd ~/Downloads/RNAsq.U2OS_G1_1.coords.gz | awk '{print $2"\t"$4"\t"$5"\t"$1"\t"$3}' >activ.genic.txstrt-end.hg19.bed
	
	#gunzip -cd RNAsq.U2OS_G1_1.coords.gz| head -n 2
	#name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	proteinID	alignID
	#uc001aaa.3	chr1	+	11873	14409	11873	11873	3	11873,12612,13220,	12227,12721,14409,		uc001aaa.3
	
	# printed chr, txStart, txEnd,#gene name, strand.
#	gunzip -cd RNAsq.U2OS_G1_1.coords.gz| head -n 2 |  awk '{print $2"\t"$4"\t"$5"\t"$1"\t"$3}'
		#chrom	txStart	txEnd	#name	strand
		#chr1	11873	14409	uc001aaa.3	+



	
	# U2OS inactive gene coordinates: bedt subtract active genes from genic.hg19.bed genome-wide gene coords.
	bedtools subtract -a <(awk 'NR>1' /data/myanco/accessory_files/regions/genic.hg19.bed | awk '{print $2"\t"$4"\t"$5"\t""0""\t""0""\t"$3}' | sort -k1,1V -k2,2n | uniq) -b <(awk 'NR>1' activ.genic.txstrt-end.hg19.bed | awk '{print $1"\t"$2"\t"$3"\t""0""\t""0""\t"$5}' | sort -k1,1V -k2,2n | uniq) -s | sort -k1,1V -k2,2n | uniq >inactive.genic.txstart-end.hg19.bed

# check there are no overlaps:
#bedtools intersect -a <(awk 'NR>1' activ.genic.txstrt-end.hg19.bed | cut -f 1-3,5 | awk '{print $1"\t"$2"\t"$3"\t""0""\t""0""\t"$4}' | sort -k1,1V -k2,2n | uniq) -b inactive.genic.txstart-end.hg19.bed -s | head  # empty

# ========================
#  TSS
# =========================
# active genes TSS: ±1.5kb from txStart on each strand, then bedt merge FORCE STRANDEDNESS!!! (need 6 cols in the bed file to force strandedness with -s...fill in 2 cols with just 0's)
awk 'NR>1' activ.genic.txstrt-end.hg19.bed | awk '{if ($5=="+") print $1"\t"$2-1500"\t"$2+1500"\t""0""\t""0""\t"$5; else print $1"\t"$3-1500"\t"$3+1500"\t""0""\t""0""\t"$5}' | sort -k1,1V -k2,2n | uniq | awk '{if ($1!~/_/) print $0}' | awk '{if ($1!~/M/) print $0}' | awk '{if ($1!~/Y/) print $0}' >3kb.around.uniqACTV-TSS.hg19.bed6

sort -k1,1V -k2,2n 3kb.around.uniqACTV-TSS.hg19.bed6 | bedtools merge -i stdin -s -o distinct -c 6 | sort -k1,1V -k2,2n | uniq >3kb.around.ACTV-TSS.hg19.mrg.bed4
	# 33011 3kb.around.uniqACTV-TSS.hg19.bed6
	# 25942 3kb.around.ACTV-TSS.hg19.mrg.bed4
	
# inactive genes TSS: ±1.5kb from txStart on each strand, then bedt merge FORCE STRANDEDNESS!!! (need 6 cols in the bed file to force strandedness with -s...fill in 2 cols with just 0's)
cat inactive.genic.txstart-end.hg19.bed | awk '{if ($6=="+") print $1"\t"$2-1500"\t"$2+1500"\t""0""\t""0""\t"$6; else print $1"\t"$3-1500"\t"$3+1500"\t""0""\t""0""\t"$6}' | sort -k1,1V -k2,2n | uniq | awk '{if ($1!~/_/) print $0}' | awk '{if ($1!~/M/) print $0}' | awk '{if ($1!~/Y/) print $0}' >3kb.around.uniqINACTV-TSS.hg19.bed6

sort -k1,1V -k2,2n 3kb.around.uniqINACTV-TSS.hg19.bed6 | bedtools merge -i stdin -s -o distinct -c 6 | sort -k1,1V -k2,2n | uniq >3kb.around.INACTV-TSS.hg19.mrg.bed4

	# 13872 3kb.around.uniqINACTV-TSS.hg19.bed6
	# 10693 3kb.around.INACTV-TSS.hg19.mrg.bed4

# check for overlaps bn INACTV-TSS and ACTV-TSS: there could be some if the backing up from the tx Start runs into another gene TSS region...
		#bedtools intersect -a 3kb.around.uniqACTV-TSS.hg19.bed6 -b 3kb.around.uniqINACTV-TSS.hg19.bed6 -s | wc -l
		# 349
		#bedtools intersect -a <(awk '{print $1"\t"$2"\t"$3"\t""0""\t""0""\t"$4}' 3kb.around.ACTV-TSS.hg19.mrg.bed4) -b <(awk '{print $1"\t"$2"\t"$3"\t""0""\t""0""\t"$4}' 3kb.around.INACTV-TSS.hg19.mrg.bed4) -s | wc -l
		# 246
# There are some overlaps bn actv and inactv unmrgd TSS and mrgdTSS, prolly bc the process involves backing up from the tx start coord 1500kb...then it could run into another TSS if the genes are closely spaced.

# =======================
# genes, nonTSS:
# =======================
# cat active TSS±1.5kb and inactive TSS±1.5kb since they sometimes overlap, then subtract that from the active and inactive genes:

cat 3kb.around.uniqACTV-TSS.hg19.bed6 3kb.around.uniqINACTV-TSS.hg19.bed6 | sort -k1,1V -k2,2n >3kb.around.allTSS.UNmrgd.bed6
		#46883 3kb.around.allTSS.UNmrgd.bed6
# active genes minus 3kbAroundAllTSS:
	# bedtools subtract -a <(awk 'NR>1' activ.genic.txstrt-end.hg19.bed | awk '{print $1"\t"$2"\t"$3"\t""0""\t""0""\t"$5}' | sort -k1,1V -k2,2n | uniq) -b 3kb.around.allTSS.UNmrgd.bed6 -s |  sort -k1,1V -k2,2n | uniq | wc -l
# 36956
	# bedtools subtract -a <(awk 'NR>1' activ.genic.txstrt-end.hg19.bed | awk '{print $1"\t"$2"\t"$3"\t""0""\t""0""\t"$5}' | sort -k1,1V -k2,2n | uniq) -b 3kb.around.allTSS.UNmrgd.bed6 -s |  sort -k1,1V -k2,2n | uniq | bedtools merge -i stdin -s -c 6 -o distinct | wc -l
# 25925


	# bedtools subtract -a <(sort -k1,1V -k2,2n inactive.genic.txstart-end.hg19.bed | uniq) -b 3kb.around.allTSS.UNmrgd.bed6 -s |  sort -k1,1V -k2,2n | uniq | wc -l
	# 9693
	# bedtools subtract -a <(sort -k1,1V -k2,2n inactive.genic.txstart-end.hg19.bed | uniq) -b 3kb.around.allTSS.UNmrgd.bed6 -s |  sort -k1,1V -k2,2n | uniq | bedtools merge -i stdin -s -c 6 -o distinct | wc -l
	# 7742

		# might as well merge...just interested in all the genic regiosn with the TSS subtracted out...dont care about gene id's anymore?
		bedtools subtract -a <(awk 'NR>1' activ.genic.txstrt-end.hg19.bed | awk '{print $1"\t"$2"\t"$3"\t""0""\t""0""\t"$5}' | sort -k1,1V -k2,2n | uniq) -b 3kb.around.allTSS.UNmrgd.bed6 -s |  sort -k1,1V -k2,2n | uniq | bedtools merge -i stdin -s -c 6 -o distinct | awk '{print $1"\t"$2"\t"$3"\t""0""\t""0""\t"$4}' | sort -k1,1V -k2,2n | uniq >nonTSS.activegenes.mrgd.bed6

	bedtools subtract -a <(sort -k1,1V -k2,2n inactive.genic.txstart-end.hg19.bed | uniq) -b 3kb.around.allTSS.UNmrgd.bed6 -s |  sort -k1,1V -k2,2n | uniq | bedtools merge -i stdin -s -c 6 -o distinct | awk '{print $1"\t"$2"\t"$3"\t""0""\t""0""\t"$4}' | sort -k1,1V -k2,2n | uniq >nonTSS.inactivegenes.mrgd.bed6


	# check no overlaps bn active genebodies and inactive genebodies
	#	bedtools intersect -a nonTSS.activegenes.mrgd.bed6 -b nonTSS.inactivegenes.mrgd.bed6 -s | head
	# nothing!
	# also the reciprocal yields nothing:
	#bedtools intersect -b nonTSS.activegenes.mrgd.bed6 -a nonTSS.inactivegenes.mrgd.bed6 -s | head

		# When don't force strandness there IS Overlap. when do, there is no overlap bn active and inactive genes!
		#chr1	1657291	1675938	0	0	-
		#grep '1657291' nonTSS.activegenes.mrgd.bed6
			# chr1	1657291	1675938	0	0	-
		#nano A.bed
			#chr1	1657291	1675938	0	0	-
		#bedtools intersect -a nonTSS.inactivegenes.mrgd.bed6 -b A.bed -wa
# chr1	1657553	1663343	0	0	+


# check that the merged active / inactive gene bodies dont overlap with any TSS;
	# bedtools intersect -a nonTSS.activegenes.mrgd.bed6 -b 3kb.around.allTSS.UNmrgd.bed6 -s | head   # nothing!
	# bedtools intersect -a nonTSS.inactivegenes.mrgd.bed6 -b 3kb.around.allTSS.UNmrgd.bed6 -s | head   # nothing!

# last step derive intergenic regions: anything in the whole genome that doesnt overlap with TSS and genes.
# don't care about strand info? since genes are ds, just anything outside of the genic region.
	# gene bodies files have weird chrs like chUn_* ...
#bedtools complement \
#-i <(cat 3kb.around.allTSS.UNmrgd.bed6 nonTSS.activegenes.mrgd.bed6 nonTSS.inactivegenes.mrgd.bed6 | sort -k1,1V -k2,2n | awk '{if ($1!~/chrM/) print $0}' | awk '{if ($1!~/_/) print $0}' | awk '{if ($1!~/chrY/) print $0}') \
#-g <(sort -k1,1V -k2,2n /data/myanco/accessory_files/hg19.chrom_chr.sizes) | awk '{if ($1!~/chrM/) print $0}' | wc -l # 18432 lines

#bedtools complement \
#-i <(cat 3kb.around.allTSS.UNmrgd.bed6 nonTSS.activegenes.mrgd.bed6 nonTSS.inactivegenes.mrgd.bed6 | sort -k1,1V -k2,2n | awk '{if ($1!~/chrM/) print $0}' | awk '{if ($1!~/_/) print $0}' | awk '{if ($1!~/chrY/) print $0}') \
#-g <(sort -k1,1V -k2,2n /data/myanco/accessory_files/hg19.chrom_chr.sizes) | awk '{if ($1!~/chrM/) print $0}' | sort -k1,1V -k2,2n | bedtools merge -i stdin | wc -l # still 18432 lines

bedtools complement \
-i <(cat 3kb.around.allTSS.UNmrgd.bed6 nonTSS.activegenes.mrgd.bed6 nonTSS.inactivegenes.mrgd.bed6 | sort -k1,1V -k2,2n | awk '{if ($1!~/chrM/) print $0}' | awk '{if ($1!~/_/) print $0}' | awk '{if ($1!~/chrY/) print $0}') \
-g <(sort -k1,1V -k2,2n /data/myanco/accessory_files/hg19.chrom_chr.sizes) | awk '{if ($1!~/chrM/) print $0}' | sort -k1,1V -k2,2n >igenic.bedtcompl.hg19.bed3


# check that all regions add up to the whole genome:
#cat 3kb.around.allTSS.UNmrgd.bed6 nonTSS.activegenes.mrgd.bed6 nonTSS.inactivegenes.mrgd.bed6 igenic.bedtcompl.hg19.bed3 | sort -k1,1V -k2,2n | wc -l  # 98982
#cat 3kb.around.allTSS.UNmrgd.bed6 nonTSS.activegenes.mrgd.bed6 nonTSS.inactivegenes.mrgd.bed6 igenic.bedtcompl.hg19.bed3 | sort -k1,1V -k2,2n | cut -f 1-3 | uniq | wc -l # 98979

#cat 3kb.around.allTSS.UNmrgd.bed6 nonTSS.activegenes.mrgd.bed6 nonTSS.inactivegenes.mrgd.bed6 igenic.bedtcompl.hg19.bed3 | sort -k1,1V -k2,2n | cut -f 1-3 | uniq | bedtools merge -i stdin | wc -l
# 1036
#cat 3kb.around.allTSS.UNmrgd.bed6 nonTSS.activegenes.mrgd.bed6 nonTSS.inactivegenes.mrgd.bed6 igenic.bedtcompl.hg19.bed3 | sort -k1,1V -k2,2n | cut -f 1-3 | uniq | bedtools merge -i stdin | awk '{if ($1!~/chrM/) print $0}' | awk '{if ($1!~/_/) print $0}' | awk '{if ($1!~/chrY/) print $0}' | wc -l
# 23; full chr's reconstituted

#cat 3kb.around.allTSS.UNmrgd.bed6 nonTSS.activegenes.mrgd.bed6 nonTSS.inactivegenes.mrgd.bed6 igenic.bedtcompl.hg19.bed3 | sort -k1,1V -k2,2n | cut -f 1-3 | uniq | bedtools merge -i stdin | awk '{if ($1!~/chrM/) print $0}' | awk '{if ($1!~/_/) print $0}' | awk '{if ($1!~/chrY/) print $0}' >test.chrs.bd

#diff -s <(sort -k1,1V -k2,2n test.chrs.bd) <(awk '{print "chr"$1"\t"$2"\t"$3}' /data/myanco/accessory_files/hg19.chrom.sizes.bed | awk '{if ($1!~/chrM/) print $0}' | sort -k1,1V -k2,2n)
# Files /dev/fd/63 and /dev/fd/62 are identical

#  eliminate chrY, chrM, chrUN from all the files as they may cause problems later on?
#awk '{if ($1!~/chrM/) print $0}' igenic.bedtcompl.hg19.bed3 | awk '{if ($1!~/_/) print $0}' | awk '{if ($1!~/chrY/) print $0}' | sort -k1,1V -k2,2n | uniq | bedtools merge -i stdin | wc -l  # cannot be merged.
	# igenic.bedtcompl.hg19.bed3 # is fine, all weird chr's already excluded.

# filter out weird chr's from:
#activ.genic.txstrt-end.hg19.bed
#inactive.genic.txstart-end.hg19.bed
#nonTSS.activegenes.mrgd.bed6
#nonTSS.inactivegenes.mrgd.bed6
#3kb.around.uniqACTV-TSS.hg19.bed6
#3kb.around.ACTV-TSS.hg19.mrg.bed4
#3kb.around.uniqINACTV-TSS.hg19.bed6
#3kb.around.INACTV-TSS.hg19.mrg.bed4

cp igenic.bedtcompl.hg19.bed3 /data/myanco/accessory_files/regions/
for f in *filt*; do cp $f /data/myanco/accessory_files/regions/$f; done


# length breakdown;
 awk '{print $3-$2}' igenic.bedtcompl.hg19.bed3 |  awk '{total += $1; count++ } END {print total/count}'
# 90718
 awk '{print $3-$2}' igenic.bedtcompl.hg19.bed3 |  awk '{total += $1; count++ } END {print total/1e6}'
 # 1672.11 Mb


# gene bodies:
awk '{print $3-$2}' nonTSS.activegenes.mrgd.filt.bed6 | awk '{total += $1; count++ } END {print total/1e6}'
# 1080.29
awk '{print $3-$2}' nonTSS.inactivegenes.mrgd.filt.bed6 | awk '{total += $1; count++ } END {print total/1e6}'
# 215.609

awk '{print $3-$2}' 3kb.around.ACTV-TSS.hg19.mrg.filt.bed4 |  awk '{total += $1; count++ } END {print total/1e6}'
# 83.4527
awk '{print $3-$2}' 3kb.around.INACTV-TSS.hg19.mrg.filt.bed4 |  awk '{total += $1; count++ } END {print total/1e6}'
# 34.3056

=1672.11+1080.29+215.609+83.4527+34.3056 =3085.7673 Mb   # haploid genome: 3.1e6 bp

# but... intergenic is both strands, other features are stranded...!
 
 # igenic regions 54%, active gene bodies 35%, inactive gene bodies 7%, active TSS±1.5kb 3%, inactive TSS±1.5kb 1%




# ==================================
# ==================================
# redo:
# ==================================
# ==================================

gunzip -cd knowngene.hg38.GENCODEv39.gz | awk 'NR>1' | awk '{print $2"\t"$4"\t"$5"\t"$3"\t"$2":"$4"-"$5}' | sort -k1,1V -k2,2n | wc -l 
# 266064
gunzip -cd knowngene.hg38.GENCODEv39.gz | awk 'NR>1' | awk '{print $2"\t"$4"\t"$5"\t"$3"\t"$2":"$4"-"$5}' | sort -k1,1V -k2,2n | bedtools merge -i stdin | wc -l
# 37630


gunzip -cd knowngene.hg38.GENCODEv39.gz | awk 'NR>1' | awk '{print $2"\t"$4"\t"$5"\t"$3"\t"$2":"$4"-"$5}' | sort -k1,1V -k2,2n | awk '{print $3-$2}' | awk '{total+=$1}END{print total/1e6}'
# 10127.2
gunzip -cd knowngene.hg38.GENCODEv39.gz | awk 'NR>1' | awk '{print $2"\t"$4"\t"$5"\t"$3"\t"$2":"$4"-"$5}' | sort -k1,1V -k2,2n | bedtools merge -i stdin | awk '{print $3-$2}' | awk '{total+=$1}END{print total/1e6}'
#1839.23


gunzip -cd knowngene.hg38.GENCODEv39.gz | awk 'NR>1' | awk '{print $2"\t"$4"\t"$5"\t"$3"\t"$2":"$4"-"$5}' | sort -k1,1V -k2,2n | bedtools merge -i stdin | grep -v 'chrM\|chrY\|chrUn\|alt\|fix\|random' >knowngene.hg38.GENCODv39.canon.mrg.bed

# 32545 knowngene.hg38.GENCODv39.canon.mrg.bed

awk '{print $3-$2}' knowngene.hg38.GENCODv39.canon.mrg.bed | awk '{total +=$1}END{print total/1e6}'
# 1759.2 Mb of the genome is taken up by genes without distinguishing strand, and with merging overlapping bins.

# which of these bins are NOT overlapped by the list of u2os-active genes: use bedt subtract.
gunzip -cd u2os-activgenes.hg38.bed.gz | awk 'NR>1' | awk '{print $2"\t"$4"\t"$5"\t"$3"\t"$2":"$4"-"$5}' | sort -k1,1V -k2,2n | bedtools merge -i stdin | grep -v 'chrM\|chrY\|chrUn\|alt\|fix\|random' >u2os.actvGenes.hg38.canon.mrg.bed


awk '{print $3-$2}' u2os.actvGenes.hg38.canon.mrg.bed | awk '{total +=$1}END{print total/1e6}'
# 1210.03 Mb covered by genes active in u2os.

bedtools subtract -a knowngene.hg38.GENCODv39.canon.mrg.bed -b u2os.actvGenes.hg38.canon.mrg.bed | sort -k1,1V -k2,2n | awk '{print $3-$2}' | awk '{total+=$1}END{print total/1e6}'
# 549.171 # genes inactive in u2os; cant be merged.

#  549.171+1210.03   =1759.201


bedtools subtract -a knowngene.hg38.GENCODv39.canon.mrg.bed -b u2os.actvGenes.hg38.canon.mrg.bed | sort -k1,1V -k2,2n >u2os.inactvGenes.hg38.canon.mrg.bed

# now find TSS: but this is dep on strand info....!



>knowngenes.hg38.GENCODv39.bed




# 4/26/22



/data/myanco/accessory_files/hg38/

# TSS: ±1`.5kb from the tx start site
# how to know whether to start from the txstart or txend coord? --> need strand info!!
# put the u2os coords of u2os activ genes back into bedt overlap w knowngenes and retain strand info?

#	-s	Require same strandedness.  That is, only report hits in B
#		that overlap A on the _same_ strand.
#		- By default, overlaps are reported without respect to strand.

#	-S	Require different strandedness.  That is, only report hits in B
#		that overlap A on the _opposite_ strand.
#		- By default, overlaps are reported without respect to strand.



#gunzip -cd u2os-activgenes.hg38.bed.gz | head
# #name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	proteinID	alignID
#ENST00000456328.2	chr1	+	11868	14409	11868	11868	3	11868,12612,13220,	12227,12721,14409,		uc286dmu.1

# find regions that do not overlap with hg38:

#gunzip -cd  knowngene.hg38.GENCODEv39.gz | head
#name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	proteinID	alignID
#ENST00000456328.2	chr1	+	11868	14409	11868	11868	3	11868,12612,13220,	#12227,12721,14409,		uc286dmu.1

#gunzip -cd knowngene.hg38.GENCODEv39.gz | cut -f 2,4-7 | head
#chrom	txStart	txEnd	cdsStart	cdsEnd
#chr1	11868	14409	11868	11868

#gunzip -cd knowngene.hg38.GENCODEv39.gz | awk '{print $2"\t"$4"\t"$5"\t"$2":"$6"-"$7}' | awk 'NR>1' | head
# chr1	11868	14409	chr1:11868-11868
# chr1	12009	13670	chr1:12009-12009

# print the full transcript coordinates of all known genes --> bedt intersect w the full txt coords of all u2os-active genes --> -v which ones dont overlap = transcript coordinates of genes NOT active in u2os.
#bedtools intersect -a <(gunzip -cd knowngene.hg38.GENCODEv39.gz | awk '{print $2"\t"$4"\t"$5"\t"$2":"$6"-"$7}' | awk 'NR>1' | sort -k1,1V -k2,2n) -b <(gunzip -cd u2os-activgenes.hg38.bed.gz | awk '{print $2"\t"$4"\t"$5"\t"$2":"$6"-"$7}' | awk 'NR>1' | sort -k1,1V -k2,2n) -v
# 68k lines
# if merge the ouptput it becomes 26k lines

#gunzip -cd knowngene.hg38.GENCODEv39.gz | wc -l
# 266k lines, when merged by txStart-End it is 37k lines.
#gunzip -cd u2os-activgenes.hg38.bed.gz | wc -l
# 166k lines, when merged by txStart-End it is 15k lines


# should the final bedt int -v output be bedt merged? diff lines correspond to diff possibl overlapping transcripts...?

#with bedt merging: total 37k , 15k u2os active genes, 26k regions that dont overlap.

# bedtools merge retain strand info from merged bins but don't use the strand info in merging: -o collapse to print strand info from merged bins..
#	-o	Specify the operation that should be applied to -c.
#		    collapse (i.e., print a delimited list (duplicates allowed)),
		    

#cat tmp1
#chr1	0	1000	+
#chr1	999	2000	-
#bedtools merge -i tmp1 -c 4 -o collapse
#chr1	0	2000	+,-




#bedtools intersect -a <(gunzip -cd knowngene.hg38.GENCODEv39.gz | awk '{print $2"\t"$4"\t"$5"\t"$3"\t"$2":"$6"-"$7}' | awk 'NR>1' | sort -k1,1V -k2,2n) -b <(gunzip -cd u2os-activgenes.hg38.bed.gz | awk '{print $2"\t"$4"\t"$5"\t"$3"\t"$2":"$6"-"$7}' | awk 'NR>1' | sort -k1,1V -k2,2n) -v
# retain strand info for the sake of getting the TSS?
# but don't factor it into the bedt overlap? if a u2os gene doesnt overlap the known genes regardless of strand..


# split up the known genome-wide genes into u2os active vs inactive.

	# done: find active genes u2os from Bostrom et al single-cell rna paper. they supplied ist of gene name abbrevs which i searched against the GenCode genome-wide gene list and pulled out the ones that overlapped.

	# step 1 find inactive genes transcript start thru end; retain strand info. --> bedtools intersect -v gencode genome-wide genes vs u2os-active genes. do NOT merge bins.
# step 2 split up active and inactive genes into TSS vs gene bodies.
			# for TSS: strand info needed. if transcript is on + strand, take ±1.5kb from the transcript START coord. if transcript is on - strand, take ±1.5kb from the transcript STOP coord. --> merge.
			# for gene bodies: strand info not needed?? subtract out the 3kb TSS from the gene bodies with bedtools subtract. --> merge.

# intergenic regions: bedtools complement hg38 vs list of all GenCode genes.
		# intergenic regions: tads vs itads. strand info not needed.


# check that all pairwise bedt int across genomic regions yields no overlaps,
# check that the full genome size is reconstituted and no one is mergable when cat all the genomic regions back together.
#			GCA_hg38.chrom_chr.nochrY.canon.bed                # use for bedt complement to get intergenic regions. 23 chr's, 0 to chr size.
      # total size of genome: 
#			awk '{total +=$3} END {print total/1e9}' GCA_hg38.chrom_chr.nochrY.canon.bed
# 3.03104 Gb
# ========================================================================================

knowngene.hg38.GENCODEv39.gz # genome-wide genes
u2os-activgenes.hg38.bed.gz  # u2os active genes
u2os.inactvGenes.hg38.canon.mrg.bed # u2os inactive genes but merged
 

#gunzip -cd knowngene.hg38.GENCODEv39.gz | cut -f 2,4-5 | wc -l
# 266065 lines unmerged, 37630 lines merged.
#gunzip -cd u2os-activgenes.hg38.bed.gz | cut -f 2,4-5 | wc -l
# 165772 lines unmerged, 15350 lines merged.

#bedtools intersect -a <(gunzip -cd knowngene.hg38.GENCODEv39.gz | awk '{print $2"\t"$4"\t"$5"\t"$3"\t"$2":"$6"-"$7}' | awk 'NR>1' | sort -k1,1V -k2,2n) -b <(gunzip -cd u2os-activgenes.hg38.bed.gz | awk '{print $2"\t"$4"\t"$5"\t"$3"\t"$2":"$6"-"$7}' | awk 'NR>1' | sort -k1,1V -k2,2n) -v
# 68224 lines unmerged, 26064 lines merged.

#165772+68224: active + inactive both unmerged adds up to less than the total unmerged lines of all genes. = 233996  

#37630+15350 = # 52980: merged active + inactive also not adding up to total merged gencode genes.

#bedtools subtract -a <(gunzip -cd knowngene.hg38.GENCODEv39.gz | awk '{print $2"\t"$4"\t"$5"\t"$3"\t"$2":"$6"-"$7}' | awk 'NR>1' | sort -k1,1V -k2,2n) -b <(gunzip -cd u2os-activgenes.hg38.bed.gz | awk '{print $2"\t"$4"\t"$5"\t"$3"\t"$2":"$6"-"$7}' | awk 'NR>1' | sort -k1,1V -k2,2n) | wc -l
# 80715
# this one prints all the INTERVALS of A that don't overlap; bedt intersect -v prints just the full a intervals that have no overlaps with b.

# do the bedtools overlap on terminal instead of on ucsc table browser?

# 18191 /data/myanco/accessory_files/hg38/genesU2OS/geneNames.PosReadCount.u2osG1_1.uniq.bed
# 18k gene id's, yet 165k lines of transcripts --> due to alternate promoters?

#grep -f /data/myanco/accessory_files/hg38/genesU2OS/geneNames.PosReadCount.u2osG1_1.uniq.bed <(gunzip -cd knowngene.hg38.GENCODEv39.gz | awk '{print $2"\t"$4"\t"$5"\t"$3"\t"$2":"$6"-"$7}' | awk 'NR>1' | sort -k1,1V -k2,2n)


# or go back to bedtools ucsc and ask for lines that do NOT overlap with u2os active genes?

#gunzip -cd knowngene.hg38.GENCODEv39.gz  | cut -f 2 | uniq | wc -l
#440 chr's/scaffolds/etc

#gunzip -cd u2os-activgenes.hg38.bed.gz | cut -f 2 | sort -k1,1V | uniq | wc -l
# 278 --> noncanon stuff captured.

#gunzip -cd u2os-activgenes.hg38.bed.gz | cut -f 2 | awk 'NR>1' | sort -k1,1V | uniq | grep -v 'alt\|Un\|EBV\|chrM\|fix\|chrY\|random'  | wc -l
# 23 chrs

# filter out noncanonical stuff from initial list of gencode genes, and from u2os-active genes.
# 

gunzip -cd knowngene.hg38.GENCODEv39.gz | awk 'NR>1' | awk '{print $2"\t"$4"\t"$5"\t"$2":"$6"-"$7":"$3}' | sort -k1,1V -k2,2n | grep -v 'alt\|Un\|EBV\|chrM\|fix\|chrY\|random' >knowngene.hg38.GENCODEv39.canon.bed
#243892 lines canon known gencode genes, 32545 merged, 

#rm knowngene.hg38.GENCODv39.canon.mrg.bed


#gunzip -cd u2os-activgenes.hg38.bed.gz | awk 'NR>1' | awk '{print $2"\t"$4"\t"$5"\t"$2":"$6"-"$7":"$3}' | sort -k1,1V -k2,2n | grep -v 'alt\|Un\|EBV\|chrM\|fix\|chrY\|random'

gunzip -cd u2os-activgenes.hg38.bed.gz | awk 'NR>1' | awk '{print $2"\t"$4"\t"$5"\t"$2":"$6"-"$7":"$3}' | sort -k1,1V -k2,2n | grep -v 'alt\|Un\|EBV\|chrM\|fix\|chrY\|random' >u2os.actv.txts.hg38.canon.bed
# 153990 unmerged, 13881 merged



#bedtools intersect -a knowngene.hg38.GENCODEv39.canon.bed -b u2os.actv.txts.hg38.canon.bed -v  | wc -l
#59764

# 59764+153990 =213754   # active u2os unmrgd + inactive u2os unmrgd = 214k; but full known txts is 244k unmrgd.
# maybe cause full known txts has some redundant lines with same tx start and end?
#cut -f 1-3 knowngene.hg38.GENCODEv39.canon.bed | uniq | wc -l
#236372
#wc -l knowngene.hg38.GENCODEv39.canon.bed
#243892 knowngene.hg38.GENCODEv39.canon.bed

#bedtools intersect -a knowngene.hg38.GENCODEv39.canon.bed -b u2os.actv.txts.hg38.canon.bed -v  | cut -f 1-3 | uniq | wc -l
#57978
#cut -f 1-3 u2os.actv.txts.hg38.canon.bed | uniq | wc -l
#148755

#57978+148755 = 206733

# maybe not all the u2os active transcript lines overlap the gencode lines:
#bedtools intersect -b knowngene.hg38.GENCODEv39.canon.bed -a u2os.actv.txts.hg38.canon.bed  | cut -f 1-3 | uniq | wc -l
# 1978313
#bedtools intersect -a u2os.actv.txts.hg38.canon.bed -b knowngene.hg38.GENCODEv39.canon.bed -wa  -v | wc -l
# 0 --> all u2os actv transcripts can be found in gencode bed file.


#bedtools subtract -a knowngene.hg38.GENCODEv39.canon.bed -b u2os.actv.txts.hg38.canon.bed | wc -l
# 71172 --> known genes get fragmented up.


# for u2os inactive transcripts use bedtools intersect -v so as to not fragment up genecode transcripts.


#awk '{print $3-$2}' u2os.actv.txts.hg38.canon.bed  | awk '{total +=$1}END{print total/1e6}'
# 7087.08  Mb
#awk '{print $3-$2}' u2os.actv.txts.hg38.canon.bed  | awk '{total +=$1; count++}END{print total/count}'
# 46023 bp avg size = 46 kb
# avg size of a human gene is about 27kb...
# longest transcript is TTN. 280kb.

# but here longest transcripts are all over 2Mb each!
#awk '{print $3-$2}' u2os.actv.txts.hg38.canon.bed | sort -k1,1n | tail
# 2,057,682
# see where they map:
#awk '{print $0"\t"$3-$2}' u2os.actv.txts.hg38.canon.bed | sort -k5,5n | tail

#chr20	13995515	16053197	chr20:13995763-16049876:+	2057682 # MACROD2, 2057kb.
#chr20	13995368	16053197	chr20:13995763-16049876:+	2057829

# double counting is at play here ^. select each TSS separately if they differ, but merge at the end! after having subtracted out TSS±1.5kb regions from the txn start (if + strand) or stop (if - strand) coordinate.

#chr8	2935360	4994914	chr8:2938584-4994416:-	2059554

#chr9	8528151	10612723	chr9:8528424-8733843:-	2084572

#chrX	31119221	33211549	chrX:31121918-33211312:-	2092328

#chr11	83457918	85627270	chr11:83459817-85626617:-	2169352
#chr11	83455172	85627344	chr11:83459817-85598696:-	2172172

#chr9	8314245	10613002	chr9:8317873-8733843:-	2298757

#chr7	146116800	148420998	chr7:146116876-148415616:+	2304198

#chr16	5239801	7711458	chr16:5239886-7709641:+	2471657



#chr16:5239801-7711458
# this gene, RbFox1, is near a telomere..


#chr7:146,116,800-148,420,998
#IGV:  LOC105375 --> a locus..
# UCSC: CNTNAP2 , ENSG00000274127

#gunzip -cd knowngene.hg38.GENCODEv39.gz | grep '148420998'
#ENST00000361727.8	chr7	+	146116800	148420998	146116876	148415616	24	146116800,146774270,146839710,147043906,147108146,147120978,147128692,147132244,147300140,147395608,147485934,147562137,147639105,147903564,147977861,148118117,148147490,148172241,148217287,148229645,148267032,148383648,148409390,148415416,	146116973,146774381,146839904,147044054,147108350,147121163,147128836,147132509,147300290,147395780,147486041,147562257,147639306,147903721,147977989,148118288,148147709,148172478,148217524,148229779,148267126,148383888,148409471,148420998,	Q9UHC6	uc003weu.4

#ENST00000463592.3	chr7	+	148339460	148420998	148383842	148415616	4	148339460,148383648,148409390,148415416,	148339633,148383888,148409471,148420998,	Q9UHC6	uc003wev.3

#grep 'Q9UHC6' genesU2OS/geneNames.PosReadCount.u2osG1_1.uniq.bed
#grep 'ENST00000463592' genesU2OS/geneNames.PosReadCount.u2osG1_1.uniq.bed
#grep 'uc003wev.3' genesU2OS/geneNames.PosReadCount.u2osG1_1.uniq.bed

#grep 'uc003weu' genesU2OS/geneNames.PosReadCount.u2osG1_1.uniq.bed

#grep 'Q9' genesU2OS/geneNames.PosReadCount.u2osG1_1.uniq.bed
# COQ9

#gunzip -cd knowngene.hg38.GENCODEv39.gz | grep 'COQ9'

# search COQ9 on UCSC, hg38:
#COQ9 (ENST00000262507.11) at chr16:57447479-57461270 - Homo sapiens coenzyme Q9 (COQ9), mRNA; nuclear gene for mitochondrial product. (from RefSeq NM_020312)
#COQ4 (ENST00000300452.8) at chr9:128322839-128334072 - Homo sapiens coenzyme Q4 (COQ4), transcript variant 1, mRNA; nuclear gene for mitochondrial product. (from RefSeq NM_016035)

# bedtools intersect -a chr16_57447479-57461270.bed -b u2os.actv.txts.hg38.canon.bed
 # lots of hits; multiple transcripts encompassed.

#gunzip -cd knowngene.hg38.GENCODEv39.gz | grep 'CNTNAP2'  # one of longest transcripts supposedly active in u2os. no hits.
#gunzip -cd knowngene.hg38.GENCODEv39.gz | grep 'ENSG00000274127'

#chr7.146116800-148420998.bed

#bedtools intersect -a chr7.146116800-148420998.bed -b knowngene.hg38.GENCODEv39.canon.bed
# lots of hits, mergable to just 1 bin.

#bedtools intersect -a chr7.146116800-148420998.bed -b u2os.actv.txts.hg38.canon.bed  | sort -k1,1V -k2,2n | bedtools merge -i stdin | wc -l
#1
# again # lots of hits, mergable to just 1 bin.

#bedtools intersect -a knowngene.hg38.GENCODEv39.canon.bed -b u2os.actv.txts.hg38.canon.bed -v | wc -l
# 59764

bedtools intersect -a knowngene.hg38.GENCODEv39.canon.bed -b u2os.actv.txts.hg38.canon.bed -v | sort -k1,1V -k2,2n >u2os.inactv.GNCD39txts.hg38.canon.bed

# split 4th col to retrieve strand info. if strand is + take the tx start coord. if strand is - take the tx stop coord.

u2os.actv.txts.hg38.canon.bed
u2os.inactv.GNCD39txts.hg38.canon.bed

# 
#awk '{split($4,a,/[:]/); print $1"\t"$2"\t"$3"\t"a[3]}' u2os.actv.txts.hg38.canon.bed | cut -f 4 | sort -k1,1 | uniq -c
#  75696 -
#  78294 +
  
#  awk '{split($4,a,/[:]/); print $1"\t"$2"\t"$3"\t"a[3]}' u2os.actv.txts.hg38.canon.bed | head
  
#   chr1	14403	29570	-   # on IGV this STOP coord, col3, does correspond to the start of a transcript on the minus strand.


awk '{split($4,a,/[:]/); print $1"\t"$2"\t"$3"\t"a[3]}' u2os.inactv.GNCD39txts.hg38.canon.bed | awk '{if ($4=="+") print $1"\t"($2-1500)"\t"($2+1500)"\t"$4; else print $1"\t"($3-1500)"\t"($3+1500)"\t"$4}' | sort -k1,1V -k2,2n >u2os.inactvTSS.hg38.canon.bed

awk '{split($4,a,/[:]/); print $1"\t"$2"\t"$3"\t"a[3]}' u2os.actv.txts.hg38.canon.bed | awk '{if ($4=="+") print $1"\t"($2-1500)"\t"($2+1500)"\t"$4; else print $1"\t"($3-1500)"\t"($3+1500)"\t"$4}' | sort -k1,1V -k2,2n >u2os.actvTSS.hg38.canon.bed


#cut -f 4 u2os.inactvTSS.hg38.canon.bed | sort -k1,1 | uniq -c
#  28200 -
#  31564 +
#cut -f 4 u2os.actvTSS.hg38.canon.bed | sort -k1,1 | uniq -c
#  75696 -
#  78294 +

# are any of the TSS mergeable within each set? are any of the active overlapping the inactive TSS?

# 153,990 u2os.actvTSS.hg38.canon.bed; 42751 mrgd.
# 59,764 u2os.inactvTSS.hg38.canon.bed; 29782 mrgd.

#bedtools intersect -a u2os.actvTSS.hg38.canon.bed -b u2os.inactvTSS.hg38.canon.bed | wc -l
# 28,812

# probably due to strand differences...and /or 3kb regions running into one another.
# which one to take? conflict...!
# override silent genes and give u2os active genes the priority?


# to get gene bodies: take each transcript, if its on the pos strand subtract out the ±1.5kb around the start coordinate, if on the negative strand subtract out the ±1.5kb around the end coordinate.
# OR add the TSS together, subtract out all the TSS from all the transcripts so as to never count things twice.

# do it on a strand basis?? match the strand --> or don't care...
# use bedt subtract, not bedt intersect,bc bedt subtract parses up the original coordinates into just the areas NOT overlapped by the b feature ie all TSS.
#bedtools subtract -a u2os.actv.txts.hg38.canon.bed -b <(cat u2os.inactvTSS.hg38.canon.bed u2os.actvTSS.hg38.canon.bed | sort -k1,1V -k2,2n ) | awk '{print $3-$2}' | awk '{total +=$1; count++} END {print total/1e6}'
#318,321 lines, 6000 Mb or 6 Gb; avg size 19 kb. --> prolly many double bins; print only unique lines:
#bedtools subtract -a u2os.actv.txts.hg38.canon.bed -b <(cat u2os.inactvTSS.hg38.canon.bed u2os.actvTSS.hg38.canon.bed | sort -k1,1V -k2,2n ) | sort -k1,1V -k2,2n | uniq | awk '{print $3-$2}' | awk '{total +=$1; count++} END {print total/1e6}'
# 245,778 lines; 4561 Mb or 4.6 Mb; avg size 18.6 kb

#bedtools subtract -a u2os.actv.txts.hg38.canon.bed -b <(cat u2os.inactvTSS.hg38.canon.bed u2os.actvTSS.hg38.canon.bed | sort -k1,1V -k2,2n ) | sort -k1,1V -k2,2n | uniq | bedtools merge -i stdin | awk '{print $3-$2}' | awk '{total +=$1; count++} END {print total/1e6}' 
# if merge it is 42k lines totaling 1070.21 Mb and avg size 25kkb.


#bedtools merge -i  u2os.actv.txts.hg38.canon.bed | wc -l
#13881, unmrgd 154k lines

#merge the transcripts file before bedt subtr?
#bedtools subtract -a <(sort -k1,1V -k2,2n u2os.actv.txts.hg38.canon.bed | bedtools merge -i stdin) -b <(cat u2os.inactvTSS.hg38.canon.bed u2os.actvTSS.hg38.canon.bed | sort -k1,1V -k2,2n ) | sort -k1,1V -k2,2n | uniq | bedtools merge -i stdin | awk '{print $3-$2}' | awk '{total +=$1; count++} END {print total/1e6}' 
# 1070.21 Mb and 25kb length, same as merging after bedt subtract.

# bedtools subtract out the tss again AFTER merging the tss-subtracted-transcripts: is it fragmented again?
#bedtools subtract -a <(sort -k1,1V -k2,2n u2os.actv.txts.hg38.canon.bed | bedtools merge -i stdin) -b <(cat u2os.inactvTSS.hg38.canon.bed u2os.actvTSS.hg38.canon.bed | sort -k1,1V -k2,2n ) | sort -k1,1V -k2,2n | uniq | bedtools merge -i stdin >tmp1
# 42333 tmp1
#bedtools subtract -a tmp1 -b <(cat u2os.inactvTSS.hg38.canon.bed u2os.actvTSS.hg38.canon.bed | sort -k1,1V -k2,2n ) | sort -k1,1V -k2,2n | wc -l
# 42333
#bedtools intersect -a tmp1 -b <(cat u2os.inactvTSS.hg38.canon.bed u2os.actvTSS.hg38.canon.bed | sort -k1,1V -k2,2n ) | wc -l
# 0

#bedtools subtract -a <(sort -k1,1V -k2,2n u2os.actv.txts.hg38.canon.bed) -b <(cat u2os.inactvTSS.hg38.canon.bed u2os.actvTSS.hg38.canon.bed | sort -k1,1V -k2,2n ) | sort -k1,1V -k2,2n | uniq >u2os.actv.txts.noTSS.hg38.canon.unmrg.bed
# 245778 u2os.actv.txts.noTSS.hg38.canon.unmrg.bed; 42333 merged.
#bedtools subtract -a <(sort -k1,1V -k2,2n u2os.inactv.GNCD39txts.hg38.canon.bed ) -b <(cat u2os.inactvTSS.hg38.canon.bed u2os.actvTSS.hg38.canon.bed | sort -k1,1V -k2,2n ) | sort -k1,1V -k2,2n | uniq >u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed
# 67908 u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed; 19639 merged.

#bedtools intersect -a u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed -b <(cat u2os.inactvTSS.hg38.canon.bed u2os.actvTSS.hg38.canon.bed | sort -k1,1V -k2,2n ) | wc -l
# 0

# are inactive and active gene bods overlap? --> no:
#bedtools intersect -a u2os.actv.txts.noTSS.hg38.canon.unmrg.bed -b u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed | wc -l
#0

# ############
# intergenic:
# #############
# bedt complement w whole genome, canon,no chrY.

# each individ gene body / tss file is mergebale. mrge each one and then add up: should not be mergeable anymore?
#cat <( bedtools merge -i u2os.actv.txts.noTSS.hg38.canon.unmrg.bed)  <(bedtools merge -i u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed) <(bedtools merge -i u2os.inactvTSS.hg38.canon.bed) <(bedtools merge -i u2os.actvTSS.hg38.canon.bed) | wc -l
#134505
#cat <( bedtools merge -i u2os.actv.txts.noTSS.hg38.canon.unmrg.bed)  <(bedtools merge -i u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed) <(bedtools merge -i u2os.inactvTSS.hg38.canon.bed) <(bedtools merge -i u2os.actvTSS.hg38.canon.bed) | sort -k1,1V -k2,2n | bedtools merge -i stdin | wc -l
#29925
# maybe fts are book-ended?
#cat <( bedtools merge -i u2os.actv.txts.noTSS.hg38.canon.unmrg.bed)  <(bedtools merge -i u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed) <(bedtools merge -i u2os.inactvTSS.hg38.canon.bed) <(bedtools merge -i u2os.actvTSS.hg38.canon.bed) | sort -k1,1V -k2,2n | awk '{print $0"\t"$1"\t"$3-1}' | bedtools merge -i stdin | wc -l
# fts are more than 1bp overlap.. again mergable to just #29925 bins.
# 
#pairwise intersect:
#bedtools intersect -a <(bedtools merge -i u2os.actv.txts.noTSS.hg38.canon.unmrg.bed) -b <(bedtools merge -i u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed) | wc -l
# 0
#bedtools intersect -a <(bedtools merge -i u2os.actv.txts.noTSS.hg38.canon.unmrg.bed) -b <(bedtools merge -i u2os.inactvTSS.hg38.canon.bed) | wc -l
# 0
#bedtools intersect -a <(bedtools merge -i u2os.actv.txts.noTSS.hg38.canon.unmrg.bed) -b <(bedtools merge -i u2os.actvTSS.hg38.canon.bed) | wc -l
# 0
#bedtools intersect -a <(bedtools merge -i u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed) -b <(bedtools merge -i u2os.inactvTSS.hg38.canon.bed) | wc -l
# 0
#bedtools intersect -a <(bedtools merge -i u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed) -b <(bedtools merge -i u2os.actvTSS.hg38.canon.bed) | wc -l
# 0
#bedtools intersect -a <(bedtools merge -i u2os.inactvTSS.hg38.canon.bed) -b <(bedtools merge -i u2os.actvTSS.hg38.canon.bed) | wc -l
# 2210 # inactive and active TSS overlap!  fix: if an active tss region overlaps with an inactive, make the active tss count.
# OR bedt subtract out the inactive tss 3kb regions from the active tss regions. that way cut short the 3kb region if it encounters an inactive pr. howeveer the machinery for the active one could still be around ... --> therefore active tss has preeminence.

#bedtools intersect -a u2os.inactvTSS.hg38.canon.bed -b u2os.actvTSS.hg38.canon.bed | wc -l
# 28,812 inactive tss overlpd by active tss.
# 59,764 u2os.inactvTSS.hg38.canon.bed
# 153,990 u2os.actvTSS.hg38.canon.bed

#bedtools intersect -a u2os.inactvTSS.hg38.canon.bed -b u2os.actvTSS.hg38.canon.bed -v 
# 54,383 of the 59.7k inactive tss are not overlapping any inactive tss.

#bedtools intersect -a u2os.inactvTSS.hg38.canon.bed -b u2os.actvTSS.hg38.canon.bed -v | sort -k1,1V -k2,2n >u2os.inactvTSS-noactvTSSovlps.hg38.canon.unmrg.bed

# redo the gene bodies: should yield same result?
#bedtools subtract -a <(sort -k1,1V -k2,2n u2os.actv.txts.hg38.canon.bed) -b <(cat u2os.inactvTSS.hg38.canon.bed u2os.actvTSS.hg38.canon.bed | sort -k1,1V -k2,2n ) | sort -k1,1V -k2,2n | uniq >u2os.actv.txts.noTSS.hg38.canon.unmrg.bed
# 245778 u2os.actv.txts.noTSS.hg38.canon.unmrg.bed; 42333 merged.
#bedtools subtract -a <(sort -k1,1V -k2,2n u2os.actv.txts.hg38.canon.bed) -b <(cat u2os.inactvTSS-noactvTSSovlps.hg38.canon.unmrg.bed u2os.actvTSS.hg38.canon.bed | sort -k1,1V -k2,2n ) | sort -k1,1V -k2,2n | uniq | wc -l
# 245895  # slightly more bins from subtracting out fewer things???
bedtools subtract -a <(sort -k1,1V -k2,2n u2os.actv.txts.hg38.canon.bed) -b <(cat u2os.inactvTSS-noactvTSSovlps.hg38.canon.unmrg.bed u2os.actvTSS.hg38.canon.bed | sort -k1,1V -k2,2n ) | sort -k1,1V -k2,2n | uniq >u2os.actv.txts.noTSS.hg38.canon.unmrg.bed

#bedtools subtract -a <(sort -k1,1V -k2,2n u2os.inactv.GNCD39txts.hg38.canon.bed ) -b <(cat u2os.inactvTSS.hg38.canon.bed u2os.actvTSS.hg38.canon.bed | sort -k1,1V -k2,2n ) | sort -k1,1V -k2,2n | uniq >u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed
# 67908 u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed; 19639 merged.
#bedtools subtract -a <(sort -k1,1V -k2,2n u2os.inactv.GNCD39txts.hg38.canon.bed ) -b <(cat u2os.inactvTSS-noactvTSSovlps.hg38.canon.unmrg.bed u2os.actvTSS.hg38.canon.bed | sort -k1,1V -k2,2n ) | sort -k1,1V -k2,2n | uniq | wc -l
# 69038
bedtools subtract -a <(sort -k1,1V -k2,2n u2os.inactv.GNCD39txts.hg38.canon.bed ) -b <(cat u2os.inactvTSS-noactvTSSovlps.hg38.canon.unmrg.bed u2os.actvTSS.hg38.canon.bed | sort -k1,1V -k2,2n ) | sort -k1,1V -k2,2n | uniq >u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed


# all pairw intersect 0 ovlp:
#bedtools intersect -a <(bedtools merge -i u2os.actv.txts.noTSS.hg38.canon.unmrg.bed) -b <(bedtools merge -i u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed) | wc -l

#bedtools intersect -a <(bedtools merge -i u2os.actv.txts.noTSS.hg38.canon.unmrg.bed) -b <(bedtools merge -i u2os.actvTSS.hg38.canon.bed) | wc -l

#bedtools intersect -a <(bedtools merge -i u2os.actv.txts.noTSS.hg38.canon.unmrg.bed) -b <(bedtools merge -i u2os.inactvTSS-noactvTSSovlps.hg38.canon.unmrg.bed) | wc -l

#bedtools intersect -a <(bedtools merge -i u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed) -b <(bedtools merge -i u2os.actvTSS.hg38.canon.bed) | wc -l

#bedtools intersect -a <(bedtools merge -i u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed) -b <(bedtools merge -i u2os.inactvTSS-noactvTSSovlps.hg38.canon.unmrg.bed) | wc -l

#bedtools intersect -a <(bedtools merge -i u2os.actvTSS.hg38.canon.bed) -b <(bedtools merge -i u2os.inactvTSS-noactvTSSovlps.hg38.canon.unmrg.bed) | wc -l

#cat <(bedtools merge -i u2os.actvTSS.hg38.canon.bed) <( bedtools merge -i  u2os.inactvTSS-noactvTSSovlps.hg38.canon.unmrg.bed ) <( bedtools merge -i u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed) <(bedtools merge -i  u2os.actv.txts.noTSS.hg38.canon.unmrg.bed) | sort -k1,1V -k2,2n | uniq | wc -l
#132960 bins; 1736.25 Mb; 13.058 kb  avg size
#cat <(bedtools merge -i u2os.actvTSS.hg38.canon.bed) <( bedtools merge -i  u2os.inactvTSS-noactvTSSovlps.hg38.canon.unmrg.bed ) <( bedtools merge -i u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed) <(bedtools merge -i  u2os.actv.txts.noTSS.hg38.canon.unmrg.bed) | sort -k1,1V -k2,2n | uniq | bedtools merge -i stdin | wc -l
# 30393; 1736.25 Mb; 57.126.8 kb avg size 

# --> total area covered is still the same when merge the bins! therefore nothing is redundant?? yet book-ended undo is still mergeable..
#cat <(bedtools merge -i u2os.actvTSS.hg38.canon.bed) <( bedtools merge -i  u2os.inactvTSS-noactvTSSovlps.hg38.canon.unmrg.bed ) <( bedtools merge -i u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed) <(bedtools merge -i  u2os.actv.txts.noTSS.hg38.canon.unmrg.bed) | sort -k1,1V -k2,2n | bedtools merge -i stdin -d -1 | wc -l
# 132960 bins. thus it IS book-ended!
# and bedt complement of the genome + all these regions merged or unmerged after cat'ing together should yield same. if all mergebale feats are book-ended, then they should not have any extra igenic gaps in between when merging or not merging the genic feats.

# 1736 Mb is 58% of 3Gb whole genome.

#bedtools complement -i <(cat u2os.actv.txts.noTSS.hg38.canon.unmrg.bed  u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed u2os.inactvTSS.hg38.canon.bed u2os.actvTSS.hg38.canon.bed | sort -k1,1V -k2,2n) -g GCA_hg38.chrom_chr.nochrY.canon.sizes | wc -l
# 29948
# bedtools complement -i <(cat u2os.actv.txts.noTSS.hg38.canon.unmrg.bed  u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed u2os.inactvTSS.hg38.canon.bed u2os.actvTSS.hg38.canon.bed | sort -k1,1V -k2,2n | bedtools merge -i stdin ) -g GCA_hg38.chrom_chr.nochrY.canon.sizes | wc -l
# 29948; not further mergeable.

bedtools complement -i <(cat u2os.actv.txts.noTSS.hg38.canon.unmrg.bed  u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed u2os.inactvTSS-noactvTSSovlps.hg38.canon.unmrg.bed u2os.actvTSS.hg38.canon.bed | sort -k1,1V -k2,2n | bedtools merge -i stdin ) -g GCA_hg38.chrom_chr.nochrY.canon.sizes | sort -k1,1V -k2,2n >igenic.noU2OSgenic.GCA.canon.nochrY.hg38.bed
# 30416 igenic.noU2OSgenic.GCA.canon.nochrY.hg38.bed, not mergable
# 1294.79 Mb total, avg size 42.5693 kbp 

# 1294.79 Mb igenic + 1737.07 Mb genic (cat u2os.actv.txts.noTSS.hg38.canon.unmrg.bed  u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed u2os.inactvTSS.hg38.canon.bed u2os.actvTSS.hg38.canon.bed | sort -k1,1V -k2,2n | bedtools merge -i stdin| awk '{print $3-$2}' | awk '{total +=$1} END {print total/1e6}' )
# 1737.07+1294.79 =  3031.8599999999997
# awk '{total+=$2} END {print total/1e6}' GCA_hg38.chrom_chr.nochrY.canon.sizes
# 3031.04 Mb total genome.
 
 
 cat igenic.noU2OSgenic.GCA.canon.nochrY.hg38.bed u2os.actv.txts.noTSS.hg38.canon.unmrg.bed u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed u2os.inactvTSS-noactvTSSovlps.hg38.canon.unmrg.bed u2os.actvTSS.hg38.canon.bed | cut -f 1-3 | sort -k1,1V -k2,2n | bedtools merge -i stdin
# 23 bins 1 per chromosome!!
cat igenic.noU2OSgenic.GCA.canon.nochrY.hg38.bed u2os.actv.txts.noTSS.hg38.canon.unmrg.bed u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed u2os.inactvTSS-noactvTSSovlps.hg38.canon.unmrg.bed u2os.actvTSS.hg38.canon.bed | cut -f 1-3 | sort -k1,1V -k2,2n | bedtools merge -i stdin | awk '{total +=$3} END {print total/1e6}'
# 3031.04 Mb



# divide igenic:
/data/myanco/accessory_files/HiC_U2OS_hg38/360min_rep12cat.TADs.mrg.hg38.bed
/data/myanco/accessory_files/HiC_U2OS_hg38/iTADcatrep12.bedtcompl.hg38-nochrY-canon.bed

# awk '{total +=($3-$2)} END {print total/1e6}' /data/myanco/accessory_files/HiC_U2OS_hg38/360min_rep12cat.TADs.mrg.hg38.bed
# 2612.52 Mb, 1000 TADs
#awk '{total +=($3-$2)} END {print total/1e6}' /data/myanco/accessory_files/HiC_U2OS_hg38/iTADcatrep12.bedtcompl.hg38-nochrY-canon.bed
# 418.522 Mb, 1021 iTADs
#2612.52+418.522 = 3031.042

# majority of itads should overlap w igenic not genic feats? --> yes, 245mb of 418mb or 51% of itads overlap w igenic regions. by contrast 1050.13/2612.52 or 40% of tads overlap igenic.
# fragment the igenic regions with bedt subtract? or just get the igenic bins which overlap majority by tads or itads?
#--> fragment, more accurate.
#  
bedtools subtract -a igenic.noU2OSgenic.GCA.canon.nochrY.hg38.bed -b /data/myanco/accessory_files/HiC_U2OS_hg38/360min_rep12cat.TADs.mrg.hg38.bed | sort -k1,1V -k2,2n >igenic.noU2OSgenic.iTADovlp.bed
# 244.654 Mb, avg 51987.7 bp, 4706 bins not mergeable --> 245 Mb of the total 418 Mb of itads overlap w igenic regions.
bedtools subtract -a igenic.noU2OSgenic.GCA.canon.nochrY.hg38.bed -b /data/myanco/accessory_files/HiC_U2OS_hg38/iTADcatrep12.bedtcompl.hg38-nochrY-canon.bed | sort -k1,1V -k2,2n >igenic.noU2OSgenic.TADovlp.bed
# 1050.13 Mb, avg 39923 bp, 26,304 bins not mergeable. --> 1050 mb of the 2612 Mb of tads overlap w igenic regions.

#1049.47+244.503 =1293.973
#awk '{print $3-$2}' igenic.noU2OSgenic.GCA.canon.nochrY.hg38.bed | awk '{total +=$1} END {print total/1e6}'
#1293.97.  igenic hetero  + euchr cover all igenic regions..



# heterochromatin igenic: less fragmented up (longer avg size), but less proportion of igenic.
# 

# check that the regions are mutually exclusive:
# book-ended igenic itad vs tads. -d -1 bedt merge = not mergeable.


# after divide igenic into tad vs itad:
# generate a bed file covering all feats, visualize on igv hg38 to see that everyone is mutually exclusive and that the full genome is covered.
# all bins together should be mergeable into 23 bins 1 per chr = whole genome.
#cat <(bedtools merge -i u2os.actvTSS.hg38.canon.bed) <(bedtools merge -i u2os.inactvTSS-noactvTSSovlps.hg38.canon.unmrg.bed) <(bedtools merge -i u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed) <(bedtools merge -i u2os.actv.txts.noTSS.hg38.canon.unmrg.bed) igenic.noU2OSgenic.iTADovlp.bed igenic.noU2OSgenic.TADovlp.bed | sort -k1,1V -k2,2n | bedtools merge -i stdin | awk '{print $3-$2}' | awk '{total +=$1} END {print total/1e6}'
#3031.04
#cat u2os.actvTSS.hg38.canon.bed u2os.inactvTSS-noactvTSSovlps.hg38.canon.unmrg.bed u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed u2os.actv.txts.noTSS.hg38.canon.unmrg.bed igenic.noU2OSgenic.iTADovlp.bed igenic.noU2OSgenic.TADovlp.bed | cut -f 1-3 | sort -k1,1V -k2,2n | bedtools merge -i stdin | awk '{print $3-$2}' | awk '{total +=$1} END {print total/1e6}'
# 3031.04
#awk '{print $3-$2}' GCA_hg38.chrom_chr.nochrY.canon.sizes | awk '{total +=$1} END {print total/1e6}'
#-3031.04




# merge the genic features as they are still mututally exclusiv when merged:
u2os.actvTSS.hg38.canon.bed
u2os.inactvTSS-noactvTSSovlps.hg38.canon.unmrg.bed
u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed
u2os.actv.txts.noTSS.hg38.canon.unmrg.bed



sort -k1,1V -k2,2n u2os.actvTSS.hg38.canon.bed | bedtools merge -i stdin | sort -k1,1V -k2,2n >u2os.actvTSS.hg38.mrg.bed

sort -k1,1V -k2,2n u2os.inactvTSS-noactvTSSovlps.hg38.canon.unmrg.bed | bedtools merge -i stdin | sort -k1,1V -k2,2n >u2os.inactvTSS.hg38.mrg.bed

sort -k1,1V -k2,2n u2os.inactv.txts.noTSS.hg38.canon.unmrg.bed | bedtools merge -i stdin | sort -k1,1V -k2,2n >u2os.inactvGenBods.hg38.mrg.bed

sort -k1,1V -k2,2n u2os.actv.txts.noTSS.hg38.canon.unmrg.bed | bedtools merge -i stdin | sort -k1,1V -k2,2n >u2os.actvGenBods.hg38.mrg.bed



-rw-r--r-- 1 myanco s3it_t_cluster_users 109K Apr 27 12:11 igenic.noU2OSgenic.iTADovlp.bed
-rw-r--r-- 1 myanco s3it_t_cluster_users 614K Apr 27 12:11 igenic.noU2OSgenic.TADovlp.bed
-rw-r--r-- 1 myanco s3it_t_cluster_users 998K Apr 27 12:21 u2os.actvTSS.hg38.mrg.bed
-rw-r--r-- 1 myanco s3it_t_cluster_users 646K Apr 27 12:21 u2os.inactvTSS.hg38.mrg.bed
-rw-r--r-- 1 myanco s3it_t_cluster_users 469K Apr 27 12:21 u2os.inactvGenBods.hg38.mrg.bed
-rw-r--r-- 1 myanco s3it_t_cluster_users 989K Apr 27 12:21 u2os.actvGenBods.hg38.mrg.bed

# to do : check if they are mutually exclusive (bedt int pairwise), add up to 23 bins merged, and visualize on igv hg38


#bedtools intersect -a igenic.noU2OSgenic.iTADovlp.bed -b igenic.noU2OSgenic.TADovlp.bed | wc -l
#bedtools intersect -a igenic.noU2OSgenic.iTADovlp.bed -b u2os.actvTSS.hg38.mrg.bed | wc -l
#bedtools intersect -a igenic.noU2OSgenic.iTADovlp.bed -b u2os.inactvTSS.hg38.mrg.bed | wc -l
#bedtools intersect -a igenic.noU2OSgenic.iTADovlp.bed -b u2os.inactvGenBods.hg38.mrg.bed | wc -l
#bedtools intersect -a igenic.noU2OSgenic.iTADovlp.bed -b u2os.actvGenBods.hg38.mrg.bed | wc -l
#bedtools intersect -a igenic.noU2OSgenic.TADovlp.bed -b u2os.actvTSS.hg38.mrg.bed | wc -l
#bedtools intersect -a igenic.noU2OSgenic.TADovlp.bed -b u2os.inactvTSS.hg38.mrg.bed | wc -l
#bedtools intersect -a igenic.noU2OSgenic.TADovlp.bed -b u2os.inactvGenBods.hg38.mrg.bed | wc -l
#bedtools intersect -a igenic.noU2OSgenic.TADovlp.bed -b u2os.actvGenBods.hg38.mrg.bed | wc -l
#bedtools intersect -a u2os.actvTSS.hg38.mrg.bed -b u2os.inactvTSS.hg38.mrg.bed | wc -l
#bedtools intersect -a u2os.actvTSS.hg38.mrg.bed -b u2os.inactvGenBods.hg38.mrg.bed | wc -l
#bedtools intersect -a u2os.actvTSS.hg38.mrg.bed -b u2os.actvGenBods.hg38.mrg.bed | wc -l
#bedtools intersect -a u2os.inactvTSS.hg38.mrg.bed -b u2os.inactvGenBods.hg38.mrg.bed | wc -l
#bedtools intersect -a u2os.inactvTSS.hg38.mrg.bed -b u2os.actvGenBods.hg38.mrg.bed | wc -l
#bedtools intersect -a u2os.inactvGenBods.hg38.mrg.bed -b u2os.actvGenBods.hg38.mrg.bed | wc -l

# always 0.
cat igenic.noU2OSgenic.iTADovlp.bed igenic.noU2OSgenic.TADovlp.bed u2os.actvTSS.hg38.mrg.bed u2os.inactvTSS.hg38.mrg.bed u2os.inactvGenBods.hg38.mrg.bed u2os.actvGenBods.hg38.mrg.bed | sort -k1,1V -k2,2n | wc -l
# 163970, merged 23

cat igenic.noU2OSgenic.iTADovlp.bed igenic.noU2OSgenic.TADovlp.bed u2os.actvTSS.hg38.mrg.bed u2os.inactvTSS.hg38.mrg.bed u2os.inactvGenBods.hg38.mrg.bed u2os.actvGenBods.hg38.mrg.bed | sort -k1,1V -k2,2n | awk '{print $3-$2}' - | awk '{total +=$1} END {print total/1e6}'
# 3031.04, same if merge.  18.5kb avg feature size


# compare to hg19 stats
hg19 stats:
 # igenic regions 54%, active gene bodies 35%, inactive gene bodies 7%, active TSS±1.5kb 3%, inactive TSS±1.5kb 1%

awk '{print $3-$2}' igenic.noU2OSgenic.TADovlp.bed | awk '{total +=$1} END {print total/1e6}'
# 1050.13 Mb
awk '{print $3-$2}' igenic.noU2OSgenic.iTADovlp.bed | awk '{total +=$1} END {print total/1e6}'
# 244.654 Mb

awk '{print $3-$2}' u2os.actvGenBods.hg38.mrg.bed | awk '{total +=$1} END {print total/1e6}'
# 1070.24 Mb
awk '{print $3-$2}' u2os.inactvGenBods.hg38.mrg.bed | awk '{total +=$1} END {print total/1e6}'
# 415.687 Mb

awk '{print $3-$2}' u2os.actvTSS.hg38.mrg.bed | awk '{total +=$1} END {print total/1e6}'
# 158.566 Mb
awk '{print $3-$2}' u2os.inactvTSS.hg38.mrg.bed | awk '{total +=$1} END {print total/1e6}'
# 91.7641 Mb



=1050.13+244.654+1070.24+415.687+158.566+91.7641
# 3031.0411 Mb total


1050.13/3031.04
244.654/3031.04
1070.24/3031.04
415.687/3031.04
158.566/3031.04
91.7641/3031.04


0.346458641258446
0.08071618982263513
0.3530933277027027
0.1371433567356419
0.05231405722128379
0.030274790171030404


#35% euchr igen, 8% het igen, active gene bodies 35%, inactive gene bodies 14%, active TSS±1.5kb 5%, inactive TSS±1.5kb 3%.
# hg19: igenic > actvgenbod>inactvgenbods>actvTSS>inactvTSS
# hg38: same

# use in bwaob, barplots.

hg38.gnmcfts.mutualexclsv.txt
TotalMb	Feat	PctgGenLen
1050.13	igenicEu	0.346458641258446
244.654	igenicHet	0.08071618982263513
1070.24	actvGB	0.3530933277027027
415.687	inactvGB	0.1371433567356419
158.566	actvTSS	0.05231405722128379
91.7641	inactvTSS	0.030274790171030404


# bwaob; multibw summ





















 







