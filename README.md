# Identification of new *Crassostrea* species: *C. zhujiangensis*
This walkthrough comprehensively examines all code and outputs related to the bioinformatic analysis presented in the paper titled "Morphological differentiation, genomic divergence, and reproductive isolation identify the southern and northern populations of C. ariakensis as two distinct species."

# 1. Genomic divergence
**Data trimming**
```bash
for i in `cat sample.list`
do
fastp -i ${i}_1.fastq.gz -I ${i}_2.fastq.gz -o ${i}_clean_1.fq.gz -O ${i}_clean_2.fq.gz —adapter_sequence auto —detect_adapter_for_pe —unpaired1 output_um_1.fastq.gz —unpaired2 output_um_2.fastq.gz —failed_out output_failed.fastq.gz —cut_front —cut_front_window_size=1 —cut_front_mean_quality=20 —cut_tail —cut_tail_window_size=1 —cut_tail_mean_quality=20 —cut_right —cut_right_window_size=4 —cut_right_mean_quality=20 —length_required=36 —thread 1 --trim_front1 5 --trim_front2 5
done
```
**Mapping reads to *C. ariakensis* reference geneome (Wu et al., 2022)**
```bash
bwa index ref.fa
for i in `cat sample.list`
do
bwa mem -M -t 4 ref.fa ${i}_clean_1.fq.gz ${i}_clean_2.fq.gz | samtools view -bS > $i.bam
# mark duplications with samtools v1.7
samtools sort -@ 8 -n -o namesort.bam $i.bam && \
rm -f $i.bam && \
samtools fixmate -@ 8 -m namesort.bam fixmate.bam && \
rm -f namesort.bam && \
samtools sort -@ 8 -o positionsort.bam fixmate.bam && \
rm -f fixmate.bam && \
samtools markdup -@ 8 -r positionsort.bam $i.dup.bam && \
rm -f positionsort.bam && \
samtools index -@ 8 $i.dup.bam
# remove dup and low quality mapped reads
samtools view -@ 8 -h -F 0x100 -F 0x400 $i.dup.bam | mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/'| samtools view -@ 8 -q 30 -bS > $i.bam && \
samtools index $i.bam && \
rm -f $i.dup.bam
done
````
**SNP calling and filtering with VCFtools v0.1.16 (Danecek et al., 2011).**
```bash
bcftools mpileup --threads 10 -a AD,DP,SP -Ou -f ref.fa *.dup.bam | bcftools call --threads 10 -f GQ,GP -mO z -o oyster.vcf.gz
VCF_IN=oyster.vcf.gz
VCF_OUT=filtered.vcf.gz
MAF=0.1
MISS=0.9
QUAL=30
MIN_DEPTH=10
MAX_DEPTH=50
vcftools --gzvcf $VCF_IN \
--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > $VCF_OUT
```
**Phylogenetic analysis based on SNP data utilizing (Lee et al., 2014).**
```bash
# installation of SNPhylo
conda create -n SNPhylo
conda activate SNPhylo
conda install python=2.7
conda install -c conda-forge r-base
conda install -c conda-forge r-igraph
BiocManager::install(c("getopt", "gdsfmt","SNPRelate","ape","phangorn"))
git clone https://github.com/thlee/SNPhylo.git
cd SNPhylo
sh setup.sh
# 
plink --vcf filterlowDP.vcf.gz --maf 0.05 --geno 0.1 --recode  vcf-iid --out filterlowDP2 --allow-extra-chr
plink --vcf filterlowDP2.vcf  --allow-extra-chr --indep-pairwise 50 10 0.2 --out filterlowDP2
plink --allow-extra-chr --extract filterlowDP2.prune.in --make-bed --out filterlowDP2.prune.in --recode vcf-iid --vcf filterlowDP2.vcf
~/tools/SNPhylo/snphylo.sh -v filterlowDP2.prune.in.vcf -A -b -B 300
```
**Principal components analysis (PCA) with PLINK v1.9(Purcell et al., 2007).**
```bash
# perform linkage pruning - i.e. identify prune sites
plink --vcf oyster.filtered.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out oyster
# prune and create pca
plink --vcf oyster.filtered.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract oyster.prune.in \
--threads 10 \
--make-bed --pca --out oyster
```
# 2. Genetic divergence
**Data trimming**
```bash
for i in `cat sample.list`
do
fastp -i ${i}_1.fastq.gz -I ${i}_2.fastq.gz -o ${i}_clean_1.fq.gz -O ${i}_clean_2.fq.gz —adapter_sequence auto —detect_adapter_for_pe —unpaired1 output_um_1.fastq.gz —unpaired2 output_um_2.fastq.gz —failed_out output_failed.fastq.gz —cut_front —cut_front_window_size=1 —cut_front_mean_quality=20 —cut_tail —cut_tail_window_size=1 —cut_tail_mean_quality=20 —cut_right —cut_right_window_size=4 —cut_right_mean_quality=20 —length_required=36 —thread 1 --trim_front1 5 --trim_front2 5
done
```
**Mapping reads to *C. ariakensis* reference geneome (Wu et al., 2022) composed of ten primary super scaffolds.**
```bash
bwa index ref.fa
for i in `cat sample.list`
do
bwa mem -M -t 4 ref.fa ${i}_clean_1.fq.gz ${i}_clean_2.fq.gz | samtools view -bS > $i.bam
# mark duplications with samtools v1.7
samtools sort -@ 8 -n -o namesort.bam $i.bam && \
rm -f $i.bam && \
samtools fixmate -@ 8 -m namesort.bam fixmate.bam && \
rm -f namesort.bam && \
samtools sort -@ 8 -o positionsort.bam fixmate.bam && \
rm -f fixmate.bam && \
samtools markdup -@ 8 -r positionsort.bam $i.dup.bam && \
rm -f positionsort.bam && \
samtools index -@ 8 $i.dup.bam
# remove dup and low quality mapped reads
samtools view -@ 8 -h -F 0x100 -F 0x400 $i.dup.bam | mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/'| samtools view -@ 8 -q 30 -bS > $i.bam && \
samtools index $i.bam && \
rm -f $i.dup.bam
done
````
**Assembling the unmapped reads**
```bash
for i in `cat sample.list`
do
mkdir $i
cd $i
samtools fastq -@ 10 -f 4 -1 unmapped.R1.fastq -2 unmapped.R2.fastq -s unmapped.RS.fastq $i.dup.bam
spades.py -1 unmapped.R1.fastq -2 unmapped.R2.fastq -s unmapped.RS.fastq --careful --cov-cutoff auto -o spades_assembly -t 30
cd ../
done
```
**Extracting *COI* sequences**
```bash
for i in `cat sample.list`
do
cd $i
makeblastdb -in scaffolds.fasta -dbtype nucl
blastn -query FJ841964.1.fa -db scaffolds.fasta -out $i.blt -outfmt 6
head -n 1 $i.blt > $i.head1.blt
coi=`sed -n -e 1p $i.head1.blt | awk '{print $2}'`
start=`sed -n -e 1p $i.head1.blt | awk '{print $9}'`
end=`sed -n -e 1p $i.head1.blt | awk '{print $10}'`
samtools faidx scaffolds.fasta && samtools faidx scaffolds.fasta ${coi}:${start}-${end} > ${i}_coi.dna.fa
python cds2prot.py ${i}_coi.dna.fa ${i}_coi.prot.fa ${i}_coi.nucl.fa
cd ../
done
```
**Multiple alignment with MUSCLE v3.8.31 (Edgar & Soc, 2004) and Kimura 2-parameter (K2P) distances calculation with the DNADIST program (Felsenstein, 1988).**
```bash
muscle -in merge.fa -out merge.afa
convertFasta2Phylip.sh merge.afa > infile
dnadist
```
**Visualization of K2P values**
```r
library(reshape2)

df <- read.table("dnadist.txt", row.names=1)
colnames(df) <- rownames(df)

colsouth <- df[,c("SR00000001","SR00000002","SR00000003","SR00000004","SR00000005","SR00000006","SR00000007","SR00000008","SR00000009","SR00000010","SR00000011","SR00000012","SR00000013","SR00000014","SR00000015","SR00000016","SR00000017","SR00000018","SR00000019","SR00000020","SR00000021","SR00000022","SR00000023","SR00000024","SR00000025","SR00000026","SR00000027","SR00000028","SR00000029","SR00000030","SR00000031","SR00000032","SR00000033","SR00000034","SR00000035","SR00000036","SR00000037","SR00000038","SR00000039","SR00000040","SR00000041","SR00000042","SR00000043","SR00000044","SR00000045","SR00000046","SR00000047","SR00000048","SR00000049","SR00000050","SR00000051","SR00000052","SR00000053","SR00000054","SR00000055","SR00000056","SR00000057","SR00000058","SR00000059","SR00000060","SR00000061","SR00000062","SR00000063","SR00000064","SR00000065","SR00000066","SR00000067","SR00000068","SR00000069")]
rowsouth_colsouth <- colsouth[c("SR00000001","SR00000002","SR00000003","SR00000004","SR00000005","SR00000006","SR00000007","SR00000008","SR00000009","SR00000010","SR00000011","SR00000012","SR00000013","SR00000014","SR00000015","SR00000016","SR00000017","SR00000018","SR00000019","SR00000020","SR00000021","SR00000022","SR00000023","SR00000024","SR00000025","SR00000026","SR00000027","SR00000028","SR00000029","SR00000030","SR00000031","SR00000032","SR00000033","SR00000034","SR00000035","SR00000036","SR00000037","SR00000038","SR00000039","SR00000040","SR00000041","SR00000042","SR00000043","SR00000044","SR00000045","SR00000046","SR00000047","SR00000048","SR00000049","SR00000050","SR00000051","SR00000052","SR00000053","SR00000054","SR00000055","SR00000056","SR00000057","SR00000058","SR00000059","SR00000060","SR00000061","SR00000062","SR00000063","SR00000064","SR00000065","SR00000066","SR00000067","SR00000068","SR00000069"),]
rowsouth_colsouth <- apply(rowsouth_colsouth, 1, median)
summary(rowsouth_colsouth)

rownorth_colsouth <- colsouth[c("SR00000070","SR00000071","SR00000072","SR00000073","SR00000074","SR00000075","SR00000076","SR00000077","SR00000078","SR00000079","SR00000080","SR00000081","SR00000082","SR00000083","SR00000084","SR00000085","SR00000086","SR00000087","SR00000088","SR00000089","SR00000090","SR00000091","SR00000092","SR00000093","SR00000094","SR00000095","SR00000096","SR00000097","SR00000098","SR00000099","SR00000100","SR00000101","SR00000102","SR00000103","SR00000104","SR00000105","SR00000106","SR00000107","SR00000108","SR00000109","SR00000110","SR00000111","SR00000112","SR00000113","SR00000114","SR00000115","SR00000116","SR00000117","SR00000118","SR00000119","SR00000120","SR00000121","SR00000122","SR00000123","SR00000124","SR00000125","SR00000126","SR00000127","SR00000128","SR00000129","SR00000130","SR00000131","SR00000132","SR00000133","SR00000134","SR00000135","SR00000136","SR00000137","SR00000138","SR00000139","SR00000140","SR00000141","SR00000142","SR00000143","SR00000144","SR00000145","SR00000146","SR00000147","SR00000148","SR00000149","SR00000150","SR00000151","SR00000152","SR00000153","SR00000154","SR00000155","SR00000156","SR00000157","SR00000158","SR00000159","SR00000160","SR00000161","SR00000162","SR00000163","SR00000164","SR00000165","SR00000166","SR00000167","SR00000168","SR00000169","SR00000170","SR00000171","SR00000172","SR00000173","SR00000174","SR00000175","SR00000176","SR00000177","SR00000178","SR00000179","SR00000180","SR00000181","SR00000182","SR00000183","SR00000184","SR00000185","SR00000186","SR00000187","SR00000188","SR00000189","SR00000190","SR00000191","SR00000192","SR00000193","SR00000194","SR00000195","SR00000196","SR00000197","SR00000198","SR00000199","SR00000200","SR00000201","SR00000202","SR00000203","SR00000204","SR00000205","SR00000206","SR00000207","SR00000208","SR00000209","SR00000210","SR00000211","SR00000212","SR00000213","SR00000214","SR00000215","SR00000216","SR00000217","SR00000218","SR00000219","SR00000220","SR00000221","SR00000222","SR00000223","SR00000224","SR00000225","SR00000226","SR00000227","SR00000228","SR00000229","SR00000230","SR00000231","SR00000232","SR00000233","SR00000234","SR00000235","SR00000236","SR00000237","SR00000238","SR00000239","SR00000240","SR00000241","SR00000242","SR00000243","SR00000244","SR00000245","SR00000246","SR00000247","SR00000248","SR00000249","SR00000250","SR00000251","SR00000252","SR00000253","SR00000254","SR00000255","SR00000256","SR00000257","SR00000258","SR00000259","SR00000260","SR00000261"),]
rownorth_colsouth <- apply(rownorth_colsouth, 1, median)
summary(rownorth_colsouth)

colnorth <- df[,c("SR00000070","SR00000071","SR00000072","SR00000073","SR00000074","SR00000075","SR00000076","SR00000077","SR00000078","SR00000079","SR00000080","SR00000081","SR00000082","SR00000083","SR00000084","SR00000085","SR00000086","SR00000087","SR00000088","SR00000089","SR00000090","SR00000091","SR00000092","SR00000093","SR00000094","SR00000095","SR00000096","SR00000097","SR00000098","SR00000099","SR00000100","SR00000101","SR00000102","SR00000103","SR00000104","SR00000105","SR00000106","SR00000107","SR00000108","SR00000109","SR00000110","SR00000111","SR00000112","SR00000113","SR00000114","SR00000115","SR00000116","SR00000117","SR00000118","SR00000119","SR00000120","SR00000121","SR00000122","SR00000123","SR00000124","SR00000125","SR00000126","SR00000127","SR00000128","SR00000129","SR00000130","SR00000131","SR00000132","SR00000133","SR00000134","SR00000135","SR00000136","SR00000137","SR00000138","SR00000139","SR00000140","SR00000141","SR00000142","SR00000143","SR00000144","SR00000145","SR00000146","SR00000147","SR00000148","SR00000149","SR00000150","SR00000151","SR00000152","SR00000153","SR00000154","SR00000155","SR00000156","SR00000157","SR00000158","SR00000159","SR00000160","SR00000161","SR00000162","SR00000163","SR00000164","SR00000165","SR00000166","SR00000167","SR00000168","SR00000169","SR00000170","SR00000171","SR00000172","SR00000173","SR00000174","SR00000175","SR00000176","SR00000177","SR00000178","SR00000179","SR00000180","SR00000181","SR00000182","SR00000183","SR00000184","SR00000185","SR00000186","SR00000187","SR00000188","SR00000189","SR00000190","SR00000191","SR00000192","SR00000193","SR00000194","SR00000195","SR00000196","SR00000197","SR00000198","SR00000199","SR00000200","SR00000201","SR00000202","SR00000203","SR00000204","SR00000205","SR00000206","SR00000207","SR00000208","SR00000209","SR00000210","SR00000211","SR00000212","SR00000213","SR00000214","SR00000215","SR00000216","SR00000217","SR00000218","SR00000219","SR00000220","SR00000221","SR00000222","SR00000223","SR00000224","SR00000225","SR00000226","SR00000227","SR00000228","SR00000229","SR00000230","SR00000231","SR00000232","SR00000233","SR00000234","SR00000235","SR00000236","SR00000237","SR00000238","SR00000239","SR00000240","SR00000241","SR00000242","SR00000243","SR00000244","SR00000245","SR00000246","SR00000247","SR00000248","SR00000249","SR00000250","SR00000251","SR00000252","SR00000253","SR00000254","SR00000255","SR00000256","SR00000257","SR00000258","SR00000259","SR00000260","SR00000261")]
rownorth_colnorth <- colnorth[c("SR00000070","SR00000071","SR00000072","SR00000073","SR00000074","SR00000075","SR00000076","SR00000077","SR00000078","SR00000079","SR00000080","SR00000081","SR00000082","SR00000083","SR00000084","SR00000085","SR00000086","SR00000087","SR00000088","SR00000089","SR00000090","SR00000091","SR00000092","SR00000093","SR00000094","SR00000095","SR00000096","SR00000097","SR00000098","SR00000099","SR00000100","SR00000101","SR00000102","SR00000103","SR00000104","SR00000105","SR00000106","SR00000107","SR00000108","SR00000109","SR00000110","SR00000111","SR00000112","SR00000113","SR00000114","SR00000115","SR00000116","SR00000117","SR00000118","SR00000119","SR00000120","SR00000121","SR00000122","SR00000123","SR00000124","SR00000125","SR00000126","SR00000127","SR00000128","SR00000129","SR00000130","SR00000131","SR00000132","SR00000133","SR00000134","SR00000135","SR00000136","SR00000137","SR00000138","SR00000139","SR00000140","SR00000141","SR00000142","SR00000143","SR00000144","SR00000145","SR00000146","SR00000147","SR00000148","SR00000149","SR00000150","SR00000151","SR00000152","SR00000153","SR00000154","SR00000155","SR00000156","SR00000157","SR00000158","SR00000159","SR00000160","SR00000161","SR00000162","SR00000163","SR00000164","SR00000165","SR00000166","SR00000167","SR00000168","SR00000169","SR00000170","SR00000171","SR00000172","SR00000173","SR00000174","SR00000175","SR00000176","SR00000177","SR00000178","SR00000179","SR00000180","SR00000181","SR00000182","SR00000183","SR00000184","SR00000185","SR00000186","SR00000187","SR00000188","SR00000189","SR00000190","SR00000191","SR00000192","SR00000193","SR00000194","SR00000195","SR00000196","SR00000197","SR00000198","SR00000199","SR00000200","SR00000201","SR00000202","SR00000203","SR00000204","SR00000205","SR00000206","SR00000207","SR00000208","SR00000209","SR00000210","SR00000211","SR00000212","SR00000213","SR00000214","SR00000215","SR00000216","SR00000217","SR00000218","SR00000219","SR00000220","SR00000221","SR00000222","SR00000223","SR00000224","SR00000225","SR00000226","SR00000227","SR00000228","SR00000229","SR00000230","SR00000231","SR00000232","SR00000233","SR00000234","SR00000235","SR00000236","SR00000237","SR00000238","SR00000239","SR00000240","SR00000241","SR00000242","SR00000243","SR00000244","SR00000245","SR00000246","SR00000247","SR00000248","SR00000249","SR00000250","SR00000251","SR00000252","SR00000253","SR00000254","SR00000255","SR00000256","SR00000257","SR00000258","SR00000259","SR00000260","SR00000261"),]
rownorth_colnorth <- apply(rownorth_colnorth, 1, median)
summary(rownorth_colnorth)

type <- c(rep("south.vs.south",length(rowsouth_colsouth)),rep("north.vs.north", length(rownorth_colnorth)),rep("north.vs.south", length(rownorth_colsouth)))
res <- data.frame(type=type, dist=c(rowsouth_colsouth, rownorth_colnorth, rownorth_colsouth))

summary(rowsouth_colsouth)
summary(rownorth_colnorth)
summary(rownorth_colsouth)

library(ggplot2)
library(Hmisc)

pdf('dnadist2.pdf', width=3, height=2.8)
p <- ggplot(res, aes(x=factor(type, level=c("south.vs.south", "north.vs.north", "north.vs.south")), y=dist)) + geom_violin(trim=FALSE)

p + geom_jitter(shape=16, position=position_jitter(0.1)) + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red") + theme_minimal()
  
dev.off()

pdf('dnadist.pdf', width=3, height=2.7)
ggplot(res, aes(x=factor(type, level=c("south.vs.south", "north.vs.north", "north.vs.south")), y=dist)) + 
    geom_boxplot(binaxis='y', stackdir='center') + theme_minimal()

dev.off()
```
**Phylogenetic analysis based on *COI* sequences utiliziong RAxML v8.2.12**
```bash
convertFasta2Phylip.sh merge.afa > merge.phy
raxmlHPC-PTHREADS -f a -# 100 -m GTRGAMMA -p 12345 -x 12345 -s merge.phy -n output.tree -T 40 -o cangulata
```
**Variability in shell characters**
```r
library(ggpubr)

df <- read.table("morphology.txt", header=T)

df$length2height <- df$length/df$height
df$width2weight <- df$width/df$weight

cmpr <- list(c("Northern_Cari","Southern_Cari"))

pdf("length2height.pdf",height=3,width=3)
# Create bar plots of means
ggbarplot(df, x = "type", y = "length2height", 
          add = c("mean_se", "jitter"),
          fill = "type", color = "black", palette = c("#D4605F", "#D4605F"),
          position = position_dodge(0.3)) +
          stat_compare_means(comparisons = cmpr) + 
          scale_y_continuous(expand = c(0,0), limits=c(0,1.3)) + 
          labs(title = "", x = "", y = "Shell length to height ratio") + 
          theme(legend.position = "none")
dev.off()
```

# References
EDGAR, R.C. 2004 MUSCLE: a multiple sequence alignment method with reduced time and space complexity. BMC Bioinformatics, 5: 1-19.

FELSENSTEIN, J. 1988 Phylogenies from molecular sequences: inference and reliability. Annual Review of Genetics, 22: 521-565.

LEE, T.-H., GUO, H., WANG, X., KIM, C. & PATERSON, A.H. 2014 SNPhylo: a pipeline to construct a phylogenetic tree from huge SNP data. Bmc Genomics, 15.

PURCELL, S., NEALE, B., TODD-BROWN, K., THOMAS, L., FERREIRA, M. A. R., BENDER, D., ... SHAM, P.C. 2007 PLINK: A tool set for whole-genome association and population-based linkage analyses. American Journal of Human Genetics, 81, 559-575.

STAMATAKIS, A., LUDWIG, T. & MEIER, H. 2005 RAxML-III: a fast program for maximum likelihood-based inference of large phylogenetic trees. Bioinformatics, 21: 456-463.

WU, B., CHEN, X., YU, M.J., REN, J.F., HU, J., SHAO, C.W., ZHOU, L.Q., SUN, X.J., YU, T., ZHENG, Y.X., WANG, Y., WANG, Z.Y., ZHANG, H., FAN, G.Y. & LIU, Z.H. 2022 Chromosome-level genome and population genomic analysis provide insights into the evolution and environmental adaptation of Jinjiang oyster Crassostrea ariakensis. Molecular Ecology Resources, 22: 1529-1544.
