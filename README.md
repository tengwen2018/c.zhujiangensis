![image](https://github.com/tengwen2018/c.zhujiangensis/assets/115065520/fd3ad4da-d76f-4024-9ada-b3dd18401beb)![image](https://github.com/tengwen2018/c.zhujiangensis/assets/115065520/8f3e1aa0-35fb-4097-abae-9083393c8fc9)# Identification of new *Crassostrea* species: *C. zhujiangensis*
This walkthrough comprehensively examines all code and outputs related to the bioinformatic analysis presented in the paper titled "Morphological differentiation, genomic divergence, and reproductive isolation identify the southern and northern populations of C. ariakensis as two distinct species."

# 1. Genomic divergence
**Data trimming**
```bash
for i in `cat sample.list`
do
fastp -i ${i}_1.fastq.gz -I ${i}_2.fastq.gz -o ${i}_clean_1.fq.gz -O ${i}_clean_2.fq.gz —adapter_sequence auto —detect_adapter_for_pe —unpaired1 output_um_1.fastq.gz —unpaired2 output_um_2.fastq.gz —failed_out output_failed.fastq.gz —cut_front —cut_front_window_size=1 —cut_front_mean_quality=20 —cut_tail —cut_tail_window_size=1 —cut_tail_mean_quality=20 —cut_right —cut_right_window_size=4 —cut_right_mean_quality=20 —length_required=36 —thread 1 --trim_front1 5 --trim_front2 5
done
```
**Map reads to *C. ariakensis* reference geneome (Wu et al., 2022)**
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



# References
LEE, T.-H., GUO, H., WANG, X., KIM, C. & PATERSON, A.H. 2014 SNPhylo: a pipeline to construct a phylogenetic tree from huge SNP data. Bmc Genomics, 15.
PURCELL, S., NEALE, B., TODD-BROWN, K., THOMAS, L., FERREIRA, M. A. R., BENDER, D., ... SHAM, P.C. 2007 PLINK: A tool set for whole-genome association and population-based linkage analyses. American Journal of Human Genetics, 81, 559-575.
WU, B., CHEN, X., YU, M.J., REN, J.F., HU, J., SHAO, C.W., ZHOU, L.Q., SUN, X.J., YU, T., ZHENG, Y.X., WANG, Y., WANG, Z.Y., ZHANG, H., FAN, G.Y. & LIU, Z.H. 2022 Chromosome-level genome and population genomic analysis provide insights into the evolution and environmental adaptation of Jinjiang oyster Crassostrea ariakensis. Molecular Ecology Resources, 22: 1529-1544.
