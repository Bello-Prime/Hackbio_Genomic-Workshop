#! bin/bash

echo –e “\n Beginning analysis…”
#download datasets from the database in fastq.qz format
mkdir –p normal_data
cd normal_data
wget $normal_data download link

echo – “\n Running fastp for trimming adapters…”
mkdir –p trimmed
for file in *.fastq.gz
do
base=$(basename $file .fastq.gz)
echo $base
fastp –i $file –o trimmed/${base}.fastq.gz
--html trimmed/${base}_fastp.html
done

echo –e “\n Running fastqc and multiqc…”
cd trimmed

mkdir –p fastqc_reports
for file in *.fastq.gz
do
fastqc $file –o ../fastqc_reports
done

cd ..

multiqc fastqc_reports –o fastqc_reports

#download reference genome
echo –e “\n getting the reference genome…”
mkdir –p ref
cd ref
echo –e “insert link to download reference genome of interest… 
#insert link to refseq
Read refgenomelink
wget $refgenomelink
#gunzip reference genome
gunzip *.fa
 cd ..
#Removing low quality sequences using Trimmomatic
ls *.fastq.gz > list.txt
cat list.txt
#create trimming.sh with nano
nano trimming.sh (copy and paste this script)
 mkdir -p trimmed_reads
for sample in `cat list.txt`
   do
trimmomatic PE -threads 8 ${sample} *.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:keepBothReads \
LEADING:3 TRAILING:10 MINLEN:25
   done

cd trimmed
bash trimming.sh

echo –e “\n indexing the reference genome…”
#Read mapping using bwa
#Index reference file
bwa index *.fa

# align by creating a aligner.sh in nano
nano aligner.sh
mkdir Mapped
#Perform alignment
bwa mem -R 4 ref/*.fa trimmed/*N*r1* trimmed/*N*r2*  | samtools view -@ 4 -S -b | samtools sort -@ 4 >  mapped/data.sorted.bam

samtools index mapped/${sample}.sorted.bam
#save and exit nano
bash aligner.sh
#Mapped reads filtering
echo –e “\n samtools filtering…”
#create filter.sh for one command filtering with nano
nano filter.sh
for sample in `list.txt`
do
samtools view -q 1 -f 0x2 -F 0x8 -b ${sample}.sorted.bam > ${sample}.filtered1.bam
done
#save and exit nano
bash filter.sh

#Duplicates removal
echo –e “\n PCR duplicate removal…”
#create duplicate.sh for one command duplicate removal with nano
  for sample in `cat list.txt`
  do
samtools sort -n -o Mapped/${sample}.namesort.bam $Mapped/{sample}.filtered.bam
samtools fixmate -m Mapped/${sample}.namesort.bam Mapped/${sample}.fixmate.bam
samtools sort -@ 4 -o Mapped${sample}.positionsort.bam Mapped${sample}.fixmate.bam
samtools markdup -@4 -r Mapped/${sample}.positionsort.bam Mapped/${sample}.clean.bam
 done  
#save and exit nano

#Left Align BAM
echo –e “\n echo left alignment…”
for sample in `cat list.txt`
do   
cat Mapped/${sample}.clean.bam  | bamleftalign -f *.fa -m 5 -c > Mapped/${sample}.leftAlign.bam 
done

Bash leftalignBAM.sh

#Recalibrate read mapping qualities
echo –e “\n recalibrating…”
for sample in `cat list.txt`
do
samtools calmd -@ 4 -b Mapped${sample}.leftAlign.bam *.fa > Mapped/${sample}.recalibrate.bam
done
bash recalibrate.sh

#Refilter read mapping qualities
echo –e “\n refiltering…”
for sample in `cat list.txt`
do
bamtools filter -in Mapped/${sample}.recalibrate.bam -mapQuality “<=254”  > Mapped/${sample}.refilter.bam
done


#convert variant to pileup
echo –e “\n  converting data to pileup for variant calling…”
mkdir Variants
for sample in `cat list.txt`
do
samtools mpileup -f *.fa Mapped/${sample}.refilter.bam --min-MQ 1 --min-BQ 28 \ > Variants/${sample}.pileup
done

#Call variants
echo –e “\n calling variants…”
VarScan somatic Variants/data.pileup \ Variants/finalresult  --sample-purity 1   --output-vcf 1

#merge vcf using bcftools
echo –e “\n Merge VCF files…”
bgzip Variants/finalresult.snp.vcf > Variants/finalresult.snp.vcf.gz
bgzip Variants/finalresult.indel.vcf > Variants/finalresult.indel.vcf.gz
tabix Variants/finalresult.snp.vcf.gz
tabix Variants/finalresult.indel.vcf.gz
bcftools merge –force-samples Variants/finalresult.snp.vcf.gz Variants/finalresult.indel.vcf.gz > Variants/finalresult.vcf

#variant annotation
snpEff  Variants/finalresult.vcf > Variants/finalresult.ann.vc
       





     






