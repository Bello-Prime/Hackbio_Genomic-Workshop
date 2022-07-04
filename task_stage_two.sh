1.	#! bin/bash
2.	#ssh into the server and enter password blindly

3.	cd einstein
4.	mkdir ridwan
5.	mkdir raw_data
6.	cd raw_data

7.	#download datasets from zenodo
8.	wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
9.	wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
10.	wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
11.	wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz     
12.	
13.	# download reference genome
14.	mkdir ref
15.	cd ref
16.	wget https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz
17.	#unzip reference
18.	unzip hg19.chr5_12_17.fa.gz
19.	
20.	cd ..
21.	#download trimmomatic
22.	wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
23.	#unzip Trimmomatic-0.39.zip
24.	unzip Trimmomatic-0.39.zip
25.	        cp Trimmomatic-039/adapters/TruSeq3-PE.fa $PWD
26.	
27.	        nano list.txt 
28.	   
29.	   SLGFSK-N_231335 
30.	   SLGFSK-T_231336
31.	        cat list.txt
32.	
33.	        #create trimming.sh with nano
34.	        nano trimming.sh 
35.	        
36.	        
37.	mkdir -p trimmed_reads
38.	        for sample in `cat list.txt`
39.	        do
40.	        trimmomatic PE -threads 8 ${sample}_r1_chr5_12_17.fastq.gz ${sample}_r2_chr5_12_17.fastq.gz \
41.	        trimmed_reads/${sample}_r1_paired.fq.gz trimmed_reads/${sample}_r1_unpaired.fq.gz \
42.	        trimmed_reads/${sample}_r2_paired.fq.gz trimmed_reads/${sample}_r2_unpaired.fq.gz \
43.	        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:keepBothReads \
44.	        LEADING:3 TRAILING:10 MINLEN:25
45.	        
46.	        done 
47.	             
48.	        cd ..
49.	        bash trimming.sh
50.	
51.	#Read mapping using bwa
52.	        #Index reference file   
53.	        bwa index hg19.chr5_12_17.fa 
54.	        
55.	
56.	        # align by creating a aligner.sh in nano
57.	        nano aligner.sh
58.	        
59.	      mkdir Mapping
60.	        
61.	        #Perform alignment
62.	        bwa mem -R '@RG\tID:231335\tSM:Normal' ref/hg19.chr5_12_17.fa trimmed_reads/SLGFSK-N_231335_r1_paired.fq.gz \
63.	        trimmed_reads/SLGFSK-N_231335_r2_paired.fq.gz > Mapping/SLGFSK-N_231335.sam
64.	        
65.	
66.	        bwa mem -R '@RG\tID:231336\tSM:Tumor' ref/hg19.chr5_12_17.fa trimmed_reads/SLGFSK-T_231336_r1_paired.fq.gz \
67.	        trimmed_reads/SLGFSK-T_231336_r2_paired.fq.gz > Mapping/SLGFSK-T_231336.sam     
68.	        #save and exit nano
69.	        
70.	
71.	        bash aligner.sh
72.	
73.	#Conversion of the SAM file to BAM file, sorting and indexing
74.	        nano bam.sh
75.	        for sample in `cat list.txt`
76.	        samtools view -@ 4 -S -b ${sample}.sam | samtools sort -@ 4 > ${sample}.sorted.bam
77.	        samtools index ${sample}.sorted.bam
78.	        done
79.	        
80.	        exit nano 
81.	        bash bam.sh
82.	      
83.	
84.	        #Mapped reads filtering
85.	        
86.	        nano filter.sh
87.	        for sample in `list.txt` 
88.	     #command filtering with nano 
89.	        do
90.	      
91.	        samtools view -q 1 -f 0x2 -F 0x8 -b ${sample}.sorted.bam > ${sample}.filtered.bam
92.	        Done
93.	
94.	        #save and exit nano
95.	      
96.	        bash filter.sh
97.	        
98.	
99.	        #Duplicates removal
100.	        #to remove duplicates from the filtered file
101.	        #create duplicate.sh for one command duplicate removal with nano
102.	        nano duplicate.sh
103.	
104.	        for sample in `cat list.txt`
105.	        do
106.	        samtools sort -n -o Mapping/${sample}.namesort.bam $Mapping/{sample}.filtered.bam
107.	        samtools fixmate -m Mapping/${sample}.namesort.bam Mapping/${sample}.fixmate.bam
108.	        samtools sort -@ 32 -o Mapping/${sample}.positionsort.bam Mapping/${sample}.fixmate.bam
109.	        samtools markdup -@32 -r Mapping/${sample}.positionsort.bam Mapping/${sample}.clean.bam
110.	        done       
111.	#save and exit nano
112.	        
113.	        #Left Align BAM  
114.	nano leftalignBAM.sh 
115.	        for sample in `cat list.txt`
116.	        do      
117.	        cat Mapping/${sample}.clean.bam  | bamleftalign -f hg19.chr5_12_17.fa -m 5 -c > Mapping/${sample}.leftAlign.bam    
118.	        
119.	        Bash leftalignBAM.sh
120.	
121.	        #Recalibrate read mapping qualities 
122.	        for sample in `cat list.txt`
123.	        do
124.	                  samtools calmd -@ 32 -b Mapping${sample}.leftAlign.bam hg19.chr5_12_17.fa > Mapping/${sample}.recalibrate.bam
125.	        done 
126.	        
127.	Bash recalibrate.sh
128.	
129.	#Refilter read mapping qualities
130.	        for sample in `cat list.txt`
131.	        done
132.	        bamtools filter -in Mapping/${sample}.recalibrate.bam -mapQuality “<=254”  > Mapping/${sample}.refilter.bam
133.	        
134.	        
135.	        #Variant calling and classification
136.	        wget  https://sourceforge.net/projects/varscan/files/VarScan.v2.3.9.jar    
137.	        
138.	
139.	        #convert variant to pileup
140.	        mkdir Variants
141.	        
142.	
143.	        for sample in `cat list.txt`
144.	        do
145.	        samtools mpileup -f hg19.chr5_12_17.fa Mapping/${sample}.refilter.bam --min-MQ 1 --min-BQ 28 \  > Variants/${sample}.pileup
146.	done
147.	        
148.	        #Call variants
149.	                  java -jar VarScan.v2.3.9.jar somatic Variants/SLGFSK-N_231335.pileup \
150.	        Variants/SLGFSK-T_231336.pileup Variants/SLGFSK \
151.	        --normal-purity 1  --tumor-purity 0.5 --output-vcf 1 
152.	        
153.	        
154.	
155.	        #merge vcf using bcftools
156.	                   bgzip Variants/SLGFSK.snp.vcf > Variants/SLGFSK.snp.vcf.gz
157.	                   bgzip Variants/SLGFSK.indel.vcf > Variants/SLGFSK.indel.vcf.gz
158.	                   tabix Variants/SLGFSK.snp.vcf.gz
159.	                   tabix Variants/SLGFSK.indel.vcf.gz
160.	                  bcftools merge –force-samples Variants/SLGFSK.snp.vcf.gz Variants/SLGFSK.indel.vcf.gz > Variants/SLGFSK.vcf
161.	        
162.	
163.	        #variant annotation
164.	        snpEff hg19 Variants/SLGFSK.vcf > Variants/SLGFSK.ann.vc
