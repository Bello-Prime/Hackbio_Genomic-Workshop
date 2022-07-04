1.	mkdir stage1
2.	cd stage1
3.	wget https://raw.githubusercontent.com/HackBio-Internship/wale-home-tasks/main/DNA.fa
4.	
5.	#Question1
6.	#bash code for counting number of DNA sequence
7.	grep -c "^>" DNA.fa 
8.	
9.	
10.	#Question 2
11.	# bash code for counting total occurrence of A,G,T,C
12.	grep -o 'A\|T\|G\|C' DNA.fa | sort | uniq -c
13.	
14.	
15.	#Question3
16.	#Setting up a conda environment
17.	wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh
18.	chmod +x Miniconda3-py38_4.12.0-Linux-x86_64.sh
19.	./Miniconda3-py38_4.12.0-Linux-x86_64.sh
20.	conda activate base
21.	conda --version
22.	
23.	
24.	#Installing fastqc, multiqc and fastp
25.	conda install -c bioconda fastqc
26.	conda install -c bioconda multiqc
27.	conda install -c bioconda fastp
28.	
29.	
30.	#Downloading Alsen, Chara and Drysdale R1 & R2 datasets
31.	wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Alsen_R1.fastq.gz?raw=true/  -O Alsen_R1.fastq.gz
32.	wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Alsen_R2.fastq.gz?raw=true/ -O Alsen_R2.fastq.gz
33.	wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Chara_R1.fastq.gz?raw=true -O Chara_R1.fastq.gz
34.	wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Chara_R2.fastq.gz?raw=true -O Chara_R2.fastq.gz
35.	wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Drysdale_R1.fastq.gz?raw=true/ -O Drysdale_R1.fastq.gz
36.	wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Drysdale_R2.fastq.gz?raw=true/ -O Drysdale_R2.fastq.gz
37.	
38.	
39.	#creating a folder 'output'
40.	mkdir output
41.	
42.	
43.	#implementing fastqc on the datasets
44.	fastqc *.fastq.gz -O output/
45.	
46.	
47.	#implementing multiqc on the qc datasets
48.	multiqc output
49.	mv multiqc_report.html output 
50.	mv multiqc_data output
51.	
52.	
53.	#implementing fastp on datasets
54.	fastp -i Alsen_R1.fastq.gz -o output/Alsen_R1.fastq.gz
55.	fastp -i Alsen_R2.fastq.gz -o output/Alsen_R2.fastq.gz
56.	fastp -i Chara_R1.fastq.gz -o output/Chara_R1.fastq.gz
57.	fastp -i Chara_R2.fastq.gz -o output/Chara_R2.fastq.gz
58.	fastp -i Drysdale_R1.fastq.gz -o output/Drysdale_R1.fastq.gz
59.	fastp -i Drysdale_R2.fastq.gz -o output/Drysdale_R2.fastq.gz
