 #Simple Bash script

  339  firstname="Ridwan"
  340  lastname="Bello"
  341  echo $firstname $lastname
  342  echo $firstname
  343  echo $lastname

#Biocomputing task

  345  mkdir Ridwan
  346  mkdir biocomputing
  347  cd biocomputing/
  348  sudo apt-get install wget
  349  wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna
  350  wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
  351  mv wildtype.fna ~/Ridwan
  352  cd ..
  353  ls
  354  cd Ridwan/
  355  ls
  356  cat Ridwan
  357  clear
  358  cd
  359  cd biocomputing/
  360  ls
  361  rm wildtype.gbk.1
  362  cd
  363  cd Ridwan/
  364  grep tatata Ridwan
  365  grep -c "tatatata" Ridwan > mutant.txt
  366  cat mutant.txt
  367  grep -o tatatata Ridwan > Mutant.txt
  368  cat Mutant.txt
  369  clear

#Biocomputing Task 2

  370  sudo apt-get install figlet
  371  figlet Ridwan
  372  figlet -f script Ridwan
  373  cd ..
  374  mkdir compare
  375  cd compare/
  376  wget https://www.bioinformatics.babraham.ac.uk/training/Introduction%20to%20Unix/unix_intro_data.tar.gz
  379  gunzip unix_intro_data.tar.gz
  380  ls
  381  tar -xvf unix_intro_data.tar
  382  cd seqmonk_genomes/Saccharomyces\ cerevisiae/
  383  cd EF4/
  384  grep rRNA Mito.dat
  385  cp -r Mito.dat ~/compare/
  386  cd
  387  ls
  388  cd compare/
  389  ls
  390  nano Mito.dat
  391  mv Mito.dat Mitochondrion.txt
  392  ls
  393  cd FastQ_Data/
  394  ls
  395  wc -l lane8_DD_P4_TTAGGC_L008_R1.fastq.gz
  396  wc -l *.fastq.gz > total.txt
  397  cat total.txt
  399  exit

