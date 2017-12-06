#!/bin/bash
#########################################################
#  Reference-guided de novo assembly - IDBA
# ====================================================
# by Heidi Lischer, 2015/2016
#########################################################

# set variables #########################################
workPathFiles=/home/hlischer/A_ly
ref=/home/hlischer/A_tal/GCF_000001735.3_TAIR10_genomic.fna
refRed=/home/hlischer/A_tal/TAIR10_10kb
primerFile=/home/hlischer/Programs/AdapterSeq_new.fa
primerFileMP=/home/hlischer/Programs/AdapterSeqMP_new.fa

NThreads=8      # set the number of threads of every parallelizable step
maxReadLength=100
kmer=61         #define best K (need to be adapted)

# paired-end libraries -------------------
name=Aly_sim           # set name of your species
lib=(150 200 400)      # set insertion libraries
insLow=(50 50 100)     # lower bound of insertion size
insHigh=(300 350 750)  # upper bound of insertion size
libsd=(34 36 87)       # sd of insertion size

# list of files with forward reads according to lib array
reads1=(Aly_150_R1.fq Aly_200_R1.fq Aly_500_R1.fq)
# list of files with rewerse reads according to lib array
reads2=(Aly_150_R2.fq Aly_200_R2.fq Aly_500_R2.fq)
# short names of libraries
shortNames=(Aly_150 Aly_200 Aly_500) 


# mate-pair libraries ---------------------
mateName=Aly_sim_Mate
mateLib=(3000 5000 7000 11000 15000)      # set insertion libraries
mateInsLow=(2000 4000 6000 9000 13000)    # lower bound of insertion size
mateInsHigh=(4000 6000 8000 13000 17000)  # upper bound of insertion size
mateLibsd=(400 400 400 400 400)           # sd of insertion size

# list of files with forward reads according to lib array
mateReads1=(Aly_3kb_1.fq Aly_5kb_1.fq Aly_7kb_1.fq Aly_11kb_1.fq Aly_15kb_1.fq)
# list of files with rewerse reads according to lib array
mateReads2=(Aly_3kb_2.fq Aly_5kb_2.fq Aly_7kb_2.fq Aly_11kb_2.fq Aly_15kb_2.fq)
# short names of libraries
mateShortNames=(Aly_3kb Aly_5kb Aly_7kb Aly_11kb Aly_15kb)


# set work path ---------------------------
workPath=${workPathFiles}/${name}_idba
# log file
log=${workPath}/log_${name}_idba.txt


# Programs --------------------------------
progPath=/home/hlischer/Programs
progFastQC=${progPath}/FastQC/fastqc
progTrimmomatic=${progPath}/Trimmomatic-0.32/trimmomatic-0.32.jar
progSamtools=samtools
progVcfutils=${progPath}/bcftools-1.3/vcfutils.pl
progBcftools=${progPath}/bcftools-1.3/bcftools
progBamtools=${progPath}/bamtools/bin/bamtools-2.3.0
progBedtools=${progPath}/bedtools2-2.19.1/bin/bedtools
progPicard=${progPath}/picard-tools-1.109
progBowtie2=${progPath}/bowtie2-2.2.1/bowtie2
progSeqtk=${progPath}/seqtk-master/seqtk
progIdba=${progPath}/idba-1.1.1/bin/
progAmos=${progPath}/amos-3.1.0/bin/
progNucmer=${progPath}/mummer-4.0.0beta2/nucmer
progGatk=${progPath}/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar
progSoapdenovo2=${progPath}/SOAPdenovo2-src-r240
progRemovShortSeq=${progPath}/RemoveShortSeq.jar
progGetBlocks=${progPath}/GetBlocks.jar
progFastaToAmos=${progPath}/FastaToAmos.jar
progWriteSoapConfig=${progPath}/WriteSoapConfig.jar
progFastaStats=${progPath}/FastaStats.jar
progSplitSeqLowCov=${progPath}/SplitSeqLowCov.jar

#########################################################



# run pipeline ##########################################
mkdir ${workPath}


# 1. Step: quality/adapter trimming and quality check:
#######################################################
  # quality check ----------
  echo "quality check of raw reads..."
  echo "quality check of raw reads..." > $log
  
  cd $workPathFiles
  fastqcOut=${workPathFiles}/FastQC_het
  mkdir -p ${fastqcOut}
    
  for i in ${!lib[*]}  #for all indexes in the array
  do 
    ${progFastQC} -t ${NThreads} -o ${fastqcOut} ${reads1[i]} ${reads2[i]}
  done  
  
  for i in ${!mateLib[*]}  #for all indexes in the array
  do 
    ${progFastQC} -t ${NThreads} -o ${fastqcOut} ${mateReads1[i]} ${mateReads2[i]}
  done  
  
  
  # quality/adapter trimming --------
  # - remove Illumina adapters provided in the primer files
  # - remove leading and trailing low quality basses (<3) or N
  # - 4 base sliding window -> remove when average quality is < 15
  # - remove reads which are shorter than 40 bp
  echo "quality/adapter trimming..."
  echo "quality/adapter trimming..." >> $log
  
  trimOut=${workPathFiles}/Trim_het
  mkdir $trimOut
  
  read1TrimPair=()
  read1TrimUnPair=()
  read2TrimPair=()
  read2TrimUnPair=()
  for i in ${!lib[*]}  #for all indexes in the array
  do
    read1TrimPair[i]=${trimOut}/${shortNames[i]}_R1_trimPair.fastq 
    read1TrimUnPair[i]=${trimOut}/${shortNames[i]}_R1_trimUnPair.fastq 
    read2TrimPair[i]=${trimOut}/${shortNames[i]}_R2_trimPair.fastq 
    read2TrimUnPair[i]=${trimOut}/${shortNames[i]}_R2_trimUnPair.fastq 
    
    java -jar ${progTrimmomatic} PE -threads ${NThreads} ${reads1[i]} ${reads2[i]} ${read1TrimPair[i]} ${read1TrimUnPair[i]} ${read2TrimPair[i]} ${read2TrimUnPair[i]} ILLUMINACLIP:${primerFile}:2:30:7:5:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40 >> $log
  done
  
  mateRead1TrimPair=()
  mateRead1TrimUnPair=()
  mateRead2TrimPair=()
  mateRead2TrimUnPair=()
  for i in ${!mateLib[*]}  #for all indexes in the array
  do
    mateRead1TrimPair[i]=${trimOut}/${mateShortNames[i]}_R1_trimPair.fastq 
    mateRead1TrimUnPair[i]=${trimOut}/${mateShortNames[i]}_R1_trimUnPair.fastq 
    mateRead2TrimPair[i]=${trimOut}/${mateShortNames[i]}_R2_trimPair.fastq 
    mateRead2TrimUnPair[i]=${trimOut}/${mateShortNames[i]}_R2_trimUnPair.fastq 
    
    java -jar ${progTrimmomatic} PE -threads ${NThreads} ${mateReads1[i]} ${mateReads2[i]} ${mateRead1TrimPair[i]} ${mateRead1TrimUnPair[i]} ${mateRead2TrimPair[i]} ${mateRead2TrimUnPair[i]} ILLUMINACLIP:${primerFileMP}:2:30:7:5:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40 >> $log
  done
  
  
  # quality check ----------
  echo "quality check of trimmed reads..."
  echo "quality check of trimmed reads..." >> $log
  
  for i in ${!lib[*]}  #for all indexes in the array
  do 
    ${progFastQC} -t ${NThreads} -o ${fastqcOut} ${read1TrimPair[i]} ${read2TrimPair[i]} ${read1TrimUnPair[i]} ${read2TrimUnPair[i]}
  done
  
  for i in ${!mateLib[*]}  #for all indexes in the array
  do 
    ${progFastQC} -t ${NThreads} -o ${fastqcOut} ${mateRead1TrimPair[i]} ${mateRead2TrimPair[i]} ${mateRead1TrimUnPair[i]} ${mateRead2TrimUnPair[i]}
  done
  
  

# 2. Step: map reads against reference 
#          and define blocks and superblocks
#######################################################
  # prepare Reference ---------
  echo "prepare reference..."
  echo "prepare reference..." >> $log
  
  #remove scaffolds shorter than 10 kb
  java -jar ${progRemovShortSeq} -i $ref -o ${refRed}.fa -length 10000 
  
  #create index files
  ${progSamtools} faidx ${refRed}.fa
  java -jar ${progPicard}/CreateSequenceDictionary.jar R=${refRed}.fa O=${refRed}.dict
  

  # map reads against reference ----------
  echo "run reference mapping..."
  echo "run reference mapping..." >> $log
  
  cd ${workPath}
  
  #merge unpaired files
  readTrimUnPair=${trimOut}/${name}_trimUnpair_mod.fastq
  cat ${read1TrimUnPair[*]} ${read2TrimUnPair[*]} > ${readTrimUnPair}
  libUnpair=Unpair
  
  # index reference file  
  echo "run reference mapping..."  
  ${progBowtie2}-build ${refRed}.fa ${refRed}
  
  mappedAll=()
  unmapped=()
  mapped=()
  mappedFiltered=()
  bowtieFailPair=()
  bowtieFailUnPair1=()
  bowtieFailUnPair2=()
  count=0
  for i in ${!lib[*]}  #for all indexes in the array
  do 
    mappedAll[i]=${workPath}/${shortNames[i]}_all.sorted.bam
    unmapped[i]=${workPath}/${shortNames[i]}_unmapped.sorted.bam
    mapped[i]=${workPath}/${shortNames[i]}.sorted.bam
    mappedFiltered[i]=${workPath}/${shortNames[i]}.sorted.filtered.bam
    bowtieFailPair[i]=${workPath}/${shortNames[i]}_failPair.fastq
    bowtieFailUnPair1[i]=${workPath}/${shortNames[i]}_failUnPairR1.fastq
    bowtieFailUnPair2[i]=${workPath}/${shortNames[i]}_failUnPairR2.fastq    
    (
      ${progBowtie2} --fast-local -p 1 -q --phred33 -I ${insLow[i]} -X ${insHigh[i]} -x ${refRed} -1 ${read1TrimPair[i]} -2 ${read2TrimPair[i]} -U ${read1TrimUnPair[i]},${read2TrimUnPair[i]} | ${progSamtools} view -bS - | ${progSamtools} sort - -T ${shortNames[i]} -o ${mappedAll[i]}
      ${progSamtools} index ${mappedAll[i]}
    
      #filter unmapped reads    
      ${progSamtools} view -b -F 4 ${mappedAll[i]} > ${mapped[i]}
      ${progSamtools} index ${mappped[i]}
    
      #get unmapped reads
      ${progSamtools} view -b -f 4 ${mappedAll[i]} > ${unmapped[i]}
      ${progSamtools} view -b -f 9 ${unmapped[i]} > ${unmapped[i]%.sorted.bam}_pair.sorted.bam
      java -jar ${progPicard}/SamToFastq.jar INPUT=${unmapped[i]%.sorted.bam}_pair.sorted.bam FASTQ=${bowtieFailPair[i]%.fastq}.1.fastq SECOND_END_FASTQ=${bowtieFailPair[i]%.fastq}.2.fastq
      rm ${unmapped[i]%.sorted.bam}_pair.sorted.bam
      ${progSamtools} view -b -F 8 -f 64 ${unmapped[i]} | ${progBamtools} convert -format fastq -out ${bowtieFailUnPair1[i]}  
      ${progSamtools} view -b -F 8 -f 128 ${unmapped[i]} | ${progBamtools} convert -format fastq -out ${bowtieFailUnPair2[i]}
    
      ${progBamtools} stats -in ${mappedAll[i]} >> $log
      echo "--> ${mappedAll[i]}" >> $log
    
      #filter for mapping quality >=10    
      ${progSamtools} view -b -q 10 ${mapped[i]} > ${mappedFiltered[i]}
      ${progBamtools} stats -in ${mappedFiltered[i]} >> $log
      echo "--> ${mappedFiltered[i]}" >> $log
    
      #check insertion size
      #java -jar ${progPicard}/CollectInsertSizeMetrics.jar R=${refRed}.fa I=${mapped[i]} O=${mapped[i]%.bam}_insertSize.txt HISTOGRAM_FILE=${mapped[i]%.bam}_insertSizeHist.pdf
    ) &
    let count+=1
    [[ $((count%${NThreads})) -eq 0 ]] && wait
  done
  wait
    
  mateMappedAll=()
  mateUnmapped=()
  mateMapped=()
  mateBowtieFailPair=()
  mateBowtieFailUnPair1=()
  mateBowtieFailUnPair2=()
  count=0
  for i in ${!mateLib[*]}  #for all indexes in the array
  do 
    mateMappedAll[i]=${workPath}/${mateShortNames[i]}_all.sorted.bam
    mateMapped[i]=${workPath}/${mateShortNames[i]}.sorted.bam
    mateUnmapped[i]=${workPath}/${mateShortNames[i]}_unmapped.sorted.bam
    mateBowtieFailPair[i]=${workPath}/${mateShortNames[i]}_failPair.fastq
    mateBowtieFailUnPair1[i]=${workPath}/${mateShortNames[i]}_failUnPairR1.fastq
    mateBowtieFailUnPair2[i]=${workPath}/${mateShortNames[i]}_failUnPairR2.fastq      
    ( 
      ${progBowtie2} --fast-local -p 1 -q --phred33 -I ${mateInsLow[i]} -X ${mateInsHigh[i]} -x ${refRed} -1 ${mateRead1TrimPair[i]} -2 ${mateRead2TrimPair[i]} --rf | ${progSamtools} view -bS - | ${progSamtools} sort - -T ${mateShortNames[i]} -o ${mateMappedAll[i]}
      ${progSamtools} index ${mateMappedAll[i]}
    
      #filter unmapped reads    
      ${progSamtools} view -b -F 4 ${mateMappedAll[i]} > ${mateMapped[i]}
      ${progSamtools} index ${mateMapped[i]}
    
      #get unmapped reads
      ${progSamtools} view -b -f 4 ${mateMappedAll[i]} > ${mateUnmapped[i]}
      ${progSamtools} view -b -f 9 ${mateUnmapped[i]} > ${mateUnmapped[i]%.sorted.bam}_pair.sorted.bam
      java -jar ${progPicard}/SamToFastq.jar INPUT=${mateUnmapped[i]%.sorted.bam}_pair.sorted.bam FASTQ=${mateBowtieFailPair[i]%.fastq}.1.fastq SECOND_END_FASTQ=${mateBowtieFailPair[i]%.fastq}.2.fastq
      rm ${mateUnmapped[i]%.sorted.bam}_pair.sorted.bam
      ${progSamtools} view -b -F 8 -f 64 ${mateUnmapped[i]} | ${progBamtools} convert -format fastq -out ${mateBowtieFailUnPair1[i]} 
      ${progSamtools} view -b -F 8 -f 128 ${mateUnmapped[i]} | ${progBamtools} convert -format fastq -out ${mateBowtieFailUnPair2[i]}
    
      ${progBamtools} stats -in ${mateMappedAll[i]} >> $log
      echo "--> ${mateMappedAll[i]}" >> $log
    
      #check insertion size
      #java -jar ${progPicard}/CollectInsertSizeMetrics.jar R=${refRed}.fa I=${mateMapped[i]} O=${mateMapped[i]%.bam}_insertSize.txt HISTOGRAM_FILE=${mateMapped[i]%.bam}_insertSizeHist.pdf
    ) &
    let count+=1
    [[ $((count%${NThreads})) -eq 0 ]] && wait
  done
  wait  
  
  #merge alignment files
  mappedMerged=${workPath}/${name}_mate.sorted_mapped
  ${progSamtools} merge ${mappedMerged}.bam ${mapped[*]} ${mateMapped[*]}
  
  
  # get blocks and superblocks ----------
  echo "get blocks and superblocks..."
  echo "get blocks and superblocks..." >> $log
  
  #get coverage along genome
  covFile=${workPath}/${name}_mate_coverage.txt
  ${progBedtools} genomecov -ibam ${mappedMerged}.bam -bga > ${covFile}
  #only proparly paired reads
  ${progSamtools} view -bf 0x2 ${mappedMerged}.bam | ${progSamtools} sort - -n | ${progBedtools} bamtobed -i - -bedpe | awk '$1 == $4' | cut -f 1,2,6 | sort -k 1,1 | ${progBedtools} genomecov -i - -bga -g ${refRed}.fa.fai > ${covFile%.txt}Paired.txt
  
  #get blocks with a minimal coverage of 10 (paired-end) reads and create superblocks of at least 12000bp length and a minimal overlap of 300bp (max. overlap = 3*300bp)
  blocks=${workPath}/blocks.txt
  superblocks=${workPath}/superblocks.txt
  java -jar ${progGetBlocks} -i ${covFile} -paired ${covFile%.txt}Paired.txt -o ${blocks} -oSuper ${superblocks} -mCov 10 -sLength 12000 -sOverlap 300 -maxLength 100000



# 3. Step: do deNovo assembly within superblocks
#######################################################
  echo "deNovo assembly within superblocks..."
  echo "deNovo assembly within superblocks..." >> $log
  
  cd ${workPath}
  idbaRes=${workPath}/idbaResults
  mkdir ${idbaRes}
  
  #split superbolcks.txt into ${NThreads} files to run it in parallel
  size=$(($(wc -l < ${superblocks})/$((NThreads))+1))
  split -l $size -d ${superblocks} ${superblocks%.txt}_
  
  count=0
  for (( j=0; j<$((NThreads)); j++ ))
  do
    file=${superblocks%.txt}_0${j}
    fileout=${superblocks%.txt}_0${j}_run.sh
    logout=${superblocks%.txt}_0${j}.log
    
    array=(${file//_/ })
    number=${array[${#array[*]}-1]}
    if [ ${number} != "00" ]
    then
      number=`echo $number|sed 's/^0*//'`
    fi
    
    if [[ $number =~ ^[0-9]+$ ]]
    then
      blockNb=$(($number*$size+1))
    else
      blockNb=1
    fi
    
    printf "#"'!'"/bin/bash\n"                 > ${fileout}
    printf "\n"                                >> ${fileout}
    printf "mkdir ${file}_temp\n"              >> ${fileout}
    printf "cd ${file}_temp\n"                 >> ${fileout}
    printf "\n" >> ${fileout}

    printf "blockNb=$blockNb\n" >> ${fileout}
    printf "start=\`date +%%s\`\n" >> ${fileout}
    printf "for block in \$(cat ${file})\n" >> ${fileout}
    printf "do\n" >> ${fileout}
    printf "  echo \$blockNb >> $logout\n" >> ${fileout}
    printf "\n" >> ${fileout}
    printf "  #extract sequence names within specified region\n" >> ${fileout}
    
    seqNames=()
    seqBam=()
    subSeq1=()
    subSeq2=()
    for i in ${!lib[*]}
    do
      seqNames[i]=sequences_${shortNames[i]}
      seqBam[i]=sequences_${shortNames[i]}.bam
      subSeq1[i]=subseq_${shortNames[i]}_R1.fastq
      subSeq2[i]=subseq_${shortNames[i]}_R2.fastq
      printf "  ${progSamtools} view -b ${mapped[i]} \$block | ${progSamtools} sort - -T ${shortNames[i]} -no ${seqBam[i]}\n" >> ${fileout}
      printf "  ${progBedtools} bamtofastq -i ${seqBam[i]} -fq 1_${subSeq1[i]} -fq2 1_${subSeq2[i]}\n" >> ${fileout}

      #extract paired reads with one pair unmapped
      printf "  ${progSamtools} view  -b -f 72 ${seqBam[i]} | ${progBamtools} convert -format fastq | paste - - - - | sort -k1,1 -t \" \" | tr \"\\\t\" \"\\\n\" > 2_${subSeq1[i]}\n" >> ${fileout}
      printf "  ${progSamtools} view  -f 72 ${seqBam[i]} | cut -f 1 | awk \'{print \$0\"/2\"}\' > ${seqNames[i]}_R2.txt\n" >> ${fileout}
      printf "  ${progSeqtk} subseq ${bowtieFailUnPair2[i]} ${seqNames[i]}_R2.txt | paste - - - - | sort -k1,1 -t \" \" | tr \"\\\t\" \"\\\n\" > 2_${subSeq2[i]}\n" >> ${fileout}
  
      printf "  ${progSamtools} view  -b -f 136 ${seqBam[i]} | ${progBamtools} convert -format fastq | paste - - - - | sort -k1,1 -t \" \" | tr \"\\\t\" \"\\\n\" > 3_${subSeq1[i]}\n" >> ${fileout}
      printf "  ${progSamtools} view  -f 136 ${seqBam[i]} | cut -f 1 | awk \'{print \$0\"/1\"}\' > ${seqNames[i]}_R1.txt\n" >> ${fileout}
      printf "  ${progSeqtk} subseq ${bowtieFailUnPair1[i]} ${seqNames[i]}_R1.txt | paste - - - - | sort -k1,1 -t \" \" | tr \"\\\t\" \"\\\n\" > 3_${subSeq2[i]}\n" >> ${fileout}
  
      printf "  cat 1_${subSeq1[i]} 2_${subSeq1[i]} 3_${subSeq1[i]} > ${subSeq1[i]}\n" >> ${fileout}
      printf "  cat 1_${subSeq2[i]} 2_${subSeq2[i]} 3_${subSeq2[i]} > ${subSeq2[i]}\n" >> ${fileout}
      printf "\n" >> ${fileout}
    done
    
    printf "  #idba --------\n" >> ${fileout}
    subSeqFasta=()
    for i in ${!lib[*]}
    do
      subSeqFasta[i]=subseq_${lib[i]}.fa
      printf "  ${progIdba}fq2fa --merge ${subSeq1[i]} ${subSeq2[i]} ${subSeqFasta[i]}\n" >> ${fileout}
    done
    subSeqFasta[${#subSeqFasta[*]}]=subseq_UnPair.fa
    printf "  ${progIdba}fq2fa ${subSeqUnPair} ${subSeqFasta[${#subSeqFasta[*]}-1]}\n" >> ${fileout}
    printf "  cat ${subSeqFasta[*]} > subseq.fa\n" >> ${fileout}
    printf "  ${progIdba}idba_ud -o idbaResults --read subseq.fa --num_threads 1 --pre_correction\n" >> ${fileout}
    
    printf "  if [ ! -f idbaResults/contig.fa ]\n" >> ${fileout}
    printf "  then\n" >> ${fileout}
    printf "    echo \"\$blockNb idba failed\" >> $logout\n" >> ${fileout}
    printf "  fi\n" >> ${fileout}
    printf "  cp idbaResults/contig.fa ${idbaRes}/sblock\${blockNb}-contigs.fa\n" >> ${fileout}
    printf "\n" >> ${fileout}    
           
    printf "  rm -rf ${file}_temp/*\n" >> ${fileout}
    printf "  ((blockNb++))\n" >> ${fileout}
    printf "done\n" >> ${fileout}
    printf "\n" >> ${fileout}   
    printf "rm -rf ${file}_temp\n" >> ${fileout}
    printf "end=\`date +%%s\`\n" >> ${fileout}
    printf "echo \$((end-start))\n" >> ${fileout}
    printf "\n" >> ${fileout}
    
    chmod +x ${fileout}
    ${fileout} &
    let count+=1
    [[ $((count%${NThreads})) -eq 0 ]] && wait
  done
  wait
  
    
  # deNovo assembly of unassembled reads ----------
  echo "deNovo assembly of unassembled reads..."
  echo "deNovo assembly of unassembled reads..." >> $log
    
  unassFolder=${workPath}/Unassembled
  mkdir ${unassFolder}
  cd ${unassFolder}
  
  #merge unassembled unpaired files
  bowtieFailUnpairMerged=${workPath}/${name}_failUnp.fastq
  cat ${bowtieFailUnpair1[*]} ${bowtieFailUnpair2[*]} > $bowtieFailUnpairMerged
  
  subSeqFasta=()
  for i in ${!lib[*]}
  do
    subSeqFasta[i]=${shortNames[i]}_failPair.fa
    ${progIdba}fq2fa --merge ${bowtieFailPair[i]%.fastq}.1.fastq ${bowtieFailPair[i]%.fastq}.2.fastq ${subSeqFasta[i]}
  done
  subSeqFasta[${#subSeqFasta[*]}]=${name}_failUnp.fa
  ${progIdba}fq2fa ${bowtieFailUnpairMerged} ${subSeqFasta[${#subSeqFasta[*]}-1]}
  cat ${subSeqFasta[*]} > subseq.fa
  ${progIdba}idba_ud -o idbaResults --read subseq.fa --num_threads ${NThreads} --pre_correction
  
  
  
# 4. Step: get non-redundant supercontigs
####################################################### 
  echo "get supercontigs..."
  echo "get supercontigs..." >> $log
  cd ${workPath}
  
  #merge deNovo assembled superblocks into one FASTA file
  idbaContigs=${idbaRes}/idbaContigs.fa 
  cat ${idbaRes}/*-contigs.fa > ${idbaContigs}
  
  #merge Unassembled files
  orgUnass=${unassFolder}/idba/idbaResults/contig.fa
  #remove short seq (<500)
  echo "remove seq < 500 in ${orgUnass}" >> $log
  unass500=${unassFolder}/Unassembled_Idba_500.fa
  java -jar ${progRemovShortSeq} -i ${orgUnass} -o ${unass500} -length 500  >> $log
  
  #merge all files
  superblockSeq=${workPath}/deNovo_Superblocks.fa
  cat ${idbaContigs} ${unass500} > ${superblockSeq}
    
  #remove short seq (<200)
  echo "remove seq < 200 in ${superblockSeq}" >> $log
  superblockSeq200=${superblockSeq%.fa}_200.fa
  java -jar ${progRemovShortSeq} -i ${superblockSeq} -o ${superblockSeq200} -length 200 -u  >> $log
 
  #remove redundency with AMOScmp
  amosFolder=${workPath}/AMOScmp
  mkdir ${amosFolder}
  cd ${amosFolder}
  
  #assemble all assembled superblocks with AMOScmp to supercontigs (with the help of reference)
  #changed parameters in AMOScmp: (casm-layout -t 1000 (maximum ignorable trim length), make-consensus -o 10 (minimum overlap base)) 
  superblockSeqAmos=${superblockSeq200%.fa}_Amos.afg 
  java -jar ${progFastaToAmos} -i ${superblockSeq200} -o ${superblockSeqAmos}
  supercontigs=Amos_supercontigs
  amosSupercontigs=${amosFolder}/${supercontigs}.fasta
  
  echo "run AMPScmp..." >> $log
  #${progAmos}AMOScmp -D TGT=${superblockSeqAmos} -D REF=${refRed}.fa ${supercontigs}
  
  # running AMPScmp step by step and use multithread nucmer to spead it up
  ## Building AMOS bank
  echo "  build AMPS bank..." >> $log
  ${progAmos}/bank-transact -c -z -b ${supercontigs}.bnk -m ${superblockSeqAmos}

  ## Collecting clear range sequences
  echo "  clear range sequences..." >> $log
  ${progAmos}/dumpreads ${supercontigs}.bnk > ${supercontigs}.seq

  ## Running nucmer
  echo "  run nucmer..." >> $log
  ${progNucmer} --maxmatch --threads=${NThreads} --prefix=${supercontigs} ${refRed}.fa ${supercontigs}.seq

  ## Running layout
  echo "  run layout..." >> $log
  ${progAmos}/casm-layout -t 1000 -U ${supercontigs}.layout -C ${supercontigs}.conflict -b ${supercontigs}.bnk ${supercontigs}.delta

  ## Running consensus
  echo "  run consensus..." >> $log
  ${progAmos}/make-consensus -o 10 -B -b ${supercontigs}.bnk

  ## Outputting contigs
  echo "  output contigs..." >> $log
  ${progAmos}/bank2contig ${supercontigs}.bnk > ${supercontigs}.contig

  ## Outputting fasta
  echo "  output fasta..." >> $log
  ${progAmos}/bank2fasta -b ${supercontigs}.bnk > ${supercontigs}.fasta
      
 
 
# 5. Step: map reads on supercontigs
#          and de novo assemble unmapped reads
####################################################### 
  echo "map reads on supercontigs and correct them..."
  echo "map reads on supercontigs and correct them..." >> $log
  
  #make seqnames unique
  amosSupercontigsUnique=${amosSupercontigs%.fasta}_unique.fa
  java -jar ${progRemovShortSeq} -i ${amosSupercontigs} -o ${amosSupercontigsUnique} -length 1 -u
   
  #get statistics
  echo ${amosSupercontigsUnique} >> $log
  java -jar ${progFastaStats} -i ${amosSupercontigsUnique} -min 200 >> $log
  
  #prepare reference
  ${progBowtie2}-build ${amosSupercontigsUnique} ${amosSupercontigsUnique%.fa}
  
  supercontMappedAll=()
  supercontUnmapped=()
  supercontFailPair=()
  supercontFailUnpair=()
  supercontMappedFiltered=()
  count=0
  for i in ${!lib[*]}  #for all indexes in the array
  do 
    supercontMappedAll[i]=${amosFolder}/${shortNames[i]}_all.sorted.bam
    supercontUnmapped[i]=${amosFolder}/${shortNames[i]}_unmapped.sorted.bam
    supercontFailPair[i]=${amosFolder}/${shortNames[i]}_failPair.fastq
    supercontFailUnpair[i]=${amosFolder}/${shortNames[i]}_failUnp.fastq
    supercontMappedFiltered[i]=${amosFolder}/${shortNames[i]}.filtered.sorted.bam
    (
      ${progBowtie2} --sensitive -p 1 -q --phred33 -I ${insLow[i]} -X ${insHigh[i]} -x ${amosSupercontigsUnique%.fa} -1 ${read1TrimPair[i]} -2 ${read2TrimPair[i]} | ${progSamtools} view -bS - | ${progSamtools} sort - -T ${shortNames[i]} -o ${supercontMappedAll[i]}
      ${progSamtools} index ${supercontMappedAll[i]}
    
      #get unmapped reads
      ${progSamtools} view -b -f 4 ${supercontMappedAll[i]} > ${supercontUnmapped[i]}
      ${progSamtools} view -b -f 9 ${supercontUnmapped[i]} > ${supercontUnmapped[i]%.sorted.bam}_pair.sorted.bam
      java -jar ${progPicard}/SamToFastq.jar INPUT=${supercontUnmapped[i]%.sorted.bam}_pair.sorted.bam FASTQ=${supercontFailPair[i]%.fastq}.1.fastq SECOND_END_FASTQ=${supercontFailPair[i]%.fastq}.2.fastq
      rm ${supercontUnmapped[i]%.sorted.bam}_pair.sorted.bam
      ${progSamtools} view -b -F 8 ${supercontUnmapped[i]} | ${progBamtools} convert -format fastq -out ${supercontFailUnpair[i]}
    
      ${progBamtools} stats -in ${supercontMappedAll[i]} >> $log
      echo "--> ${supercontMappedAll[i]}" >> $log
     
      #filter for mapping quality >=10    
      ${progSamtools} view -b -F 4 -q 10 ${supercontMappedAll[i]} > ${supercontMappedFiltered[i]}
      ${progBamtools} stats -in ${supercontMappedFiltered[i]} >> $log
      echo "--> ${supercontMappedFiltered[i]}" >> $log
    ) &
    let count+=1
    [[ $((count%${NThreads})) -eq 0 ]] && wait
  done
  wait
       

  # deNovo assemble unassembled reads ----------
  echo "deNovo assemble unassembled reads..."
  echo "deNovo assemble unassembled reads..." >> $log
  
  supercontFailUnpairMerged=${amosFolder}/${name}_failUnp.fastq
  cat ${supercontFailUnpair[*]} > ${supercontFailUnpairMerged}
  
  supercontUnassFolder=${amosFolder}/Unassembled
  mkdir ${supercontUnassFolder}
  cd ${supercontUnassFolder}
  
  subSeqFasta=()
  for i in ${!lib[*]}
  do
    subSeqFasta[i]=${shortNames[i]}_failPair.fa
    ${progIdba}fq2fa --merge ${supercontFailPair[i]%.fastq}.1.fastq ${supercontFailPair[i]%.fastq}.2.fastq ${subSeqFasta[i]}
  done
  subSeqFasta[${#subSeqFasta[*]}]=${name}_failUnp.fa
  ${progIdba}fq2fa ${supercontFailUnpairMerged} ${subSeqFasta[${#subSeqFasta[*]}-1]}
  cat ${subSeqFasta[*]} > subseq.fa
  ${progIdba}idba_ud -o idbaResults --read subseq.fa --num_threads ${NThreads} --pre_correction
  
  #remove contigs shorter than 200 bp
  supercontSeqUnass=${supercontUnassFolder}/Unass-contigs_200.fa
  java -jar ${progRemovShortSeq} -i idbaResults/contig.fa -o ${supercontSeqUnass} -length 100 >> $log
  echo ${supercontSeqUnass} >> $log
  java -jar ${progFastaStats} -i ${supercontSeqUnass} -min 200 >> $log
  
 
 
# 6. Step: map reads to all supercontics and correct them 
########################################################## 
  echo "merge contigs..." 
  echo "merge contigs..." >> $log
  
  mergedFolder=${workPath}/merged_corr
  mkdir $mergedFolder
  cd $mergedFolder
  
  merged=${mergedFolder}/${name}_supercontSeq_Unass.fa
  cat ${amosSupercontigsUnique} ${supercontSeqUnass} > $merged
  
  echo "${merged} >> $log
  java -jar ${progFastaStats} -i ${merged} -min 200 >> $log
  
    ${progBowtie2}-build ${merged} ${merged%.fa}
  
  mergedMappedAll=()
  mergedMappedFiltered=()
  count=0
  for i in ${!lib[*]}  #for all indexes in the array
  do 
    mergedMappedAll[i]=${mergedFolder}/${shortNames[i]}_all.sorted.bam
    mergedMappedFiltered[i]=${mergedFolder}/${shortNames[i]}.filtered.sorted.bam
    (
      ${progBowtie2} --sensitive -p 1 -q --phred33 -I ${insLow[i]} -X ${insHigh[i]} -x ${merged%.fa} -1 ${read1TrimPair[i]} -2 ${read2TrimPair[i]} | ${progSamtools} view -bS - | ${progSamtools} sort - -T ${shortNames[i]} -o ${mergedMappedAll[i]}
      ${progSamtools} index ${mergedMappedAll[i]}
        
      ${progBamtools} stats -in ${mergedMappedAll[i]} >> $log
      echo "--> ${mergedMappedAll[i]}" >> $log
     
      #filter for mapping quality >=10
      ${progSamtools} view -b -F 4 -q 10 ${mergedMappedAll[i]} > ${mergedMappedFiltered[i]}
    
      ${progBamtools} stats -in ${mergedMappedFiltered[i]} >> $log
      echo "--> ${mergedMappedFiltered[i]}" >> $log
    ) &
    let count+=1
    [[ $((count%${NThreads})) -eq 0 ]] && wait
  done
  wait
  
  
  # error correction ----------  
  # add RG header of second file (lost while merging)
  for i in ${!lib[*]}  #for all indexes in the array
  do
    echo -e "@RG\tID:${shortNames[i]}.filtered.sorted\tPL:illumina\tPU:${lib[i]}\tLB:${lib[i]}\tSM:${shortNames[i]}" >> rg
  done
  ${progSamtools} view -H ${mergedMappedFiltered[1]} | cat - rg > header
      
  #realign reads
  mergedMappedMerged=${mergedFolder}/${name}.filtered_RG.sorted.bam
  ${progSamtools} merge -r -h header ${mergedMappedMerged} ${mergedMappedFiltered[*]}
  rm rg header
 
  ${progSamtools} index ${mergedMappedMerged}
  ${progSamtools} faidx ${merged}
  java -jar ${progPicard}/CreateSequenceDictionary.jar R=${merged} O=${merged%.fa}.dict
  java -jar ${progGatk} -T RealignerTargetCreator -R ${merged} -I ${mergedMappedMerged} -o target_intervals.list
  mergedMappedMergedReal=${mergedMappedMerged%.bam}_realigned.bam
  java -jar ${progGatk} -T IndelRealigner -R ${merged} -I ${mergedMappedMerged} -targetIntervals target_intervals.list -o ${mergedMappedMergedReal}

  #get alternative seq
  mergedCorr=${merged%.fa}_corr.fq
  ${progSamtools} mpileup -uf ${merged} ${mergedMappedMergedReal} | ${progBcftools} call -c - | ${progVcfutils} vcf2fq -d 1 > ${mergedCorr}    
  
  #remove start and end N
  mergedCorrWN=${mergedCorr%.fq}WN.fa 
  echo ${mergedCorrWN} >> $log
  java -jar ${progRemovShortSeq} -i ${mergedCorr} -o ${mergedCorrWN} -length 100 -n -fq >> $log
  
  #get statistics
  echo ${mergedCorrWN} >> $log
  java -jar ${progFastaStats} -i ${mergedCorrWN} -min 200 >> $log
  

  # split sequences at places with no coverage ----------
  ${progBowtie2}-build ${mergedCorrWN} ${mergedCorrWN%.fa}
  
  mergedCorrMappedAll=()
  mergedCorrMappedFiltered=()
  count=0
  for i in ${!lib[*]}  #for all indexes in the array
  do 
    mergedCorrMappedAll[i]=${mergedFolder}/${shortNames[i]}_corrWN_all.sorted.bam
    mergedCorrMappedFiltered[i]=${mergedFolder}/${shortNames[i]}_corrWN.filtered.sorted.bam
    (
      ${progBowtie2} --sensitive -p 1 -q --phred33 -I ${insLow[i]} -X ${insHigh[i]} -x ${mergedCorrWN%.fa} -1 ${read1TrimPair[i]} -2 ${read2TrimPair[i]} | ${progSamtools} view -bS - | ${progSamtools} sort - -T ${shortNames[i]} -o ${mergedCorrMappedAll[i]}
      ${progSamtools} index ${mergedCorrMappedAll[i]}
   
      ${progBamtools} stats -in ${mergedCorrMappedAll[i]} >> $log
      echo "--> ${mergedCorrMappedAll[i]}" >> $log
    
      #filter for mapping quality >=10
      ${progSamtools} view -b -F 4 -q 10 ${mergedCorrMappedAll[i]} > ${mergedCorrMappedFiltered[i]}
    
      ${progBamtools} stats -in ${mergedCorrMappedFiltered[i]} >> $log
      echo "--> ${mergedCorrMappedFiltered[i]}" >> $log
    ) &
    let count+=1
    [[ $((count%${NThreads})) -eq 0 ]] && wait
  done
  wait

  mergedCorrMappedFilteredMerged=${mergedFolder}/${name}_corrWN.filtered.sorted.bam
  ${progSamtools} merge ${mergedCorrMappedFilteredMerged} ${mergedCorrMappedFiltered[*]}
  ${progBedtools} genomecov -ibam ${mergedCorrMappedFilteredMerged} -bga > ${mergedCorrWN%.fa}_filteredCov.txt
  #only proparly paired reads
  ${progSamtools} faidx $mergedCorrWN
  ${progSamtools} view -bf 0x2 ${mergedCorrMappedFilteredMerged} | ${progSamtools} sort - -n | ${progBedtools} bamtobed -i - -bedpe | awk '$1 == $4' | cut -f 1,2,6 | sort -k 1,1 | ${progBedtools} genomecov -i - -bga -g ${mergedCorrWN}.fai  > ${mergedCorrWN%.fa}_filteredPairedCov.txt
  
  java -jar ${progSplitSeqLowCov} -i ${mergedCorrWN%.fa}_filteredCov.txt -paired ${mergedCorrWN%.fa}_filteredPairedCov.txt -o ${mergedCorrWN%.fa}_filteredNotCov.txt -mCov 1 -fasta ${mergedCorrWN} -fastaOut ${mergedCorrWN%.fa}_splitFiltered.fa >> $log
  echo ${mergedCorrWN%.fa}_splitFiltered.fa >> $log
  java -jar ${progFastaStats} -i ${mergedCorrWN%.fa}_splitFiltered.fa -min 200 >> $log
  


# 7. Step: scaffolding and gap closing
#######################################################
  echo "scaffolding..." 
  echo "scaffolding..." >> $log
  
  scafFolder=${mergedFolder}/scaffold_gapClosed
  mkdir $scafFolder
  cd $scafFolder
  
  for i in ${!lib[*]}
  do
    if [ $i == 0 ]
    then
      libList=${lib[i]}
      forwardReads=${read1TrimPair[i]}
      reverseReads=${read2TrimPair[i]}
    else
      libList=${libList},${lib[i]}
      forwardReads=${forwardReads},${read1TrimPair[i]}
      reverseReads=${reverseReads},${read2TrimPair[i]}
    fi      
  done
  
  for i in ${!mateLib[*]}
  do
    if [ $i == 0 ]
    then
      mateLibList=${mateLib[i]}
      mateForwardReads=${mateRead1TrimPair[i]}
      mateReverseReads=${mateRead2TrimPair[i]}
    else
      mateLibList=${mateLibList},${mateLib[i]}
      mateForwardReads=${mateForwardReads},${mateRead1TrimPair[i]}
      mateReverseReads=${mateReverseReads},${mateRead2TrimPair[i]}
    fi      
  done
  
  #write config file
  soapConf=${scafFolder}/soap.config
  java -jar ${progWriteSoapConfig} -insLength ${libList} -r1 ${forwardReads} -r2 ${reverseReads} -max ${maxReadLength} -ru 2 -mateInsLength ${mateLibList} -mateR1 ${mateForwardReads} -mateR2 ${mateReverseReads} -mateRu 2 -rank -o ${soapConf}
  scafFile=${name}_${kmer}
  ${progSoapdenovo2}/prepare/finalFusion -D -c ${mergedCorrWN%.fa}_splitFiltered.fa -K ${kmer} -g ${scafFile} -p ${NThreads}
  ${progSoapdenovo2}/SOAPdenovo-127mer map -s ${soapConf} -g ${scafFile} -p ${NThreads}
  ${progSoapdenovo2}/SOAPdenovo-127mer scaff -g ${scafFile} -p ${NThreads} -F
   
   
  #remove scaffolds < 200 bp ----------
  scafSeq=${scafFolder}/${name}_scafSeq.fa
  echo ${scafFile}.scafSeq >> $log
  java -jar ${progRemovShortSeq} -i ${scafFile}.scafSeq -o ${scafSeq} -length 200 >> $log
  java -jar ${progRemovShortSeq} -i ${scafSeq} -o ${scafSeq%.fa}_500.fa -length 500 >> $log
  java -jar ${progRemovShortSeq} -i ${scafSeq} -o ${scafSeq%.fa}_1000.fa -length 1000 >> $log
 
  #get statistics
  echo ${scafSeq} >> $log
  java -jar ${progFastaStats} -i ${scafSeq} -min 200 >> $log
  java -jar ${progFastaStats} -i ${scafSeq} -min 500 >> $log
  java -jar ${progFastaStats} -i ${scafSeq} -min 1000 >> $log

  #map reads against scaffolds
  ${progBowtie2}-build ${scafSeq} ${scafSeq%.fa}
   
  scafMappedAll=()
  scafMappedFiltered=()
  count=0
  for i in ${!lib[*]}  #for all indexes in the array
  do 
    scafMappedAll[i]=${scafFolder}/${shortNames[i]}_all.sorted.bam
    scafMappedFiltered[i]=${scafFolder}/${shortNames[i]}.filtered.sorted.bam
    (
      ${progBowtie2} --sensitive -p 1 -q --phred33 -I ${insLow[i]} -X ${insHigh[i]} -x ${scafSeq%.fa} -1 ${read1TrimPair[i]} -2 ${read2TrimPair[i]} | ${progSamtools} view -bS - | ${progSamtools} sort - -T ${shortNames[i]} -o ${scafMappedAll[i]}
      ${progSamtools} index ${scafMappedAll[i]}
    
      ${progBamtools} stats -in ${scafMappedAll[i]} >> $log
      echo "--> ${scafMappedAll[i]}" >> $log
    
      #filter for mapping quality >=10    
      ${progSamtools} view -b -F 4 -q 10 ${scafMappedAll[i]} > ${scafMappedFiltered[i]}
    
      ${progBamtools} stats -in ${scafMappedFiltered[i]} >> $log
      echo "--> ${scafMappedFiltered[i]}" >> $log
    ) &
    let count+=1
    [[ $((count%${NThreads})) -eq 0 ]] && wait
  done
  wait
  