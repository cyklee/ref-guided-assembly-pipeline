#########################################################
#  README -
#  Reference-guided de novo assembly
# ====================================================
# by Heidi Lischer (heidi.lischer@ieu.uzh.ch), 2015/15
#########################################################

We adapted and extended the reference-guided assembly approach from Schneeberger et al. [1]. 
The main idea of this approach is to first map reads against a reference genome of a 
related species to reduce the complexity of de novo assembly within continuous covered 
regions. In a further step, reads with no similarity to the related genome are integrated. 

Reference-guided de novo assembly pipeline:
1. Step: quality/adapter trimming and quality check
2. Step: map reads against reference and define blocks and superblocks
3. Step: do deNovo assembly within superblocks
4. Step: get non-redundant supercontigs
5. Step: map reads on supercontigs and de novo assemble unmapped reads
6. Step: map reads to all supercontics and correct them 
7. Step: scaffolding and gap closing


[1] http://www.pnas.org/content/early/2011/06/01/1107739108.full.pdf+html?with-ds=yes


 Dependencies
--------------
(in brackets: used version of softwares)

Third party programs:
- fastqc (v 0.10.1): Fastq quality checking
- trimmomatic-0.32 (v 0.32): quality and adapter trimming
- samtools (v 1.3): Tools for alignments in the SAM format
- bcftools (v 1.3): Tools for variant calling and manipulating VCFs and BCFs
- bamtools (v 2.3.0): Tools for alignments in the BAM format
- bedtools (v 2.19.0): Tool for genomic arithmetics
- picardtools (v 1.109): Processing alignment files
- bowtie2 (v 2.2.1): NGS alignment tool
- seqtk (v 1.0-r45): Fast Fasta/Fastq manipulation tool
- AMOScmp (v 3.1.0): Comprarative genome assembly
- MUMmer (v 3.23): Comprarative genome assembly
- GenomeAnalysisTK (v 3.1-1): Genome analysis toolkit
- SOAPdenovo2 (v r240): De novo genome assembler and scaffolder

and one of these de novo assemblers:
- AllPaths-LG (v 51279)
- idba (v 1.1.1)
- abyss-pe (v 1.5.2)
- SOAPdenovo2 (v r240)


 Additional in-house scripts
-----------------------------
- RemoveShortSeq.jar: Remove short sequences and make unique identifiers from FASTA/FASTQ 
- GetBlocks_new.jar: Create blocks and superblocks
- FastaToAmos.jar: Transform FASTA to Amos format
- WriteSoapConfig.jar: writes a SOAP config file
- FastaStats.jar: outputs FASTA statistics
- SplitSeqLowCov.jar: splits sequences at low coverage


 Running
---------
In first lines of the main scripts you have to adapt the paths to your system.
(Everything between: 
 # set variables #########################################
 ...
 #########################################################)


If the parameters are set you can run the pipeline as follows:
bash refGuidedDeNovoAssembly_ALLPATHS.sh
or
bash refGuidedDeNovoAssembly_IDBA.sh
or
bash refGuidedDeNovoAssembly_ABYSS.sh
or
bash refGuidedDeNovoAssembly_SOAP.sh
