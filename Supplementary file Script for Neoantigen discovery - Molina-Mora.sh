# In silico pipeline to identify tumor-specific antigens for cancer immunotherapy using exome sequencing data
Pipeline developed by Ph.D. Jose A. Molina Mora

# Samples: exome sequencing data. Here we provide the pipeline for a one tumor sample and one normal sample (somatic variant calling approach).
# The working directory must contain samples (sequencing data in fastq format, usually zipped .gz), and the Human reference genome (fasta format). 

# STEP 1. Quality control: visualization and trimming (removal of adapters and low quality bases)

# visualize quality before trimming
fastqc *.fastq*

# trimming step: run for each sample!
trimmomatic PE -threads 4 reads-R1.fastq.gz reads-R2.fastq.gz -baseout reads-trim.fastq.gz SLIDINGWINDOW:4:30

# visualize quality after trimmming
fastqc *P.fastq*


# STEP 2. Variant Calling analysis

# Index the reference genome
bwa index human.fasta

# Map the read using a reference genome: repeat for the tumor sample. 
bwa mem -t 4 human.fasta reads-normal-trim_1P.fastq reads-normal-trim_2P.fastq | samtools sort > normal.bam
bwa mem -t 4 human.fasta reads-tumor-trim_1P.fastq reads-tumor-trim_2P.fastq | samtools sort > tumor.bam

# Quality control of alignment
qualimap bamqc -bam normal.bam -outfile result_qualimap.pdf

# Variant calling: Prepare mpileup files (required by VarScan)
samtools mpileup -q 20 -Q 25 -B -d 1000 -f human.fasta normal.bam > normal.mpileup
samtools mpileup -q 20 -Q 25 -B -d 1000 -f human.fasta tumor.bam > tumor.mpileup

# Run variant calling with all the steps including filtering
java -jar VarScan.jar somatic normal.mpileup tumor.mpileup --min-coverage 8 --min-coverage-normal 8 --min-coverage-tumor 6 --min-var-freq 0.10 --output-vcf
java -jar VarScan.jar somaticFilter VCFfile.snp.vcf --indel-file VCFfile.indel.vcf --output-file VCFfile.filtered1.vcf

# Select somatic mutations with desired filters: using vcffilter within VCFlib package
vcffilt -f "QUAL > 20" VCFfile.filtered1.vcf > VCFfile.filtered2.vcf

# Annotations and final filtering
snpEff.jar ann -v -stats HTML_human.html GRCh37 VCFfile.filtered2.vcf  | SnpSift.jar filter "! exists ID"> VCFfile.annotated1.vcf 
table_annovar.pl VCFfile.annotated1.vcf humandb/ -buildver hg19 -out VCFfile.final.vcf -vcfinput


# STEP 3. Neoantigens predictions: due to a substantial reduction in data is expected (only tumor-exclusive mutations are reported in the VCFfile.final.vcf), the subsequent analyses can be easily run in Webservers. 

# Run MuPeXi and NetMHCpan (within MuPeXi): two HLA alleles are shown as example (a config.ini file needs to be prepared according to MuPeXi or run in the webserver)
MuPeXI.py -v VCFfile.final.vcf -a HLA-A01:01,HLA-B08:01 -c path_to_config.ini --netmhc-full-anal

# The candidate peptides were analyzed using customs databases. See maintext for details.






