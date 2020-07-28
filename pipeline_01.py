#Raw data download
/BiO/Install/sratoolkit.2.9.6-ubuntu64/bin/ fastq-dump --gzip --split-3 SRR1002940

#Short reads QC
/BiO/Install/FastQC_0.10.1/fastqc -t 4 SRR1002940.r1.temp.fq
/BiO/Install/FastQC_0.10.1/fastqc -t 4 SRR1002940.r2.temp.fq

#Short reads trimming
java -jar /BiO/Install/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 -phred33 
SRR1002940.r1.temp.fq SRR1002940.r2.temp.fq SRR1002940.r1.trim.fq SRR1002940.r1.unpair.fq SRR1002940.r2.trim.fq SRR1002940.r2.unpair.fq 
ILLUMINACLIP:/BiO/Install/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:151:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#Mappint to REF gene (BWA mem)
bwa mem -t 4 -R '@RG\tPL:Illumina\tID:YUHL\tSM:SRR1002940\tLB:HiSeq' /BiO/Education/WGS/REF/hg19.fa SRR1002940.r1.trim.fq SRR1002940.r2.trim.fq > SRR1002940.sam

#Mark Duplicate(PICARD)
##process one
java -jar /BiO/Install/picard-tools-2.22.3/picard.jar AddOrReplaceReadGroups TMP_DIR=TEMP_PICARD VALIDATION_STRINGENCY=LENIENT SO=coordinate I= SRR1002940.sam O= SRR1002940_sorted.bam RGID=SRR1002940 RGLB=HiSeq RGPL=Illumona RGPU=unit1 RGSM=SRR1002940 CREATE_INDEX=true

##process two
java -jar /BiO/Install/picard-tools-2.22.3/picard.jar MarkDuplicates TMP_DIR=TEMP_PICARD VALIDATION_STRINGENCY=LENIENT I= SRR1002940_sorted.bam O= SRR1002940_dedup.sam M= SRR1002940.duplicate_metrics REMOVE_DUPLICATES=true AS=true

##process three
java -jar /BiO/Install/picard-tools-2.22.3/picard.jar SortSam TMP_DIR=TEMP_PICARD VALIDATION_STRINGENCY=LENIENT SO=coordinate I= SRR1002940_dedup.sam O=SRR1002940_dedup.bam CREATE_INDEX=true

#Base Quality Score Recalibration (GATK)
##Analyze pattern
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar BaseRecalibrator -R /BiO/Education/WGS/REF/hg19.fa -I SRR1002940_dedup.bam --known-sites /BiO/Education/WGS/REF/dbsnp_138.hg19.vcf --known-sites /BiO/Education/WGS/REF/1000GENOMES-phase_3_indel.vcf -O SRR1002940_recal_pass1.table

##Recalibration
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar ApplyBQSR -R /BiO/Education/WGS/REF/hg19.fa -I SRR1002940_dedup.bam --bqsr-recal-file SRR1002940_recal_pass1.table -O SRR1002940_recal_pass1.bam

##Do it again
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar BaseRecalibrator -R /BiO/Education/WGS/REF/hg19.fa -I SRR1002940_recal_pass1.bam --known-sites /BiO/Education/WGS/REF/dbsnp_138.hg19.vcf --known-sites /BiO/Education/WGS/REF/1000GENOMES-phase_3_indel.vcf -O SRR1002940_recal_pass2.table

##Do it again
java -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar ApplyBQSR -R /BiO/Education/WGS/REF/hg19.fa -I SRR1002940_recal_pass1.bam -bqsr SRR1002940_recal_pass2.table -O SRR1002940_recal_pass2.bam

#Calling Variants - HaplotypeCaller
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar HaplotypeCaller -R /BiO/Education/WGS/REF/hg19.fa -I SRR1002940_recal_pass2.bam -O SRR1002940.rawVariants.g.vcf -ERC GVCF --standard-min-confidence-threshold-for-calling 20

#Applying Genotype - GenotypeGVCFs
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar GenotypeGVCFs -R /BiO/Education/WGS/REF/hg19.fa -V SRR1002940.rawVariants.g.vcf -O SRR1002940_genotype.vcf

#Extracting Variants (SNP.INDEL) - SelectVariants
##SNP
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar SelectVariants -R /BiO/Education/WGS/REF/hg19.fa -V SRR1002940_genotype.vcf --select-type-to-include SNP -O SRR1002940.rawSNPs.vcf

##INDEL
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar SelectVariants -R /BiO/Education/WGS/REF/hg19.fa -V SRR1002940_genotype.vcf --select-type-to-include INDEL -O SRR1002940.rawINDELs.vcf

#Filtration - VariantFiltration
##SNP
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar VariantFiltration -R /BiO/Education/WGS/REF/hg19.fa -V SRR1002940.rawSNPs.vcf -O SRR1002940.rawSNPs.filtered.vcf --filter-name "." --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || haplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0"

##INDEL
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar VariantFiltration -R /BiO/Education/WGS/REF/hg19.fa -V SRR1002940.rawINDELs.vcf -O SRR1002940.rawINDELs.filtered.vcf --filter-name "." --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"

#Merge Variants Information 
##SortVcf - gatk
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar SortVcf -I SRR1002940.rawSNPs.filtered.vcf -I SRR1002940.rawINDELs.filtered.vcf -O SRR1002940.Filtered_gatk.variant.vcf

##MergeVcfs - picard
java -Xmx8g -jar /BiO/Install/picard-tools-2.22.3/picard.jar MergeVcfs I= SRR1002940.rawSNPs.filtered.vcf I= SRR1002940.rawINDELs.filtered.vcf O= SRR1002940.Filtered_picard.variant.vcf

#Annotation - Annovar
egrep "^#|PASS" SRR1002940.Filtered_gatk.variant.vcf > SRR1002940.Filtered.variant.PASS.vcf


perl /BiO/Install/annovar/table_annovar.pl SRR1002940.Filtered.variant.PASS.vcf /BiO/Education/WGS/humandb/ -buildver hg19 -out SRR1002940 -remove -protocol refGene,cytoBand,avsnp138,clinvar_20190305 -operation g,r,f,f -nastring . -vcfinput








