#!/bin/bash

## MATE - Genetics an Biotechnology Institute - Gödöllő - Hungary

#Pipeline for RNA-seq alignemnts and variant calling in a family model (many samples)
#Memory and thread usage is configured for 10 samples running parallel, Adjust for your run (Careful!!)
#Progressing from fastq to combined

## It uses STAR for alignemnt and GATK for variants handling . Genome should be indexed prior to running this pipeline, gtf option is added

## Adjust the following parameters for the paths of your files

genome="OryCun_3.0.fasta"
knownsites="hoax.vcf.gz" ## mandatory file, could be empty
path="./"
vault="" #location of fastq files to process
Mem = 50000000000 ## Change this according to your server/cluster


date
echo "files used:"
ls $genome
ls $knownsites


if [ ! -d ${path}/output ]; then
    mkdir ${path}/output
fi

if [ ! -d ${path}/bam ]; then
    mkdir ${path}/bam
fi

if [ ! -d ${path}/log ]; then
    mkdir ${path}/log
fi

if [ ! -d ${path}/vcf ]; then
    mkdir ${path}/vcf
fi

#Make a list of IDs to process
IDlist=$(
    while read ID; do
        echo $ID
    done < ${path}/sample.list
)

for ID in ${IDlist}; do



          if [ ! -d ${path}/bam/${ID} ]; then
                mkdir ${path}/bam/${ID}
                fi


        if [ ! -f ${path}/bam/${ID}/bamAligned.sortedByCoord.out.bam ]; then

cd ${path}/${ID}
        STAR --genomeDir index/Orycun3.0 --readFilesCommand cat --readFilesIn ${ID}.sra_1.fastq ${ID}.sra_2.fastq --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM ${Mem} --outSAMunmapped None --twopassMode Basic --outFilterMultimapNmax 1 --outSAMreadID Standard --runThreadN 24 --quantMode GeneCounts --outFileNamePrefix ${path}/bam/${ID}/bam > ${path}/log/STAR_align_oryc3_${ID}.log 2> ${path}/log/STAR_align_oryc3_${ID}.err  #STAR Aligner is running
        fi

cd ${path}
done

for ID in ${IDlist}; do
    if [ ! -f ${path}/bam/${ID}/bamAligned.sortedByCoord.out.bai ]; then
        samtools index ${path}/bam/${ID}/bamAligned.sortedByCoord.out.bam 2> ${path}/log/${ID}.index.log &
    else
        echo "skip indexing bam file for ${ID} - file already exists"
    fi
done
wait


#Mark duplicates
for ID in ${IDlist}; do
    if [ ! -f ${path}/output/${ID}.dedup.bam ]; then
	java -Xmx300G -jar /molbio/bin/gatk-package-4.1.8.1-local.jar MarkDuplicates -I ${path}/bam/${ID}//bamAligned.sortedByCoord.out.bam -M ${path}/output/${ID}.dupmetrics.txt -O ${path}/output/${ID}.dedup.bam 2> ${path}/log/${ID}.dedup.log
    else
	echo "skip deduplication for ${ID} - file already exists"
    fi
done
wait
## SplitNCigar

for ID in ${IDlist}; do 
         if [ ! -f ${path}/output/${ID}.SNC.bam ]; then
	java -Xmx300G -jar /molbio/bin/gatk-package-4.1.8.1-local.jar SplitNCigarReads -R $genome -I ${path}/output/${ID}.dedup.bam -O ${path}/output/${ID}.SNC.bam 2> ${path}/log/${ID}.SplitNCi.bam.log
         fi
done

#BQSR and apply
for ID in ${IDlist}; do
    if [ ! -f ${path}/output/${ID}_BQSR.txt ]; then
	java -Xmx300G -jar /molbio/bin/gatk-package-4.1.8.1-local.jar BaseRecalibrator -I ${path}/output/${ID}.SNC.bam --known-sites $knownsites -O ${path}/output/${ID}_BQSR.txt -R $genome 2> ${path}/log/${ID}.bqsr.bam.log
    else
	echo "skip BQSR for ${ID} - file already exists"
    fi
done
wait

for ID in ${IDlist}; do
    if [ ! -f ${path}/output/${ID}.bqsr.bam ]; then
	java -Xmx400G -jar /molbio/bin/gatk-package-4.1.8.1-local.jar ApplyBQSR -bqsr ${path}/output/${ID}_BQSR.txt -I ${path}/output/${ID}.dedup.bam -O ${path}/output/${ID}.bqsr.bam 2> ${path}/log/${ID}.bqsr.log &
    else
	echo "skip applyBQSR for ${ID} - file already exists"
    fi
done
wait

#HaplotypeCaller
for ID in ${IDlist}; do
    if [ ! -f ${path}/vcf/${ID}.gvcf ];
    then
	java -Xmx40G -jar /molbio/bin/gatk-package-4.3.0.0-local.jar HaplotypeCaller -ERC GVCF -I ${path}/output/${ID}.bqsr.bam -O ${path}/vcf/${ID}.gvcf -R $genome  2>> ${path}/log/${ID}.gvcf.log &
    else
	echo "skip variantcall ${ID} - file already exists"
    fi
done
wait

#Combining and filtering
#Create working directory
if [ ! -d ${path}/combined ]; then
    mkdir ${path}/combined
fi

#CombineGVCFs
for ID in ${IDlist}; do
    echo --variant ${path}/vcf/${ID}.gvcf
done > tmp.txt
cat tmp.txt | tr "\n" " " > combine.list
rm tmp.txt

if [ ! -f ${path}/combined/allsamples.combined.gvcf ];
then
    java -Xmx300G -jar /molbio/bin/gatk-package-4.1.8.1-local.jar CombineGVCFs -R $genome --arguments_file combine.list -O ${path}/combined/allsamples.combined.gvcf.gz 2> ${path}/log/allsamples.combine.log
    tabix -p vcf ${path}/combined/allsamples.combined.gvcf.gz
else
    echo "allsamples.combined.gvcf already exists! Remove it if you are sure about running CombineGVCFs again"
fi

#GenotypeGVCFs
if [ ! -f ${path}/combined/allsamples.genotyped.vcf ];
then 
    java -Xmx300G -jar /molbio/bin/gatk-package-4.1.8.1-local.jar GenotypeGVCFs -R $genome -V ${path}/combined/allsamples.combined.gvcf.gz -O ${path}/combined/allsamples.genotyped.vcf.gz 2> ${path}/log/allsamples.genotype.log
    tabix -p vcf ${path}/combined/allsamples.genotyped.vcf.gz
else
    echo "allsamples.genotyped.vcf already exists! Remove it if you are sure about running GenotypeGVCFs again"
fi

#VariantFiltration - preparing and filtering
#Select Variants
#SNP
if [ ! -f ${path}/combined/snp.allsamples.vcf.gz ]; then
    java -jar -Xmx300G /molbio/bin/gatk-package-4.1.8.1-local.jar SelectVariants -V ${path}/combined/allsamples.genotyped.vcf.gz -select-type SNP -O ${path}/combined/snp.allsamples.vcf.gz 2> ${path}/log/select-snp.log
else
    echo "${path}/combined/snp.allsamples.vcf.gz already exists! Remove it if you are sure about running SelectVariants again"
fi
#INDEL
if [ ! -f ${path}/combined/indel.allsamples.vcf.gz ]; then
    java -jar -Xmx300G /molbio/bin/gatk-package-4.1.8.1-local.jar SelectVariants -V ${path}/combined/allsamples.genotyped.vcf.gz -select-type INDEL -select-type MIXED -O ${path}/combined/indel.allsamples.vcf.gz 2> ${path}/log/select-indel.log
else
    echo "${path}/combined/indel.allsamples.vcf.gz already exists! Remove it if you are sure about running SelectVariants again"
fi
#Filter
#SNP
if [ ! -f ${path}/combined/snp.filtered.vcf.gz ]; then
    java -jar -Xmx300G /molbio/bin/gatk-package-4.1.8.1-local.jar VariantFiltration -V ${path}/combined/snp.allsamples.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ${path}/combined/snp.filtered.vcf.gz 2> ${path}/log/filter-snp.log
else
    echo "${path}/combined/snp.filtered.vcf.gz already exists! Remove it if you are sure about running VariantFiltration again"
fi

#INDEL
if [ ! -f ${path}/combined/indel.filtered.vcf.gz ]; then
    java -jar -Xmx300G /molbio/bin/gatk-package-4.1.8.1-local.jar VariantFiltration -V ${path}/combined/indel.allsamples.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" -O ${path}/combined/indel.filtered.vcf.gz 2> ${path}/log/filter-indel.log
else
    echo "${path}/combined/indel.filtered.vcf.gz already exists! Remove it if you are sure about running VariantFiltration again"
fi
#Merge filtered
if [ ! -f  ${path}/combined/allsamples.merged.filtered.vcf.gz ]; then
    java -jar -Xmx300G /molbio/bin/gatk-package-4.1.8.1-local.jar MergeVcfs -I ${path}/combined/snp.filtered.vcf.gz -I ${path}/combined/indel.filtered.vcf.gz -O ${path}/combined/allsamples.merged.filtered.vcf.gz 2> ${path}/log/merge-filtered.log
else
    echo " ${path}/combined/allsamples.merged.filtered.vcf.gz already exists! Remove it if you are sure about running MergeVcfs again"
fi
#Select variants that have passed filters
if [ ! -f  ${path}/combined/allsamples.merged.filtered.pass.vcf.gz ]; then
    bcftools view ${path}/combined/allsamples.merged.filtered.vcf.gz -f PASS -o ${path}/combined/allsamples.merged.filtered.pass.vcf
    bgzip allsamples.merged.filtered.pass.vcf
    tabix -p vcf allsamples.merged.filtered.pass.vcf.gz
else
    echo "Final file already exists"
fi
echo "Processing fq files into combined, finished vcf finished"
echo "Please make sure all your output files are in order before proceeding further"

