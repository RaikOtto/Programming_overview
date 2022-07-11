#!/bin/bash

source activate first-steps
##sbatch --job-name=MASTER_variant_caller --output=/fast/users/ottor_c/work/MAPTor_NET/logs/%x-%j.MAaster.Variants.log --ntasks=4 --nodes=1 --mem=64G  --partition=long --time=1440:00 /fast/users/ottor_c/work/MAPTor_NET/Scripts/Align_Variant_Calling_RNA_seq.sh

#find ./* -type f -name '*.sra' -print0 | xargs -0 -Imysongs mv -i mysongs .

PREFIX=$HOME"/work/"
KALLISTO_INDEX=$PREFIX"/Data/Homo_sapiens.GRCh38.kallisto.idx"
REF_GEN_V="hg38"
INPUT_FOLDER="/data/projects/Sers-DKTKMasterTranscriptomics/work/2021-03-30_download/download/"
TMP_FOLDER="/fast/users/ottor_c/scratch/MAPTOR/"
OUTPUT_FOLDER=$PREFIX"/MAPTor_NET/Results/VCFs/Master/"

STAR_INDICES_FOLDER="/fast/users/ottor_c/work/Data/STAR"
REF_GEN_DIR="/fast/users/ottor_c/work/Data/GRCh38"
REF_GEN=$REF_GEN_DIR"/GRCh38.primary_assembly.genome.fa"
PICARD_CALL="picard "
GATK="gatk "
SJDBG_FILE="/fast/users/ottor_c/work/Data/GRCh38/gencode.v29.annotation.gtf"
DBSNP_PATH="/fast/projects/cubit/18.12/static_data/db/dbSNP/b150/GRCh38/All_20170710.vcf.gz"
KNOWN_INDELS="/fast/users/ottor_c/work/Data/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz"
KNOWN_INDELS_MILLS="/fast/users/ottor_c/work/Data/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

STUDY="Master"

NR_CPUS=4

if [ ! -d $OUTPUT_FOLDER ]; then
	mkdir $OUTPUT_FOLDER
fi

#FWD_ENDING="_1.fq.gz"
#FWD_ENDING="_R1.fastq.gz"
FWD_ENDING="_R1_complete_filtered.fastq.gz"
BKW_ENDING=${FWD_ENDING//"1"/"2"}

FASTQ_FILES=(`find $INPUT_FOLDER -type f | grep $FWD_ENDING'$'`)

function align_reads {

	FWD_READS=$1
	BKW_READS=${FWD_READS//$FWD_ENDING/$BKW_ENDING}

    SAMPLE_INPUT_FOLDER=${FWD_READS_UNZIP%'/'*}'/'
    SAMPLE_NAME_FOLDER=${FWD_READS//'_R1.fastq.gz'/""}
    SAMPLE_NAME_FOLDER=(${SAMPLE_NAME_FOLDER//\// })
    SAMPLE_NAME_FOLDER=${SAMPLE_NAME_FOLDER[-1]}
    SAMPLE_TMP_FOLDER=$TMP_FOLDER"/"$SAMPLE_NAME_FOLDER"/"
    OUTPUT_PATH_SAMPLE=$OUTPUT_FOLDER"/"$SAMPLE_NAME_FOLDER"/"
    
    FWD_READS_UNZIP=$SAMPLE_TMP_FOLDER"/"$SAMPLE_NAME_FOLDER"_R1.fastq"
    BKW_READS_UNZIP=$SAMPLE_TMP_FOLDER"/"$SAMPLE_NAME_FOLDER"_R2.fastq"
    FWD_READS_FILT=${FWD_READS_UNZIP//'_R1.fastq'/"_R1_val_1.fq"}
    BKW_READS_FILT=${BKW_READS_UNZIP//'_R2.fastq'/"_R2_val_2.fq"}
    
    BAM_ALIGNED=$SAMPLE_TMP_FOLDER"/"$SAMPLE_NAME_FOLDER".Aligned.out.bam"
    SA_INDEX=$SAMPLE_TMP_FOLDER"/"$SAMPLE_NAME_FOLDER".SAindex"
    BAM_PICARD=$SAMPLE_TMP_FOLDER"/"$SAMPLE_NAME_FOLDER".picard.bam"
    BAM_PICARD_SORTED=$SAMPLE_TMP_FOLDER"/"$SAMPLE_NAME_FOLDER".picard.sorted.bam"
    BAM_MARKED_DUPLICATES=${BAM_ALIGNED/".bam"/".marked_dups.bam"}
    BAM_MARKED_DUPLICATES_TXT=${BAM_ALIGNED/".bam"/".marked_dups.txt"}
    BAM_REALIGNED_INTERVALS=${BAM_ALIGNED/.bam/.realigned_intervals.interval_list}
    BAM_REALIGNED=${BAM_ALIGNED/.bam/.realigned.bam}
    BAM_RECALIBRATED=${BAM_ALIGNED/.bam/.recalibrated.bam}
    BAM_RECALLED=${BAM_ALIGNED/.bam/.recalled.bam}
    BAM_RECALLED_RECALIBRATED_LIST=${BAM_ALIGNED/.bam/.recalled.recallibrated.list}
    BAM_SPLIT_BAM=${BAM_ALIGNED/.bam/.split.bam}
    BAM_BQSR_RECAL_TABLE=${BAM_SPLIT_BAM/.bam/.BQSR.table}
    BAM_BQSR_RECAL_BAM=${BAM_SPLIT_BAM/.bam/.BQSR_recal.bam}
    
    VCF_FILE_PRE_FILTERING=$SAMPLE_TMP_FOLDER"/"$SAMPLE_NAME_FOLDER".GATK.RNA_seq.pre_filtering.hg38.vcf"
    VCF_FILE_POST_FILTERING=$SAMPLE_TMP_FOLDER"/"$SAMPLE_NAME_FOLDER".GATK.RNA_seq.filtered.hg38.vcf"
    
    RESULT_VCF_FILE=$HOME"/work/Results/VCFs/"$SAMPLE_NAME_FOLDER".GATK.RNA_seq.hg38.vcf"
    RESULT_FILE_ANNOT=$HOME"/work/Results/Annotations/"$SAMPLE_NAME_FOLDER".annot.csv.avinput"

    if [ ! -d $SAMPLE_TMP_FOLDER ]; then
        mkdir -p $SAMPLE_TMP_FOLDER
    fi;
    
    if [ ! -f $SAMPLE_TMP_FOLDER"/Under_construction.txt" ]; then
        touch $SAMPLE_TMP_FOLDER"/Under_construction.txt"
    else
        continue
    fi;
    
    if [ -f $RESULT_FILE_ANNOT ]; then
        echo "Results file found, skipping. " $RESULT_FILE_ANNOT
        continue
    fi;

    echo "Processing "$SAMPLE_NAME_FOLDER

    if [ ! -f $FWD_READS_UNZIP ]; then
        echo "Unzipping "$FWD_READS" "
        gunzip -c $FWD_READS > $FWD_READS_UNZIP;
    fi;
    
    if [ ! -f $BKW_READS_UNZIP ]; then
        echo "Unzipping "$BKW_READS" "
        gunzip -c $BKW_READS > $BKW_READS_UNZIP;
    fi;
    
    if [ ! -d $OUTPUT_PATH_SAMPLE ]; then
        mkdir $OUTPUT_PATH_SAMPLE
    fi;
    
    if [ ! -f $FWD_READS_FILT ]; then
        echo "Adapter trimming of sample "$SAMPLE_NAME_FOLDER
        trim_galore --paired $FWD_READS_UNZIP $BKW_READS_UNZIP -o $SAMPLE_TMP_FOLDER --no_report_file
    fi;
    
    cd $SAMPLE_TMP_FOLDER
    
    if [ ! -f $BAM_ALIGNED ]; then
    
        STAR \
            --genomeDir $STAR_INDICES_FOLDER \
            --readFilesIn $FWD_READS_FILT $BKW_READS_FILT  \
            --outFileNamePrefix $SAMPLE_TMP_FOLDER"/"$SAMPLE_NAME"." \
            --runThreadN $NR_CPUS \
            --quantMode TranscriptomeSAM GeneCounts \
            --outFilterMultimapNmax 20 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --sjdbGTFfile $SJDBG_FILE \
            --outFilterMismatchNoverLmax 0.04 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --outSAMtype BAM Unsorted \
            --bamRemoveDuplicatesType UniqueIdentical
    
    fi;
    
    ## 2 pass
    
    if [ !  -f $SAMPLE_TMP_FOLDER"/Genome" ]; then
        echo "2nd pass genome generate"
        STAR \
            --runMode genomeGenerate \
            --genomeDir $SAMPLE_TMP_FOLDER \
            --genomeFastaFiles $REF_GEN \
            --sjdbFileChrStartEnd $SAMPLE_TMP_FOLDER"/"".SJ.out.tab" \
            --runThreadN $NR_CPUS
    fi
    
    if [ !  -f $SAMPLE_TMP_FOLDER"/"$SAMPLE_NAME_FOLDER".Aligned.out.bam" ]; then
        echo "2nd pass align"
        STAR \
            --genomeDir $SAMPLE_TMP_FOLDER \
            --outFileNamePrefix $SAMPLE_TMP_FOLDER"/"$SAMPLE_NAME_FOLDER"."\
            --readFilesIn $FWD_READS_FILT $BKW_READS_FILT  \
            --runThreadN $NR_CPUS \
            --outFilterMultimapNmax 19 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverLmax 0.04 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --outSAMtype BAM Unsorted
    fi
    
    if [ ! -f $BAM_PICARD ]; then
        echo $BAM_PICARD
        $PICARD_CALL AddOrReplaceReadGroups \
        RGSM= default \
        RGLB= library \
        RGPL= ILLUMINA \
        RGPU= default \
        INPUT=$BAM_ALIGNED \
        OUTPUT=$BAM_PICARD \
        CREATE_INDEX= True \
        SORT_ORDER= coordinate \
        VALIDATION_STRINGENCY= LENIENT
    fi
    
    if [ ! -f $BAM_PICARD_SORTED ]; then
        echo "Sorting before duplication detection "
        samtools sort -T $SAMPLE_TMP_FOLDER"/"$SAMPLE_NAME_FOLDER"_sorting_file_" -o $BAM_PICARD_SORTED $BAM_PICARD
        samtools index $BAM_PICARD_SORTED
    fi
    
    # picard map duplicates
    
    if [ ! -f $BAM_MARKED_DUPLICATES ]; then
    
        echo "Marking duplicates"
        $PICARD_CALL MarkDuplicates -Xmx64g \
        I=$BAM_PICARD_SORTED \
        O=$BAM_MARKED_DUPLICATES \
        M=$BAM_MARKED_DUPLICATES_TXT \
        AS=TRUE
    fi
    
    # SplitNCigarReads
    
    if [ ! -f $BAM_SPLIT_BAM ]; then
        echo "Indexing Split Read step"
        samtools index $BAM_MARKED_DUPLICATES
        
        echo "Splitting SplitNCigarReads"
        
        $GATK SplitNCigarReads\
        -R $REF_GEN \
        -I $BAM_MARKED_DUPLICATES \
        -O $BAM_SPLIT_BAM
        
        echo "Indexing "$BAM_SPLIT_BAM
        samtools index $BAM_SPLIT_BAM
    fi
    
    
    echo "BQSR determination"
    
    if [ ! -f $BAM_BQSR_RECAL_TABLE ]; then
        
        $GATK BaseRecalibrator \
        -R $REF_GEN \
        -I $BAM_MARKED_DUPLICATES \
        -O $BAM_BQSR_RECAL_TABLE \
        -known-sites $DBSNP_PATH \
        -known-sites $KNOWN_INDELS \
        -known-sites $KNOWN_INDELS_MILLS
    fi
    
    echo "BQSR recalculation"
    
    if [ ! -f $BAM_BQSR_RECAL_BAM ]; then
        
        $GATK ApplyBQSR \
        -R $REF_GEN \
        -I $BAM_MARKED_DUPLICATES \
        --bqsr-recal-file $BAM_BQSR_RECAL_TABLE \
        --use-original-qualities \
        -O $BAM_BQSR_RECAL_BAM
    fi
    
    # variant calling
    echo "GATK haplotype caller"
    
    if [ ! -f $VCF_FILE_PRE_FILTERING ]; then
    
        echo "Calling Variants"
        
        $GATK HaplotypeCaller \
            -R $REF_GEN \
            -I $BAM_SPLIT_BAM \
            -O $VCF_FILE_PRE_FILTERING \
            --dbsnp $DBSNP_PATH \
            --native-pair-hmm-threads $NR_CPUS
            #-I $BAM_BQSR_RECAL_BAM \
            #--stand_call_conf 20 \
    fi
    
    if [ ! -f $VCF_FILE_POST_FILTERING ]; then
        
        #post-filtering
        echo "GATK filtering"
        
        $GATK VariantFiltration \
            -filter "QD < 2.0" --filter-name "QD2" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -filter "SOR > 3.0" --filter-name "SOR3" \
            -filter "FS > 60.0" --filter-name "FS60" \
            -filter "MQ < 40.0" --filter-name "MQ40" \
            -R $REF_GEN \
            -V $VCF_FILE_PRE_FILTERING \
            -O $VCF_FILE_POST_FILTERING
    fi

    echo "Annotating"

    if [ ! -f $RESULT_FILE_ANNOT ]; then
        perl $HOME/work/annovar/table_annovar.pl \
            $VCF_FILE_POST_FILTERING \
            --buildver hg38 \
            --out $RESULT_FILE_ANNOT \
            -remove \
            --protocol refGene,exac03,avsnp150 \
            -operation gx,f,f \
            -nastring . \
            -polish \
            -vcfinput $HOME/work/annovar/humandb
    fi

    echo "Finished"
    break
}

for index in ${!FASTQ_FILES[@]}; do
    
    #if [[ "$index" -eq 0 ]]; then
    #    continue;
    #fi
    
    #if [[ "$index" -eq 1 ]]; then
    #    continue;
    #fi
    
    SAMPLE=${FASTQ_FILES[index]}
    echo "$index ${SAMPLE}"
    align_reads $SAMPLE
    #break
done