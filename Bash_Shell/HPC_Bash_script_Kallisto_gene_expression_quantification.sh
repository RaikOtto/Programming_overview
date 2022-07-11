#!/bin/bash

source activate first-steps
##sbatch --job-name=BAMs_new --output=/fast/users/ottor_c/work/MAPTor_NET/logs/BAMs_new_%x-%j.Kallisto.log --ntasks=4 --nodes=1 --mem=40G  --partition=long --time=1440:00 /fast/users/ottor_c/work/MAPTor_NET/Scripts/Kallisto_quantification.sh

#find ./* -type f -name '*.sra' -print0 | xargs -0 -Imysongs mv -i mysongs .

PREFIX=$HOME"/work/"
KALLISTO_INDEX=$PREFIX"/Data/Homo_sapiens.GRCh38.kallisto.idx"
REF_GEN_V="hg38"
#STUDY="Groetzinger"
NR_CUPS=4

INPUT_FOLDER="/data/projects/Sers-DKTKMasterTranscriptomics/work/2021-03-30_download/download/"
#INPUT_FOLDER="/fast/users/ottor_c/scratch/JGA/"
TMP_FOLDER="/fast/users/ottor_c/scratch/MAPTOR/"
#TMP_FOLDER="/fast/users/ottor_c/scratch/JGA/"
OUTPUT_FOLDER=$PREFIX"/MAPTor_NET/BAMs_new//"
#OUTPUT_FOLDER=$PREFIX"/Deko_Projekt/Results/JGA/"

if [ ! -d $OUTPUT_FOLDER ]; then
	mkdir $OUTPUT_FOLDER
fi

#FWD_ENDING="_1.fq.gz"
#FWD_ENDING="_R1.fastq.gz"
FWD_ENDING="_R1_complete_filtered.fastq.gz"
BKW_ENDING=${FWD_ENDING//"1"/"2"}

FASTQ_FILES=(`find $INPUT_FOLDER -type f | grep $FWD_ENDING'$'`)

function kallisto_quant {

	FWD_READS=$1
	BKW_READS=${FWD_READS//$FWD_ENDING/$BKW_ENDING}
    
    SAMPLE_INPUT_FOLDER=${FWD_READS_UNZIP%'/'*}'/'
    #SAMPLE_NAME_FOLDER=${FWD_READS//FWD_ENDING/""}
    SAMPLE_NAME_FOLDER=${FWD_READS//$FWD_ENDING/""}
    SAMPLE_NAME_FOLDER=(${SAMPLE_NAME_FOLDER//\// })
    SAMPLE_NAME_FOLDER=${SAMPLE_NAME_FOLDER[-1]}
    
    SAMPLE_NAME_FOLDER=${SAMPLE_NAME_FOLDER[6]}
    SAMPLE_NAME_FOLDER=${SAMPLE_NAME_FOLDER//"H021-"/""}
    SAMPLE_NAME_FOLDER=${SAMPLE_NAME_FOLDER//"P021-"/""}

    SAMPLE_TMP_FOLDER=$TMP_FOLDER"/"$SAMPLE_NAME_FOLDER
    OUTPUT_PATH_SAMPLE=$OUTPUT_FOLDER"/"$SAMPLE_NAME_FOLDER"/"
    
    FWD_READS_UNZIP=$SAMPLE_TMP_FOLDER"/"$SAMPLE_NAME_FOLDER"_R1.fastq"
    BKW_READS_UNZIP=$SAMPLE_TMP_FOLDER"/"$SAMPLE_NAME_FOLDER"_R2.fastq"
    FWD_READS_FILT=${FWD_READS_UNZIP//'_R1.fastq'/"_R1_val_1.fq"}
    BKW_READS_FILT=${BKW_READS_UNZIP//'_R2.fastq'/"_R2_val_2.fq"}
    
    #FWD_READS_UNZIP=$SAMPLE_TMP_FOLDER"/"$SAMPLE_NAME_FOLDER"_1.fq"
    #BKW_READS_UNZIP=$SAMPLE_TMP_FOLDER"/"$SAMPLE_NAME_FOLDER"_2.fq"
    #FWD_READS_FILT=${FWD_READS_UNZIP//'_1.fq'/"_1_val_1.fq"}
    #BKW_READS_FILT=${BKW_READS_UNZIP//'_2.fq'/"_2_val_2.fq"}
    echo $SAMPLE_TMP_FOLDER
 
    
    if [ -d $SAMPLE_TMP_FOLDER ]; then
        continue
    else
        mkdir -p $SAMPLE_TMP_FOLDER
    fi;
    
    if [ ! -d $OUTPUT_PATH_SAMPLE ]; then
        mkdir $OUTPUT_PATH_SAMPLE
    fi;
    
    echo "Processing "$SAMPLE_NAME_FOLDER
    
    if [ ! -f $FWD_READS_UNZIP ]; then
        echo "Unzipping "$FWD_READS" "
        gunzip -c $FWD_READS > $FWD_READS_UNZIP;
        gunzip -c $BKW_READS > $BKW_READS_UNZIP;
    fi;

    if [ ! -f $FWD_READS_FILT ]; then
        echo "Adapter trimming of sample "$SAMPLE_NAME_FOLDER
        trim_galore --paired $FWD_READS_UNZIP $BKW_READS_UNZIP -o $SAMPLE_TMP_FOLDER --no_report_file
    fi;
    
    fastqc $FWD_READS_FILT $BKW_READS_FILT
    
    if [ ! -f $OUTPUT_PATH_SAMPLE"/abundance.tsv" ]; then

        kallisto quant \
            --index=$KALLISTO_INDEX \
            --output-dir=$OUTPUT_PATH_SAMPLE \
            --threads=$NR_CUPS \
            -b $NR_CUPS \
            $FWD_READS_FILT \
            $BKW_READS_FILT
    fi;

    echo "Finished processing sample "$SAMPLE_NAME_FOLDER
}

for INDEX in ${!FASTQ_FILES[@]}; do

    SAMPLE=${FASTQ_FILES[INDEX]}
    echo $INDEX $SAMPLE
    kallisto_quant $SAMPLE 

done
