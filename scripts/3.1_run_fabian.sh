# Description: This script runs the FABIAN variant tool on the input VCF file for a given cancer type.
# The input VCF file is split into multiple files if there are more than 10000 mutations.
# The output is saved in the data/fabian_output directory.

CANCER_TYPE=$1
echo $CANCER_TYPE

# create output directory
OUTPUT_DIR=data/fabian_output
mkdir -p $OUTPUT_DIR

# split VCF files if there are more than 10000 mutations
# as this exceeds the input limit set by FABIAN
vcf_file=data/fabian_input/${CANCER_TYPE}_input.vcf
num_lines=$(wc -l $vcf_file | awk '{print $1}')
echo $num_lines
if [ $num_lines -gt 10000 ]; then
    output_prefix=data/fabian_input/${CANCER_TYPE}_input_
    split -l 10000 ${vcf_file} ${output_prefix}
    # rename files to have .vcf extension
    for file in ${output_prefix}*; do
        mv $file $file.vcf
    done
fi

# get FABIAN output for each VCF file 
vcf_files=$(ls data/fabian_input/${CANCER_TYPE}_input_*.vcf)
if [ -z "$vcf_files" ]; then
    vcf_files=$vcf_file
fi
for vcf_file in ${vcf_files[@]}; do
    echo $vcf_file
    
    # generate output file name
    filename=$(basename $vcf_file)
    filename=${filename%.vcf}
    filename=${filename/input/output}
    echo $filename

    if [ -f ${OUTPUT_DIR}/${filename}.tsv ]; then
        echo "FABIAN output file already exists. Skipping..."
        continue
    fi

    # submit VCF file to FABIAN
    printf "($(date +%T)) Submitting " && \
    FABIANID=$( curl -sLD - -o /dev/null \
    -F "mode=vcf" \
    -F "filename=@${vcf_file}" \
    -F "genome=hg19" \
    -F "tfs_filter=all" \
    -F "models_filter=tffm_d" \
    -F "models_filter=tffm_fo" \
    -F "models_filter=pwm" \
    -F "dbs_filter=jaspar2022" \
    -F "dbs_filter=cisbp_1.02" \
    -F "dbs_filter=HOCOMOCOv11" \
    -F "dbs_filter=hPDI" \
    -F "dbs_filter=jolma2013" \
    -F "dbs_filter=SwissRegulon" \
    -F "dbs_filter=UniPROBE" \
    https://www.genecascade.org/fabian/analyse.cgi \
    | grep -m 1 "Location: " | grep -o "\([0-9]\+_[0-9]\+\)" ) && \
    i=1; until curl -sfo fabian.data_${FABIANID}.zip \
    https://www.genecascade.org/temp/QE/FABIAN/${FABIANID}/fabian.data.zip; \
    do printf "\r($(date +%T)) Waiting for $FABIANID"; \
    [ $i == 30 ] && sleep $i || sleep $((i++)); done && \
    printf "\r($(date +%T)) Saved file fabian.data_${FABIANID}.zip\n"

    # move result to output folder
    unzip fabian.data_${FABIANID}.zip -d ${OUTPUT_DIR}/${filename}
    mv ${OUTPUT_DIR}/${filename}/fabian.data ${OUTPUT_DIR}/${filename}.tsv

    # remove temporary files
    rm fabian.data_${FABIANID}.zip
    rm -rf ${OUTPUT_DIR}/${filename}
done

# remove split VCF files
vcf_files=$(ls data/fabian_input/${CANCER_TYPE}_input_*.vcf)
if [ -n "$vcf_files" ]; then
    rm -f $vcf_files
fi
