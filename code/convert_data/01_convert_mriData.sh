root_dir=/cerebro/cerebro1/dataset/spine_7T/

cd  $root_dir"/spine_7T_analysis/config/"
source spine_7T_env_032024.sh

ID=101


# Create participant directory
cd /cerebro/cerebro1/dataset/spine_7T/sourcedata/

mkdir "sub-"$ID  # create folder
mkdir "sub-"$ID"/pmu"
mkdir "sub-"$ID"/mri"
mkdir "sub-"$ID"/behav"

# Claim and download the file
#1: login to the BIC, usually at login.bic.mni.mcgill.ca
cd /data/dicom/ 
find_mri "acdc_spine_7T_"$ID

cd $root_dir"/sourcedata/"
file=/data/dicom/acdc_spine_7T_101_20251106_091631571
rsync -a "/data/dicom/"$file $root_dir"/sourcedata/sub-"$ID"/" # download the data

#rsync -a /data/dicom/acdc_spine_7T_90_20250728_163015298 $main_dir"/sourcedata/sub-"$ID"/" # download the data

# sort dicom files
# PMU and behavioral data should be copyed manually into the pmu and behav folder
cd $root_dir"/spine_7T_analysis/code/convert_data/"
python sortDCM.py -d $root_dir"/sourcedata/sub-"$ID"/"$file -o $root_dir"/sourcedata/sub-"$ID"/mri/"


#Convert in BIDS
cd $root_dir"/spine_7T_analysis/code/convert_data/"
dcm2bids -d $root_dir"/sourcedata/sub-$ID/mri/" -p $ID -c $root_dir"/spine_7T_analysis/config/config_bids_6Nov25.txt" -o $root_dir"/rawdata/"

#Compress physio files ------------------------------------------------
cd $root_dir"/sourcedata/sub-"$ID"/pmu/"
EXTENSIONS=("ext" "puls" "resp")
# Collect basenames from the known extensions only
basenames=$(for ext in "${EXTENSIONS[@]}"; do
    find . -maxdepth 1 -type f -name "*.${ext}" \
        | sed 's|^\./||' | sed "s/\.${ext}$//"
done | sort -u)

for base in $basenames; do
    files=()

    # Collect existing files for this basename
    for ext in "${EXTENSIONS[@]}"; do
        f="${base}.${ext}"
        [[ -f "$f" ]] && files+=("$f")
    done

    # Must have exactly 3 files (.ext, .puls, .resp)
    if [[ ${#files[@]} -ne ${#EXTENSIONS[@]} ]]; then
        echo "Skipping $base — missing .ext, .puls, or .resp"
        continue
    fi

    # Create archive
    tar -czf "${base}.tar.gz" "${files[@]}"
    echo "Created: ${base}.tar.gz"
done

#Convert physio to BIDS
cd $root_dir"/spine_7T_analysis/code/convert_data/"
python physio2bids.py -t $root_dir"sourcedata/"$ID"/pmu/"$ID"_rest"$func_run".tar.gz" -s $ID -o $root_dir"/rawdata/" -v True

