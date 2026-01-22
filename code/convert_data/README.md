# Project: Spinal cord fMRI analysis at 7T


```bash
git clone https://github.com/CarolineLndl/spine_7t_fmri_analysis.git
export PATH_CODE="${PATH_PROJECT}/spine_7t_fmri_analysis"
```

### 1. Dependencies 🔗

- [dcm2niix](https://github.com/rordenlab/dcm2niix) (installed automatically with the following code)


### 2. Get data into BIDS format 📑
#### 2.1 Convert mri data
The dataset downloaded from OpenNeuro is already organized in the BIDS format.
If you have DICOM data that you wish to add to the dataset,
use `dcm2bids` to convert the DICOMS into a BIDS dataset:

```bash
cd ${PATH_CODE}/code/

dcm2bids -d "${PATH_DATA}/sourcedata/sub-$ID/mri" \
          -p "${ID}" \
          -c "${PATH_CODE}/config/config_bids.txt" \
          -o "${PATH_DATA}/spine_7t_fmri_data"
```

- `${ID}` is the subject ID (e.g., 095 103)
- For full data conversion instructions, see: `${PATH_CODE}/code/convert_data/01_convert_mriData.sh`

#### 2.2 Convert physio data 
Use `${PATH_CODE}/code/convert_data/02_convert_physioData.sh` to convert raw physio data into BIDS format.

```bash
cd "${PATH_CODE}/code/convert_data"
bash 02_convert_physioData.sh
```

#### 2.3 Convert behavioral data 
Use `${PATH_CODE}/code/convert_data/03_convert_behData.sh` to convert raw behavioral data into BIDS format.
Participant ID and task names should be modified in the .sh script.

```bash
cd "${PATH_CODE}/code/convert_data"
bash 03_convert_behData.sh
```

