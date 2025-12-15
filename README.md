# Project: Spinal cord fMRI analysis at 7T

## Overview
Processing of spinal cord functional data acquired at 7T.

---

## 1. Getting Started

### 1.1 Dependencies 🔗
To run the processing pipeline, external dependencies are required. They are listed below but more details on how to install them is described further in a), b) and c).
- Conda (environment file: `spine_7t_fmri_analysis/config/requirements.txt`)
- Spinal Cord Toolbox 7.1
- FSL
- MATLAB (for denoising step only)

#### a. Clone this repository and download the dataset
We recommend creating a folder for this project (`PATH_PROJECT`: "spine_7t"). In this folder, clone this repository (`PATH_CODE`: "spine_7t_fmri_analysis") and download the dataset (`PATH_DATA`: "spine_7t_fmri_data").

```bash
mkdir "spine_7t"
cd "spine_7t"
git clone https://github.com/CarolineLndl/spine_7t_fmri_analysis.git
git clone https://github.com/OpenNeuroDatasets/ds007067.git spine_7t_fmri_data
export PATH_PROJECT="$(pwd)"
export PATH_DATA=${PATH_PROJECT}/spine_7t_fmri_data/
export PATH_CODE=${PATH_PROJECT}/spine_7t_fmri_analysis/
```

The dataset is also accessible to download manually here: https://openneuro.org/datasets/ds007067/download

#### b. Set up the toolbox paths
<details>
<summary>👉 How to install dependencies</summary>

**Toolboxes for preprocessing**
- Spinal Cord Toolbox 7.1: [Installation instructions](https://spinalcordtoolbox.com/en/latest/user_section/installation.html)
- FSL: see here [Installation instructions](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation)

**Toolboxes for denoising:**
We use some Matlab toolboxes for denoising. Specific MATLAB versions are compatible with specific Python versions. Depending on your version of MATLAB, edit the conda environment file (`spine_7t_fmri_analysis/config/requirements.txt`) with your version of Python and your version of the MATLAB engine.
- Verify which version of MATLAB is compatible with your Python version (*vis versa*): see here [Compatibility table](https://www.mathworks.com/support/requirements/python-compatibility.html)
- Install MATLAB: see here [Installation instructions](https://www.mathworks.com/help/install/)
- Install MATLAB engine for Python: see here [Installation instructions](https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html).

```bash
# Example for MATLAB R2023b
LD_LIBRARY_PATH="/export01/local/matlab23b/sys/os/glnxa64:$toolbox_home/libraries"
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:"/export01/local/matlab23b/bin/" # The LD_LIBRARY_PATH environment variable tells the system where to find shared libraries
cd /export/local/matlab23b/extern/engines #Navigate to MATLAB Folder, engines subfolder
python -m pip install matlabengine==23.2.10 #Install MATLAB engine for Python
```
</details>

If you have a single version of SCT and FSL, `SCT_DIR` and `FSLDIR` should already be set. You can manually check them with `echo "${SCT_DIR}"` and `echo "${FSLDIR}"`. Their binaries should also be in your `PATH` variable. If you have multiple versions installed, make sure the correct versions are set. An example of how to do that is shown below.

```bash
SCT_DIR=${PATH_CODE}/toolboxes/spinalcordtoolbox
FSLDIR=/cerebro/cerebro1/dataset/bmpd/derivatives/thibault_test/code/toolbox/fsl
export PATH="${SCT_DIR}/bin:${PATH}"   # spinalcordtoolbox
export PATH=${FSLDIR}/bin:${PATH} # FSL
export FSLDIR PATH
source ${FSLDIR}/etc/fslconf/fsl.sh
```

#### c. Setup the conda environment
<details>
<summary>👉 How to create the conda environment </summary>

Make sure conda is installed: see here [Installation instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
Create the appropriate conda environment:
*If you just what to run the preprocessing you can remove matlabengine from the requirements.txt file.*

```bash
conda create --name spine_7T_env_py10 python=3.10
conda activate spine_7T_env_py10
pip install -r config/requirements.txt
```
</details>

Load conda environment:

```bash
anaconda_dir=$(conda info --base)
source ${anaconda_dir}/etc/profile.d/conda.sh
source activate spine_7T_env_py10
```

### 1.2 Data organization 📑
Files are organized according to the BIDS standard:
<details>
<summary>Click to expand folder tree</summary>

```
├── spine_7t_fmri_analysis  # GitHub repository
│   ├── code
│     ├── convert_data
│     │   ├── 00_convert_mriData.sh
│     │   └── ...
│     ├── ...
│   ├── config
│     ├── config_spine_7t_fmri.json
│     ├── participants.tsv
│     └── ...
│   ├── template
│     ├── ...
│   └── log
│       ├── ...
├── spine_7t_fmri_data # Data directory
│   ├── derivatives
│   │   │   ├── ...
│   │   ├── manual  # Manually corrected files
│   │   │   └── sub-100
│   │   │       ├── anat
│   │   │       │   ├── sub-100_T2star_space-orig_label-ivd_mask.nii.gz
│   │   │       │   └── sub-100_T2star_space-orig_label-ivd_mask.nii.gz
│   │   │       └── func
│   │   │           ├── task-motor_acq-shimBase+3mm
│   │   │           │   ├── sub-100_task-motor_acq-shimBase+3mm_bold_moco_mean_seg.nii.gz
│   │   │           │   ├── sub-100_task-motor_acq-shimBase+3mm_bold_tmean_centerline.csv
│   │   │           │   └── sub-100_task-motor_acq-shimBase+3mm_bold_tmean_centerline.nii.gz
│   │   │           ├── ...
│   │   └── processing
│   │       ├── preprocessing
│   │       │   ├── QC  # QC reports
│   │       │   │   ├── ...
│   │       │   └── sub-100
│   │       │       ├── anat
│   │       │       │   ├── sct_deepseg
│   │       │       │   │   ├── sub-100_T2star_seg.json
│   │       │       │   │   └── sub-100_T2star_seg.nii.gz
│   │       │       │   ├── sct_label_vertebrae
│   │       │       │   │   ...
│   │       │       │   ├── sct_register_to_template
│   │       │       │   │   ...
│   │       │       │   └── sub-100_T2star.nii.gz
│   │       │       └── func
│   │       │           ├── task-motor_acq-shimBase+3mm
│   │       │           │   ├── sct_deepseg
│   │       │           │   │   ├── sub-100_task-motor_acq-shimBase+3mm_bold_moco_mean_seg.json
│   │       │           │   │   └── sub-100_task-motor_acq-shimBase+3mm_bold_moco_mean_seg.nii.gz
│   │       │           │   ├── sct_fmri_moco
│   │       │           │   │   ...
│   │       │           │   ├── sct_get_centerline
│   │       │           │   │   ...
│   │       │           │   ├── sct_propseg
│   │       │           │   │   ...
│   │       │           └── task-motor_acq-shimSlice+3mm
│   │       │               ...
│   │       └── ...  # Other processing steps (denoising, first-level analysis, etc)
│   ├── dataset_description.json
│   ├── sub-100
│   │   ├── anat
│   │   │   ├── sub-100_T2star.json
│   │   │   └── sub-100_T2star.nii.gz
│   │   └── func
│   │       ├── sub-100_task-motor_acq-shimBase+3mm_bold.json
│   │       ├── sub-100_task-motor_acq-shimBase+3mm_bold.nii.gz
│   │       └── ...
│   ├── sourcedata  # Original DICOM and behavioral data
│   │   ├── sub-100
│   │   │   ├── beh
│   │   │   │   ├── *.csv
│   │   │   │   ├── *.log
│   │   │   │   ├── *.psydat
│   │   │   │   └── ...
│   │   │   ├── mri
│   │   │   │   ├── 01-localizer_iso_ND
│   │   │   │   ├── *.dcm
│   │   │   │   └── ...
│   │   │   ├── ...
│   │   └── pmu
│   │       ├── ...
```

</details>

### 1.3 Get data into BIDS format 🗂️
#### Convert mri data
Use `dcm2bids` to convert raw mri data:

```bash
cd ${PATH_CODE}/code/

dcm2bids -d ${PATH_DATA}/sourcedata/sub-$ID/mri/ \
          -p $ID \
          -c ${PATH_CODE}/config/config_bids.txt \
          -o ${PATH_DATA}/spine_7t_fmri_data/
```

- `$ID` is the subject ID (e.g., 095 103)
- For full data conversion instructions, see: `${PATH_CODE}/code/convert_data/01_convert_mriData.sh`

#### Convert physio data (need to be update)
Use `${PATH_CODE}/code/convert_data/02_convert_physioData.sh` to convert raw physio data into BIDS format.

```bash
cd ${PATH_CODE}/code/convert_data/
bash 02_convert_physioData.sh
```

---

## 2. Analysis Pipelines 📊
Files for preprocessing are in this repository.

- **code/**: Functions and code to run the analyses. Do not modify the file.
  - `preprocessing.py` > library of preprocessing functions
  - `preprocessing_workflow.py` > orchestrates preprocessing steps using the functions
  - `run_all_processing.sh` > shell script to launch any combination of workflows (so far only one workflow)
  - **convert_data/**: Scripts to convert raw mri and physio data into BIDS format.
- **config/**: Configuration files for paths and parameters.
  - `config_spine_7t_fmri.json` is used by `preprocessing_workflow.py`
  - `participants.tsv` contains demographical information and important info for preprocessing (*e.g.,* slice number for vertebrae labeling initiation)
- **template images**: Used for analyses; do not modify.
- **log**: Log files generated during processing run from bash script (the folder is not tracked by git).

### 2.1 Preprocessing 🤯
▸ runs preprocessing steps automatically with output log from STDOUT   
▸ By default all the steps are rerun even if some outputs already exist. If manual corrections were made, these files will be used as input for subsequent steps.  
▸ if you already setup the PATH_CODE and PATH_DATA you don't need to specify --path_data --path_code  
▸ Specify individuals to process (--ids XXX), the default option run preprocessing on all participants in the `participants.tsv`

```bash
bash ${PATH_CODE}/code/run_all_processing.sh --path_data ${PATH_DATA} --path_code ${PATH_CODE} --ids 090 101 106

```

⚠️ *Each step manually modified will imply that all subsequent steps need to be re-run. </span>* <br><br>

##### Visual check and manual corrections ✏️
<details>
<summary>For more details, click to expand </summary>

  - **I.a Motion correction (mask)** : ✏️
  check the automatic centerline detection and the mask used for motion correction, if needed, manually correct the centerline you can modify the line 43 of the run_all_processing.sh:
  ```
  nohup python -u ../code/preprocessing_workflow.py --ids "${IDs[@]}" --redo True --manual_centerline True \
  ```

  The output files can be found in:
  ```
  /spine_7t_fmri_analysis/derivatives/manual/sub-<ID>/func/
      └── <task*_acq*>/
          ├── sub-<ID>_<task_acq>_bold_tmean_centerline.csv
          └── sub-<ID>_<task_acq>_bold_tmean_centerline.nii.gz

  ```

  - **II Segmentation** ✏️
  Check the segmentation results, if needed, manually correct the segmentation in fsleyes using the anatomical image or mean functional image as background.
 When saving the corrected segmentation, make sure to keep the same name as the original segmentation file but save it in the `manual` folder:
  ```
  /spine_7t_fmri_analysis/derivatives/manual/sub-<ID>/func
      └── <task*_acq*>/
          └── sub-<ID>_<task_acq>_bold_moco_mean_seg.nii.gz
  ```

  - **III Labeling of inter vertebral disk** ✏️
  Check the automatic labeling of the inter vertebral disks on the anatomical image, if needed (now default is manual), you can modify the line 43 of the run_all_processing.sh :
  ```
  nohup python -u ../code/preprocessing_workflow.py --ids "${IDs[@]}" --redo True --auto_vert_labels False \
  ```
  The output files can be found in:
  ```
  /spine_7t_fmri_analysis/derivatives/manual/sub-<ID>/anat
      └── sub-<ID>_T2star_space-orig_label-ivd.nii.gz
  ```
</details>


##### ‼️ What we want to try to improve
> - **I. Motion correction:** try different parameters for the mask size, or different reference images (mean functional, middle volume, etc).
> - **IV. Registration to template:** check if the parameters for the registration are ok.

### 2.2 Denoising (work in progress) 🧹

Should be run after preprocessing.
- ⚠️ csf segmentation should be checked and manually corrected if needed before running the denoising.
- Details on the different steps are in the .py script and will be added in the Readme later.

#### Two options to run preprocessing:

▸ runs steps automatically: recommanded to run all steps at once
▸ By default all the steps are rerun even if some outputs already exist.
```bash
bash ${PATH_CODE}/code/run_all_processing.sh
```

### 2.3 First-level Analysis (TBD) 📈
