# Project: Spinal cord fMRI analysis at 7T

## Overview
Processing of spinal cord functional data acquired at 7T.

---

## 1. Getting Started

### 1.1 Dependencies 🔗
Your environment should include:
- Python (3.10.14 was used)
- Spinal Cord Toolbox 7.1
- Conda environment: `spine_7t_fmri_analysis/config/requirements.txt`
- FSL
- dcm2niix
- MATLAB (for denoising step only)

#### a. Set up your project paths
```bash
PATH_PROJECT=/cerebro/cerebro1/dataset/spine_7t/
PATH_DATA=$PATH_PROJECT/spine_7t_fmri_data/
PATH_CODE=$PATH_PROJECT/spine_7t_fmri_analysis/
```

#### b. Set up the toolbox paths
<details>
<summary>👉 How to install dependencies</summary>

**Toolboxes for preprocessing**
- Spinal Cord Toolbox 7.1: [Installation instructions](https://spinalcordtoolbox.com/en/latest/user_section/installation.html)
- FSL: see here [Installation instructions](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation)

**Toolboxes for denoising:**
- Verify wich version of MATLAB is compatible with your Python version (*vis versa*): see here [Compatibility table](https://www.mathworks.com/support/requirements/python-compatibility.html)
- Install MATLAB: see here [Installation instructions](https://www.mathworks.com/help/install/)
- Install MATLAB engine for Python: see here [Installation instructions](https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html)

```bash
# Example for MATLAB R2023b
LD_LIBRARY_PATH="/export01/local/matlab23b/sys/os/glnxa64:$toolbox_home/libraries"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/export01/local/matlab23b/bin/" # The LD_LIBRARY_PATH environment variable tells the system where to find shared libraries
cd /export/local/matlab23b/extern/engines #Navigate to MATLAB Folder, engines subfolder
python -m pip install matlabengine==23.2.10 #Install MATLAB engine for Python
```
</details>

```bash
SCT_DIR=$PATH_CODE/toolboxes/spinalcordtoolbox
FSLDIR=/cerebro/cerebro1/dataset/bmpd/derivatives/thibault_test/code/toolbox/fsl
export PATH="$SCT_DIR/bin:$PATH"   # spinalcordtoolbox
export PATH=${FSLDIR}/bin:${PATH} # FSL
export FSLDIR PATH
. $FSLDIR/etc/fslconf/fsl.sh

#MATLAB_DIR=/export01/local/matlab23b
#LD_PREFIX="${MATLAB_DIR}/sys/os/glnxa64:/cerebro/cerebro1/dataset/bmpd/derivatives/thibault_test/code/toolbox/libraries"
#export  LD_LIBRARY_PATH=/export01/local/matlab23b/bin/glnxa64/
```
> ⚠️ *Qt version conflict: SCT and MATLAB may use incompatible Qt libraries. If you don’t need MATLAB, consider commenting out the MATLAB path in the setup script to avoid errors. If you need MATLAB for denoising, uncomment the MATLAB path, but be aware that Qt-related errors may appear when using SCT manually.*

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
│   │       ├── sub-100_task-rest_acq-shimBase+3mm_physio.json
│   │       ├── sub-100_task-rest_acq-shimBase+3mm_physio.tsv.gz
│   │       └── ...
│   ├── sourcedata  # Original DICOM and behavioral data
│   │   ├── sub-100
│   │   │   ├── behav
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
          -o $root_dir/spine_7t_fmri_data/
```

- `$ID` is the subject ID (e.g., 103)
- For full data conversion instructions, see: `${PATH_CODE}/code/convert_data/01_convert_mriData.sh`

#### Convert physio data
Use `${PATH_CODE}/code/convert_data/02_convert_physioData.sh` to convert raw physio data into BIDS format.

```bash
cd ${PATH_CODE}/code/convert_data/
bash 02_convert_physioData.sh
``` 

---

## 2. Analysis Pipelines 📊
Files for preprocessing are in this repository.

- **code/**: Functions and code to run the analyses. Do not modify the file.
  - **convert_data/**: Scripts to convert raw mri and physio data into BIDS format.
- **config/**: Configuration files for paths and parameters.
  - `config_spine_7t_fmri.json` is used by `01_spine7T_preprocessing.ipynb`
  - `participants.tsv` contains demographical information and important info for preprocessing (*e.g.,* slice number for vertebrae labeling initiation)
- **template images**: Used for analyses; do not modify.
- **log**: Log files generated during processing run from bash script (the folder is not tracked by git).

### 2.1 Preprocessing 🤯
▸ runs preprocessing steps automatically with with output log from STDOUT
▸ By default all the steps are rerun even if some outputs already exist. If manual corrections were made, these files will be used as input for subsequent steps.

```bash
bash ${PATH_CODE}/code/run_batch_preprocessing.sh
```

⚠️ *Each step manually modified will imply that all subsequent steps need to be re-run. </span>* <br><br>
  
##### Visual check and manual corrections ✏️ 
<details>
<summary>For more details, click to expand </summary>

  - **I.a Motion correction (mask)** : ✏️
  check the automatic centerline detection and the mask used for motion correction, if needed, manually correct the centerline with:
  ```
  ctrl_sc_files_, mask_sc_files_=preprocess_Sc.moco_mask(ID=ID,i_img=mean_func_f[ID][tag][run_nb],
                                                                       radius_size=25,task_name=tag,
                                                                       manual=True,
                                                                       redo_ctrl=True,
                                                                       redo_mask=True,
                                                                       verbose=verbose)
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
  Check the automatic labeling of the inter vertebral disks on the anatomical image, if needed (now default is manual), manually correct the labeling with:
  ```
  vert_labels_files.append(preprocess_Sc.label_vertebrae(ID=ID,
                                                               i_img=raw_anat[ID_nb],
                                                               seg_img=seg_anat_sc_files[ID_nb],
                                                               c="t2",
                                                               initz=f"{z_value},{vert}",auto=False,
                                                               redo=False,
                                                               verbose=verbose))
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

### 2.2 Denoising 🧹

Should be run after preprocessing.
- ⚠️ csf segmentation should be checked and manually corrected if needed before running the denoising.
- Details on the different steps are in the .py script and will be added in the Readme later.

#### Two options to run preprocessing:

▸ runs steps automatically: recommanded to run all steps at once 
▸ By default all the steps are rerun even if some outputs already exist.
```bash
bash ${PATH_CODE}/code/run_batch_denoising.sh
```

### 2.3 First-level Analysis (TBD) 📈
