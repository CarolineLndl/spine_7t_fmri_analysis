# Project: Spinal cord fMRI analysis at 7T

## Overview
Processing of spinal cord functional data acquired at 7T.

---

## Getting Started

### Set up your project paths

Create a folder that will contain the code of this repository as well as the source and processed data, then define the variable in SHELL:

```bash
export PATH_PROJECT=<PATH_TO_PROJECT>
```

### Download data 📀

See: https://openneuro.org/datasets/ds007067/download

<details>
<summary>Files are organized according to the BIDS standard.</summary>

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

Define variable:
```bash
export PATH_DATA="${PATH_PROJECT}/ds007067"
```

### Clone repository

```bash
git clone https://github.com/CarolineLndl/spine_7t_fmri_analysis.git
export PATH_CODE="${PATH_PROJECT}/spine_7t_fmri_analysis"
```

### Dependencies 🔗

#### External dependencies

- [Spinal Cord Toolbox v7.1](https://spinalcordtoolbox.com/en/latest/user_section/installation.html)
- [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation)
- [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

#### Setup the conda environment

Your environment should include:
- Python (tested with 3.10.14, but other versions could work)
- [dcm2niix](https://github.com/rordenlab/dcm2niix) (installed automatically with the following code)

Create the appropriate conda environment:

```bash
conda create --name spine_7T_env_py10 python=3.10
conda activate spine_7T_env_py10
pip install -r "${PATH_CODE}/config/requirements.txt"
```

### Data organization 📑

### 1.3 Get data into BIDS format 🗂️
#### Convert mri data
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

#### Convert physio data (need to be update)
Use `${PATH_CODE}/code/convert_data/02_convert_physioData.sh` to convert raw physio data into BIDS format.

```bash
cd "${PATH_CODE}/code/convert_data"
bash 02_convert_physioData.sh
```

---

## Analysis Pipeline ⚙️

<details><summary>Here is a brief description of the files used for data analysis.</summary>

- **`code/`**: Functions and code to run the analyses. Do not modify the file.
  - `preprocessing.py` > library of preprocessing functions
  - `preprocessing_workflow.py` > orchestrates preprocessing steps using the functions
  - `run_all_processing.sh` > shell script to launch any combination of workflows (so far only one workflow)
  - **`convert_data/`**: Scripts to convert raw mri and physio data into BIDS format.
- **`config/`**: Configuration files for paths and parameters.
  - `config_spine_7t_fmri.json` is used by `preprocessing_workflow.py`
  - `participants.tsv` contains demographical information and important info for preprocessing (*e.g.,* slice number for vertebrae labeling initiation)
- **`template`**: Used for analyses; do not modify.
- **`log`**: Log files generated during processing run from bash script (the folder is not tracked by git).

</details>

Run the pipeline:

```bash
bash "${PATH_CODE}/code/run_all_processing.sh" --path-data "${PATH_DATA}" --path-code "${PATH_CODE}" --ids 090 101 106 --tasks motor --no-denoising
```

- Runs preprocessing steps automatically with output log from STDOUT.
- By default all the steps are rerun even if some outputs already exist. If manual corrections were made, these files will be used as input for subsequent steps.
- If you have already setup `PATH_CODE` and `PATH_DATA`, you don't need to specify `--path-data` and `--path-code`.
- Specify individuals to process (`--ids XXX`), the default option run preprocessing on all participants in the `participants.tsv`.
- Specify task to process (`--tasks` `motor` or `rest`), the default option runs preprocessing on all tasks defined in the `config_file_7t_fmri.json`
- You can remove `--no-denoising` to run denoising after preprocessing (see section 2.2).

> [!WARNING]  
> Each step manually modified will imply that all subsequent steps need to be re-run.

### Visual check and manual corrections ✏️

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
  Check the segmentation results, if needed, manually correct the segmentation in FSLeyes using the anatomical image or mean functional image as background.
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

### 2.2 Denoising  🧹

Should be run after preprocessing.
- ⚠️ csf segmentation should be checked and manually corrected if needed before running the denoising.

#### Description of the denoising steps
- **I. Extract motion parameters:** to regressed out the residual motion effects, we extracted slice-wise motion parameters from the moco files generated during the motion correction step.
- **II. Compute outliers calculation:** to identify volumes with excessive motionas and included as nuisance regressors in the denoising step.
- **III. CompCor calculation:** to eliminate the non-neural aspects of the signal, we used the CompCor (Behzadi et al., 2007) approach by extracting the mean signal and the first five principal components of the unsmoothed signal recorded from the CSF (CSF-mask in functional space).
- **IV. Signal cleaning:** we regressed out the nuisance regressors (motion parameters, outliers, mean CSF signal and first five principal components from the CSF) and applied a high-pass temporal filter (cut-off frequency = 0.01 Hz) to the functional data.

#### Run denoising
- Runs preprocessing steps automatically with output log from STDOUT.
- By default all the steps are rerun even if some outputs already exist.
- If you already have setup `PATH_CODE` and `PATH_DATA`, you don't need to specify `--path-data` and `--path-code`.
- Specify individuals to denoise (`--ids XXX`), the default option run denoising on all participants in the `participants.tsv`.
- Specify task to denoise (`--tasks` `motor` or `rest`), the default option run denoising on all tasks defined in the `config_file_7t_fmri.json`.
- You can remove `--no-preprocess` to run preprocessing before denoising (see section 2.1).

```bash
bash "${PATH_CODE}/code/run_all_processing.sh" --path-data "${PATH_DATA}" --path-code "${PATH_CODE}" --ids 090 101 106 --tasks motor --no-preprocess

```

### 2.3 First-level Analysis (TBD) 📈
