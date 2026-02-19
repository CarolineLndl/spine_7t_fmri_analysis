import os
import glob
import json
import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib.pyplot as plt
from nilearn.plotting import plot_design_matrix
from nilearn.glm.first_level import FirstLevelModel
from nilearn.glm.second_level import SecondLevelModel

from preprocess import Preprocess_main, Preprocess_Sc

#####################################################
class Postprocess_main:
    '''
    The Postprocess_main class is used to setup the Post-processing path and execute the Post-processing steps.

    Attributes
    ----------
    config : dict
        Defining all the parameters of the analysis including the path to the raw data, the participants to analyze, the design of the experiment, and the preprocessing parameters
    IDs : list
        List of participant IDs to process (e.g., ['A001', 'A002'])
    verbose : bool
        Whether to print information during the each step (default: True)
    '''

    def __init__(self, config, IDs=None,verbose=True):
        if IDs==None:
            raise ValueError("Please provide the participant ID (e.g., _.stc(ID='A001')).")
        
        # Class attributes -------------------------------------------------------------------------------------
        self.config = config # load config info
        self.participant_IDs= IDs # list of the participants to analyze
        self.raw_dir = os.path.join(self.config["raw_dir"])  # directory of the raw data
        self.derivatives_dir = os.path.join(self.config["raw_dir"], self.config["derivatives_dir"])  # directory of the derivatives data
        self.first_level_dir = os.path.join(self.config["raw_dir"], self.config["first_level"]["dir"])  # directory of the derivatives data
        self.second_level_dir = os.path.join(self.config["raw_dir"], self.config["second_level"]["dir"])  # directory of the second-level analysis data
        self.manual_dir = os.path.join(self.config["raw_dir"], self.config["manual_dir"])  # directory of the manual corrections

        # Create directories -------------------------------------------------------------------------------------
        for ID in self.participant_IDs:
            ID_first_level_dir=self.first_level_dir.format(ID)
            os.makedirs(ID_first_level_dir, exist_ok=True)

            # Create a folder for each task in participant folder
            if "design_exp" in self.config.keys():
                for ses_name in self.config["design_exp"]['ses_names']:
                    ses_dir=ses_name if int(self.config["design_exp"]["ses_nb"])>1 else ""
                    if "acq_names" in self.config["design_exp"].keys():
                        for task_name in self.config["design_exp"]['task_names']:
                            for acq_name in self.config["design_exp"]['acq_names']:
                                tag="task-" + task_name + "_acq-" + acq_name
                                os.makedirs(ID_first_level_dir + tag ,exist_ok=True)
        



    def run_first_level_glm(self, ID=None, i_fname=None,events_file=None,mask_file=None,task_name=None,run_name=None,contrasts = ["trial_RH-rest", "trial_RH", "rest"],smoothing_fwhm=1.5,verbose=True,redo=False):
        """
        Run first-level GLM for a specific subject and task.

        Parameters
        ----------
        ID : str
            Participant ID (e.g., "093")
        i_fname : str
            Filename of the input fMRI image (4D NIfTI file)
        events_file : str
            Filename of the events TSV file
        mask_file : str
            Filename of the mask NIfTI file where to restrict the analysis
        task_name : str
            Task name (e.g., "motor_acq-shimBase+3mm")
        contrasts : list of str, optional
            List of contrasts to compute (default is ["trial_RH-rest", "trial_RH", "rest"])
        smoothing_fwhm : float, optional
            Full-width at half-maximum for spatial smoothing (default is 1.5 mm)
        verbose : bool, optional
            Whether to print information during processing (default is True)
        redo : bool, optional
            Whether to redo the analysis even if results already exist (default is False)

        Returns
        -------
        None
        """
        # --- Input validation -------------------------------------------------------------
        if ID is None:
            raise ValueError("Please provide the participant ID (e.g., _.stc(ID='A001')).")
        if i_fname is None:
            raise ValueError("Please provide the filename of the input image.")
        if events_file is None:
            raise ValueError("Please provide the filename of the events TSV file.")
        if run_name is None or run_name=="":
            run_tag=""
        else:
            run_tag="_" + run_name
        # --- Define directories and load files -----------------------------------------------------------
        first_level_dir = self.first_level_dir.format(ID) + task_name + "/"
        os.makedirs(first_level_dir, exist_ok=True)

        df_events = pd.read_csv(events_file, sep="\t") # Load event file
        df_events=df_events.iloc[1:-1] #remove the first raw

        # Load json file
        json_file = os.path.join(self.raw_dir, f"sub-{ID}/func/sub-{ID}_{task_name}{run_tag}_bold.json")
        with open(json_file, "r") as f:
            json_data = json.load(f)
        tr = json_data.get("RepetitionTime")

        # Load fMRI image
        img = nib.load(i_fname)
        n_scans = img.shape[3]
        frame_times = np.arange(n_scans) * tr

        # --- Fit first-level model -----------------------------------------------------------
        design_mat_file = os.path.join(first_level_dir, f"sub-{ID}_{task_name}{run_tag}_design_matrix.png")
        if not os.path.exists(design_mat_file) or redo:
            model = FirstLevelModel(
                t_r=tr,
                noise_model="ar1",
                min_onset=0,
                standardize=False,
                hrf_model="glover + derivative + dispersion",
                drift_model=None,
                signal_scaling=False,
                high_pass=None,
                smoothing_fwhm=smoothing_fwhm,
                mask_img=mask_file
            )

            fmri_glm = model.fit(i_fname, events=df_events)

            # Plot design matrix 
            design_mat = fmri_glm.design_matrices_[0]
            
            
            fig, ax1 = plt.subplots(1, 1, figsize=(6, 4), constrained_layout=True)
            plot_design_matrix(design_mat, ax=ax1)
            ax1.set_title(f"Design Matrix: sub-{ID}, {task_name}", fontsize=12)
            plt.savefig(design_mat_file)

        else:
            if verbose:
                print(f"First-level results already exist for sub-{ID} {task_name} {run_name}. Skipping computation.")
        
        # --- Compute contrasts and save -----------------------------------------------------------
        stat_maps=[]
        for i, contrast in enumerate(contrasts):
            if smoothing_fwhm is not None:
                tag="_s"
            else:
                tag=""
            stat_maps.append(os.path.join(first_level_dir, f"sub-{ID}_{task_name}{run_tag}_{contrast}{tag}.nii.gz"))
            
            if not os.path.exists(stat_maps[i]) or redo:
                results = fmri_glm.compute_contrast(contrast, output_type="z_score")
                results.to_filename(stat_maps[i])
        
        return stat_maps

    def run_second_level_glm(self,i_fnames=None,design_matrix=None,mask_fname=None,smoothing_fwhm=None,task_name=None,run_name=None,contrasts = ["trial_RH-rest", "trial_RH", "rest"],verbose=True,redo=False):

        # ongoing test nilearn: https://nilearn.github.io/stable/modules/generated/nilearn.glm.second_level.SecondLevelModel.html
        # --- Input validation -------------------------------------------------------------
        if i_fnames is None:
            raise ValueError("Please provide the list of filenames of the input contrast images.")
        
        # Load design matrix file if provided, otherwise create a default design matrix with an intercept only
        if design_matrix is None:
            design_matrix = pd.DataFrame([1] * len(i_fnames),columns=["intercept"])
        
        # --- Define directories  -----------------------------------------------------------
        second_level_dir = self.second_level_dir.format(task_name) + "/"
        os.makedirs(second_level_dir, exist_ok=True)
        
        # --- Estimate and Fit second-level model -----------------------------------------------------------
        second_level_model = SecondLevelModel(mask_img=mask_fname,smoothing_fwhm=smoothing_fwhm, n_jobs=2, verbose=1) # define the model to the contrast images and the design matrix
        second_level_model.fit(i_fnames, design_matrix=design_matrix)  # fit the model to the contrast images and the design matrix
        
        # --- Compute contrasts and save -----------------------------------------------------------
        z_map = second_level_model.compute_contrast(second_level_contrast="intercept",output_type="z_score")
        z_map_file = os.path.join(second_level_dir, f"n{len(i_fnames)}_{task_name}_intercept_z_map.nii.gz")
        z_map.to_filename(z_map_file)
        
        return z_map_file