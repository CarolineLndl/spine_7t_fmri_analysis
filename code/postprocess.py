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
from nilearn.glm.second_level import non_parametric_inference

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
    
    def plot_first_level_maps(self, i_fnames_pair=None, output_dir=None,stat_min=2.3, stat_max=5,background_fname=None, underlay_fname=None,task_name=None, verbose=True, redo=False):
    
        if output_dir is None:
            output_dir = os.path.join(self.first_level_dir)
        if i_fnames_pair is None or len(i_fnames_pair) == 0:
            raise ValueError("i_fnames_pair is empty")

        # --- Load PAM50 template and optional underlay --------------------------
        template_data = nib.as_closest_canonical(nib.load(background_fname)).get_fdata()
        underlay_data = None
        if underlay_fname is not None:
            underlay_data = nib.as_closest_canonical(nib.load(underlay_fname)).get_fdata()

        for subject_idx, maps_pair in enumerate(i_fnames_pair):
            if len(maps_pair) != 2:
                raise ValueError("Each subject should have exactly 2 statistical maps")

            # --- Create a 2x2 grid for this subject -------------------------------
            fig, axs = plt.subplots(2, 2, figsize=(8, 8), constrained_layout=True)
            axs = axs.flatten()  # 0,1 = coronal; 2,3 = axial

            for i, i_fname in enumerate(maps_pair):
                statmap_data = nib.as_closest_canonical(nib.load(i_fname)).get_fdata()
                stat_thresh = np.where(statmap_data > stat_min, statmap_data, 0)

                # Coronal (top row)
                mip_cor = np.max(stat_thresh, axis=1).T
                mip_cor = np.where(mip_cor > stat_min, mip_cor, np.nan)
                y_slice = statmap_data.shape[1] // 2
                template_cor = template_data[:, y_slice, :].T

                ax_cor = axs[i]
                ax_cor.imshow(template_cor, cmap="gray", origin="lower")
                if underlay_data is not None:
                    ax_cor.imshow(underlay_data[:, y_slice, :].T, cmap="gray", origin="lower")
                ax_cor.imshow(mip_cor, cmap="hot", origin="lower", vmin=stat_min, vmax=stat_max)
                ax_cor.axvline(x=y_slice, color="white", linestyle="--", linewidth=0.8, alpha=0.6)
                ax_cor.axis("off")
                if subject_idx == 0:
                    ax_cor.text(0.05, 0.05, "L", transform=ax_cor.transAxes, color="white", fontsize=10, ha="left", va="bottom")
                    ax_cor.text(0.95, 0.05, "R", transform=ax_cor.transAxes, color="white", fontsize=10, ha="right", va="bottom")

                # Axial (bottom row)
                mip_axi = np.max(stat_thresh, axis=2).T
                mip_axi = np.where(mip_axi > stat_min, mip_axi, np.nan)
                z_slice = statmap_data.shape[2] // 2
                template_axi = template_data[:, :, z_slice].T

                ax_axi = axs[i+2]
                ax_axi.imshow(template_axi, cmap="gray", origin="lower")
                if underlay_data is not None:
                    ax_axi.imshow(underlay_data[:, :, z_slice].T, cmap="gray", origin="lower")
                ax_axi.imshow(mip_axi, cmap="hot", origin="lower", vmin=stat_min, vmax=stat_max)
                ax_axi.axis("off")
                if subject_idx == 0:
                    ax_axi.text(0.02, 0.5, "L", transform=ax_axi.transAxes, color="white", fontsize=10, ha="left", va="center")
                    ax_axi.text(0.98, 0.5, "R", transform=ax_axi.transAxes, color="white", fontsize=10, ha="right", va="center")
                    ax_axi.text(0.5, 0.98, "A", transform=ax_axi.transAxes, color="white", fontsize=10, ha="center", va="top")
                    ax_axi.text(0.5, 0.02, "P", transform=ax_axi.transAxes, color="white", fontsize=10, ha="center", va="bottom")

            # --- Save figure for this subject -------------------------------------
            out_file = os.path.join(output_dir, f"first_level_maps_{task_name}_sub-{subject_idx+1}.png")
            fig.savefig(out_file, dpi=300)
            plt.close(fig)
            

        # Plot subject number and task name on the figure
        # add one scale bar for all the maps
        # plot R/L and S/I orientation on the first plot
        # save in the right folder


    def plot_first_level_maps1(self,i_fnames=None, output_dir=None,stat_min=2.3,stat_max=5, background_fname=None,underlay_fname=None,task_name=None,verbose=True,redo=False):
        if output_dir==None:
            output_dir = os.path.join(self.first_level_dir, task_name + "/")
        if i_fnames is None or len(i_fnames) == 0:
            raise ValueError("i_fnames is empty")
        
        # --- Figure layout -----------------------------------------------------------
        maps_per_row=6 # number of columns to plot the maps
        col_nb_per_map = 3  # sagittal + coronal + axial
        gap_cols = 1   # one empty column between maps
        total_cols = maps_per_row * (col_nb_per_map + gap_cols)  # total axes per row
        row_nb=int(np.ceil(len(i_fnames)/maps_per_row)) # number of rows needed to plot all the maps with 6 columns
        fig, axes = plt.subplots(row_nb, total_cols, figsize=(total_cols, row_nb*3), constrained_layout=True)
        plt.subplots_adjust(wspace=0.1, hspace=0.1)  # smaller spacing
        axes = np.atleast_1d(axes).flatten()# Make axes iterable in all cases

        # --- Load files -----------------------------------------------------------
        template_data = nib.as_closest_canonical(nib.load(background_fname)).get_fdata()
        underlay_data = None
        if underlay_fname is not None:
            underlay_data = nib.as_closest_canonical(nib.load(underlay_fname)).get_fdata()
        
       # --- Plot maps -------------------------------------------------------------
        for i, i_fname in enumerate(i_fnames):
            idx = i * (col_nb_per_map + gap_cols)  # skip over gap
            ax_sag = axes[idx]
            ax_cor = axes[idx + 1]
            ax_axi = axes[idx + 2]

            statmap_data = nib.as_closest_canonical(nib.load(i_fname)).get_fdata()
            stat_thresh = np.where(statmap_data > stat_min, statmap_data, 0)

            # Maximal Ibtensity Projection along the x-axis → sagittal view
            mip_sag = np.max(stat_thresh, axis=0).T  # shape (Y, Z) → transpose if needed
            mip_sag = np.where(mip_sag > stat_min, mip_sag, np.nan)  # mask low values
            mip_cor = np.max(stat_thresh, axis=1).T  # shape (X, Z) → transpose if needed
            mip_cor = np.where(mip_cor > stat_min, mip_cor, np.nan)  # mask low values
            mip_axi = np.max(stat_thresh, axis=2).T  # shape (X, Y) → transpose if needed
            mip_axi = np.where(mip_axi > stat_min, mip_axi, np.nan)  # mask low values

            # Sagittal view (X fixed, Y–Z plane)
            x_slice = statmap_data.shape[0] // 2
            template_sag = template_data[x_slice, :, :].T
            ax_sag.imshow(template_sag, cmap="gray", origin="lower")
            if underlay_data is not None:
                ax_sag.imshow(underlay_data[x_slice, :, :].T, cmap="gray", alpha=0.3, origin="lower")
            ax_sag.imshow(mip_sag, cmap="hot", origin="lower", vmin=stat_min, vmax=stat_max)
            ax_sag.axis("off")
            # Add R/L S/I labels only for the first map
            if i == 0:
                ax_sag.text(0.05, 0.05, "P",transform=ax_sag.transAxes,color="white", fontsize=12,ha="left", va="bottom")
                ax_sag.text(0.95, 0.05, "A",transform=ax_sag.transAxes,color="white", fontsize=12,ha="right", va="bottom")


            # Coronal view (Y fixed, X–Z plane)
            y_slice = statmap_data.shape[1] // 2
            template_cor = template_data[:, y_slice, :].T
            ax_cor.imshow(template_cor, cmap="gray", origin="lower")
            if underlay_data is not None:
                ax_cor.imshow(underlay_data[:, y_slice, :].T, cmap="gray", alpha=0.3, origin="lower")

            ax_cor.axvline(x=y_slice,color="white",linestyle="--",linewidth=0.8,alpha=0.6)
            ax_cor.imshow(mip_cor, cmap="hot", origin="lower", vmin=stat_min, vmax=stat_max)
            ax_cor.axis("off")
            if i == 0:
                ax_cor.text(0.05, 0.05, "L",transform=ax_cor.transAxes,color="white", fontsize=12,ha="left", va="bottom")
                ax_cor.text(0.95, 0.05, "R",transform=ax_cor.transAxes,color="white", fontsize=12,ha="right", va="bottom")

            # Axial view (Z fixed, X–Y plane)
            z_slice = statmap_data.shape[2] // 2
            axi_slice = statmap_data[:, :, z_slice].T
            template_axi = template_data[:, :, z_slice].T
            ax_axi.imshow(template_axi, cmap="gray", origin="lower")
            if underlay_data is not None:
                ax_axi.imshow(underlay_data[:, :, z_slice].T, cmap="gray", alpha=0.3, origin="lower")
            
            ax_axi.imshow(mip_axi, cmap="hot", origin="lower", vmin=stat_min, vmax=stat_max)
            ax_axi.axis("off")
            if i == 0:
                ax_axi.text(0.02, 0.50, "L",transform=ax_axi.transAxes,color="white", fontsize=12,ha="left", va="bottom")
                ax_axi.text(0.98, 0.50, "R",transform=ax_axi.transAxes,color="white", fontsize=12,ha="right", va="bottom")
                ax_axi.text(0.5, 0.98, "A", transform=ax_axi.transAxes,color="white", fontsize=12,ha="left", va="top")
                ax_axi.text(0.5, 0.02, "P",transform=ax_axi.transAxes,color="white", fontsize=12,ha="left", va="bottom")

            # --- Hide unused axes (including gap columns) ---
            for g in range(1, gap_cols+1):
                gap_idx = idx + col_nb_per_map + (g-1)
                if gap_idx < len(axes):
                    axes[gap_idx].axis("off")
            for j in range(len(i_fnames) * (col_nb_per_map + gap_cols), len(axes)):
                axes[j].axis("off")


        out_file = os.path.join(output_dir,f"first_level_maps_{task_name}.png")
        fig.savefig(out_file, dpi=300)
        plt.close(fig)
            


        # Plot subject number and task name on the figure
        # add one scale bar for all the maps
        # plot R/L and S/I orientation on the first plot
        # save in the right folder


    def run_second_level_glm(self,i_fnames=None,design_matrix=None,mask_fname=None,smoothing_fwhm=None,task_name=None,parametric=False,run_name=None,verbose=True,redo=False):
        '''
        Run second-level GLM for a specific task.
        # ongoing test nilearn: https://nilearn.github.io/stable/modules/generated/nilearn.glm.second_level.SecondLevelModel.html

        Parameters
        ----------
        i_fnames : list of str
            List of filenames of the input contrast images in the same space (e.g., ["sub-A001_task-motor_contrast-trial_RH-rest_inTemplate.nii.gz", "sub-A002_task-motor_contrast-trial_RH-rest_inTemplate.nii.gz"])
        design_matrix : pandas DataFrame, optional
            Design matrix for the second-level analysis (default is None, which will create a design matrix
            with an intercept only)
        mask_fname : str, optional
            Filename of the mask NIfTI file where to restrict the analysis (default is None)
        smoothing_fwhm : float, optional
            Full-width at half-maximum for spatial smoothing (default is None, which means no smoothing)
        task_name : str
            Task name (e.g., "motor_acq-shimBase+3mm")
        run_name : str, optional
            Run name (e.g., "run-1") (default is None, which means no run name will be added to the output filename)
        verbose : bool, optional
            Whether to print information during processing (default is True)
        redo : bool, optional
            Whether to redo the analysis even if results already exist (default is False)
        Returns
        -------
        z_map_file : str
            Filename of the output z-map NIfTI file (e.g., "n20_motor_acq-shimBase+3mm_intercept_z_map.nii.gz")
        '''
        


        # --- Input validation -------------------------------------------------------------
        if i_fnames is None:
            raise ValueError("Please provide the list of filenames of the input contrast images.")
        
        # --- Define directories  -----------------------------------------------------------
        second_level_dir = self.second_level_dir.format(task_name) + "/"
        os.makedirs(second_level_dir, exist_ok=True)

        # Load design matrix file if provided, otherwise create a default design matrix with an intercept only
        if design_matrix is None:
            design_matrix = pd.DataFrame([1] * len(i_fnames),columns=["intercept"])

        if parametric ==True:
            stat_map_file = os.path.join(second_level_dir, f"n{len(i_fnames)}_{task_name}_intercept_z_map.nii.gz")
            if not os.path.exists(stat_map_file) or redo:
                if verbose:
                    print(f"Computing parametric second-level analysis for task {task_name}.")
                
                
                
                # --- Estimate and Fit second-level model -----------------------------------------------------------
                second_level_model = SecondLevelModel(mask_img=mask_fname,smoothing_fwhm=smoothing_fwhm, n_jobs=2, verbose=1) # define the model to the contrast images and the design matrix
                second_level_model.fit(i_fnames, design_matrix=design_matrix)  # fit the model to the contrast images and the design matrix
                
                # --- Compute contrasts and save -----------------------------------------------------------
                z_map = second_level_model.compute_contrast(second_level_contrast="intercept",output_type="z_score")
                z_map.to_filename(stat_map_file)
        
        else:
            stat_map_file = os.path.join(second_level_dir, f"n{len(i_fnames)}_{task_name}_")
            if not os.path.exists(stat_map_file + 'logp_max_t.nii.gz') or redo:
                if verbose:
                    print(f"Computing non-parametric second-level analysis for task {task_name}.")
            
                out_dict = non_parametric_inference(
                    i_fnames,
                    design_matrix=design_matrix,
                    mask=mask_fname,
                    model_intercept=True,
                    n_perm=10000,  # 500 for the sake of time. Ideally, this should be 10,000.
                    two_sided_test=False,
                    smoothing_fwhm=1,
                    n_jobs=2,
                    threshold=0.01, # voxel level threshold for cluster definition (uncorrected p-value)
                    #tfce=True,
                    verbose=1,
                    )
                
                print(out_dict)
                out_dict["t"].to_filename(stat_map_file+ 't.nii.gz')
                out_dict["logp_max_t"].to_filename(stat_map_file+ 'logp_max_t.nii.gz')
                out_dict["logp_max_size"].to_filename(stat_map_file+ 'logp_max_size.nii.gz')
                out_dict["logp_max_mass"].to_filename(stat_map_file+ 'logp_max_mass.nii.gz')
                #out_dict["tfce"].to_filename(stat_map_file+ 'tfce.nii.gz')
                #out_dict["logp_max_tfce"].to_filename(stat_map_file+ 'logp_max_tfce.nii.gz')

                #mask the t-map with the significant cluster in the logp_max_size map
                logp_max_size_img = nib.load(stat_map_file+ 'logp_max_size.nii.gz')
                logp_max_size_data = logp_max_size_img.get_fdata()
                #threshold the logp_max_size map at p<0.05
                logp_max_size_data_thresholded = logp_max_size_data > -np.log10(0.05)
                #mask the t-map with the thresholded logp_max_size map
                t_img = nib.load(stat_map_file+ 't.nii.gz')
                t_data = t_img.get_fdata()
                t_data_masked = t_data * logp_max_size_data_thresholded
                t_masked_img = nib.Nifti1Image(t_data_masked, t_img.affine, t_img.header)
                t_masked_img.to_filename(stat_map_file+ 't_clustercorrected.nii.gz')

        
        return stat_map_file
    
  