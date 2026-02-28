    #!/usr/bin/env python
# coding: utf-8

# # Spinal cord fMRI denoising 

# @ author of the script: Caroline Landelle, caroline.landelle@mcgill.ca // landelle.caroline@gmail.com
#
# Description: This workflow provides code for first level analyses 
# I. Run first level analysis for each subject and task
# II. Normalize the resulting stat maps to PAM50 template space
#
#------------------------------------------------------------------
#------ Initialization
#------------------------------------------------------------------
# Main imports ------------------------------------------------------------
import re, json, sys, os, glob, argparse
import pandas as pd
from nilearn.glm import threshold_stats_img
import nibabel as nib
import numpy as np
from collections import defaultdict

# Get the environment variable PATH_CODE
path_code = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

with open(path_code + '/config/config_spine_7t_fmri.json') as config_file: # the notebook should be in 'xx/notebook/' folder #config_proprio
    config = json.load(config_file) # load config file should be open first and the path inside modified

parser = argparse.ArgumentParser()
parser.add_argument("--ids", nargs='+', default=[""])
parser.add_argument("--tasks", nargs='+', default=[""])
parser.add_argument("--verbose", default="False")
parser.add_argument("--redo", default="True")
parser.add_argument("--path-data", required=True)
args = parser.parse_args()

IDs = args.ids
tasks = args.tasks
verbose = args.verbose.lower() == "true"
redo = args.redo.lower() == "true"
path_data = os.path.abspath(args.path_data)

config["raw_dir"]=path_data
config["code_dir"]=path_code

participants_tsv = pd.read_csv(path_code + '/config/participants.tsv', sep='\t',dtype={'participant_id': str})

new_IDs=[]
if IDs == [""]:
    for ID in participants_tsv["participant_id"]:
        new_IDs.append(ID)

    IDs = new_IDs

if tasks != [""]:
    config["design_exp"]["task_names"] = tasks

#Import scripts
sys.path.append(path_code + "/code/") # Change this line according to your directory
import postprocess, preprocess

postprocess=postprocess.Postprocess_main(config,IDs=IDs)
preprocess_Sc=preprocess.Preprocess_Sc(config,IDs=IDs)

# initialize directories
preprocessing_dir = os.path.join(config["raw_dir"], config["preprocess_dir"]["main_dir"])
denoising_dir= os.path.join(config["raw_dir"], config["denoising"]["dir"])
manual_dir = os.path.join(config["raw_dir"], config["manual_dir"])
main_fig_dir = os.path.join(config["raw_dir"], "derivatives/processing/figures/")#os.path.join(config["raw_dir"], config["figures_dir"]["main_dir"])
fig_task_dir = os.path.join(main_fig_dir, "task")
os.makedirs(main_fig_dir, exist_ok=True)
os.makedirs(fig_task_dir, exist_ok=True)

#------------------------------------------------------------------
#------ First level
#------------------------------------------------------------------
print("")
print("=== First level analysis script Start ===", flush=True)
print("Participant(s) included : ", IDs, flush=True)
print("===================================", flush=True)
print("")

df_task = pd.DataFrame(columns=["ID", "task", "acq", "run", "contrast", "max_zscore", "active_voxels"])
fname_task_metrics = os.path.join(fig_task_dir, "first_level_metrics_n" + str(len(IDs)) + ".csv")
norm_mask=[]
for ID_nb, ID in enumerate(IDs):
    print("", flush=True)
    print(f'=== First level start for :  {ID} ===', flush=True)

    for task_name in config["design_exp"]["task_names"]:
        for acq_name in config["design_exp"]["acq_names"]:
            tag="task-" + task_name + "_acq-" + acq_name
            raw_func=glob.glob(os.path.join(config["raw_dir"], f'sub-{ID}', 'func', f'sub-{ID}_{tag}_*bold.nii.gz'))
            for func_file in raw_func:
                # Check run number if multiple run exists
                match = re.search(r"_?(run-\d+)", func_file)
                if match:
                    run_name=match.group(1)
                    print(run_name)
                else:
                    run_name=""

                denoised_fmri=glob.glob(os.path.join(denoising_dir.format(ID ), tag, config["denoising"]["denoised_dir"],"*"+run_name+"*_nostd_s.nii.gz"))[0]
                # Select manual seg if exists
                mask_file_list = glob.glob(os.path.join(preprocessing_dir.format(ID), 'func',tag, config["preprocess_f"]["func_seg"].format(ID,tag,"")))
                mask_file = mask_file_list[0] if len(mask_file_list) > 0 else None
                manual_seg_file_list = glob.glob(os.path.join(f"{manual_dir}", f"sub-{ID}", "func", config["preprocess_f"]["func_seg"].format(ID,tag,run_name)))
                manual_seg_file = manual_seg_file_list[0] if len(manual_seg_file_list) > 0 else ""
                mask_final_file = manual_seg_file if os.path.exists(manual_seg_file) else mask_file
                if mask_final_file is None:
                    raise RuntimeError(f"No mask file found for subject {ID}, task {tag}, run {run_name}. Please check the preprocessing outputs and manual corrections.")

                # Select segmentation
                mask_final_file = os.path.join(preprocessing_dir.format(ID), 'func', tag, f"sub-{ID}_{tag}_bold_moco_mean_seg.nii.gz")
                if not os.path.exists(mask_final_file):
                    raise RuntimeError(f"No mask file found for subject {ID}, task {tag}. Please check the preprocessing outputs and manual corrections.")

                # Select warp file
                warp_file = os.path.join(preprocessing_dir.format(ID), 'func', tag, f"sub-{ID}_{tag}_from-func_to_PAM50_mode-image_xfm.nii.gz")
                
                
                if not os.path.exists(warp_file):
                    raise RuntimeError(f"No warp file found for subject {ID}, task {tag}. Please check the preprocessing outputs and manual corrections.")

                events_file=glob.glob(os.path.join(config["raw_dir"], f'sub-{ID}', 'func', f'sub-{ID}_{tag}_*events.tsv'))[0]

                #------------------------------------------------------------------
                #------ I. Run first level GLM
                stat_maps=postprocess.run_first_level_glm(ID=ID,
                                                          i_fname=denoised_fmri,
                                                          events_file=events_file,
                                                          mask_file=mask_final_file,
                                                          task_name=tag,
                                                          run_name=run_name,
                                                          redo=redo,
                                                          verbose=verbose)

                #------ II. Apply correction and extract metrics
                for i, contrast_fname in enumerate(stat_maps):
                    # Apply correction
                    corr_type="fpr";alpha=0.01;cluster=0
                    
                    fname_thr_img=stat_maps[i].split(".")[0] +f"_{corr_type}_{str(alpha)[2:]}_{str(cluster)}cluster.nii.gz"
                    
                    if not os.path.exists(fname_thr_img) or redo:
                        thresholded_map, threshold = threshold_stats_img(stat_maps[i],
                                                                        alpha=alpha,
                                                                        height_control=corr_type,
                                                                        cluster_threshold=cluster,
                                                                            two_sided=False)
                        thresholded_map.to_filename(fname_thr_img)

                    # extract max z-score and total number of active voxels in the spinal cord mask
                    if "RH-rest" in contrast_fname:
                        if run_name != "run-02":
                            if not os.path.exists(fname_task_metrics) or redo:
                                img = nib.load(fname_thr_img)
                                mask = nib.load(mask_final_file)
                                masked_data = img.get_fdata() * mask.get_fdata()# Extract values within the mask
                                max_zscore = np.max(masked_data[masked_data > 0]) if np.any(masked_data > 0) else 0 # Find max z-score in the mask
                                active_voxels = np.sum(masked_data > 0)# Count active voxels in the mask

                                #store in a dataframe for each individual, task_name,acq_name, run_name and contrast
                                if len(df_task) == 0:
                                    df_task = pd.DataFrame([[ID,task_name,acq_name,run_name,"RH-rest",max_zscore,active_voxels]], columns=df_task.columns)
                                
                                df_task = pd.concat([df_task, pd.DataFrame([[ID,task_name,acq_name,run_name,"RH-rest",max_zscore,active_voxels]], columns=df_task.columns)], ignore_index=True)
                                
                
                #------ III. Normalize the resulting stat maps to PAM50 template space

                for i, contrast_fname in enumerate(stat_maps):
                    norm_stat_maps=preprocess_Sc.apply_warp(
                            i_img=[stat_maps[i]], # input clean image
                            ID=[ID],
                            o_folder=[os.path.dirname(stat_maps[i]) + "/"], # output folder
                            dest_img=os.path.join(path_code, "template", config["PAM50_t2"]), # PAM50 template
                            warping_field=warp_file,
                            tag="_inTemplate",
                            mean=False,
                            n_jobs=1,
                            verbose=False,
                            redo=redo)
                
                # Normalize the individual masks to template space
                norm_mask.append(preprocess_Sc.apply_warp(
                            i_img=[mask_final_file], # input clean image
                            ID=[ID],
                            o_folder=[os.path.dirname(stat_maps[i]) + "/"], # output folder
                            dest_img=os.path.join(path_code, "template", config["PAM50_t2"]), # PAM50 template
                            warping_field=warp_file,
                            tag="_inTemplate",
                            mean=False,
                            n_jobs=1,
                            threshold=0.1,
                            verbose=False,
                            redo=redo)[0])
    
    print(f'=== First level done for : {ID} ===', flush=True)
    print("=========================================", flush=True)

#------ IV. Extract the commun mask for all participants and tasks
first_level_dir = os.path.join(config["raw_dir"], config["first_level"]["dir"])
common_mask_fname = os.path.join(first_level_dir.split("sub")[0], "common_mask_PAM50.nii.gz")
cropped_PAM50_fname = os.path.join(first_level_dir.split("sub")[0], "PAM50_cord_cropped.nii.gz")

if not os.path.exists(cropped_PAM50_fname) or redo:
    norm_mask_data = [nib.as_closest_canonical(nib.load(f)).get_fdata() for f in norm_mask]
    n_files = len(norm_mask_data)

    # Compute common mask (n-1)---
    sum_mask = np.sum(norm_mask_data, axis=0)
    common_mask_data = (sum_mask >= n_files-3).astype(np.uint8)
    common_mask_fname = os.path.join(first_level_dir.split("sub")[0], "common_mask_PAM50.nii.gz")
    common_mask_img = nib.Nifti1Image(common_mask_data, affine=nib.load(norm_mask[0]).affine)
    common_mask_img.to_filename(common_mask_fname)
    common_mask_data = common_mask_img.get_fdata()

    # ---  Extract the z-slices that contain the common mask ---
    z_indices = np.where(np.any(common_mask_data > 0, axis=(0,1)))[0]
    z_min, z_max = z_indices[[0, -1]]
    z_size = z_max - z_min + 1

    #Crop the PAM50 template to the common mask z-slices
    
    pam50_fname = os.path.join(path_code, "template", config["PAM50_cord"])
    cmd = f"fslroi {pam50_fname} {cropped_PAM50_fname} 0 -1 0 -1 {z_min} {z_size}"
    os.system(cmd)

    print(f"Cropped PAM50 saved: {cropped_PAM50_fname}, slices z={z_min}-{z_max}")

#------ V. Plot first level results for each task and participant

# --- Select stat map files ---
i_fnames=[]
for task_name in config["design_exp"]["task_names"]:
    for acq_name in config["design_exp"]["acq_names"]:
        tag="task-" + task_name + "_acq-" + acq_name
        for ID in IDs:
            # define the run name if multiple runs exist    
            raw_func=sorted(glob.glob(os.path.join(config["raw_dir"], f'sub-{ID}', 'func', f'sub-{ID}_{tag}_*bold.nii.gz')))
            func_file = raw_func[0]# take only the first run
            match = re.search(r"_?(run-\d+)", func_file)
            run_name = match.group(1) if match else ""  # extract run number if exists
            i_fnames.append(glob.glob(os.path.join(first_level_dir.format(ID), f"{tag}", f"*{tag}*{run_name}*trial_RH-rest*inTemplate.nii.gz"))[0])

# --- Group by participant ID ---
subject_files = defaultdict(list)
for f in i_fnames:
    sub_id = os.path.basename(f).split("_")[0]
    subject_files[sub_id].append(f)

# --- Keep only first two maps per subject as pairs ---
i_fnames_pairs = []
for sub_id in sorted(subject_files.keys()):
    pair = subject_files[sub_id][:2]  # take first two files
    if len(pair) == 2:
        i_fnames_pairs.append(pair)

print(i_fnames_pairs)
output_dir=os.path.join(config["raw_dir"], config["figures_dir"]["main_dir"], "task")
postprocess.plot_first_level_maps(i_fnames_pair=i_fnames_pairs,
                                          output_dir=output_dir,
                                          background_fname=os.path.join(path_code, "template", config["PAM50_t2"]),
                                          #underlay_fname=os.path.join(path_code, "template", config["PAM50_cord"]),
                                          task_name=tag,
                                          verbose=True,
                                          redo=True)

#if not os.path.exists(fname_task_metrics) or redo:
 #   df_task.to_csv(fname_task_metrics, index=False)
  #  print(f"Task metrics saved to: {fname_task_metrics}")


                

#------------------------------------------------------------------
#------ Second level
#------------------------------------------------------------------
#TODO: add a commun mask for secon level analyses

print("")
print("=== Second level analysis script Start ===", flush=True)
print("Participant(s) included : ", IDs, flush=True)
print("===================================", flush=True)
print("")

# list of first evel contrast images in template space for each participant and task
second_level=False
if second_level==True:
    for task_name in config["design_exp"]["task_names"]:
        for acq_name in config["design_exp"]["acq_names"]:
            i_fnames=[]
            tag="task-" + task_name + "_acq-" + acq_name
            for ID in IDs:
                # define the run name if multiple runs exist
                raw_func=sorted(glob.glob(os.path.join(config["raw_dir"], f'sub-{ID}', 'func', f'sub-{ID}_{tag}_*bold.nii.gz')))

                # take only the first run
                func_file = raw_func[0]

                # extract run number if exists
                match = re.search(r"_?(run-\d+)", func_file)
                run_name = match.group(1) if match else ""
                

                # find the corresponding first-level file
                
                i_fnames.append(glob.glob(os.path.join(first_level_dir.format(ID), f"{tag}", f"*{tag}*{run_name}*trial_RH-rest*inTemplate.nii.gz"))[0])

            print(i_fnames)
            
            z_map_file=postprocess.run_second_level_glm(i_fnames=i_fnames,
                                                            mask_fname=f"{path_code}/template/{config['PAM50_cord']}",
                                                            task_name=tag,
                                                            run_name="",
                                                            redo=True,
                                                            verbose=verbose)
        
        #corr_type="fpr";alpha=0.001;cluster=50
                    
        #fname_thr_img=z_map_file.split(".")[0] +f"_{corr_type}_{str(alpha)[2:]}_{str(cluster)}cluster.nii.gz"
        #print(alpha)      
        #if not os.path.exists(fname_thr_img) or redo:
         #   thresholded_map, threshold = threshold_stats_img(z_map_file,
               #                                              alpha=alpha,
              #                                               height_control=corr_type,
             #                                                cluster_threshold=cluster,
            #                                                 parametric=False,
           #                                                  two_sided=redo)
          #  thresholded_map.to_filename(fname_thr_img)