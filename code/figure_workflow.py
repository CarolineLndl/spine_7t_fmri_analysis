#!/usr/bin/env python
# coding: utf-8

#  Spinal Cord fMRI figures
# ____________________________________________________

# ### Project: acdc_spine_7T
# ____________________________________________________

#------------------------------------------------------------------
#------ Initialization
#------------------------------------------------------------------
# Imports
import sys, json, glob, os, re, shutil, argparse

import nibabel as nib
import pandas as pd
from utils import compute_tsnr_map, extract_tsnr_metric
from IPython.display import Image, display
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ids", nargs='+', default=[""])
    parser.add_argument("--tasks", nargs='+', default=[""])
    parser.add_argument("--verbose", default="False")
    parser.add_argument("--redo", default="True")
    parser.add_argument("--path-data", required=True)
    args = parser.parse_args()

    # get path of the parent location of this file, and go up one level
    path_code = os.path.dirname(os.path.abspath(__file__)).rsplit('/', 1)[0]
    fname_config = os.path.join(path_code, 'config', 'config_spine_7t_fmri.json')
    with open(fname_config) as config_file:
        config = json.load(config_file)  # load config file should be open first and the path inside modified

    IDs = args.ids
    tasks = args.tasks
    verbose = args.verbose.lower() == "true"
    redo = args.redo.lower() == "true"
    path_data = args.path_data

    config["raw_dir"] = path_data
    config["code_dir"] = path_code

    # Load participants info
    participants_tsv = pd.read_csv(os.path.join(path_code, 'config', 'participants.tsv'), sep='\t',dtype={'participant_id': str})

    new_IDs=[]
    if IDs == [""]:
        for ID in participants_tsv["participant_id"]:
            new_IDs.append(ID)
        IDs = new_IDs

    if tasks != [""]:
        config["design_exp"]["task_names"] = tasks


    #------------------------------------------------------------------
    #------ Figures
    #------------------------------------------------------------------
    print("=== Figure script Start ===", flush=True)
    print("Participant(s) included : ", IDs, flush=True)
    print("===================================", flush=True)
    print("")

    path_main_fig = os.path.join(config["raw_dir"], config["figures_dir"]["main_dir"])
    path_fig_tsnr = os.path.join(path_main_fig, "tsnr")
    path_fig_data = os.path.join(path_fig_tsnr, "data")

    #------------------------------------------------------------------
    #------ Compute tSNR
    # On tSNR map in PAM50 space : sub-{}_task-{}_acq-{}_bold_moco_mean_coreg_in_PAM50
    # On tSNR map in Original space : sub-{}_task-{}_acq-{}_bold_moco
    # Todo: Use moco with nearest neighbour
    # Use the run with the most volumes
    # Use the same number of volumes for each tsnr calculation
    # Todo: refactor into figure.py file
    #------------------------------------------------------------------
    print("=== Compute tSNR map on longest moco nearest neighbour run ===", flush=True)
    df_tsnr = pd.DataFrame(columns=["ID", "task", "acq", "tsnr_mean"])

    # Find the minimum number of volumes across all runs to standardize tSNR calculation
    min_vols_for_tsnr = 1000
    for ID in IDs:
        for task in config["design_exp"]["task_names"]:
            for acq_name in config["design_exp"]["acq_names"]:
                selected_file = find_moco_for_tsnr_calculation(config, ID, task, acq_name)
                if selected_file is None:
                    continue
                n_vols = nib.load(selected_file).shape[3]
                if n_vols < min_vols_for_tsnr:
                    min_vols_for_tsnr = n_vols

    print(f"Minimum number of volumes across all runs: {min_vols_for_tsnr}", flush=True)
    # Minimum number of volumes across all runs: 30 (2026-01-28)
    # Compute_tsnr
    for ID in IDs:
        for task in config["design_exp"]["task_names"]:
            for acq_name in config["design_exp"]["acq_names"]:
                selected_file = find_moco_for_tsnr_calculation(config, ID, task, acq_name)
                if selected_file is None:
                    continue

                path_tsnr_sub_folder = os.path.join(path_fig_data, f"sub-{ID}", f"task-{task}_acq-{acq_name}")
                print("=== Compute tSNR map ===", flush=True)
                fname_tsnr = compute_tsnr_map(selected_file, path_tsnr_sub_folder, redo, min_vols_for_tsnr)

                # Extract metrics
                print("=== Extract tSNR metrics ===", flush=True)
                if fname_tsnr is not None:
                    # Todo use nn mask moco
                    fname_mask = os.path.join(
                        config["raw_dir"],
                        config["preprocess_dir"]["main_dir"].format(ID),
                        "func",
                        f"task-{task}_acq-{acq_name}",
                        "sct_deepseg",
                        os.path.basename(selected_file).replace("_bold_moco.nii.gz", "_bold_moco_mean_seg.nii.gz")
                    )

                    if not os.path.exists(fname_mask):
                        raise RuntimeError(f"Mask file not found: {fname_mask}")

                    tsnr_mean = extract_tsnr_metric(fname_tsnr, fname_mask)
                    print(f"sub-{ID} task-{task} acq-{acq_name} tsnr_mean: {tsnr_mean}", flush=True)
                    df_tsnr = pd.concat([pd.DataFrame([[ID, task, acq_name, tsnr_mean]], columns=df_tsnr.columns), df_tsnr], ignore_index=True)

    df_tsnr.to_csv(os.path.join(path_fig_tsnr, "tsnr_metrics.csv"), index=False)

    # Todo: Compute tSNR map in PAM50 space

    # Todo: Generate figures

    print("=== Figure script End ===", flush=True)


def find_moco_for_tsnr_calculation(config, ID, task, acq_name):
    files = glob.glob(os.path.join(
        config["raw_dir"],
        config["preprocess_dir"]["main_dir"].format(ID),
        "func",
        f"task-{task}_acq-{acq_name}",
        "sct_fmri_moco",
        f"sub-{ID}_task-{task}_acq-{acq_name}*_bold_moco.nii.gz"
    ))
    if len(files) == 0:
        return None
    elif len(files) == 1:
        selected_file = files[0]
    else:
        max_volumes = 0
        selected_file = None
        for f in files:
            img = nib.load(f)
            n_volumes = img.shape[3]
            if n_volumes > max_volumes:
                max_volumes = n_volumes
                selected_file = f
    return selected_file


if __name__ == "__main__":
    main()
