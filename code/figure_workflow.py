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
import matplotlib.pyplot as plt
from IPython.display import Image, display


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
    # ------------------------------------------------------------------

    # On tSNR map in PAM50 space : sub-{}_task-{}_acq-{}_bold_moco_mean_coreg_in_PAM50
    # On tSNR map in Original space : sub-{}_task-{}_acq-{}_bold_moco
    # Todo: Use nn for moco
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

                # Compute tSNR map in native space
                path_tsnr_sub_folder = os.path.join(path_fig_data, f"sub-{ID}", f"task-{task}_acq-{acq_name}")
                print("=== Compute tSNR map ===", flush=True)
                fname_tsnr = compute_tsnr_map(selected_file, path_tsnr_sub_folder, redo, min_vols_for_tsnr)

                # Warp tSNR in PAM50 space
                fname_tsnr_in_template = fname_tsnr.replace("_bold_moco_tSNR.nii.gz", "_bold_moco_tsnr_in_PAM50.nii.gz")
                if not os.path.exists(fname_tsnr_in_template) or redo:
                    print("=== Warp tSNR map to PAM50 space ===", flush=True)

                    fname_warp_from_func_to_template = os.path.join(
                        config["raw_dir"],
                        config["preprocess_dir"]["main_dir"].format(ID),
                        "func",
                        f"task-{task}_acq-{acq_name}",
                        "sct_register_multimodal",
                        os.path.basename(selected_file).replace("_bold_moco.nii.gz", "_from-func_to_PAM50_mode-image_xfm.nii.gz")
                    )

                    if not os.path.exists(fname_warp_from_func_to_template):
                        raise RuntimeError(f"Warp file not found: {fname_warp_from_func_to_template}")
                    fname_template = os.path.join(config["code_dir"], "template", config["PAM50_t2"])
                    cmd_coreg = f"sct_apply_transfo -i {fname_tsnr} -d {fname_template} -w {fname_warp_from_func_to_template} -o {fname_tsnr_in_template} -x nn"
                    os.system(cmd_coreg)

                # Extract metrics from native space
                print("=== Extract tSNR metrics ===", flush=True)
                if fname_tsnr is not None:
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

    fname_tsnr_metrics = os.path.join(path_fig_tsnr, "tsnr_metrics.csv")
    df_tsnr.to_csv(fname_tsnr_metrics, index=False)

    # Todo: Generate figures
    # Create metric violin plot figure
    print("=== Generate tSNR violin plot figure ===", flush=True)
    df_tsnr = pd.read_csv(fname_tsnr_metrics)
    list_baseline_tsnr = []
    list_slicewise_tsnr = []

    name_baseline = [a for a in config["design_exp"]["acq_names"] if a.find("Base") != -1][0]
    name_slicewise = [a for a in config["design_exp"]["acq_names"] if a.find("Slice") != -1][0]

    for ID in IDs:
        df_sub = df_tsnr[df_tsnr["ID"] == int(ID)]
        print(df_sub, flush=True)
        done = False
        # Ordering matters here
        tasks = ["rest", "motor"]
        for task in tasks:
            if "rest" in task and len(df_sub[df_sub["task"] == task]) >= 2:
                done = True
                df_task = df_sub[df_sub["task"] == task]
                if len(df_task) != 2:
                    raise RuntimeError(f"We don't have 2 tSNR metric for sub-{ID} task-{task}")

                tsnr_baseline = df_task[df_task["acq"] == name_baseline]["tsnr_mean"].values
                tsnr_slicewise = df_task[df_task["acq"] == name_slicewise]["tsnr_mean"].values
                list_baseline_tsnr.append(tsnr_baseline[0])
                list_slicewise_tsnr.append(tsnr_slicewise[0])
                continue

            if not done:
                # If no rest task, do something with the motor tasks
                print(f"No rest task found for sub-{ID}, using motor task instead", flush=True)

    fname_tsnr_violin_plot = os.path.join(path_fig_tsnr, "snr_violin_plot.png")
    create_snr_violin_plot(list_baseline_tsnr, list_slicewise_tsnr, fname_tsnr_violin_plot)

    print("=== Figure script End ===", flush=True)


def create_snr_violin_plot(snr_baseline: list, snr_shim: list, fname_output: str):
    # Create a violin plot with matplotolib showing the distribution of SNR values for baseline and shim
    # Link corresponding values for each subject with a line
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(1, 1, 1)
    data = [snr_baseline, snr_shim]
    parts = ax.violinplot(data, showmeans=True, showmedians=True)
    for pc in parts['bodies']:
        pc.set_facecolor('lightblue')
        pc.set_edgecolor('black')
        pc.set_alpha(0.7)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(['Baseline 2nd order shim', 'Slice-wise f0xyz shim'])
    ax.set_ylabel('Mean tSNR for each subject')
    # Link corresponding values for each subject with a line
    for i in range(len(snr_baseline)):
        ax.plot([1, 2], [snr_baseline[i], snr_shim[i]], color='gray', linestyle='--', marker='o')
    fig.tight_layout()
    fig.savefig(fname_output, dpi=1000)



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
