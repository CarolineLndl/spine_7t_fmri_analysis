#!/usr/bin/env python
# coding: utf-8

#  Spinal Cord fMRI figures
# ____________________________________________________

# ### Project: acdc_spine_7T
# ____________________________________________________

# ------------------------------------------------------------------
# ------ Initialization
# ------------------------------------------------------------------
# Imports
import sys, json, glob, os, re, shutil, argparse

import nibabel as nib
import numpy as np
import pandas as pd
from utils import compute_tsnr_map, extract_mean_within_mask
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from nilearn.image import smooth_img
import warnings
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
    path_code = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
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

    for task_available in config['design_exp']['task_names']:
        if task_available not in config['design_exp']['task_names']:
            raise ValueError(f"All tasks need to be run to generate the figure")

    # Load participants info
    participants_tsv = pd.read_csv(os.path.join(path_code, 'config', 'participants.tsv'), sep='\t',
                                   dtype={'participant_id': str})

    new_IDs = []
    if IDs == [""]:
        for ID in participants_tsv["participant_id"]:
            new_IDs.append(ID)
        IDs = new_IDs

    if tasks != [""]:
        config["design_exp"]["task_names"] = tasks

    # ------------------------------------------------------------------
    # ------ Figures
    # ------------------------------------------------------------------
    print("=== Figure script Start ===", flush=True)
    print("Participant(s) included : ", IDs, flush=True)
    print("===================================", flush=True)
    print("")

    path_main_fig = os.path.join(config["raw_dir"], config["figures_dir"]["main_dir"])
    path_fig_tsnr = os.path.join(path_main_fig, "tsnr")
    path_fig_data = os.path.join(path_fig_tsnr, "data")

    # ------------------------------------------------------------------
    # ------ Compute tSNR
    # ------------------------------------------------------------------

    # On tSNR map in PAM50 space : sub-{}_task-{}_acq-{}_bold_moco_mean_coreg_in_PAM50
    # On tSNR map in Original space : sub-{}_task-{}_acq-{}_bold_moco
    # Todo: Use nn for moco
    # Use the run with the most volumes
    # Use the same number of volumes for each tsnr calculation
    # ------------------------------------------------------------------
    print("=== Compute tSNR map on longest moco neighbour run ===", flush=True)
    fname_tsnr_metrics = os.path.join(path_fig_tsnr, "tsnr_metrics.csv")
    fname_template = os.path.join(config["code_dir"], "template", config["PAM50_t2"])
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
                tag = "task-" + task + "_acq-" + acq_name

                selected_file = find_moco_for_tsnr_calculation(config, ID, task, acq_name)
                if selected_file is None:
                    continue

                # Compute tSNR map in native space
                path_tsnr_sub_folder = os.path.join(path_fig_data, f"sub-{ID}", f"task-{task}_acq-{acq_name}")
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
                        os.path.basename(selected_file).replace("_bold_moco.nii.gz",
                                                                "_from-func_to_PAM50_mode-image_xfm.nii.gz")
                    )

                    if not os.path.exists(fname_warp_from_func_to_template):
                        raise RuntimeError(f"Warp file not found: {fname_warp_from_func_to_template}")

                    cmd_coreg = f"sct_apply_transfo -i {fname_tsnr} -d {fname_template} -w {fname_warp_from_func_to_template} -o {fname_tsnr_in_template} -x nn"
                    os.system(cmd_coreg)

                # Extract metrics from native space
                if fname_tsnr is not None:

                    fname_mask = os.path.join(config["raw_dir"],
                                              config["preprocess_dir"]["main_dir"].format(ID),
                                              'func',
                                              tag,
                                              f"sub-{ID}_{tag}_bold_moco_mean_seg.nii.gz")

                    if not os.path.exists(fname_mask):
                        raise RuntimeError(f"Mask file not found: {fname_mask}")

                    tsnr_mean = extract_mean_within_mask(fname_tsnr, fname_mask)
                    if len(df_tsnr) == 0:
                        df_tsnr = pd.DataFrame([[ID, task, acq_name, tsnr_mean]], columns=df_tsnr.columns)
                    df_tsnr = pd.concat(
                        [pd.DataFrame([[ID, task, acq_name, tsnr_mean]], columns=df_tsnr.columns), df_tsnr],
                        ignore_index=True)

    df_tsnr.to_csv(fname_tsnr_metrics, index=False)

    print("=== Generate tSNR violin plot figure ===", flush=True)
    name_baseline = [a for a in config["design_exp"]["acq_names"] if a.find("Base") != -1][0]
    name_slicewise = [a for a in config["design_exp"]["acq_names"] if a.find("Slice") != -1][0]
    df_tsnr = pd.read_csv(fname_tsnr_metrics)
    list_baseline_tsnr = []
    list_slicewise_tsnr = []
    for ID in IDs:
        df_sub = df_tsnr[df_tsnr["ID"] == int(ID)]
        done = False
        # Try rest task
        if len(df_sub[df_sub["task"] == "rest"]) >= 2:
            done = True
            df_task = df_sub[df_sub["task"] == "rest"]
            if len(df_task) != 2:
                raise RuntimeError(f"We don't have 2 tSNR metric for sub-{ID} task-rest")

            tsnr_baseline = df_task[df_task["acq"] == name_baseline]["tsnr_mean"].values
            tsnr_slicewise = df_task[df_task["acq"] == name_slicewise]["tsnr_mean"].values
            list_baseline_tsnr.append(tsnr_baseline[0])
            list_slicewise_tsnr.append(tsnr_slicewise[0])

        # If rest task not found, use motor task
        if not done:
            # Todo: If no rest task, use the motor task, we could use the volumes at rest during the motor task
            print(f"No rest task found for sub-{ID}, using motor task instead", flush=True)
            df_task = df_sub[df_sub["task"] == "motor"]
            if len(df_task) != 2:
                warnings.warn(f"We don't have 2 tSNR metric for sub-{ID} task-motor")
                continue

            tsnr_baseline = df_task[df_task["acq"] == name_baseline]["tsnr_mean"].values
            tsnr_slicewise = df_task[df_task["acq"] == name_slicewise]["tsnr_mean"].values
            list_baseline_tsnr.append(tsnr_baseline[0])
            list_slicewise_tsnr.append(tsnr_slicewise[0])

    # Create a figure
    fig = plt.figure(constrained_layout=True, figsize=(7, 5))
    gs_main = gridspec.GridSpec(1, 3, figure=fig, width_ratios=[0.08, 1, 1])

    create_tsnr_violin_plot(fig, gs_main[2], list_baseline_tsnr, list_slicewise_tsnr)

    print("=== Generate tSNR in PAM50 figure ===", flush=True)

    nii_template = nib.load(fname_template)
    data_tsnr_baseline = np.zeros_like(nii_template.get_fdata(), dtype=float)
    data_tsnr_slicewise = np.zeros_like(nii_template.get_fdata(), dtype=float)
    data_count_subjects_baseline = np.zeros_like(nii_template.get_fdata(), dtype=int)
    data_count_subjects_slicewise = np.zeros_like(nii_template.get_fdata(), dtype=int)

    fname_tsnr_baseline_avg = os.path.join(path_fig_tsnr, "data", "tsnr_baseline_avg_in_PAM50.nii.gz")
    fname_tsnr_slicewise_avg = os.path.join(path_fig_tsnr, "data", "tsnr_slicewise_avg_in_PAM50.nii.gz")

    for ID in IDs:
        tasks = ["rest", "motor"]
        done = False
        for task in tasks:
            if task == "rest":
                path_task_baseline = os.path.join(
                    path_fig_data,
                    f"sub-{ID}",
                    f"task-{task}_acq-{name_baseline}")
                if not os.path.exists(path_task_baseline):
                    continue
                done = True
                fname_tsnr_in_template_baseline = os.path.join(
                    path_task_baseline,
                    f"sub-{ID}_task-{task}_acq-{name_baseline}*_bold_moco_tsnr_in_PAM50.nii.gz"
                )
                if len(glob.glob(fname_tsnr_in_template_baseline)) != 1:
                    raise RuntimeError(
                        f"0 or more than 1 tSNR in template files found: {glob.glob(fname_tsnr_in_template_baseline)}")
                fname_tsnr_in_template_baseline = glob.glob(fname_tsnr_in_template_baseline)[0]

                path_task_slicewise = os.path.join(
                    path_fig_data,
                    f"sub-{ID}",
                    f"task-{task}_acq-{name_slicewise}")
                fname_tsnr_in_template_slicewise = os.path.join(
                    path_task_slicewise,
                    f"sub-{ID}_task-{task}_acq-{name_slicewise}*_bold_moco_tsnr_in_PAM50.nii.gz"
                )
                if len(glob.glob(fname_tsnr_in_template_slicewise)) != 1:
                    raise RuntimeError(
                        f"0 or more than 1 tSNR in template files found: {glob.glob(fname_tsnr_in_template_slicewise)}")
                fname_tsnr_in_template_slicewise = glob.glob(fname_tsnr_in_template_slicewise)[0]

                nii_baseline = nib.load(fname_tsnr_in_template_baseline)
                data_tsnr_baseline += nii_baseline.get_fdata()
                nii_slicewise = nib.load(fname_tsnr_in_template_slicewise)
                data_tsnr_slicewise += nii_slicewise.get_fdata()
                data_count_subjects_baseline += np.array(nii_baseline.get_fdata() > 0).astype(int)
                data_count_subjects_slicewise += np.array(nii_slicewise.get_fdata() > 0).astype(int)

            if not done:
                print(f"No rest task found for sub-{ID}, using motor task instead", flush=True)
                # Todo: If no rest task, use the motor task, we could use the volumes at rest during the motor task
                # This is only relevant for acqs 93, 94, 95, 96
                path_task_baseline = os.path.join(
                    path_fig_data,
                    f"sub-{ID}",
                    f"task-{task}_acq-{name_baseline}")
                if not os.path.exists(path_task_baseline):
                    warnings.warn(f"No motor task found for sub-{ID}, we need it to compute the tSNR figure")
                    continue
                done = True
                fname_tsnr_in_template_baseline = os.path.join(
                    path_task_baseline,
                    f"sub-{ID}_task-{task}_acq-{name_baseline}*_bold_moco_tsnr_in_PAM50.nii.gz"
                )
                if len(glob.glob(fname_tsnr_in_template_baseline)) != 1:
                    raise RuntimeError(
                        f"0 or more than 1 tSNR in template files found: {glob.glob(fname_tsnr_in_template_baseline)}")
                fname_tsnr_in_template_baseline = glob.glob(fname_tsnr_in_template_baseline)[0]

                path_task_slicewise = os.path.join(
                    path_fig_data,
                    f"sub-{ID}",
                    f"task-{task}_acq-{name_slicewise}")
                fname_tsnr_in_template_slicewise = os.path.join(
                    path_task_slicewise,
                    f"sub-{ID}_task-{task}_acq-{name_slicewise}*_bold_moco_tsnr_in_PAM50.nii.gz"
                )
                if len(glob.glob(fname_tsnr_in_template_slicewise)) != 1:
                    raise RuntimeError(
                        f"0 or more than 1 tSNR in template files found: {glob.glob(fname_tsnr_in_template_slicewise)}")
                fname_tsnr_in_template_slicewise = glob.glob(fname_tsnr_in_template_slicewise)[0]

                nii_baseline = nib.load(fname_tsnr_in_template_baseline)
                data_tsnr_baseline += nii_baseline.get_fdata()
                nii_slicewise = nib.load(fname_tsnr_in_template_slicewise)
                data_tsnr_slicewise += nii_slicewise.get_fdata()
                data_count_subjects_baseline += np.array(nii_baseline.get_fdata() > 0).astype(int)
                data_count_subjects_slicewise += np.array(nii_slicewise.get_fdata() > 0).astype(int)
        # Average
        data_tsnr_baseline_avg = np.divide(data_tsnr_baseline, data_count_subjects_baseline,
                                           where=data_count_subjects_baseline != 0)
        data_tsnr_slicewise_avg = np.divide(data_tsnr_slicewise, data_count_subjects_slicewise,
                                            where=data_count_subjects_slicewise != 0)

        nii_tsnr_baseline_avg = nib.Nifti1Image(data_tsnr_baseline_avg, affine=nii_template.affine,
                                                header=nii_template.header)
        nib.save(nii_tsnr_baseline_avg, fname_tsnr_baseline_avg)

        nii_tsnr_slicewise_avg = nib.Nifti1Image(data_tsnr_slicewise_avg, affine=nii_template.affine,
                                                 header=nii_template.header)
        nib.save(nii_tsnr_slicewise_avg, fname_tsnr_slicewise_avg)

    nii_tsnr_baseline_avg = nib.load(fname_tsnr_baseline_avg)
    nii_tsnr_slicewise_avg = nib.load(fname_tsnr_slicewise_avg)

    fname_tsnr_baseline_avg_smooth = os.path.join(path_fig_tsnr, "data", "tsnr_baseline_avg_smooth_in_PAM50.nii.gz")
    nii_tsnr_baseline_avg_smooth = smooth_img(nii_tsnr_baseline_avg, fwhm=[2, 2, 4])
    nib.save(nii_tsnr_baseline_avg_smooth, fname_tsnr_baseline_avg_smooth)

    fname_tsnr_slicewise_avg_smooth = os.path.join(path_fig_tsnr, "data", "tsnr_slicewise_avg_smooth_in_PAM50.nii.gz")
    nii_tsnr_slicewise_avg_smooth = smooth_img(nii_tsnr_slicewise_avg, fwhm=[2, 2, 4])
    nib.save(nii_tsnr_slicewise_avg_smooth, fname_tsnr_slicewise_avg_smooth)

    nii_tsnr_baseline_avg_smooth = nib.load(fname_tsnr_baseline_avg_smooth)
    nii_tsnr_slicewise_avg_smooth = nib.load(fname_tsnr_slicewise_avg_smooth)

    fname_template_seg = os.path.join(path_code, 'template', config["PAM50_cord"])
    fname_template_seg_dil = os.path.join(path_fig_tsnr, "data", "template_cord_dil.nii.gz")
    if not os.path.exists(fname_template_seg_dil) or redo:
        cmd = f"sct_maths -i {fname_template_seg} -dilate 1 -o {fname_template_seg_dil}"
        os.system(cmd)
    fname_template_seg_outline = os.path.join(path_fig_tsnr, "data", "template_cord_outline.nii.gz")
    if not os.path.exists(fname_template_seg_outline) or redo:
        cmd = f"sct_maths -i {fname_template_seg_dil} -sub {fname_template_seg} -o {fname_template_seg_outline}"
        os.system(cmd)
    nii_outline = nib.load(fname_template_seg_outline)

    create_tsnr_template_plot(fig, gs_main[1], nii_tsnr_baseline_avg_smooth, nii_tsnr_slicewise_avg_smooth,
                              nii_template, nii_outline)
    fig.suptitle(f"A) Average tSNR map (n={len(IDs)})              B) tSNR accross participants", fontsize=12, fontweight='bold',
                 x=0.52)
    fname_output = os.path.join(path_fig_tsnr, "tsnr_plot.png")
    fig.savefig(fname_output, dpi=1000)

    print("=== Figure script End ===", flush=True)


def create_tsnr_template_plot(fig, gs, nii_baseline, nii_slicewise, nii_template, nii_cord_outline):
    axial_slice = 244
    axial_left_bound = 47
    axial_right_bound = 95
    axial_bot_bound = 53
    axial_top_bound = 89

    axial_outline = nii_cord_outline.get_fdata()[axial_left_bound:axial_right_bound, axial_bot_bound:axial_top_bound,
                    axial_slice]
    axial_cord_outline_nan = np.full_like(axial_outline, np.nan)
    axial_cord_outline_nan[axial_outline > 0.5] = 1

    cor_slice = 72
    cor_left_bound = 47
    cor_right_bound = 95
    cor_bot_bound = 158
    cor_top_bound = 322
    cor_templ_bot = 125
    cor_templ_top = 360
    hline = cor_templ_top - axial_slice

    cor_template_slice = nii_template.get_fdata()[cor_left_bound:cor_right_bound, cor_slice, :]
    cor_baseline_slice = np.full_like(cor_template_slice, np.nan)
    cor_baseline_slice[:, cor_bot_bound:cor_top_bound] = nii_baseline.get_fdata()[cor_left_bound:cor_right_bound,
                                                         cor_slice, cor_bot_bound:cor_top_bound]
    cor_baseline_slice = cor_baseline_slice[:, cor_templ_bot:cor_templ_top]

    cor_slicewise_slice = np.full_like(cor_template_slice, np.nan)
    cor_slicewise_slice[:, cor_bot_bound:cor_top_bound] = nii_slicewise.get_fdata()[cor_left_bound:cor_right_bound,
                                                          cor_slice, cor_bot_bound:cor_top_bound]
    cor_slicewise_slice = cor_slicewise_slice[:, cor_templ_bot:cor_templ_top]

    cor_template_slice = cor_template_slice[:, cor_templ_bot:cor_templ_top]

    spinal_levels = {5: range(300, 333),  # C5
                     6: range(269, 300),  # C6
                     7: range(238, 269),  # C7
                     8: range(206, 238),  # C8
                     9: range(172, 206),  # T1
                     10: range(135, 172)}  # T2
    data_spinal_levels = np.full_like(cor_template_slice, 0, dtype=int)
    for level, range_ in spinal_levels.items():
        for r in range_:
            data_spinal_levels[:, r - cor_templ_bot] = level

    # Create a figure showing the average tSNR map in PAM50 space for baseline and slicewise shim
    vmin = 5
    vmax = 14
    gs_1 = gs.subgridspec(ncols=2, nrows=2, width_ratios=[1, 1], height_ratios=[4, 1])
    axs = gs_1.subplots()
    axs[0, 0].imshow(np.rot90(cor_template_slice), cmap='gray')
    im1 = axs[0, 0].imshow(np.rot90(cor_baseline_slice), vmin=vmin, vmax=vmax, cmap='turbo')
    axs[0, 0].hlines(hline, xmin=0, xmax=cor_right_bound - cor_left_bound - 1, color='black', linestyle='--',
                     linewidth=1)
    axs[0, 0].axis('off')
    axs[0, 0].set_aspect('equal', adjustable='box')
    axs[0, 0].set_title('Baseline\n2nd order shim', fontsize=10)
    axs[0, 1].imshow(np.rot90(cor_template_slice), cmap='gray')
    axs[0, 1].imshow(np.rot90(cor_slicewise_slice), vmin=vmin, vmax=vmax, cmap='turbo')
    axs[0, 1].hlines(hline, xmin=0, xmax=cor_right_bound - cor_left_bound - 1, color='black', linestyle='--',
                     linewidth=1)
    axs[0, 1].axis('off')
    axs[0, 1].set_aspect('equal', adjustable='box')
    axs[0, 1].set_title('Slice-wise\nf0xyz shim', fontsize=10)
    axs[1, 0].imshow(np.rot90(
        nii_baseline.get_fdata()[axial_left_bound:axial_right_bound, axial_bot_bound:axial_top_bound, axial_slice]),
                     vmin=vmin, vmax=vmax, cmap='turbo')
    axs[1, 0].imshow(np.rot90(axial_cord_outline_nan), vmin=0, vmax=1, cmap='grey', alpha=0.75)
    axs[1, 0].axis('off')
    axs[1, 0].set_aspect('equal', adjustable='box')
    axs[1, 1].imshow(np.rot90(
        nii_slicewise.get_fdata()[axial_left_bound:axial_right_bound, axial_bot_bound:axial_top_bound, axial_slice]),
                     vmin=vmin, vmax=vmax, cmap='turbo')
    axs[1, 1].imshow(np.rot90(axial_cord_outline_nan), vmin=0, vmax=1, cmap='grey', alpha=0.75)
    axs[1, 1].axis('off')
    axs[1, 1].set_aspect('equal', adjustable='box')

    # Add a colorbar under the 4 subplots
    cbar = fig.colorbar(im1, ax=axs, orientation='horizontal', shrink=0.4, ticks=[vmin, vmax])
    cbar.set_label('tSNR', fontsize=10)

    # Add orientation labels on the left axial slice and on the left coronal
    text_fontsize = 10
    axs[1, 0].text(0.08, 0.5, 'L', color='white', fontsize=text_fontsize, fontweight='bold', ha='center', va='center',
                   transform=axs[1, 0].transAxes)
    axs[1, 0].text(0.92, 0.5, 'R', color='white', fontsize=text_fontsize, fontweight='bold', ha='center', va='center',
                   transform=axs[1, 0].transAxes)
    axs[1, 0].text(0.5, 0.1, 'A', color='white', fontsize=text_fontsize, fontweight='bold', ha='center', va='center',
                   transform=axs[1, 0].transAxes)
    axs[1, 0].text(0.5, 0.9, 'P', color='white', fontsize=text_fontsize, fontweight='bold', ha='center', va='center',
                   transform=axs[1, 0].transAxes)

    axs[0, 0].text(-0.1, 0, 'L', color='black', fontsize=text_fontsize, fontweight='bold', ha='center', va='bottom',
                   transform=axs[0, 0].transAxes)
    axs[0, 0].text(1.1, 0, 'R', color='black', fontsize=text_fontsize, fontweight='bold', ha='center', va='bottom',
                   transform=axs[0, 0].transAxes)

    # Add spinal level labels on the left of the coronal slices
    l, b, w, h = axs[0, 0].get_position(False).bounds
    width = 0.015
    left = l - w - 0.02 - 0.08
    bottom = b - 0.128
    height = h + 0.122
    ax_spl_coords = (left, bottom, width, height)
    ax_spl = fig.add_axes(ax_spl_coords, anchor='SW')
    data_spinal_alpha = np.zeros_like(data_spinal_levels, dtype=float)
    data_spinal_alpha[data_spinal_levels > 0] = 1
    data_spinal_levels_2 = np.copy(data_spinal_levels).astype(float)
    data_spinal_levels_2[data_spinal_levels % 2 == 0] = 0.5
    data_spinal_levels_2[data_spinal_levels % 2 == 1] = 0.75
    ax_spl.imshow(np.rot90(data_spinal_levels_2), vmin=0, vmax=1, cmap='gray', alpha=np.rot90(data_spinal_alpha))
    ax_spl.set_xticks([])
    ax_spl.spines['top'].set_visible(False)
    ax_spl.spines['bottom'].set_visible(False)
    ax_spl.spines['left'].set_visible(False)
    ax_spl.spines['right'].set_visible(False)
    ax_spl.set_yticks([cor_template_slice.shape[1] - (np.mean(spinal_levels[i]) - cor_templ_bot) for i in range(5, 11)])
    ax_spl.set_yticklabels(['C5', 'C6', 'C7', 'C8', 'T1', 'T2'])
    ax_spl.tick_params(axis='y', which='both', bottom=False, top=False, left=False, right=False, labelcolor='grey',
                       labelsize=8)
    ax_spl.set_ylabel('Spinal levels', fontsize=8, color='grey')
    ax_spl.set_aspect('auto', adjustable='datalim')
    ax_spl.set_label('Spinal levels')


def create_tsnr_violin_plot(fig, gs, snr_baseline: list, snr_shim: list):
    # Create a violin plot with matplotolib showing the distribution of SNR values for baseline and shim
    # Link corresponding values for each subject with a line

    gs1 = gs.subgridspec(ncols=1, nrows=2, height_ratios=[9, 0.6])
    axs = gs1.subplots()
    axs[1].axis('off')

    ax = axs[0]
    data = [snr_baseline, snr_shim]
    parts = ax.violinplot(data, showmeans=True, showmedians=True)
    for pc in parts['bodies']:
        pc.set_facecolor('lightblue')
        pc.set_edgecolor('black')
        pc.set_alpha(0.7)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(['Baseline\n2nd order shim', 'Slice-wise\nf0xyz shim'])
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
    ax.tick_params(axis='x', which='both', bottom=False, top=False)
    ax.set_ylabel('Mean tSNR for each subject')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    # Link corresponding values for each subject with a line
    for i in range(len(snr_baseline)):
        ax.plot([1, 2], [snr_baseline[i], snr_shim[i]], color='gray', linestyle='--', marker='o', alpha=0.7)


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
