
import json
import nibabel as nib
import pandas as pd
import numpy as np
import argparse


	
def get_files_info(json_file, fmri_nifti, timing_file):
    """
    Load files

    Parameters
    ----------
    json_file : str
        Path to the fMRI JSON sidecar file.
    fmri_nifti : str
        Path to the fMRI NIfTI file (.nii or .nii.gz).
    timing_file : str
        Path to the timing file containing trigger times.

    Returns
    -------
    tr : float
        Repetition time (TR) in seconds.
    nb_slices : int
        Number of slices in the acquisition.
    nb_volumes : int
        Number of fMRI volumes (timepoints).
    df_triggers : pandas.DataFrame
        DataFrame containing trigger times
    """

    # --- Load JSON sidecar ---
    with open(json_file, 'r') as f:
        json_data = json.load(f)
        
    tr = json_data.get("RepetitionTime")
    slice_timing = json_data.get("SliceTiming", [])
    nb_slices = len(slice_timing) if slice_timing else None

    # --- Load fMRI NIfTI to get number of volumes ---
    img = nib.load(fmri_nifti)
    nb_volumes = img.shape[3]

    # --- Load timing file to get trigger time ---
    df = pd.read_csv(timing_file, sep="\t", header=None,
                     names=["time", "type", "message"])
    df_triggers = df[df['type'] == 'DATA '].reset_index(drop=True)
    df_triggers['time'] = df_triggers['time'].astype(float)
	
    if len(df_triggers) < nb_volumes:
        raise ValueError("Number of triggers is smaller than number of fMRI volumes")

    return tr, nb_slices, nb_volumes, df_triggers

def calculate_volume_onset(df_triggers, tr, nb_dummy_scans):
	"""
	Calculate the onset times of fMRI volumes based on trigger times.

	Parameters
	----------
	df_triggers : pandas.DataFrame
		DataFrame containing trigger times.
	tr : float
		Repetition time (TR) in seconds.
	nb_dummy_scans : int
		Number of dummy scans to exclude from the beginning.

	Returns
	-------
	onset_times : list of float
		List of onset times for each fMRI volume.
	"""  
	df_nodummy = df_triggers.iloc[nb_dummy_scans:].reset_index(drop=True) # remove nb_dummy_scans rows
	df_nodummy["onset_time_real"] = df_nodummy["time"] # Onset times
	df_nodummy["onset_time_relative"] = df_nodummy["time"] - df_nodummy["time"].iloc[0] # Onset times relative to the first volume
	df_nodummy['volume'] = np.arange(1, len(df_nodummy) + 1) # calculate the volume number starting from 1
	df_nodummy['onset_time_estimate'] = (df_nodummy['volume'] - 1) * tr # estimate the volume time based on the volume number
	df_volume = df_nodummy[['volume', 'onset_time_real','onset_time_relative', 'onset_time_estimate']].copy()

	return df_volume

def task_design(design_file, df_volume, events=["start", "rest", "trial_RH"]):
    '''
    Create a dataframe with the timing information of the different events of the task
    Parameters
    ----------
    design_file : str
        Path to the design file (CSV).
    df_volume : pandas.DataFrame
        DataFrame containing volume onset times.
    events : list of str
        List of event names to extract from the design file.
    Returns
    -------
    df_events : pandas.DataFrame
        DataFrame containing the timing information of the events.
    '''
    df_design = pd.read_csv(design_file, sep=",")
    event_all, onset_all, duration_all = [], [], []

    for event in events:
        start_col = event + ".started"
        stop_col = event + ".stopped"

        if start_col not in df_design.columns:
            raise ValueError(f"Column {start_col} missing in {design_file}")
        if stop_col not in df_design.columns:
            raise ValueError(f"Column {stop_col} missing in {design_file}")

        starts = df_design[start_col].dropna().values
        stops = df_design[stop_col].dropna().values

        # calculate the number of time there is a value in the column event + '.started'
        nb_event = df_design[event + '.started'].notna().sum()
        for s, e in zip(starts, stops):
            event_all.append(event)
            onset_all.append(float(s))
            duration_all.append(float(e - s))

    df_events = pd.DataFrame({
        "onset": onset_all,
        "duration": duration_all,
        "trial_type": event_all  
    })

    df_events = df_events.sort_values(by='onset').reset_index(drop=True) # reorder the columns based on the starting time

    # Onset relative to the first volume
    first_volume_time = df_volume['onset_time_real'].iloc[0]
    df_events['onset'] -= first_volume_time

    return df_events


def main():

    parser = argparse.ArgumentParser(description="Convert behavioral logs to BIDS events.tsv")

    parser.add_argument("--json", required=True, help="fMRI JSON sidecar")
    parser.add_argument("--nifti", required=True, help="fMRI NIfTI file")
    parser.add_argument("--triggers", required=True, help="Trigger timing TSV file")
    parser.add_argument("--design", required=True, help="Task design CSV file")
    parser.add_argument("--output", required=True, help="Output events.tsv file")

    parser.add_argument("--dummy", type=int, default=0, help="Number of dummy scans")
    parser.add_argument("--events", nargs="+", default=["start","rest","trial_RH"])

    args = parser.parse_args()

    tr, nb_slices, nb_volumes, df_triggers = get_files_info(
        args.json, args.nifti, args.triggers
    )

    df_volume = calculate_volume_onset(df_triggers, tr, args.dummy)

    df_events = task_design(
        args.design,
        df_volume,
        args.events
    )

    df_events.to_csv(args.output, sep="\t", index=False)
    
    # Create JSON sidecar for events.tsv
    events_json_dict =  {"onset": {"Description": "Onset of the event in seconds relative to the first volume of the scan."},
                         "duration": {"Description": "Duration of the event in seconds."},
                         "trial_type": {"Description": "Name of the experimental condition."}}
    json_output = args.output.replace('.tsv', '.json')
    with open(json_output, 'w') as f:
        json.dump(events_json_dict, f, indent=4)
    


if __name__ == "__main__":
    main()