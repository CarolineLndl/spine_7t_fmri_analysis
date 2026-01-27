
IDs=(107)
TASKS=(motor_acq-shimBase+3mm motor_acq-shimSlice+3mm)
#TASKS=(motor_acq-shimBase+3mm motor_acq-shimSlice+3mm_run-01 motor_acq-shimSlice+3mm_run-02)
DUMMY=0

# Loop over subjects and tasks
for ID in "${IDs[@]}"; do
    for TASK in "${TASKS[@]}"; do
        DUMMY=0  # set to 0 if you don't want to skip dummy scans
        JSON=${PATH_DATA}"/sub-"$ID"/func/sub-"$ID"_task-"$TASK"_bold.json"
        NIFTI=${PATH_DATA}"/sub-"$ID"/func/sub-"$ID"_task-"$TASK"_bold.nii.gz"
        TRIGGERS=${PATH_DATA}"/sourcedata/sub-"$ID"/beh/sub-"$ID"_task-"$TASK"_beh.log"
        DESIGN=${PATH_DATA}"/sourcedata/sub-"$ID"/beh/sub-"$ID"_task-"$TASK"_beh.csv"
        OUTPUT=${PATH_DATA}"/sub-"$ID"/func/sub-"$ID"_task-"$TASK"_events.tsv"

        # Print variables
        
        python3 beh2bids.py \
            --json "$JSON" \
            --nifti "$NIFTI" \
            --triggers "$TRIGGERS" \
            --design "$DESIGN" \
            --output "$OUTPUT" \
            --dummy "$DUMMY" \
            --events start rest trial_RH
        
        echo "Converted behavioral data for subject "$ID", task "$TASK
    done
done

