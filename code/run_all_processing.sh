#!/bin/bash
# --------------------------
# User parameters
# --------------------------
# Default values
PATH_DATA="$PATH_DATA" #Defaults from environment
PATH_CODE="$PATH_CODE" #Defaults from environment
IDs=() # empty  → process all participants 

RUN_PREPROSS=true
RUN_DENOISING=true

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --path_data) PATH_DATA="$2"; shift 2 ;;
        --path_code) PATH_CODE="$2"; shift 2 ;;
        --ids) shift; while [[ $# -gt 0 && "$1" != --* ]]; do IDs+=("$1"); shift; done ;;
        --no-preprocess) RUN_PREPROSS=false; shift;;
        --no-denoising) RUN_DENOISING=false; shift;;
        *) echo "Unknown argument $1"; exit 1 ;;
    esac
done

# Show participants
[ ${#IDs[@]} -eq 0 ] && echo "No specific IDs provided: processing all participants" \
                     || echo "Processing participants: ${IDs[@]}"
[ ${#IDs[@]} -eq 0 ] && IDs=("") # If no IDs were provided, set to empty string 

# --------------------------
# Prepare log folder
# --------------------------
cd "$PATH_CODE" || exit 1
mkdir -p log
cd log || exit 1

timestamp=$(date +"%Y%m%d_%H%M%S")

# --------------------------
# Run preprocessing
# --------------------------

if [ "$RUN_PREPROSS" = true ]; then
    echo "Starting preprocessing..."
    nohup python -u ../code/preprocessing_workflow.py --ids "${IDs[@]}" --redo True \
    > "nohup_preprocessing_${timestamp}.out" 2>&1 &

    echo "Preprocessing launched in background."
    echo "Log file: log/nohup_preprocessing_${timestamp}.out"
fi

# --------------------------
# Run denoising
# --------------------------

if [ "$RUN_DENOISING" = true ]; then
    echo "Starting denoising..."
    nohup python -u ../code/denoising_workflow.py --ids "${IDs[@]}" --redo True \
    > "nohup_denoising_${timestamp}.out" 2>&1 &

    echo "Denoising launched in background."
    echo "Log file: log/nohup_denoising_${timestamp}.out"
fi