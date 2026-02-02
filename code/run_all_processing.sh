#!/bin/bash
# --------------------------
# User parameters
# --------------------------
# Default values
PATH_DATA="$PATH_DATA" #Defaults from environment
PATH_CODE="$PATH_CODE" #Defaults from environment
IDs=() # empty  → process all participants
TASKS=() # empty → process all tasks
RUN_PREPROSS=true
RUN_DENOISING=true
RUN_FIGURES=true
REDO=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --path-data) PATH_DATA="$2"; shift 2 ;;
        --path-code) PATH_CODE="$2"; shift 2 ;;
        --ids) shift; while [[ $# -gt 0 && "$1" != --* ]]; do IDs+=("$1"); shift; done ;;
        --tasks) shift; while [[ $# -gt 0 && "$1" != --* ]]; do TASKS+=("$1"); shift; done ;;
        --no-preprocess) RUN_PREPROSS=false; shift;;
        --no-denoising) RUN_DENOISING=false; shift;;
        --no-figures) RUN_FIGURES=false; shift;;
        --redo) REDO=true; shift;;
      *) echo "Unknown argument $1"; exit 1 ;;
    esac
done

# Show participants
[ ${#IDs[@]} -eq 0 ] && echo "No specific IDs provided: processing all participants" \
                     || echo "Processing participants: ${IDs[*]}"
[ ${#IDs[@]} -eq 0 ] && IDs=("") # If no IDs were provided, set to empty string


if [ ${#TASKS[@]} -eq 0 ]; then
    echo "No task specified: processing all tasks"
    TASKS_ARG=()        # do not pass --task
else
    echo "Processing tasks: ${TASKS[*]}"
    TASKS_ARG=(--tasks "${TASKS[@]}")
fi


# --------------------------
# Prepare log folder
# --------------------------
cd "${PATH_CODE}" || exit 1
mkdir -p log
cd log || exit 1

timestamp=$(date +"%Y%m%d_%H%M%S")

# --------------------------
# Run preprocessing
# --------------------------

if [ "${RUN_PREPROSS}" = true ]; then
    echo "Starting preprocessing..."
    nohup python -u ../code/preprocessing_workflow.py --path-data "${PATH_DATA}"  --ids "${IDs[@]}" "${TASKS_ARG[@]}" --redo "${REDO}" \
    > "nohup_preprocessing_${timestamp}.out" 2>&1 &

    PID=$!
    echo "Preprocessing launched in background."
    echo "Log file: log/nohup_preprocessing_${timestamp}.out"
    echo "To stop the process, run:"
    echo "kill ${PID}"
fi

# --------------------------
# Run denoising
# --------------------------

if [ "${RUN_DENOISING}" = true ]; then
    echo "Starting denoising..."
    nohup python -u ../code/denoising_workflow.py --path-data "${PATH_DATA}" --ids "${IDs[@]}" "${TASKS_ARG[@]}" --redo "${REDO}" \
    > "nohup_denoising_${timestamp}.out" 2>&1 &
    
    PID=$!
    echo "Denoising launched in background."
    echo "Log file: log/nohup_denoising_${timestamp}.out"
    echo "To stop the process, run:"
    echo "kill ${PID}"
fi


# --------------------------
# Run figures
# --------------------------

if [ "${RUN_FIGURES}" = true ]; then
    echo "Starting figure generation..."
    nohup python -u ../code/figure_workflow.py --path-data "${PATH_DATA}" --ids "${IDs[@]}" "${TASKS_ARG[@]}" --redo "${REDO}" \
    > "nohup_figures_${timestamp}.out" 2>&1 &

    PID=$!
    echo "Figure generation launched in background."
    echo "Log file: log/nohup_figures_${timestamp}.out"
    echo "To stop the process, run:"
    echo "kill ${PID}"
fi
