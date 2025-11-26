mkdir -p log
cd log
timestamp=$(date +"%Y%m%d_%H%M%S")

nohup python ../code/01_spine7T_preprocessing.py ../log > "nohup_${timestamp}.out" 2>&1 &