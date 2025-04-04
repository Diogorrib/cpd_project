#!/bin/bash
TEST_DIR="samples/"
EXEC_DIR="mpi/"
OUT_DIR="out/"
SCRIPT_DIR="scripts/"
RESULTS_DIR="results/"

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <number_of_tasks> <number_of_runs>"
    exit 1
fi

n_threads=2
num_processes=$1
num_runs=$2
results_file="${RESULTS_DIR}${num_processes}_procs.txt"

# Track jobs in the cluster
total_jobs=0
job_ids=()

# Make sure the needed directories exist
mkdir -p "$OUT_DIR"
mkdir -p "$SCRIPT_DIR"

compile_mpi() {
    echo "Compiling code in $EXEC_DIR..."

    # (Re)compile the code
    cd "$EXEC_DIR" || exit
    make clean > /dev/null 2>&1
    make > /dev/null 2>&1

    if [ $? -ne 0 ]; then
        echo "Compilation failed. Exiting..."
        exit 1
    fi

    # Return to the original directory
    cd - || exit

    echo "Compilation successful!"
}

run_script() {
    local base_name="$1"
    local index="$2"
    local args="$3"

    local job_name="${num_processes}_${base_name}_${index}"
    local script_name="${SCRIPT_DIR}${job_name}.sh"

    # Create the script for the job
    echo "#!/usr/bin/env bash"                          > "$script_name"
    echo "#SBATCH --job-name=$job_name"                 >> "$script_name"
    echo "#SBATCH --output=${OUT_DIR}${job_name}.out"   >> "$script_name"
    echo "#SBATCH --error=${OUT_DIR}${job_name}.err"    >> "$script_name"
    echo "#SBATCH --ntasks=$num_processes"              >> "$script_name"
    echo "#SBATCH --cpus-per-task=$n_threads"           >> "$script_name"
    echo "#SBATCH --export=OMP_NUM_THREADS=$n_threads"  >> "$script_name"
    echo ""                                             >> "$script_name"
    echo "srun ${EXEC_DIR}parsim $args"                 >> "$script_name"

    # Submit job and store job ID
    job_id=$(sbatch "$script_name" | awk '{print $4}')
    job_ids+=("$job_id")
}

check_result_and_get_time() {
    local base_name="$1"
    local index="$2"

    local expected_output="${TEST_DIR}${base_name}.out"
    local job_name="${num_processes}_${base_name}_${index}"
    local output_file="${OUT_DIR}${job_name}.out"
    local error_file="${OUT_DIR}${job_name}.err"

    if diff -wB <(echo "$actual_output") "$expected_output_file" >/dev/null; then
        echo "[PASS] $job_name"

        # Extract the execution time from the error file and return it
        elapsed_time=$(grep -oP '\d+(\.\d+)?s' "$error_file" | sed 's/s$//')
        echo "$elapsed_time"
    else
        echo "[FAIL] $job_name"
        echo "Expected:"
        echo "===================="
        cat "$expected_output_file"
        echo "Actual:"
        echo "===================="
        cat "$output_file"
        echo "Error log:"
        cat "$error_file"
        echo "===================="
        exit 1
    fi
}

# Delete old files
rm -f "${OUT_DIR}${num_processes}_*.out" "${OUT_DIR}${num_processes}_*.err" "${SCRIPT_DIR}${num_processes}_*.sh"
rm -f "${RESULTS_DIR}${num_processes}_procs.txt"

compile_mpi

# Get input files and start jobs
for input_file in "$TEST_DIR"*.in; do
    n_rows=$(awk '{print $3}' "$input_file")
    if [ -n "$n_rows" ] && [ "$n_rows" -ge "$num_processes" ]; then
        base_name=$(basename "$input_file" .in)
        args=$(<"$input_file")

        for ((i=0; i<$num_runs; i++)); do
            run_script "$base_name" "$((i+1))" "$args"
            ((total_jobs++))
        done
    fi
done

echo "$total_jobs jobs have been submitted. Now waiting for them to finish..."
echo ""
echo "You can also Ctrl+C to stop this script at any time. The jobs will continue running in the background. Then use collect_results_and_time.sh to check the results when they are done."
echo ""

# Wait for all jobs to finish
completed_jobs=()
count=0
while [ ${#completed_jobs[@]} -lt $total_jobs ]; do

    for job_id in "${job_ids[@]}"; do
        # Check if the job ID is in the squeue output (using grep)
        if ! squeue --job "$job_id" | grep -q "$job_id"; then
            # Job is done, add it to completed_jobs array if not already there
            if [[ ! " ${completed_jobs[@]} " =~ " ${job_id} " ]]; then
                completed_jobs+=("$job_id")
                echo "Job $job_id is complete."
            fi
        fi
    done

    remaining_jobs=$((total_jobs - ${#completed_jobs[@]}))
    echo "Waiting for $remaining_jobs jobs to finish... ($((count++)) checks so far)"

    sleep 10  # Sleep for 10 seconds before checking again
done

echo "All jobs finished. Collecting results and execution times..."

bash collect_results_and_time.sh "$num_processes" "$num_runs"

echo "Done!"
