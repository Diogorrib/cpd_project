#!/bin/bash
TEST_DIR="samples_reduced/"
EXEC_DIR="mpi/"
OUT_DIR="out/"
SCRIPT_DIR="scripts/"
RESULTS_DIR="results/"

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <number_of_tasks> <number_of_threads> <number_of_runs>"
    exit 1
fi

num_processes=$1
num_threads=$2
num_runs=$3
test_id="${num_processes}_${num_threads}"
results_file="${RESULTS_DIR}${test_id}_procs.txt"

# Track tests to get mean and results
base_names=()

check_result_and_get_time() {
    local base_name="$1"
    local index="$2"

    local expected_output="${TEST_DIR}${base_name}.out"
    local job_name="${test_id}_${base_name}_${index}"
    local output_file="${OUT_DIR}${job_name}.out"
    local error_file="${OUT_DIR}${job_name}.err"

    if diff -wB "$output_file" "$expected_output" >/dev/null; then
        echo "[PASS] $job_name" >> "$results_file"
    else
        echo "[FAIL] $job_name"     >> "$results_file"
        echo "Expected:"            >> "$results_file"
        echo "====================" >> "$results_file"
        cat "$expected_output"      >> "$results_file"
        echo "Actual:"              >> "$results_file"
        echo "====================" >> "$results_file"
        cat "$output_file"          >> "$results_file"
        echo "Error log:"           >> "$results_file"
        cat "$error_file"           >> "$results_file"
        echo "====================" >> "$results_file"
    fi

    # Extract the execution time from the error file and return it
    elapsed_time=$(grep -oP '\d+(\.\d+)?s' "$error_file" | sed 's/s$//')
    echo "$elapsed_time"
}

for output_file in ${OUT_DIR}${test_id}_*.out; do
    # Get the base name of the output file
    base_name=$(basename "$output_file" .out | sed 's/^[0-9]*_[0-9]*_//; s/_*[0-9]*$//')

    # Add the base name to the array if it's not already there
    if [[ ! " ${base_names[@]} " =~ " ${base_name} " ]]; then
        base_names+=("$base_name")
    fi
done

# Make sure the needed directories exist
mkdir -p "$RESULTS_DIR"

echo "Results for $num_processes tasks, $num_threads threads & $num_runs runs:" > "$results_file"
echo "========================================="                                >> "$results_file"

for base_name in "${base_names[@]}"; do
    acc_time=0
    for ((i=0; i<num_runs; i++)); do
        # Verify results and accumulate execution time
        elapsed_time=$(check_result_and_get_time "$base_name" "$((i+1))")
        acc_time=$(echo "$acc_time + $elapsed_time" | bc)
    done

    # Calculate mean execution time
    mean_time=$(echo "scale=2; $acc_time / $num_runs" | bc)
    echo "Mean time for $base_name over $num_runs runs: $mean_time seconds" >> "$results_file"
done
echo "=========================================" >> "$results_file"

cat "$results_file"
echo "Results above saved to $results_file"
