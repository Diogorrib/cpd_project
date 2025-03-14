#!/bin/bash

TEST_DIR="samples"
BASE_SCRIPT="test.sh"
PASSED=0
FAILED=0

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <serial|omp|mpi>"
    exit 1
fi

export OMP_NUM_THREADS=$(nproc)  # For Linux
# export OMP_NUM_THREADS=$(sysctl -n hw.ncpu)  # For macOS
echo "Using $OMP_NUM_THREADS threads for OpenMP execution..."

PROJ_DIR="$1"

for input_file in "$TEST_DIR"/*.in; do
    base_name=$(basename "$input_file" .in)
    expected_output_file="$base_name.out"

    bash "$BASE_SCRIPT" "$PROJ_DIR" "$base_name"

    # Capture the exit status of the original script
    result=$?

    # Increment counters based on the result
    if [ $result -eq 0 ]; then
        ((PASSED++))
    else
        ((FAILED++))
    fi
done

echo "Tests Passed: $PASSED"
echo "Tests Failed: $FAILED"
