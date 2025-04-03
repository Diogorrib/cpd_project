#!/bin/bash

TEST_DIR="samples"
BASE_SCRIPT="test.sh"
PASSED=0
FAILED=0

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <serial|omp|mpi> [<optional-threads>]"
    exit 1
fi

mpi_function() {
    echo "localhost slots=$OMP_NUM_THREADS" > hostfile.txt
    for input_file in "$TEST_DIR"/*.in; do
        third_arg=$(awk '{print $3}' $input_file)
        if [ "$third_arg" -ge "$OMP_NUM_THREADS" ]; then
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
        fi
    done
}

other_function() {
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
}

if [ -z "$2" ]; then
    export OMP_NUM_THREADS=$(nproc)  # For Linux
    # export OMP_NUM_THREADS=$(sysctl -n hw.ncpu)  # For macOS
else
    export OMP_NUM_THREADS=$2
fi

echo "Using $OMP_NUM_THREADS threads (or processes) for OpenMP (or MPI) execution..."

PROJ_DIR="$1"

if [[ "$PROJ_DIR" == mpi* ]]; then
    mpi_function
else
    other_function
fi

echo "Tests Passed: $PASSED"
echo "Tests Failed: $FAILED"
