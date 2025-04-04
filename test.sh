#!/bin/bash

check_file_exists() {
    if [ ! -f "$1" ]; then
        echo "[ERROR] $1 not found"
        exit 1
    fi
}

check_executable_exists() {
    if [ ! -x "$1" ]; then
        echo "[ERROR] $1 not found or not executable: try running 'make' first"
        exit 1
    fi
}

mpi_function() {
    args=$(<"$input_file")
    actual_output=$(mpirun --oversubscribe --hostfile hostfile.txt $PROJ_DIR/parsim $args 2>/tmp/parsim_stderr)
}

other_function() {
    args=$(<"$input_file")
    actual_output=$(./$PROJ_DIR/parsim $args 2>/tmp/parsim_stderr)
}

TEST_DIR="samples"

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <serial|omp|mpi> <base_name>"
    exit 1
fi

PROJ_DIR="$1"
base_name="$2"

check_executable_exists "$PROJ_DIR/parsim"

input_file="$TEST_DIR/${base_name}.in"
expected_output_file="$TEST_DIR/${base_name}.out"

check_file_exists "$input_file"
check_file_exists "$expected_output_file"

if [[ "$PROJ_DIR" == mpi* ]]; then
    mpi_function
else
    other_function
fi

# Compare actual output with expected output
if diff -wB <(echo "$actual_output") "$expected_output_file" >/dev/null; then
    echo "[PASS] $input_file"
    cat /tmp/parsim_stderr
    exit 0
else
    echo "[FAIL] $input_file"
    echo "Expected:"
    echo "===================="
    cat "$expected_output_file"
    echo "===================="
    echo "Got:"
    echo "===================="
    echo "$actual_output"
    echo "===================="
    exit 1
fi
