#!/bin/bash

PROJ_DIR="serial"
TEST_DIR="../samples"

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <base_name>"
    exit 1
fi

cd "$PROJ_DIR"

base_name="$1"
input_file="$TEST_DIR/${base_name}.in"
expected_output_file="$TEST_DIR/${base_name}.out"

base_name="${input_file%.in}"
expected_output_file="$base_name.out"

args=$(cat "$input_file")
actual_output=$(./parsim $args 2>/tmp/parsim_stderr)

if diff -wB <(echo "$actual_output") "$expected_output_file" >/dev/null; then
    echo "[PASS] $input_file"
    cat /tmp/parsim_stderr
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
fi
