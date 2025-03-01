#!/bin/bash

PROJ_DIR="omp"
TEST_DIR="../samples"
PASSED=0
FAILED=0

cd "$PROJ_DIR"

for input_file in "$TEST_DIR"/*.in; do
    base_name="${input_file%.in}"
    expected_output_file="$base_name.out"

    args=$(cat "$input_file")
    actual_output=$(./parsim $args 2>/tmp/parsim_stderr)

    if diff -wB <(echo "$actual_output") "$expected_output_file" >/dev/null; then
        echo "[PASS] $input_file"
        cat /tmp/parsim_stderr
        ((PASSED++))
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
        ((FAILED++))
    fi
done

echo "Tests Passed: $PASSED"
echo "Tests Failed: $FAILED"
