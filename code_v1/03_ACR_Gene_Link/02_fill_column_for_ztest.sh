#!/bin/bash

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input_file> <output_file> <total_columns> <fill_value>"
    echo "fill_value can be 'mean', '0', or 'NA'"
    exit 1
fi

input_file="$1"
output_file="$2"
total_columns="$3"
fill_value="$4"
#fill_value can be mean, 0 or NA.

awk -v fill_value="$fill_value" -v total_columns="$total_columns" '{
    n = NF
    sum = 0
    for (i = 1; i <= n; i++) {
        sum += $i
    }
    mean = sum / n
    printf "%s", $1
    for (i = 2; i <= n; i++) {
        printf " %s", $i
    }
    for (i = n + 1; i <= total_columns; i++) {
        if (fill_value == "mean") {
            printf " %.6f", mean
        } else if (fill_value == "0") {
            printf " 0.000000"
        } else if (fill_value == "NA") {
            printf " NA"
        } else {
            printf " Invalid fill_value"
            exit 1
        }
    }
    printf "\n"
}' "$input_file" > "$output_file"

echo "File saved: $output_file"
