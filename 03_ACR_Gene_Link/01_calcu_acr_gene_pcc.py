import numpy as np
import sys

# read input filename
if len(sys.argv) < 2:
    print("Usage: python script.py <matrix>")
    sys.exit(1)

filename = sys.argv[1]

# read matrix
matrix = np.loadtxt(filename)

# extract head line
first_row = matrix[0, :]

# Calculate the correlation coefficient between the first row and each of the other rows line by line, retain four decimal places and output to the screen.
for i in range(1, matrix.shape[0]):
    correlation = np.corrcoef(first_row, matrix[i, :])[0, 1]
    print(f"{correlation:.4f}")
