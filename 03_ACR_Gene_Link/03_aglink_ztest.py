import pandas as pd
from scipy.stats import norm
import argparse

parser = argparse.ArgumentParser(description='Perform Z-test on two datasets.')
parser.add_argument('file1_path', type=str, help='Path to the first file (control data).')
parser.add_argument('file2_path', type=str, help='Path to the second file (ACR data).')
parser.add_argument('output_path', type=str, help='Path to save the output results.')

args = parser.parse_args()

# Read file 1, with no limit on the number of columns
file1_data = pd.read_csv(args.file1_path, delim_whitespace=True, header=None)

# Read file 2, assuming the first column is the name, the second column is the paired name, and the third column is the value.
file2_data = pd.read_csv(args.file2_path, delim_whitespace=True, header=None)
file2_data.columns = ['name', 'pair_name', 'value']

output_data = []

for _, row in file1_data.iterrows():
    name = row[0] 
    values = row[1:].values.astype(float)

    if values.size == 0 or pd.isnull(values).all():
        print(f"Warning: No valid values for name '{name}'. Skipping.")
        continue
    mean_file1 = values.mean()
    std_file1 = values.std(ddof=1) 
    if std_file1 == 0:
        print(f"Warning: Standard deviation is zero for name '{name}'. Skipping.")
        continue
    matching_rows = file2_data[file2_data['name'] == name]
    
    if matching_rows.empty:
        print(f"Warning: No matching rows found in file2 for name '{name}'.")
        continue
    for _, match in matching_rows.iterrows():
        pair_name = match['pair_name']
        value_file2 = match['value']

        if pd.isnull(value_file2):
            print(f"Warning: Missing value in file2 for name '{name}', pair_name='{pair_name}'. Skipping.")
            continue
        abs_mean_file1 = abs(mean_file1)
        abs_value_file2 = abs(value_file2)
        if abs_mean_file1 > 0:
            z = (abs_value_file2 - abs_mean_file1) / std_file1
            p_value = 1 - norm.cdf(z)
        else:
            print(f"Warning: Mean is zero for name '{name}'. Skipping.")
            continue
        output_data.append([name, pair_name, mean_file1, std_file1, z, value_file2, p_value])

output_df = pd.DataFrame(output_data, columns=['name', 'pair_name', 'mean_file1', 'std_file1', 'z', 'value_file2', 'p_value'])

# Save file
output_df.to_csv(args.output_path, sep='\t', index=False)

print(f"Results have been saved to {args.output_path}")


