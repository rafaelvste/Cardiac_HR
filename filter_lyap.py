import numpy as np

# Function to load data and apply filter
def process_lyapunov_data(input_file, output_file):
    data = np.loadtxt(input_file)

    # Create new column based on the condition
    new_third_column = []
    for row in data:
        if row[2] > 1e-3:
            new_third_column.append(row[2])
        else:
            new_third_column.append(row[3])

    # Convert to array for saving and manipulating
    new_third_column = np.array(new_third_column)
    new_data = np.column_stack((data[:, :2], new_third_column))

    # Save new data to .dat file
    np.savetxt(output_file, new_data, fmt="%.16e")

    return new_data

# Input and output files
input_file = "lyapunov_spectrum_vs_b_Ie_A.dat"
output_file = "lyap_vs_b_Ie_A_filt.dat"

# Process the data
data_filtered = process_lyapunov_data(input_file, output_file)