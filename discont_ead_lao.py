# Import the libraries
import numpy as np
import matplotlib.pyplot as plt

# Import the relevant ISI_max data
data = np.loadtxt('diag_ISI_max_E.dat')
data = np.loadtxt('diag_ISI_max_B.dat')

# Transition between EAO and LAO
threshold = 35 #arbitrary threshold (defined by looking at the data),
                # all discontinuity with Delta_ISI_max > 35 is saved

# Obtain the indexes of the rows with a difference > threshold
discontinuities_idx = np.where(abs(np.diff(data[:,2]))>threshold)[0]

# Uses the indexes to obtain the corresponding (Ie,b) values of the transition 
ead_lao_points = np.zeros((discontinuities_idx.shape[0],2))
ead_lao_points = data[discontinuities_idx,:-1]

# Because each line of the diagram stops at a certain b value, a discontinuty
# in the data is always identified for the rows with b = b_max
idx_delete = np.where(ead_lao_points[:, 0] == max(ead_lao_points[:,0]))[0]
ead_lao_points_clean_E = np.delete(ead_lao_points, idx_delete, axis=0)
ead_lao_points_clean_B = np.delete(ead_lao_points, idx_delete, axis=0)

# Concatenate data from the two files
ead_lao_points_clean = np.concatenate((ead_lao_points_clean_E,
                                       ead_lao_points_clean_B), axis=0)

# Plot the results
plt.plot(ead_lao_points_clean[:,0], ead_lao_points_clean[:,1], 'o')
plt.grid()
plt.show()

# Save to files
np.savetxt("ead_lao_points.dat", ead_lao_points_clean, delimiter=" ", fmt = '%.18e %.18e')
