import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Load data from the file with high precision
data = np.loadtxt("superstability_points.dat")

x = data[:, 0]
y = data[:, 1]

# Define the fitting functions
def quadratic_model(x, a, b, c):
    return a * x**2 + b * x + c

def cubic_model(x, a, b, c, d):
    return a * x**3 + b * x**2 + c * x + d

# Perform quadratic fit
quad_params, _ = curve_fit(quadratic_model, x, y)
y_quad_fit = quadratic_model(x, *quad_params)

# Perform cubic fit
cubic_params, _ = curve_fit(cubic_model, x, y)
y_cubic_fit = cubic_model(x, *cubic_params)

# Calculate R^2 for both models
def calculate_r2(y_true, y_pred):
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    return 1 - (ss_res / ss_tot)

r2_quad = calculate_r2(y, y_quad_fit)
r2_cubic = calculate_r2(y, y_cubic_fit)

# Generate fitted values for both models
x_fit = np.linspace(1.3, 2.0, 100)
y_quad_fit_points = quadratic_model(x_fit, *quad_params)
y_cubic_fit_points = cubic_model(x_fit, *cubic_params)

# Save the fitted points for both models to .dat files
quad_output_data = np.column_stack((x_fit, y_quad_fit_points))
cubic_output_data = np.column_stack((x_fit, y_cubic_fit_points))
np.savetxt("quadratic_fit_curve.dat", quad_output_data, fmt="%.6e", header="x y", comments="")
np.savetxt("cubic_fit_curve.dat", cubic_output_data, fmt="%.6e", header="x y", comments="")

# Plot original data and fitted curves
plt.figure(figsize=(8, 6))
plt.scatter(x, y, label="Original Data", color="red", zorder=5)
plt.plot(x_fit, y_quad_fit_points, label=f"Quadratic Fit (R^2 = {r2_quad:.3f})", color="green", linestyle="--")
plt.plot(x_fit, y_cubic_fit_points, label=f"Cubic Fit (R^2 = {r2_cubic:.3f})", color="purple", linestyle="--")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Curve Fitting")
plt.legend()
plt.grid()
plt.show()
