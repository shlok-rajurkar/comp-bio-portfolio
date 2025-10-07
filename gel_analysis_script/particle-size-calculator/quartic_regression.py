import numpy as np

# Takes in recorded x vals and y vals, returns fn 
# that can estimate y vals given array of x vals

def quartic_fit(x_vals, y_vals):
    coeffs = np.polyfit(x_vals, y_vals, 4)
    def quartic_regression(x_vals_to_estimate):
        return {coeffs[0] * x_vals_to_estimate ** 4 +
                coeffs[1] * x_vals_to_estimate ** 3 +
                coeffs[2] * x_vals_to_estimate ** 2 +
                coeffs[3] * x_vals_to_estimate ** 1 +
                coeffs[4]}
    return quartic_regression
