"""
Generates design points on which you should evaluate the model.

This requires installing pyDOE which you can do by:

    pip install pyDOE

The design we are going to use is called a Latin Hypercube Design.

"""

import pyDOE as pd
import numpy as np

# Specify how many simulations you can afford:
num_simulations = 15


# Specify the lower and upper bounds of your model inputs
# I am not using exactly your numbers. I am rounding up or down.
# The purpose here is to cover the range of possible inputs as well
# as possible.
# [[0.1, 1.],  # [left, right] bounds, mu [MPa]
bounds = np.array([
    [14.5, 25],  # fixed skin tolerance (15mm to 30mm)
    # Location of the expander (2mm to 10mm from inside skin)
    [1, 6],
    [0, 0],  # theta (-15 to 15 degrees CCW)
    [1, 1], # Sx (expander profile variable, 0 to 1)
    [1.1567, 1.1567], # Critical Stretch Tcrit
    [0.150, 0.150]])  # Shear Modulus mu

# The number of inputs that you have:
num_inputs = bounds.shape[0]

# Now I sample designs in [0,1]^4 and I am going to scale them back to the
# original bounds
X_scaled = pd.lhs(num_inputs, num_simulations)

# Let's get to the original space
X = X_scaled * (bounds[:, 1] - bounds[:, 0]) + bounds[:, 0]

# add jobnumber to the output
jobnum = np.linspace(1, 15, 15).reshape(15, 1)

# Save to a text file
np.savetxt('P14_tol_h_theta_Sx_Tcrit_mu_15T_V1_TESTBATCH.txt',
           np.append(jobnum, X, axis=1), fmt='%.3f')
# You must do all these simulations...
