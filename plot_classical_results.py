import numpy as np
import matplotlib.pyplot as plt

# Load data
data = np.loadtxt("output.dat")

# Columns
it   = data[:,0].astype(int)
t    = data[:,1]
ix   = data[:,2].astype(int)
x    = data[:,3]
f_p  = data[:,4]
f_m  = data[:,5]
J    = data[:,6]
H    = data[:,7]
K    = data[:,8]

# Get unique timesteps
timesteps = np.unique(it)

# Pick the LAST timestep
final_it = timesteps[-1]

mask = (it == final_it)

x_plot = x[mask]
f_p_plot = f_p[mask]
f_m_plot = f_m[mask]
J_plot = J[mask]

# Make a movie!
plt.figure()

for step in timesteps:  
    mask = (it == step)

    plt.clf()
    plt.plot(x[mask], f_p[mask], label="f_plus")
    plt.plot(x[mask], f_m[mask], label="f_minus")

    plt.title(f"it={step}")
    plt.ylim(0,2.0)
    plt.legend()
    plt.pause(0.05)

plt.show()