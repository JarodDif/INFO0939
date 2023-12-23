import numpy as np
import matplotlib.pyplot as plt

s_tab, rms_tab = np.loadtxt("empirical_study")

plt.plot(s_tab, rms_tab)
plt.axvline(1/np.sqrt(3), color="k", linestyle="--", alpha=0.5, label="$s = 3^{-0.5}$")
plt.grid()
plt.legend()
plt.xlabel(r"Convergence factor $s = c (\Delta x) / (\Delta t)$ [-]")
plt.ylabel(r"RMS error estimator $\epsilon$ [-]")
plt.xscale("log")
plt.yscale("log")
plt.show()