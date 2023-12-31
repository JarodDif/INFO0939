import numpy as np
import matplotlib.pyplot as plt

h = [5e-3, 7.5e-3, 10e-3]

for _h in h:
    s_tab, rms_tab = np.loadtxt("empirical_study_"+str(_h).replace(".", "_"))
    plt.plot(s_tab, rms_tab, label=r"$\delta$"+f" = {_h:.2e}")

plt.axvline(3**-0.5, color="k", linestyle="--", alpha=0.8)
plt.grid()
plt.legend()
plt.xlim((1e-3, 9))
plt.xlabel(r"Courant number $Co = c (\Delta t) / (\Delta x)$ [-]")
plt.ylabel(r"RMS error estimator $\epsilon$ [-]")
plt.xscale("log")
plt.yscale("log")
x_bounds = plt.gca().get_xlim()

# Annotate the vertical line
plt.gca().annotate(text="$Co = 3^{-0.5}$", xy=((3**-0.5+0.01, 1e-1)), xytext=((1, 1e2)), alpha=0.8, arrowprops=dict(arrowstyle="->", alpha=0.8))

plt.show()
