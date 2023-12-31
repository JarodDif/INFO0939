import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data_ss_111 = pd.read_csv('./strong_scaling/data_111.csv')
data_ss_222 = pd.read_csv('./strong_scaling/data_222.csv')
data_ws     = pd.read_csv('./weak_scaling/data.csv')

ref_111     = data_ss_111["time"][0]
ref_222     = data_ss_222["time"][0]
ref_ws      = data_ss_222["time"][0]

amdahl      = lambda np, p: 1/((1-p) + p/np )

popt_ss_111 = curve_fit(amdahl, data_ss_111["processes"], ref_111/data_ss_111["time"])[0][0]
popt_ss_222 = curve_fit(amdahl, data_ss_222["processes"], ref_222/data_ss_222["time"])[0][0]
np_tab      = np.linspace(1, 64, 100)

plt.plot(data_ss_111["processes"], ref_111/data_ss_111["time"], "-ob", label=r"$100^3$")
plt.plot(np_tab, amdahl(np_tab, popt_ss_111), "--b", alpha=0.6, label=r"Amdahl's law with $p={:.2f}$".format(popt_ss_111))
plt.plot(data_ss_222["processes"], ref_222/data_ss_222["time"], "-og", label=r"$200^3$")
plt.plot(np_tab, amdahl(np_tab, popt_ss_222), "--g", alpha=0.6, label=r"Amdahl's law with $p={:.2f}$".format(popt_ss_222))

plt.xlabel('Number of processes')
plt.ylabel('Speedup')
plt.grid()
plt.legend()
plt.show()

plt.plot(data_ws["processes"], ref_ws/data_ws["time"])
plt.xlabel('Number of processes')
plt.ylabel('Speedup')
plt.grid()
plt.show()