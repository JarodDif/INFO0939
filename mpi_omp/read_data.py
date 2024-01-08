import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data_ss_111 = pd.read_csv('./strong_scaling/data_111.csv')
data_ss_222 = pd.read_csv('./strong_scaling/data_222.csv')
data_ws_1e6 = pd.read_csv('./weak_scaling/data_1e6.csv')
data_ws_8e6 = pd.read_csv('./weak_scaling/data_8e6.csv')

ref_111     = data_ss_111["time"][0]
ref_222     = data_ss_222["time"][0]
ref_1e6     = data_ws_1e6["time"][0]
ref_8e6     = data_ws_8e6["time"][0]

amdahl      = lambda np, p: 1/((1-p) + p/np )
gustafson   = lambda np, p: (1-p) + p * np

popt_ss_111 = curve_fit(amdahl, data_ss_111["processes"], ref_111/data_ss_111["time"])[0][0]
popt_ss_222 = curve_fit(amdahl, data_ss_222["processes"], ref_222/data_ss_222["time"])[0][0]

popt_ws_1e6 = curve_fit(gustafson, data_ws_1e6["processes"], ref_1e6/data_ws_1e6["time"]*data_ws_1e6["processes"])[0][0]
popt_ws_8e6 = curve_fit(gustafson, data_ws_8e6["processes"], ref_8e6/data_ws_8e6["time"]*data_ws_8e6["processes"])[0][0]

np_tab      = np.linspace(1, 64, 100)

plt.plot(data_ss_111["processes"], ref_111/data_ss_111["time"], "-ob", label=r"$100^3$")
plt.plot(np_tab, amdahl(np_tab, popt_ss_111), "--b", alpha=0.6, label=r"Amdahl's law with $p={:.2f}$".format(popt_ss_111))
plt.plot(data_ss_222["processes"], ref_222/data_ss_222["time"], "-og", label=r"$200^3$")
plt.plot(np_tab, amdahl(np_tab, popt_ss_222), "--g", alpha=0.6, label=r"Amdahl's law with $p={:.2f}$".format(popt_ss_222))

plt.xlabel('Number of cores used')
plt.ylabel('Speedup')
plt.grid()
plt.legend()
plt.show()

plt.plot(data_ws_1e6["processes"], ref_1e6/data_ws_1e6["time"]*data_ws_1e6["processes"], "-ob", label="$10^6$ elements per process")
plt.plot(np_tab, gustafson(np_tab, popt_ws_1e6), "--b", label=f"Gustafson's law with $p ={popt_ws_1e6:.2f}$")
plt.plot(data_ws_8e6["processes"], ref_8e6/data_ws_8e6["time"]*data_ws_8e6["processes"], "-og", label="$8 \\ 10^6$ elements per process")
plt.plot(np_tab, gustafson(np_tab, popt_ws_8e6), "--g", label=f"Gustafson's law with $p ={popt_ws_8e6:.2f}$")
plt.xlabel('Number of cores used')
plt.ylabel('Scaled Speedup')
plt.grid()
plt.legend()
plt.show()