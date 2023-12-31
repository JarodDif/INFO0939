import pandas as pd
import matplotlib.pyplot as plt

data_ss_111 = pd.read_csv('./strong_scaling/data_111.csv')
data_ss_222 = pd.read_csv('./strong_scaling/data_222.csv')
data_ws     = pd.read_csv('./weak_scaling/data.csv')

plt.plot(data_ss_111["processes"], data_ss_111["updates"], label=r"$100^3$")
plt.plot(data_ss_222["processes"], data_ss_222["updates"], label=r"$200^3$")
plt.xlabel('Processes')
plt.ylabel('MUpdates/s')
plt.grid()
plt.legend()
plt.show()

plt.plot(data_ws["processes"], data_ws["updates"])
plt.xlabel('Processes')
plt.ylabel('MUpdates/s')
plt.grid()
plt.show()