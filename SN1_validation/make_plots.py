import matplotlib.pyplot as plt
import pandas as pd
import glob
import seaborn as sns

sns.set_theme(
    context="paper",
    style="ticks",
    palette="deep",
    font="serif",
    font_scale=1,
    color_codes=True,
    rc=None,
)
sns.despine()

plt.figure()
for file in glob.glob("*.csv"):
    if file not in ("MMS_L2errors.csv","MMS_true.csv"):
        data = pd.read_csv(file,skiprows=3)
        sns.lineplot(x=data.iloc[:, 0],y=data.iloc[:, 1], label=file)

plt.legend()
plt.xlabel("z")
plt.ylabel("Scalar Flux")
plt.title("MMS Convergence Plot")
plt.savefig("MMS_plot.png",dpi=300)
plt.close()

data_errors = pd.read_csv("MMS_L2errors.csv",skiprows=2)
plt.figure()
sns.lineplot(x=data_errors.iloc[:,0],y=data_errors.iloc[:,1], marker="o",markersize=8,linewidth=2)
plt.xlabel("Number of Cells")
plt.title("L2 Error of MMS Solution")
plt.yscale("log")
plt.xscale("log")
plt.grid(True,which="both")
plt.savefig("L2error_plot.png",dpi=300)
plt.close()

print("Plots exported.")
print(" ")
