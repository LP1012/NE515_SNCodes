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
# for file in glob.glob("*.csv"):
#     data = pd.read_csv(file,skiprows=3)
#     sns.lineplot(x=data.iloc[:, 0],y=data.iloc[:, 1], label=file)

file = "100_out.csv"
data = pd.read_csv(file,skiprows=3)
sns.lineplot(x=data.iloc[:, 0],y=data.iloc[:, 1],label=file,linewidth=4)

plt.legend()
plt.grid("on")
plt.xlabel("z")
plt.ylabel("Scalar Flux")
plt.title("SN1 Flux Output")
plt.savefig("plots/SN1_plot.png",dpi=300)
plt.close()

print("")
print("Plot exported.")
print("")
