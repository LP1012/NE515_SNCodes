import matplotlib.pyplot as plt
import pandas as pd
import glob
import seaborn as sns
import os

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
    data = pd.read_csv(file,skiprows=3)

    base_name, ext = os.path.splitext(file)

    plt.figure()
    sns.lineplot(x=data.iloc[:, 0],y=data.iloc[:, 1], label=base_name,linewidth=4)
    plt.legend()
    plt.grid("on")
    plt.xlabel("z")
    plt.ylabel("Scalar Flux")
    plt.title(base_name + " Flux Output")
    plt.savefig("plots/" + base_name + "_plot.png",dpi=300)
    plt.close()

# file = "100_out.csv"
# data = pd.read_csv(file,skiprows=3)
# sns.lineplot(x=data.iloc[:, 0],y=data.iloc[:, 1],label=file,linewidth=4)



print("")
print("Plot exported.")
print("")
