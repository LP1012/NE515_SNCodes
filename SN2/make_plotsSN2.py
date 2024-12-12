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
    data = pd.read_csv(file,skiprows=3)
    sns.lineplot(x=data.iloc[:, 0],y=data.iloc[:, 1], label=file,linewidth=3)

    plt.grid("on")
    plt.legend()
    plt.xlabel("z")
    plt.ylabel("Scalar Flux")
    plt.title("SN2 Output Plots")
    plt.savefig("plots/" + f"{file}.png",dpi=300)
    plt.close()

print("")
print("Plots exported.")
print("")
print("")
