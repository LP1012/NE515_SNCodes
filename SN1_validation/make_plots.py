import matplotlib.pyplot as plt
import pandas as pd
import glob
import seaborn as sns
import numpy as np

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
plt.savefig("plots/MMS_plot.png",dpi=300)
plt.close()

# New dataframe for L2 errors

data_errors = pd.read_csv("MMS_L2errors.csv",skiprows=2)

# Create reference curve
x_refs = data_errors.iloc[:,0].to_numpy()
xs = np.array([1/x for x in x_refs])

y_quads = np.array([x**2 for x in xs])
y_linears = xs

# Create dz values
ys = data_errors.iloc[:,1].to_numpy()


plt.figure()
sns.lineplot(x=xs,y=ys, marker="o",markersize=8,linewidth=2,label='Experimental Error')
sns.lineplot(x=xs, y=y_quads, linewidth=2,label=r"Reference Curve $y=x^2$",linestyle="dashed")
sns.lineplot(x=xs, y=y_linears,label=r"Reference Curve $y=x$",linewidth=2,linestyle='dashed')
plt.legend()
plt.xlabel(r"$\Delta z$")
plt.ylabel("L2-norm of Error")
plt.title("L2-norm of Error of MMS Solution Against Reference Curve")
plt.yscale("log")
plt.xscale("log")
plt.grid(True,which="both")
plt.savefig("plots/L2error_plot.png",dpi=300)
plt.close()

print(" ")
print("Plots exported.")
print(" ")
