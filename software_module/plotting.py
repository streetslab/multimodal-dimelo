import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

_base_palette = sns.color_palette()

default_palette_map = {
    "A,0": _base_palette[0], # blue
    "WCG,1": _base_palette[3], # red
    "GCH,1": _base_palette[2] # green
}

def encode_narrowpeak_bed_qc(table: pd.DataFrame) -> None:
    _, axes = plt.subplots(1, 3,
                           subplot_kw=dict(box_aspect=0.85),
                           layout="constrained")
    """
    TODO: Right now, only using qValue, because it's always there in the files I have and correlates perfectly with pValue when that column is present.
    """

    sns.ecdfplot(table, x="signalValue", ax=axes[0])
    sns.ecdfplot(table, x="qValue", ax=axes[1])
    sns.scatterplot(table, x="signalValue", y="qValue", ax=axes[2])

    plt.show()
