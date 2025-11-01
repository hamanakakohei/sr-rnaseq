#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# フォント全体設定
plt.rcParams.update({
    "font.family": "DejaVu Sans", # "Arial",
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "xtick.major.pad": 1,
    "ytick.major.pad": 1,
    #"axes.titleweight": 'bold',
})

def adjust_axes(ax):
    ax.set_ylabel(ax.get_ylabel(), labelpad=2)
    ax.set_xlabel(ax.get_xlabel(), labelpad=2)
    ax.tick_params(axis="y", pad=1)
    ax.tick_params(axis="x", pad=1)


highlight_gene = "SERPINB7"

datasets = {
    "A": {
        "DEG_csv": "data/GSE140990_signatureData.csv",
        "TMM_csv": "data/GSE140990_GeneLevel_Normalized(CPM.and.TMM)_data.csv",
        "ctrl_samples": ['GSM4192039', 'GSM4192040'],
        "kd_samples": ['GSM4192043', 'GSM4192044']
    },
    "B": {
        "DEG_csv": "data/GSE137410_signatureData.csv",
        "TMM_csv": "data/GSE137410_GeneLevel_Normalized(CPM.and.TMM)_data.csv",
        "ctrl_samples": ['GSM4078731', 'GSM4078732', 'GSM4078733'],
        "kd_samples": ['GSM4078734', 'GSM4078735', 'GSM4078736', 'GSM4078737']
    }
}

fig, axes = plt.subplots(
    1, 4, figsize=(7, 3.5),
    gridspec_kw={"width_ratios": [4, 1, 4, 1]},  # MAプロットは横幅3倍
    #constrained_layout=True
)
#plt.subplots_adjust(wspace=0.0)


for i, panel_label in enumerate(["A", "B"]):
    data = datasets[panel_label]

    # DEGとTMM読み込み
    df = pd.read_csv(data["DEG_csv"])
    exp_df = pd.read_csv(data["TMM_csv"])
    gene_exp = exp_df[exp_df["gene_symbol"] == highlight_gene]
    control_vals = gene_exp[data["ctrl_samples"]].values.flatten()
    kd_vals = gene_exp[data["kd_samples"]].values.flatten()

    # MAプロット
    ax_ma = axes[i*2]  # 0, 2
    ax_ma.scatter(df["logCPM"], df["Log_FoldChange"], color="grey", alpha=0.3, s=10)

    highlight_row = df[df["Gene_symbol"] == highlight_gene]
    x_val = highlight_row["logCPM"].values[0]
    y_val = highlight_row["Log_FoldChange"].values[0]
    ax_ma.scatter(x_val, y_val, color="black", s=50, zorder=5)

    # 点から右下に線
    line_dx = 1.0
    line_dy = -4.0
    #ax_ma.plot([x_val, x_val + line_dx], [y_val, y_val + line_dy],
    #           color="black", linewidth=1.5)
    ax_ma.text(x_val + line_dx, y_val + line_dy, r"$\it{" + highlight_gene + "}$",
               color="black", fontsize=7, va="center", ha="left")

    ax_ma.axhline(0, linestyle="--", color="black", linewidth=0.8)
    ax_ma.set_xlabel("logCPM")
    ax_ma.set_ylabel("Log_FoldChange")

    # 軸目盛り
    y_lim = np.max(np.abs(df["Log_FoldChange"]))
    ax_ma.set_ylim(-y_lim, y_lim)
    ax_ma.xaxis.set_major_formatter('{:.1f}'.format)
    ax_ma.yaxis.set_major_formatter('{:.1f}'.format)

    # パネルラベル
    ax_ma.text(-0.15, 1.12, panel_label, transform=ax_ma.transAxes,
               fontsize=14, fontweight="bold", va="top", ha="right")

    # ここでタイトル追加
    ax_ma.set_title(f"{data['DEG_csv'].split('/')[-1].replace('_signatureData.csv','')}", fontsize=12)

    # 棒グラフ
    ax_bar = axes[i*2 + 1]  # 1, 3
    means = [np.mean(control_vals), np.mean(kd_vals)]
    #labels = [f"Control\n(n = {len(control_vals)})", f"KLF4 KD\n(n = {len(kd_vals)})"]
    labels = [f"Control", f"KLF4 KD"]
    bars = ax_bar.bar(labels, means, alpha=0.8) #, width=0.5)
    plt.setp(ax_bar.get_xticklabels(), rotation=315, ha="left")

    # レプリケート
    x_positions = [0, 1]
    ax_bar.scatter([x_positions[0]]*len(control_vals), control_vals,
                   edgecolor="black", facecolors='none', zorder=10)
    ax_bar.scatter([x_positions[1]]*len(kd_vals), kd_vals,
                   edgecolor="black", facecolors='none', zorder=10)

    ax_bar.set_ylabel(r"$\it{" + highlight_gene + "}$ expression level")

    # P値
    pval = highlight_row["PValue"].values[0]
    y_max = max(means)
    bar_top = y_max * 1.1
    ax_bar.plot([0, 1], [bar_top, bar_top], color="black")
    ax_bar.text(0.5, bar_top * 1.02, f"P =\n{pval:.0e}", ha="center", va="bottom", fontsize=7)

    # 枠線
    ax_bar.spines['top'].set_visible(False)
    ax_bar.spines['right'].set_visible(False)

for ax in axes.ravel():
    adjust_axes(ax)

plt.tight_layout(pad=0.2, w_pad=0.2)
# pad: 図全体と外枠の余白
# w_pad: 横方向の隙間
# h_pad: 縦方向の隙間


plt.savefig("figure1.pdf", dpi=300, bbox_inches="tight")
plt.close()
