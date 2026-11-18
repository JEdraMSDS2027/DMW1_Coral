"""Report figures: return Matplotlib figures without displaying or saving them."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

GROUP_ORDER = ["Low", "Moderate", "High"]
GROUP_COLORS = {"Low": "#a9c4e3", "Moderate": "#2a78d6", "High": "#0b2f66"}
C_POS, C_NEG, C_MUTED = "#2a78d6", "#d9622b", "#8a8984"


def configure_plotting():
    """Apply the original blue group palette and restrained report styling."""
    plt.rcParams.update(
        {
            "figure.dpi": 100,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.edgecolor": "#c9c8c2",
            "axes.grid": True,
            "grid.color": "#ecebe7",
            "axes.axisbelow": True,
            "axes.titlesize": 10,
            "axes.titleweight": "bold",
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.frameon": False,
            "legend.fontsize": 8,
        }
    )
    pd.set_option(
        "display.max_columns",
        12,
        "display.max_rows",
        70,
        "display.width",
        160,
        "display.float_format",
        "{:,.4g}".format,
    )


def effect_sizes(table, column="cohens_d", title="High vs Low bleaching", top=None):
    """Show signed differences; nominal q-values are used only as a legacy visual cue."""
    table = table.copy()
    if top:
        table = table.loc[table[column].abs().nlargest(top).index]
    table = table.iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.5, max(3, 0.3 * len(table) + 1.3)))
    colours = [C_POS if value >= 0 else C_NEG for value in table[column]]
    ax.barh(table.index, table[column], color=colours)
    ax.axvline(0, color="#444", lw=0.8)
    ax.set(xlabel=column.replace("_", " ") + " (High − Low)", title=title)
    ax.grid(axis="y", visible=False)
    fig.tight_layout()
    return fig


def time_score_plot(samples, scores, variance, first="PC1", second="PC2"):
    """Retain year colouring and cyclic month colouring within each hemisphere."""
    from matplotlib.colors import LinearSegmentedColormap

    cycle = LinearSegmentedColormap.from_list(
        "month_cycle", ["#0b2f66", "#2a78d6", "#d9a21b", "#d9622b", "#8a3a8f", "#0b2f66"]
    )
    years, months = samples.Date_Year.to_numpy(), samples.Date_Month.to_numpy()
    north = samples.Latitude_Degrees.ge(0).to_numpy()
    order = np.argsort(years, kind="stable")
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    sc = axes[0].scatter(
        scores[first].to_numpy()[order],
        scores[second].to_numpy()[order],
        c=years[order],
        s=4,
        alpha=0.5,
        cmap="viridis",
        edgecolors="none",
    )
    fig.colorbar(sc, ax=axes[0], label="year")
    axes[0].set_title("Coloured by sampling year")
    for ax, mask, label in [
        (axes[1], north, "northern hemisphere"),
        (axes[2], ~north, "southern hemisphere"),
    ]:
        sc = ax.scatter(
            scores.loc[mask, first],
            scores.loc[mask, second],
            c=months[mask],
            s=4,
            alpha=0.5,
            cmap=cycle,
            vmin=1,
            vmax=12,
            edgecolors="none",
        )
        fig.colorbar(sc, ax=ax, label="month", ticks=[1, 4, 7, 10, 12])
        ax.set_title(f"Coloured by month — {label} (n={mask.sum():,})")
    for ax in axes:
        ax.set(
            xlabel=f"{first} ({variance.loc[first,'variance_%']:.1f}%)",
            ylabel=f"{second} ({variance.loc[second,'variance_%']:.1f}%)",
        )
        ax.axhline(0, color="#bbb", lw=0.8)
        ax.axvline(0, color="#bbb", lw=0.8)
    fig.suptitle(f"{first} vs {second}", x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    return fig


def annual_patterns(samples, scores, min_samples=30):
    """Return original annual score median/IQR and bleaching-composition figures."""
    counts = samples.groupby("Date_Year").size()
    years = counts.index[counts >= min_samples]
    fig, axes = plt.subplots(1, 4, figsize=(16, 3.2))
    for ax, name in zip(axes, list(scores)[:4]):
        values = (
            scores[name].groupby(samples.Date_Year).quantile([0.25, 0.5, 0.75]).unstack().loc[years]
        )
        ax.fill_between(years, values[0.25], values[0.75], color=C_POS, alpha=0.2, lw=0)
        ax.plot(years, values[0.5], color=C_POS, lw=2, marker="o", ms=3)
        ax.axhline(0, color="#bbb", lw=0.8)
        ax.set_title(f"{name} by year (median, IQR band)")
    fig.tight_layout()
    composition, axes = plt.subplots(
        2, 1, figsize=(11, 5.2), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
    )
    shares = (
        pd.crosstab(samples.Date_Year, samples.Bleaching_Group, normalize="index").reindex(
            index=years, columns=GROUP_ORDER
        )
        * 100
    )
    bottom = np.zeros(len(shares))
    for name in GROUP_ORDER:
        axes[0].bar(
            years,
            shares[name],
            bottom=bottom,
            color=GROUP_COLORS[name],
            width=0.8,
            label=name,
            edgecolor="white",
            linewidth=0.6,
        )
        bottom += shares[name].to_numpy()
    axes[0].set(
        ylabel="% of samples",
        title=f"Bleaching groups by year (years with at least {min_samples} samples)",
    )
    axes[0].legend(ncol=3, loc="upper left", bbox_to_anchor=(0, 1.18))
    axes[1].bar(years, counts.loc[years], color=C_MUTED, width=0.8)
    axes[1].set(ylabel="samples", xlabel="year")
    composition.tight_layout()
    return fig, composition


def scree_plot(result):
    """Render the existing scree plot and return its figure(s)."""
    eigenvalues = result["model"].explained_variance_
    parallel_p95 = result["variance"].parallel_95th.to_numpy()
    var_ratio = result["model"].explained_variance_ratio_
    cumulative_ratio = np.cumsum(var_ratio)
    K_RETAIN = result["n_retained"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 3.6))
    ax = axes[0]
    k = np.arange(1, 21)
    ax.plot(k, eigenvalues[:20], "o-", color=C_POS, lw=2, ms=4, label="data")
    ax.plot(
        k,
        parallel_p95[:20],
        "--",
        color=C_MUTED,
        lw=1.5,
        label="parallel analysis (95th pct)",
    )
    ax.axhline(1, color=C_NEG, lw=1, ls=":")
    ax.text(20, 1.05, "eigenvalue = 1", ha="right", color=C_NEG, fontsize=8)
    ax.axvline(K_RETAIN + 0.5, color=C_MUTED, lw=1)
    ax.set_xticks(k)
    ax.set_xlabel("component")
    ax.set_ylabel("eigenvalue")
    ax.set_title("Scree plot (first 20 components)")
    ax.legend()
    ax = axes[1]
    ax.bar(
        np.arange(1, 21),
        100 * var_ratio[:20],
        color=C_POS,
        width=0.7,
        label="per component",
    )
    ax.plot(
        np.arange(1, 21),
        100 * cumulative_ratio[:20],
        "o-",
        color=C_NEG,
        ms=4,
        lw=2,
        label="cumulative",
    )
    ax.axhline(80, color=C_MUTED, lw=1, ls="--")
    ax.set_xticks(np.arange(1, 21))
    ax.set_xlabel("component")
    ax.set_ylabel("% of variance")
    ax.set_title("Variance explained")
    ax.legend(loc="center right")
    plt.tight_layout()

    return fig


def loading_heatmap(result):
    """Render the existing loading heatmap and return its figure(s)."""
    loadings = result["loadings"]
    var_ratio = result["model"].explained_variance_ratio_
    pc_names = list(loadings)
    K_RETAIN = result["n_retained"]
    FEATURE_ORDER = result["feature_order"]
    fam = result["families"]
    FAMILY_ORDER = result["family_order"]
    # ---- heatmap of loadings for the retained components ----
    show = pc_names[:K_RETAIN]
    L = loadings.loc[FEATURE_ORDER, show]
    fig, ax = plt.subplots(figsize=(1.0 + 0.9 * K_RETAIN, 0.27 * len(L) + 1.4))
    im = ax.imshow(L, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.grid(False)
    ax.set_xticks(range(K_RETAIN))
    ax.set_xticklabels([f"{c}\n{100*var_ratio[i]:.1f}%" for i, c in enumerate(show)])
    ax.set_yticks(range(len(L)))
    ax.set_yticklabels(L.index, fontsize=7)
    for i in range(L.shape[0]):
        for j in range(L.shape[1]):
            v = L.iloc[i, j]
            if abs(v) >= 0.4:
                ax.text(
                    j,
                    i,
                    f"{v:.2f}",
                    ha="center",
                    va="center",
                    fontsize=6.5,
                    color="white" if abs(v) > 0.7 else "black",
                )
    edges = (
        np.cumsum([sum(fam[c] == f for c in FEATURE_ORDER) for f in FAMILY_ORDER if f in set(fam)])[
            :-1
        ]
        - 0.5
    )
    for e in edges:
        ax.axhline(e, color="white", lw=2)
    fig.colorbar(im, ax=ax, shrink=0.6, label="loading (correlation with the component)")
    ax.set_title(
        "Loadings of the retained components (values shown where |loading| ≥ 0.4; "
        "rows grouped by family)"
    )
    plt.tight_layout()

    return fig


def strongest_loadings(result, n_interpret=5):
    """Render the existing strongest loadings and return its figure(s)."""
    N_INTERPRET = n_interpret
    loadings = result["loadings"]
    pc_names = list(loadings)
    var_ratio = result["model"].explained_variance_ratio_
    # ---- strongest variables per interpreted component ----
    fig, axes = plt.subplots(1, N_INTERPRET, figsize=(3.7 * N_INTERPRET, 3.8))
    for ax, c in zip(axes, pc_names[:N_INTERPRET]):
        top = (
            loadings[c]
            .reindex(loadings[c].abs().sort_values(ascending=False).index[:10])
            .iloc[::-1]
        )
        ax.barh(
            top.index,
            top.values,
            color=[C_POS if v > 0 else C_NEG for v in top],
            height=0.7,
        )
        ax.axvline(0, color="#999", lw=0.8)
        ax.set_xlim(-1, 1)
        ax.grid(axis="y", visible=False)
        ax.set_title(f"{c} — {100*var_ratio[pc_names.index(c)]:.1f}% of variance")
        ax.set_xlabel("loading")
        ax.tick_params(axis="y", labelsize=7)
    plt.tight_layout()

    return fig


def family_contributions(result):
    """Render the existing family contributions and return its figure(s)."""
    K_RETAIN = result["n_retained"]
    fam_contrib = result["family_contributions"]
    show = list(result["scores"])[:K_RETAIN]
    # ---- which environmental families does each component describe? ----
    fig, ax = plt.subplots(figsize=(1.0 + 0.85 * K_RETAIN, 0.42 * len(fam_contrib) + 1.2))
    fc = fam_contrib[show]
    ax.imshow(fc, cmap="Blues", vmin=0, vmax=max(60, fc.to_numpy().max()))
    ax.grid(False)
    ax.set_xticks(range(K_RETAIN))
    ax.set_xticklabels(show)
    ax.set_yticks(range(len(fc)))
    ax.set_yticklabels(fc.index)
    for i in range(fc.shape[0]):
        for j in range(fc.shape[1]):
            v = fc.iloc[i, j]
            ax.text(
                j,
                i,
                f"{v:.0f}",
                ha="center",
                va="center",
                fontsize=8,
                color="white" if v > 40 else "black",
            )
    ax.set_title("Share of each component (%) made up by each family of variables")
    plt.tight_layout()

    return fig


def score_plots(result, group_table, group_labels):
    """Render the existing score plots and return its figure(s)."""
    scores = result["scores"]
    var_ratio = result["model"].explained_variance_ratio_
    # ---- score plots coloured by bleaching group ----
    top2 = sorted(
        group_table.cliffs_delta.abs().sort_values(ascending=False).index[:2],
        key=lambda c: int(c[2:]),
    )  # the two best-separating components
    best1 = group_table.cliffs_delta.abs().idxmax()
    pairs = [("PC1", "PC2"), (top2[0], top2[1]), ("PC1", best1), ("PC2", best1)]
    pairs = [pr for i, pr in enumerate(pairs) if pr[0] != pr[1] and pr not in pairs[:i]]
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    axes = axes.ravel()
    for ax, (a, b) in zip(axes, pairs):
        for k in GROUP_ORDER:  # Low drawn first, High last so the rarer samples stay visible
            m = group_labels == k
            ax.scatter(
                scores.loc[m, a],
                scores.loc[m, b],
                s=4,
                alpha=0.35,
                color=GROUP_COLORS[k],
                edgecolors="none",
                label=k,
            )
        for k in GROUP_ORDER:
            m = group_labels == k
            ax.scatter(
                scores.loc[m, a].median(),
                scores.loc[m, b].median(),
                s=110,
                color=GROUP_COLORS[k],
                edgecolors="white",
                linewidths=2,
                zorder=5,
            )
        ax.set_xlabel(f"{a} ({100*var_ratio[int(a[2:])-1]:.1f}%)")
        ax.set_ylabel(f"{b} ({100*var_ratio[int(b[2:])-1]:.1f}%)")
        ax.axhline(0, color="#bbb", lw=0.8)
        ax.axvline(0, color="#bbb", lw=0.8)
    first_unused = len(pairs)
    for ax in axes[first_unused:]:
        ax.set_visible(False)
    axes[0].legend(markerscale=3, title="bleaching", loc="best")
    fig.suptitle(
        "PCA scores coloured by bleaching group (large dots = group medians)",
        x=0.01,
        ha="left",
        fontweight="bold",
        fontsize=11,
    )
    plt.tight_layout()

    return fig


def score_densities(result, group_table, group_labels, n_interpret=5):
    """Render the existing score densities and return its figure(s)."""
    N_INTERPRET = n_interpret
    scores = result["scores"]
    pc_names = list(scores)
    # ---- score densities (KDE) by bleaching group ----
    fig, axes = plt.subplots(1, N_INTERPRET, figsize=(3.4 * N_INTERPRET, 3.3))
    for ax, c in zip(axes, pc_names[:N_INTERPRET]):
        low, high = np.percentile(scores[c], [0.5, 99.5])
        grid = np.linspace(low, high, 300)
        for k in GROUP_ORDER:
            values = scores.loc[group_labels == k, c].to_numpy()
            density = stats.gaussian_kde(values)(grid)
            ax.plot(grid, density, color=GROUP_COLORS[k], lw=1.8, label=k)
            ax.fill_between(grid, density, color=GROUP_COLORS[k], alpha=0.15)
            ax.axvline(np.median(values), color=GROUP_COLORS[k], lw=0.8, ls="--")
        ax.set_title(
            f"{c}\nd = {group_table.cohens_d[c]:+.2f}," f"  δ = {group_table.cliffs_delta[c]:+.2f}"
        )
        ax.set_xlabel("score")
        ax.set_yticks([])
    axes[0].set_ylabel("density (each group integrates to 1)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, title="bleaching group", loc="lower center", ncol=3)
    plt.tight_layout(rect=(0, 0.14, 1, 1))

    return fig


def component_delta(group_table):
    """Render the existing component delta and return its figure(s)."""
    # ---- which components separate the groups? ----
    fig, ax = plt.subplots(figsize=(8, 3.2))
    ax.bar(
        group_table.index,
        group_table.cliffs_delta,
        color=[C_POS if v > 0 else C_NEG for v in group_table.cliffs_delta],
        width=0.65,
    )
    for y in (-0.33, -0.147, 0.147, 0.33):
        ax.axhline(y, color="#ccc", lw=0.8, ls="--")
    ax.axhline(0, color="#999", lw=0.8)
    ax.set_ylabel("Cliff's δ  (High vs Low)")
    ax.set_ylim(-0.7, 0.7)
    ax.set_title("How strongly each component separates High from Low bleaching")
    ax.grid(axis="x", visible=False)
    plt.tight_layout()

    return fig


def categorical_heatmap(result, omega2_table):
    """Render the existing categorical heatmap and return its figure(s)."""
    PCS = list(result["scores"])[: result["n_retained"]]
    var_ratio = result["model"].explained_variance_ratio_
    # Step 8 · How much of each component is explained by each categorical
    # feature?
    fig, ax = plt.subplots(figsize=(1.0 * len(PCS) + 3.2, 0.5 * len(omega2_table) + 1.4))
    ax.imshow(
        omega2_table.to_numpy(),
        cmap="Blues",
        vmin=0,
        vmax=max(40, omega2_table.to_numpy().max()),
        aspect="auto",
    )
    ax.set_xticks(range(len(PCS)))
    ax.set_xticklabels([f"{c}\n{100*var_ratio[i]:.0f}%" for i, c in enumerate(PCS)])
    ax.set_yticks(range(len(omega2_table)))
    ax.set_yticklabels(omega2_table.index)
    ax.grid(False)
    for i in range(omega2_table.shape[0]):
        for j in range(omega2_table.shape[1]):
            v = omega2_table.iat[i, j]
            ax.text(
                j,
                i,
                f"{v:.0f}" if v >= 1 else "<1",
                ha="center",
                va="center",
                fontsize=8,
                color="white" if v > 22 else "#222",
            )
    ax.set_title("% of each component's variance explained by each categorical feature (ω²)")
    plt.tight_layout()

    return fig


def source_delta_heatmap(pc_delta_table, within_source_pc, n_interpret=5):
    """Render the existing source delta heatmap and return its figure(s)."""
    N_INTERPRET = n_interpret
    comp_cols = [f"PC{i+1}" for i in range(n_interpret)]
    fig, ax = plt.subplots(figsize=(1.2 * N_INTERPRET + 3, 0.5 * len(pc_delta_table) + 1.3))
    M = pc_delta_table[list(comp_cols)].to_numpy(float)
    im = ax.imshow(M, cmap="RdBu_r", vmin=-0.8, vmax=0.8, aspect="auto")
    ax.grid(False)
    ax.set_xticks(range(len(comp_cols)))
    ax.set_xticklabels(list(comp_cols))
    ax.set_yticks(range(len(pc_delta_table)))
    ax.set_yticklabels(
        [
            (f"{i} (H={int(r.n_high)}, L={int(r.n_low)})" if i in within_source_pc.index else i)
            for i, r in pc_delta_table.fillna(0).iterrows()
        ]
    )
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(
                j,
                i,
                f"{M[i, j]:+.2f}",
                ha="center",
                va="center",
                fontsize=8,
                color="white" if abs(M[i, j]) > 0.45 else "#222",
            )
    fig.colorbar(im, ax=ax, label="Cliff's δ (High vs Low)")
    ax.set_title("High-vs-Low separation of each component, within each data source")
    plt.tight_layout()

    return fig


def source_distributions(df, scores, within_source_pc):
    """Render the existing source distributions and return its figure(s)."""
    group_labels = df.Bleaching_Group.to_numpy()
    # ---- distributions of PC3 / PC4 by bleaching group inside each source ----
    fig, axes = plt.subplots(
        2,
        len(within_source_pc),
        figsize=(3.0 * len(within_source_pc), 5.4),
        sharey="row",
        squeeze=False,
    )
    for j, s in enumerate(within_source_pc.index):
        for i, c in enumerate(["PC3", "PC4"]):
            data = [scores.loc[(df.Data_Source == s) & (group_labels == k), c] for k in GROUP_ORDER]
            bp = axes[i, j].boxplot(
                data,
                widths=0.6,
                showfliers=False,
                patch_artist=True,
                medianprops=dict(color="white", lw=2),
            )
            for patch, k in zip(bp["boxes"], GROUP_ORDER):
                patch.set_facecolor(GROUP_COLORS[k])
                patch.set_edgecolor("none")
            axes[i, j].set_xticks(range(1, 4))
            axes[i, j].set_xticklabels(GROUP_ORDER)
            axes[i, j].grid(axis="x", visible=False)
            axes[i, j].set_title(f"{s} — {c}  (δ={within_source_pc.loc[s, c]:+.2f})", fontsize=8.5)
    plt.tight_layout()

    return fig
