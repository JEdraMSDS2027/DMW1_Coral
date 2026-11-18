"""Effect sizes and statistical comparisons consolidated from the original notebooks.

Sample-level p-values and legacy independent-site Welch results are descriptive
comparators. The site-cluster bootstrap provides the corrected uncertainty audit.
"""

import numpy as np
import pandas as pd
from scipy import stats


def auc_high_low(values, groups):
    """Probability that a High sample scores above a Low sample.

    ``values`` is a score per sample and ``groups`` the bleaching group of
    each sample ("Low", "Moderate" or "High"). 0.5 means no separation.
    """
    values, groups = np.asarray(values, float), np.asarray(groups)
    high, low = values[groups == "High"], values[groups == "Low"]
    if not len(high) or not len(low):
        raise ValueError("AUC requires both High and Low observations.")
    u_stat = stats.mannwhitneyu(high, low, alternative="two-sided", method="asymptotic").statistic
    return u_stat / (len(high) * len(low))


def omega_squared(values, labels, min_n=50):
    """Bias-corrected share of variance in ``values`` between categories.

    Categories with fewer than ``min_n`` samples are pooled into "Other".
    The result is never negative.
    """
    labels = pd.Series(labels).astype(str)
    counts = labels.value_counts()
    labels = labels.where(labels.map(counts) >= min_n, "Other")
    data = pd.DataFrame({"x": np.asarray(values, float), "group": labels.to_numpy()})
    grand_mean = data.x.mean()
    n_obs, n_groups = len(data), data.group.nunique()
    if n_groups < 2 or n_obs <= n_groups or data.x.var() == 0:
        return 0.0
    between = data.groupby("group").x.agg(lambda v: len(v) * (v.mean() - grand_mean) ** 2).sum()
    total = ((data.x - grand_mean) ** 2).sum()
    within_mean_sq = (total - between) / (n_obs - n_groups)
    numerator = between - (n_groups - 1) * within_mean_sq
    omega = numerator / (total + within_mean_sq)
    return max(0.0, omega)


def high_low_u(values, mask, groups):
    """Mann-Whitney U statistic and pair count for High vs Low in ``mask``.

    Returns ``(u_statistic, n_pairs)``; ``(0, 0)`` if either group is empty.
    """
    high = values[mask & (groups == "High")]
    low = values[mask & (groups == "Low")]
    if len(high) == 0 or len(low) == 0:
        return 0, 0
    u_stat = stats.mannwhitneyu(high, low, alternative="two-sided", method="asymptotic").statistic
    return u_stat, len(high) * len(low)


def within_stratum_delta(frame, groups, strata_col, columns, min_n=30):
    """High-vs-Low Cliff's delta inside each level of a categorical column.

    ``columns`` maps a name to an array of scores (one per row of
    ``frame``). Only levels with at least ``min_n`` High and ``min_n`` Low
    samples are kept. Returns the per-level table and a dictionary with the
    stratified (source-adjusted) delta per column, in which each level is
    weighted by its number of High x Low pairs.
    """
    labels = frame[strata_col]
    levels = [
        level
        for level in labels.value_counts().index
        if ((labels == level) & (groups == "High")).sum() >= min_n
        and ((labels == level) & (groups == "Low")).sum() >= min_n
    ]
    rows, u_total, pair_total = [], {name: 0.0 for name in columns}, 0
    for level in levels:
        mask = (labels == level).to_numpy()
        n_high = int((mask & (groups == "High")).sum())
        n_low = int((mask & (groups == "Low")).sum())
        row = {strata_col: level, "n_high": n_high, "n_low": n_low}
        for name, scores in columns.items():
            u_stat, n_pairs = high_low_u(scores, mask, groups)
            row[name] = 2 * u_stat / n_pairs - 1
            u_total[name] += u_stat
        pair_total += n_high * n_low
        rows.append(row)
    if not rows:
        return pd.DataFrame(columns=["n_high", "n_low", *columns]).rename_axis(strata_col), {
            name: np.nan for name in columns
        }
    table = pd.DataFrame(rows).set_index(strata_col)
    adjusted = {name: 2 * u_total[name] / pair_total - 1 for name in columns}
    return table, adjusted


def benjamini_hochberg(p_values):
    """Benjamini-Hochberg adjusted p-values (q-values) controlling the FDR.

    ``p_values`` is a one-dimensional array of raw p-values. The result has
    the same order; a q-value below 0.05 means the finding stays significant
    while keeping the expected share of false discoveries at 5 %.
    """
    p_values = np.asarray(p_values, dtype=float)
    n_tests = len(p_values)
    order = np.argsort(p_values)
    scaled = p_values[order] * n_tests / np.arange(1, n_tests + 1)
    adjusted = np.minimum.accumulate(scaled[::-1])[::-1]
    result = np.empty(n_tests)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def welch_high_low(columns, groups, sites=None):
    """Welch t-test of High versus Low bleaching for several columns.

    ``columns`` maps a name to an array with one value per sample and
    ``groups`` gives the bleaching group of each sample. Returns one row per
    column with the mean difference (High - Low), its 95 % confidence
    interval, Cohen's d, the Welch p-value and the Benjamini-Hochberg q-value
    over all columns. If ``sites`` (a site label per sample) is given, the
    test is repeated on site means (the mean of the High samples of each
    site against the mean of the Low samples of each site), a legacy calculation that assumes the two sets of site means are independent.
    Shared sites violate that assumption; use the cluster-bootstrap audit for
    corrected uncertainty, and do not treat this comparator as independent-site inference.
    """
    is_high, is_low = groups == "High", groups == "Low"
    rows = []
    for name, values in columns.items():
        values = np.asarray(values, dtype=float)
        high, low = values[is_high], values[is_low]
        test = stats.ttest_ind(high, low, equal_var=False)
        interval = test.confidence_interval()
        pooled_sd = np.sqrt(
            ((len(high) - 1) * high.var(ddof=1) + (len(low) - 1) * low.var(ddof=1))
            / (len(high) + len(low) - 2)
        )
        row = {
            "variable": name,
            "mean_high": high.mean(),
            "mean_low": low.mean(),
            "difference": high.mean() - low.mean(),
            "ci95_low": interval.low,
            "ci95_high": interval.high,
            "cohens_d": (high.mean() - low.mean()) / pooled_sd,
            "t": test.statistic,
            "df": test.df,
            "p": test.pvalue,
        }
        if sites is not None:
            frame = pd.DataFrame({"value": values, "site": sites, "group": groups})
            site_means = frame.groupby(["group", "site"]).value.mean()
            site_high = site_means.loc["High"].to_numpy()
            site_low = site_means.loc["Low"].to_numpy()
            site_test = stats.ttest_ind(site_high, site_low, equal_var=False)
            site_sd = np.sqrt(
                (
                    (len(site_high) - 1) * site_high.var(ddof=1)
                    + (len(site_low) - 1) * site_low.var(ddof=1)
                )
                / (len(site_high) + len(site_low) - 2)
            )
            row["site_level_d"] = (site_high.mean() - site_low.mean()) / site_sd
            row["site_level_p"] = site_test.pvalue
        rows.append(row)
    table = pd.DataFrame(rows).set_index("variable")
    table["q_bh"] = benjamini_hochberg(table["p"].to_numpy())
    if sites is not None:
        table["site_level_q_bh"] = benjamini_hochberg(table["site_level_p"].to_numpy())
    return table


def hotelling_two_sample(high, low):
    """Hotelling's two-sample T-squared test for High versus Low.

    ``high`` and ``low`` are arrays (samples x variables). Returns the
    squared Mahalanobis distance between the two group means (a multivariate
    effect size that does not depend on the sample size), the F statistic and
    its p-value. A generalised inverse is used so that strongly correlated
    variables do not break the calculation.
    """
    n_high, n_low, n_vars = len(high), len(low), high.shape[1]
    diff = high.mean(axis=0) - low.mean(axis=0)
    pooled = (
        (n_high - 1) * np.cov(high, rowvar=False) + (n_low - 1) * np.cov(low, rowvar=False)
    ) / (n_high + n_low - 2)
    mahalanobis_sq = float(diff @ np.linalg.pinv(pooled) @ diff)
    t_squared = mahalanobis_sq * n_high * n_low / (n_high + n_low)
    df_2 = n_high + n_low - n_vars - 1
    f_stat = t_squared * df_2 / ((n_high + n_low - 2) * n_vars)
    return {
        "mahalanobis_sq": mahalanobis_sq,
        "F": f_stat,
        "df1": n_vars,
        "df2": df_2,
        "p": float(stats.f.sf(f_stat, n_vars, df_2)),
    }


def site_cluster_intervals(columns: dict, groups, sites, resamples=1000, seed=42):
    """Correct uncertainty for overlapping sites using a site-cluster bootstrap.

    Estimand: difference of High and Low site means, giving each represented site
    equal weight in each group, as in the legacy site-level comparison. Resample
    the union of sites with common weights in both groups, retaining their overlap.
    Percentile 95% intervals are pointwise, not multiplicity-adjusted. PCA and
    feature transforms remain fixed; these intervals do not include fitting uncertainty.
    """
    if resamples < 100:
        raise ValueError("Use at least 100 bootstrap resamples.")
    values = pd.DataFrame(columns)
    values["site"], values["group"] = np.asarray(sites), np.asarray(groups)
    site_means = values.groupby(["site", "group"])[list(columns)].mean()
    high = site_means.xs("High", level="group")
    low = site_means.xs("Low", level="group")
    union = high.index.union(low.index)
    high = high.reindex(union)
    low = low.reindex(union)
    has_high, has_low = high.notna().to_numpy(float), low.notna().to_numpy(float)
    high_values, low_values = high.fillna(0).to_numpy(), low.fillna(0).to_numpy()
    rng = np.random.default_rng(seed)
    draws = np.empty((resamples, len(columns)))
    for repeat in range(resamples):
        weights = np.bincount(rng.integers(len(union), size=len(union)), minlength=len(union))
        draws[repeat] = (weights @ high_values) / (weights @ has_high) - (weights @ low_values) / (
            weights @ has_low
        )
    if not np.isfinite(draws).all():
        raise ValueError("Some bootstrap samples lack a comparison group; review cluster counts.")
    bounds = np.quantile(draws, [0.025, 0.975], axis=0)
    return pd.DataFrame(
        {
            "difference_of_site_means": high.mean() - low.mean(),
            "cluster_ci95_low": bounds[0],
            "cluster_ci95_high": bounds[1],
            "sites_high": has_high[:, 0].sum().astype(int),
            "sites_low": has_low[:, 0].sum().astype(int),
            "sites_both": int((has_high[:, 0] * has_low[:, 0]).sum()),
            "bootstrap_resamples": resamples,
        }
    ).rename_axis("variable")


def component_associations(scores, samples, variance, groups_order):
    """Combine existing group medians, weighted Cohen d, AUC/delta, and epsilon squared."""
    groups = samples.Bleaching_Group.to_numpy()
    effects = welch_high_low({name: scores[name].to_numpy() for name in scores}, groups)
    rows = []
    for name in scores:
        values = scores[name].to_numpy()
        subsets = [values[groups == group] for group in groups_order]
        auc = auc_high_low(values, groups)
        rows.append(
            {
                "component": name,
                "variance_%": variance.loc[name, "variance_%"],
                **{
                    f"median_{group}": np.median(part) for group, part in zip(groups_order, subsets)
                },
                "AUC_high_vs_low": auc,
                "cliffs_delta": 2 * auc - 1,
                "cohens_d": effects.loc[name, "cohens_d"],
                "epsilon2": (stats.kruskal(*subsets).statistic - len(subsets) + 1)
                / (len(values) - len(subsets)),
                "spearman_with_bleaching_%": stats.spearmanr(
                    values, samples.Percent_Bleaching
                ).statistic,
            }
        )
    result = pd.DataFrame(rows).set_index("component")
    result["strength"] = pd.cut(
        result.cliffs_delta.abs(),
        [-1, 0.147, 0.33, 0.474, 1],
        labels=["negligible", "small", "medium", "large"],
    )
    return result


def variable_associations(frame, samples, families):
    """Preserve older variable-wise AUC and equal-weight-SD d; label its convention."""
    groups = samples.Bleaching_Group.to_numpy()
    rows = []
    for name in frame:
        x = frame[name].to_numpy()
        high, low = x[groups == "High"], x[groups == "Low"]
        rows.append(
            {
                "variable": name,
                "family": families[name],
                "legacy_unweighted_d": (high.mean() - low.mean())
                / np.sqrt((high.var(ddof=1) + low.var(ddof=1)) / 2),
                "AUC_high_vs_low": auc_high_low(x, groups),
                **{
                    f"mean_{group.lower()}": x[groups == group].mean()
                    for group in ["Low", "Moderate", "High"]
                },
            }
        )
    return (
        pd.DataFrame(rows)
        .set_index("variable")
        .sort_values("legacy_unweighted_d", key=abs, ascending=False)
    )


def hotelling_comparison(frame, scores, groups):
    """Preserve nominal sample-level Hotelling tests and D² compression ratios."""
    rows = []
    matrices = {
        "PC1–PC2": scores.iloc[:, :2],
        f"PC1–PC{scores.shape[1]} (retained)": scores,
        f"all {frame.shape[1]} original variables": frame,
    }
    for name, matrix in matrices.items():
        values = matrix.to_numpy()
        rows.append(
            {
                "inputs": name,
                "n_features": matrix.shape[1],
                **hotelling_two_sample(values[groups == "High"], values[groups == "Low"]),
            }
        )
    result = pd.DataFrame(rows).set_index("inputs")
    result["share_of_all_variables_D2"] = result.mahalanobis_sq / result.mahalanobis_sq.iloc[-1]
    return result


def depth_sensitivity(samples, scores):
    """Original imputed-versus-measured score delta and measured-depth-only bleaching delta."""
    missing = samples.Depth_m_was_missing.eq(1).to_numpy()
    groups = samples.Bleaching_Group.to_numpy()
    rows = []
    for name in scores:
        x = scores[name].to_numpy()
        u = stats.mannwhitneyu(x[missing], x[~missing], method="asymptotic").statistic
        rows.append(
            {
                "component": name,
                "imputed_vs_measured_delta": 2 * u / (missing.sum() * (~missing).sum()) - 1,
                "all_events_delta": 2 * auc_high_low(x, groups) - 1,
                "measured_depth_delta": 2 * auc_high_low(x[~missing], groups[~missing]) - 1,
            }
        )
    return pd.DataFrame(rows).set_index("component")


def temporal_associations(samples, scores, periods, min_group=30):
    """Preserve year correlation, hemisphere seasonal R², and within-period deltas."""
    rows, period_rows = [], []
    for name in scores:
        x = scores[name].to_numpy()
        row = {
            "component": name,
            "spearman_with_year": stats.spearmanr(x, samples.Date_Year).statistic,
        }
        north = samples.Latitude_Degrees.ge(0).to_numpy()
        for label, mask in [("north", north), ("south", ~north)]:
            angle = 2 * np.pi * samples.Date_Month.to_numpy()[mask] / 12
            design = np.column_stack([np.ones(mask.sum()), np.sin(angle), np.cos(angle)])
            residual = x[mask] - design @ np.linalg.lstsq(design, x[mask], rcond=None)[0]
            row[f"seasonal_R²_{label}"] = 1 - np.var(residual) / np.var(x[mask])
        rows.append(row)
    groups = samples.Bleaching_Group.to_numpy()
    for lo, hi in periods:
        mask = samples.Date_Year.between(lo, hi).to_numpy()
        nh, nl = int((mask & (groups == "High")).sum()), int((mask & (groups == "Low")).sum())
        period_rows.append(
            {
                "period": f"{lo}–{hi}",
                "samples": int(mask.sum()),
                "High": nh,
                "Low": nl,
                **{
                    name: (
                        2 * auc_high_low(scores.loc[mask, name].to_numpy(), groups[mask]) - 1
                        if min(nh, nl) >= min_group
                        else np.nan
                    )
                    for name in scores
                },
            }
        )
    return pd.DataFrame(rows).set_index("component"), pd.DataFrame(period_rows).set_index("period")


def grouped_prediction(frame, scores, samples, raw_measurements=None, folds=5, seed=42):
    """Compare original diagnostic CV with fold-fitted preprocessing/PCA correction.

    The cohort, named study features, and requested retained dimension are fixed.
    Corrected evaluation refits medians, skew selection, power transforms, scales,
    and PCA inside every training fold; labels never fit unsupervised steps.
    """
    from sklearn.model_selection import GroupKFold
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import balanced_accuracy_score, roc_auc_score
    from sklearn.decomposition import PCA
    from .preprocessing import fit_preprocessing, transform_features

    order = ["Low", "Moderate", "High"]
    labels, sites = samples.Bleaching_Group.to_numpy(), samples.Site_ID.to_numpy()
    modes = ["legacy_global_fit"] + (["fold_fitted"] if raw_measurements is not None else [])
    rows, fold_rows = [], []
    for mode in modes:
        predictions = {
            name: np.empty(len(labels), object)
            for name in ["PC1–PC2", "retained PCs", "all variables"]
        }
        probabilities = {name: np.zeros((len(labels), 3)) for name in predictions}
        for fold, (train, test) in enumerate(GroupKFold(folds).split(frame, labels, groups=sites)):
            if set(sites[train]) & set(sites[test]):
                raise AssertionError("A site appears in both training and test data.")
            if mode == "fold_fitted":
                _, tr, model, _ = fit_preprocessing(raw_measurements.iloc[train], list(frame))
                te = transform_features(raw_measurements.iloc[test], model)
                pca = PCA(n_components=scores.shape[1], random_state=seed).fit(tr)
                ts, vs = pca.transform(tr), pca.transform(te)
            else:
                tr, te = frame.iloc[train], frame.iloc[test]
                ts, vs = scores.iloc[train].to_numpy(), scores.iloc[test].to_numpy()
            matrices = {
                "PC1–PC2": (ts[:, :2], vs[:, :2]),
                "retained PCs": (ts, vs),
                "all variables": (tr, te),
            }
            for name, (a, b) in matrices.items():
                fitted = LogisticRegression(max_iter=2000, class_weight="balanced").fit(
                    a, labels[train]
                )
                if set(fitted.classes_) != set(order):
                    raise ValueError("A training fold lacks a bleaching group.")
                probabilities[name][test] = fitted.predict_proba(b)[
                    :, [list(fitted.classes_).index(group) for group in order]
                ]
                predictions[name][test] = fitted.predict(b)
            fold_rows.append(
                {
                    "mode": mode,
                    "fold": fold + 1,
                    "train_samples": len(train),
                    "test_samples": len(test),
                    "shared_sites": 0,
                }
            )
        for name in predictions:
            proba = probabilities[name]
            mask = np.isin(labels, ["High", "Low"])
            rows.append(
                {
                    "mode": mode,
                    "inputs": name,
                    "balanced_accuracy": balanced_accuracy_score(labels, predictions[name]),
                    "macro_AUC": roc_auc_score(
                        pd.Categorical(labels, order).codes,
                        proba,
                        multi_class="ovr",
                        average="macro",
                        labels=[0, 1, 2],
                    ),
                    "AUC_high_vs_low": roc_auc_score(
                        labels[mask] == "High", proba[mask, 2] / (proba[mask, 2] + proba[mask, 0])
                    ),
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(fold_rows)
