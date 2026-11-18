"""PCA extraction, retention, loadings, source fits, and robustness checks."""

import re
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

DEFAULT_COMPONENTS = 5
DEFAULT_SEED = 42


def fit_loadings(frame, n_components=DEFAULT_COMPONENTS, seed=DEFAULT_SEED):
    """Fit a PCA on ``frame`` and return its loadings and variance ratios.

    The loadings are the eigenvectors (variables x components); the second
    value is the share of variance explained by each component.
    """
    model = PCA(n_components=n_components, random_state=seed)
    model.fit(frame.to_numpy())
    names = [f"PC{i + 1}" for i in range(n_components)]
    loadings = pd.DataFrame(model.components_.T, index=frame.columns, columns=names)
    return loadings, model.explained_variance_ratio_


def loading_congruence(reference, other, n_components=DEFAULT_COMPONENTS):
    """Match each reference component to its best component in ``other``.

    Congruence is the absolute cosine between loading vectors, so it ignores
    the arbitrary sign of a component. Returns the best congruence per
    reference component and the name of the matching component.
    """
    common = reference.index.intersection(other.index)
    ref = reference.loc[common].to_numpy()[:, :n_components].copy()
    oth = other.loc[common].to_numpy().copy()
    ref /= np.linalg.norm(ref, axis=0)
    oth /= np.linalg.norm(oth, axis=0)
    cosines = np.abs(ref.T @ oth)
    matches = [f"PC{i + 1}" for i in cosines.argmax(axis=1)]
    return cosines.max(axis=1), matches


FAMILY_RULES = [  # first match wins
    ("Turbidity", r"^Turbidity$"),
    ("Depth", r"^Depth_m$"),
    ("Cyclone frequency", r"^Cyclone_Frequency$"),
    ("Distance to shore", r"^Distance_to_Shore$"),
    ("Wind speed", r"^Windspeed$"),
    ("Degree heating weeks (DHW)", r"DHW"),
    ("Anomaly frequency", r"Frequency"),
    ("TSA (thermal stress anomaly)", r"^TSA"),
    ("SSTA (SST anomaly)", r"^SSTA"),
    ("Temperature variability", r"^Temperature.*Standard_Deviation"),
    ("Sea surface temperature (SST)", r"^(ClimSST|Temperature)"),
]

FAMILY_ORDER = [
    "Sea surface temperature (SST)",
    "Temperature variability",
    "SSTA (SST anomaly)",
    "TSA (thermal stress anomaly)",
    "Anomaly frequency",
    "Degree heating weeks (DHW)",
    "Turbidity",
    "Depth",
    "Distance to shore",
    "Cyclone frequency",
    "Wind speed",
    "Other",
]


def family(name):
    return next((f for f, pat in FAMILY_RULES if re.search(pat, str(name))), "Other")


def fit_pca(
    frame: pd.DataFrame, permutations: int = 50, seed: int = 42, retention_rule: str = "parallel"
) -> dict:
    """Fit the existing full PCA with consistent signs and explicit retention.

    `parallel` follows the original narrative; `legacy_min` reproduces the old
    min(Kaiser, parallel) rule. Both counts are always returned for comparison.
    Loadings retain the original eigenvector * sqrt(sample eigenvalue) definition;
    with population-standardised inputs they approximate correlations within n/(n-1).
    """
    values = frame.to_numpy(dtype=float)
    if values.shape[0] <= values.shape[1] or not np.isfinite(values).all():
        raise ValueError("Use a finite feature matrix with more rows than features.")
    if permutations < 1 or retention_rule not in {"parallel", "legacy_min"}:
        raise ValueError("Choose positive permutations and parallel or legacy_min retention.")
    model = PCA(random_state=seed).fit(values)
    signs = np.sign(
        model.components_[np.arange(len(frame.columns)), np.abs(model.components_).argmax(axis=1)]
    )
    signs[signs == 0] = 1
    # Orient the model itself so reloaded model.transform reproduces exported scores.
    model.components_ *= signs[:, None]
    names = [f"PC{i+1}" for i in range(len(frame.columns))]
    scores = pd.DataFrame(model.transform(values), index=frame.index, columns=names)
    eigenvalues = model.explained_variance_
    rng = np.random.default_rng(seed)
    null = np.empty((permutations, len(names)))
    for iteration in range(permutations):
        shuffled = np.column_stack([rng.permutation(values[:, j]) for j in range(len(names))])
        null[iteration] = np.linalg.eigvalsh(np.corrcoef(shuffled, rowvar=False))[::-1]
    parallel95 = np.percentile(null, 95, axis=0)
    k_kaiser = int((eigenvalues > 1).sum())
    k_parallel = int((eigenvalues > parallel95).sum())
    retained = k_parallel if retention_rule == "parallel" else min(k_kaiser, k_parallel)
    if retained < 1:
        raise ValueError("The selected retention criterion retained no components.")
    families = pd.Series({name: family(name) for name in frame}, name="family")
    contributions = pd.DataFrame(model.components_.T**2 * 100, index=frame.columns, columns=names)
    order = [name for name in FAMILY_ORDER if name in set(families)]
    return {
        "model": model,
        "scores": scores,
        "loadings": pd.DataFrame(
            model.components_.T * np.sqrt(eigenvalues), index=frame.columns, columns=names
        ),
        "contributions": contributions,
        "families": families,
        "family_order": order,
        "feature_order": sorted(
            frame.columns, key=lambda name: (FAMILY_ORDER.index(families[name]), name)
        ),
        "family_contributions": contributions.groupby(families).sum().reindex(order),
        "variance": pd.DataFrame(
            {
                "eigenvalue": eigenvalues,
                "variance_%": 100 * model.explained_variance_ratio_,
                "cumulative_%": 100 * np.cumsum(model.explained_variance_ratio_),
                "parallel_95th": parallel95,
            },
            index=names,
        ),
        "n_retained": retained,
        "n_kaiser": k_kaiser,
        "n_parallel": k_parallel,
        "n_legacy_min": min(k_kaiser, k_parallel),
        "retention_rule": retention_rule,
    }


def compare_sources(frame, samples, sources, n_components=5, seed=42, min_samples=800):
    """Refit original globally transformed inputs per source; compare eigenvector directions.

    Best matches are independent, as in the original; multiple references may
    match the same component. This is a descriptive stability diagnostic.
    """
    reference, _ = fit_loadings(frame, n_components, seed)
    rows = []
    for source in sources:
        mask = samples.Data_Source.eq(source).to_numpy()
        if mask.sum() < min_samples:
            continue
        other, _ = fit_loadings(frame.loc[mask], n_components, seed)
        values, matches = loading_congruence(reference, other, n_components)
        for i, (value, match) in enumerate(zip(values, matches)):
            rows.append(
                {
                    "source": source,
                    "samples": int(mask.sum()),
                    "component": f"PC{i+1}",
                    "congruence": value,
                    "matching_component": match,
                }
            )
    return pd.DataFrame(rows)


def robustness_checks(frame, samples, result, n_components=5, seed=42, extreme=10):
    """Preserve random site-half, correlation pruning, and extreme-score refits."""
    from .pca_statistics import auc_high_low

    reference, reference_variance = fit_loadings(frame, n_components, seed)
    rng = np.random.default_rng(seed)
    sites = samples.Site_ID.unique()
    rows, stability = [], []
    for repeat in range(6):
        half = rng.choice(sites, len(sites) // 2, replace=False)
        mask = samples.Site_ID.isin(half).to_numpy()
        for part in (mask, ~mask):
            alternative, _ = fit_loadings(frame.loc[part], n_components, seed)
            values, _ = loading_congruence(reference, alternative, n_components)
            stability.append(values)
    stability = np.array(stability)
    for i in range(n_components):
        rows.append(
            {
                "check": "random site halves (12 refits)",
                "component": f"PC{i+1}",
                "congruence": stability[:, i].mean(),
                "minimum_congruence": stability[:, i].min(),
            }
        )
    correlation = frame.corr().abs()
    seen, keep = set(), []
    for name in result["feature_order"]:
        if name not in seen:
            members = [
                other
                for other in result["feature_order"]
                if other not in seen and correlation.loc[name, other] >= 0.90
            ]
            seen.update(members)
            keep.append(name)
    alternative, pruned_variance = fit_loadings(frame[keep], n_components, seed)
    values, matches = loading_congruence(reference, alternative, n_components)
    groups = samples.Bleaching_Group.to_numpy()
    pruned_scores = PCA(n_components=n_components, random_state=seed).fit_transform(
        frame[keep].to_numpy()
    )
    for i in range(n_components):
        rows.append(
            {
                "check": f"correlation pruning ({len(keep)} inputs)",
                "component": f"PC{i+1}",
                "congruence": values[i],
                "matching_component": matches[i],
                "absolute_delta": abs(2 * auc_high_low(pruned_scores[:, i], groups) - 1),
            }
        )
    extreme_mask = result["scores"].iloc[:, : result["n_retained"]].abs().gt(extreme).any(axis=1)
    if extreme_mask.any():
        alternative, _ = fit_loadings(frame.loc[~extreme_mask], n_components, seed)
        values, matches = loading_congruence(reference, alternative, n_components)
        for i in range(n_components):
            rows.append(
                {
                    "check": f"exclude |score| > {extreme}",
                    "component": f"PC{i+1}",
                    "congruence": values[i],
                    "matching_component": matches[i],
                }
            )
    summary = {
        "pruned_columns": keep,
        "pruned_variance_percent": float(100 * pruned_variance.sum()),
        "full_variance_percent": float(100 * reference_variance.sum()),
        "extreme_samples": int(extreme_mask.sum()),
        "extreme_by_group": samples.loc[extreme_mask].Bleaching_Group.value_counts().to_dict(),
    }
    return pd.DataFrame(rows), summary
