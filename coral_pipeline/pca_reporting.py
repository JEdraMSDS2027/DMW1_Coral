"""PCA summary tables, computed narrative, and one organised export step."""

import json
from pathlib import Path
from importlib.metadata import version

import joblib
import numpy as np
import pandas as pd

from .data import CATEGORICAL_COLUMNS, file_sha256
from .pca_statistics import (
    component_associations,
    welch_high_low,
    hotelling_comparison,
    variable_associations,
    omega_squared,
    within_stratum_delta,
)


def analysis_tables(samples, frame, result, n_interpret=5):
    """Assemble the original descriptive and nominal-inference tables without plotting."""
    categories = samples.copy()
    category_names = [name for name in CATEGORICAL_COLUMNS if name in samples]
    for name in category_names:
        categories[name] = categories[name].fillna("Unknown").astype(str)
    groups = samples.Bleaching_Group.to_numpy()
    scores = result["scores"].iloc[:, : result["n_retained"]]
    component_tests = welch_high_low(
        {name: scores[name].to_numpy() for name in scores}, groups, samples.Site_ID.to_numpy()
    )
    variable_tests = welch_high_low(
        {name: frame[name].to_numpy() for name in frame}, groups, samples.Site_ID.to_numpy()
    )
    group_table = component_associations(
        scores, samples, result["variance"], ["Low", "Moderate", "High"]
    )
    source = (
        categories.groupby("Data_Source")
        .agg(
            samples=("Sample_ID", "size"),
            sites=("Site_ID", "nunique"),
            pct_high=("Bleaching_Group", lambda s: 100 * s.eq("High").mean()),
            pct_low=("Bleaching_Group", lambda s: 100 * s.eq("Low").mean()),
            median_bleaching=("Percent_Bleaching", "median"),
            mean_bleaching=("Percent_Bleaching", "mean"),
            pct_depth_imputed=("Depth_m_was_missing", lambda s: 100 * s.mean()),
            first_year=("Date_Year", "min"),
            last_year=("Date_Year", "max"),
        )
        .sort_values("samples", ascending=False)
    )
    for label in ["High", "Low"]:
        source[f"n_{label.lower()}"] = (
            categories.loc[groups == label]
            .groupby("Data_Source")
            .size()
            .reindex(source.index, fill_value=0)
        )
    analyse = [
        name
        for name in [
            "Data_Source",
            "Ocean_Name",
            "Realm_Name",
            "Country_Name",
            "Ecoregion_Name",
            "Exposure",
        ]
        if name in categories
    ]
    omega = pd.DataFrame(
        {
            name: [100 * omega_squared(scores[name], categories[cat]) for cat in analyse]
            + [100 * omega_squared(scores[name], groups)]
            for name in scores
        },
        index=analyse + ["Bleaching_Group (reference)"],
    )
    columns = {name: scores[name].to_numpy() for name in list(scores)[:n_interpret]}
    within, adjusted = within_stratum_delta(categories, groups, "Data_Source", columns)
    pooled = group_table.cliffs_delta.loc[list(columns)].to_dict()
    delta = pd.concat(
        [
            within,
            pd.DataFrame([pooled, adjusted], index=["POOLED (all samples)", "SOURCE-ADJUSTED"]),
        ]
    )
    variable = variable_associations(frame, samples, result["families"])
    return {
        "pca_explained_variance": result["variance"],
        "pca_loadings": result["loadings"],
        "feature_contributions": result["contributions"],
        "family_contributions": result["family_contributions"],
        "pca_group_separation": group_table,
        "welch_tests_components": component_tests,
        "welch_tests_variables": variable_tests,
        "hotelling_high_vs_low": hotelling_comparison(frame, scores, groups),
        "variable_high_vs_low": variable,
        "variable_family_effects": variable.groupby("family").legacy_unweighted_d.agg(
            ["mean", "min", "max", "count"]
        ),
        "data_source_summary": source,
        "categorical_omega2_by_component": omega,
        "within_source_component_delta": delta,
        "within_source_only": within,
    }


def interpret_components(result, n_interpret=5):
    """Produce a readable table of dominant families and variables from current loadings."""
    rows = []
    for name in list(result["scores"])[:n_interpret]:
        loading = result["loadings"][name]
        top = loading.abs().nlargest(5).index
        families = result["family_contributions"][name].nlargest(3)
        rows.append(
            {
                "component": name,
                "variance_percent": result["variance"].loc[name, "variance_%"],
                "leading_families": "; ".join(
                    f"{label}: {value:.1f}%" for label, value in families.items()
                ),
                "strongest_variables": "; ".join(
                    f"{label} ({loading[label]:+.2f})" for label in top
                ),
            }
        )
    return pd.DataFrame(rows).set_index("component")


def results_markdown(samples, result, tables):
    """Generate numerical conclusions from computed tables rather than copied Markdown."""
    group = tables["pca_group_separation"]
    strongest = group.cliffs_delta.abs().nlargest(2).index.tolist()
    variance = result["variance"]
    k = result["n_retained"]
    lines = [
        f"The {len(samples):,} events cover {samples.Date_Year.min()}–{samples.Date_Year.max()}. "
        f"{k} components retain {variance.iloc[k-1]['cumulative_%']:.1f}% of environmental variance; "
        f"PC1–PC2 retain {variance.iloc[1]['cumulative_%']:.1f}%.",
        "The strongest pooled High-versus-Low score differences are:",
    ]
    for name in strongest:
        row = group.loc[name]
        adjusted = tables["within_source_component_delta"].loc["SOURCE-ADJUSTED", name]
        lines.append(
            f"- **{name}**: Cliff's δ = {row.cliffs_delta:+.3f}, AUC = {row.AUC_high_vs_low:.3f}, "
            f"weighted Cohen's d = {row.cohens_d:+.3f}; within-source combined δ = {adjusted:+.3f}."
        )
    ratio = tables["hotelling_high_vs_low"].share_of_all_variables_D2.iloc[1]
    lines.append(
        f"\nThe retained components preserve {100*ratio:.1f}% of the sample-level squared Mahalanobis "
        "separation between High and Low, relative to all 35 variables. This is a descriptive "
        "mean-separation ratio, not predictive accuracy or the fraction of all bleaching information."
    )
    lines.append(
        "Source-stratified and pooled estimates use different eligible populations and weights. "
        "Their difference cannot be interpreted as the causal fraction explained by survey source."
    )
    if "site_cluster_intervals" in tables:
        cluster = tables["site_cluster_intervals"]
        lines.append(
            f"\n{int(cluster.sites_both.iloc[0])} sites occur in both comparison groups. "
            "The original independent-site Welch p-values are retained as legacy results. "
            "The following site-cluster intervals provide pointwise uncertainty conditional on "
            "the fitted preprocessing/PCA, assuming sites are independent of one another."
        )
        for name in strongest:
            row = cluster.loc[name]
            lines.append(
                f"- {name}: difference of site means {row.difference_of_site_means:+.3f}, "
                f"cluster 95% CI [{row.cluster_ci95_low:+.3f}, {row.cluster_ci95_high:+.3f}]."
            )
    return "\n\n".join(lines)


def compare_reference(root, result, tables):
    """Check refactoring against the captured original calculations for this input hash."""
    root = Path(root)
    path = root / "tests" / "pca_reference.json"
    if not path.exists():
        return pd.DataFrame(
            [{"check": "Original reference fixture available", "status": "not checked"}]
        )
    reference = json.loads(path.read_text(encoding="utf-8"))
    if reference["input_sha256"] != file_sha256(root / "processed" / "coral_pca_ready.csv"):
        return pd.DataFrame(
            [
                {
                    "check": "Reference input fingerprint matches",
                    "status": "different dataset; not comparable",
                }
            ]
        )
    comparisons = {
        "Eigenvalues": (reference["eigenvalues"], result["model"].explained_variance_),
        "Oriented loadings": (reference["loadings"], result["loadings"].to_numpy()),
        "Component Cliff delta": (
            reference["component_deltas"],
            tables["pca_group_separation"].cliffs_delta,
        ),
        "Component Cohen d": (reference["component_d"], tables["pca_group_separation"].cohens_d),
        "Hotelling squared distances": (
            reference["hotelling_D2"],
            tables["hotelling_high_vs_low"].mahalanobis_sq,
        ),
        "Within-source PC3/PC4 delta": (
            [reference["within_source_pc3"], reference["within_source_pc4"]],
            tables["within_source_component_delta"].loc["SOURCE-ADJUSTED", ["PC3", "PC4"]],
        ),
    }
    rows = []
    for name, (old, new) in comparisons.items():
        np.testing.assert_allclose(old, new, rtol=1e-8, atol=1e-10)
        rows.append(
            {
                "check": name,
                "status": "pass",
                "max_absolute_difference": float(np.max(np.abs(np.asarray(old) - np.asarray(new)))),
            }
        )
    return pd.DataFrame(rows)


def export_analysis(root, samples, result, tables, figures, settings, input_metadata, output_dir=None):
    """Export the analysis to a chosen directory (default: outputs/pca).

    A separate report variant can supply its own directory without replacing
    another version's tables, figures, fitted model, or manifest.
    """
    root = Path(root)
    output = Path(output_dir) if output_dir is not None else root / "outputs" / "pca"
    output.mkdir(parents=True, exist_ok=True)
    (output / "figures").mkdir(exist_ok=True)
    exported = []
    all_tables = {
        **tables,
        "component_interpretation": interpret_components(result, settings["n_interpret"]),
        "pca_scores": pd.concat(
            [samples, result["scores"].iloc[:, : result["n_retained"]]], axis=1
        ),
    }
    for name, table in all_tables.items():
        path = output / f"{name}.csv"
        keep_index = not isinstance(table.index, pd.RangeIndex) or table.index.name is not None
        index_name = table.index.name or (
            "feature" if name in {"pca_loadings", "feature_contributions"} else "item"
        )
        table.to_csv(path, index=keep_index, index_label=index_name, date_format="%Y-%m-%d")
        exported.append(path)
    for name, figure in figures.items():
        path = output / "figures" / f"{name}.png"
        figure.savefig(path, dpi=140, bbox_inches="tight")
        exported.append(path)
    model_path = output / "pca_model.joblib"
    joblib.dump(result["model"], model_path)
    exported.append(model_path)
    (output / "findings.md").write_text(
        results_markdown(samples, result, tables) + "\n", encoding="utf-8"
    )
    manifest = {
        "settings": settings,
        "samples": len(samples),
        "features": (
            list(result["model"].feature_names_in_)
            if hasattr(result["model"], "feature_names_in_")
            else input_metadata["pca_input_columns"]
        ),
        "input_sha256": file_sha256(root / "processed" / "coral_pca_ready.csv"),
        "preprocessing_metadata_sha256": file_sha256(
            root / "processed" / "coral_preprocessing_metadata.json"
        ),
        "preprocessing_source_sha256": input_metadata["source_sha256"],
        "code_sha256": {
            path.name: file_sha256(path)
            for path in (root / "coral_pipeline").glob("*.py")
        },
        "retention_counts": {
            "kaiser": result["n_kaiser"],
            "parallel": result["n_parallel"],
            "legacy_min": result["n_legacy_min"],
            "selected": result["n_retained"],
        },
        "pca_fit_scope": "Full descriptive cohort; CV uses fold-specific fits separately.",
        "site_inference": "Original sample-level Welch/Hotelling and independent-site Welch are legacy nominal results. Cluster intervals are pointwise and conditional on the fitted axes; they are not multiplicity-adjusted.",
        "congruence_matching": "Original independent best match; not necessarily one-to-one.",
        "versions": {
            name: version(name)
            for name in ["numpy", "pandas", "scipy", "scikit-learn", "matplotlib", "joblib"]
        },
        "method_sources": [
            "https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.bootstrap.html",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC7148287/",
        ],
        "output_sha256": {str(path.relative_to(output)): file_sha256(path) for path in exported},
    }
    (output / "run_metadata.json").write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return pd.DataFrame({"file": [str(path.relative_to(output)) for path in exported]})
