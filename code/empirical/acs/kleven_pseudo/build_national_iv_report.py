# Reusable report/figure builder for the national ACS Twin1/SameSex2 housing
# IV run. Regenerates output/acs_fertility_iv/RESULTS.md and the main
# figure (PNG+PDF) from the SMALL receipts already collected under
# output/acs_fertility_iv/national_128g_results/ (fit_receipts/*.json,
# primary_receipts/*.json, national_18row_table.csv, counts*.json,
# sample_gate_receipt.json, analytic_frames_identity.json). Never reads raw
# source partitions or analytic_frames.rds. Run:
#   python3 code/empirical/acs/kleven_pseudo/build_national_iv_report.py
import csv, json, os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
RES_DIR = os.path.join(ROOT, "output", "acs_fertility_iv")
NAT_DIR = os.path.join(RES_DIR, "national_128g_results")


def load():
    rows = list(csv.DictReader(open(os.path.join(NAT_DIR, "national_18row_table.csv"))))
    counts = json.load(open(os.path.join(NAT_DIR, "counts.json")))
    sample_gate = json.load(open(os.path.join(NAT_DIR, "sample_gate_receipt.json")))
    identity = json.load(open(os.path.join(NAT_DIR, "analytic_frames_identity.json")))
    review = json.load(open(os.path.join(NAT_DIR, "review_summary.json")))
    return rows, counts, sample_gate, identity, review


def f(x):
    return float(x)


def pct(x):
    """probability-units coefficient -> percentage points, 4dp."""
    return f(x) * 100


def enhance_csv(rows):
    """Add pp-converted columns for ownership RF/IV/AR and FS (all designs,
    since FS is always a probability), WITHOUT removing or altering any
    original native-probability-unit column. Written as a *new* file so the
    original national_18row_table.csv units are never touched."""
    out_path = os.path.join(NAT_DIR, "national_18row_table_with_pp.csv")
    fieldnames = list(rows[0].keys()) + [
        "fs_coef_pp", "fs_se_pp",
        "rf_coef_pp_if_ownership", "rf_se_pp_if_ownership", "rf_ci_lower_pp_if_ownership", "rf_ci_upper_pp_if_ownership",
        "iv_coef_pp_if_ownership", "iv_se_pp_if_ownership", "iv_ci_lower_pp_if_ownership", "iv_ci_upper_pp_if_ownership",
        "ar_summary_lower_pp_if_ownership", "ar_summary_upper_pp_if_ownership",
        "ar_honesty_note",
    ]
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            is_own = r["outcome"] == "OWNERSHP_out"
            row = dict(r)
            row["fs_coef_pp"] = round(pct(r["fs_coef"]), 4)
            row["fs_se_pp"] = round(pct(r["fs_se"]), 4)
            for k_out, k_in in [("rf_coef_pp_if_ownership", "rf_coef"), ("rf_se_pp_if_ownership", "rf_se"),
                                 ("rf_ci_lower_pp_if_ownership", "rf_ci_lower"), ("rf_ci_upper_pp_if_ownership", "rf_ci_upper"),
                                 ("iv_coef_pp_if_ownership", "iv_coef"), ("iv_se_pp_if_ownership", "iv_se"),
                                 ("iv_ci_lower_pp_if_ownership", "iv_ci_lower"), ("iv_ci_upper_pp_if_ownership", "iv_ci_upper"),
                                 ("ar_summary_lower_pp_if_ownership", "ar_summary_lower"),
                                 ("ar_summary_upper_pp_if_ownership", "ar_summary_upper")]:
                row[k_out] = round(pct(r[k_in]), 4) if (is_own and r[k_in] not in ("", "NA", "None")) else ""
            n_comp = r.get("ar_n_components", "")
            bounded = r.get("ar_fully_interior_bounded", "")
            all_acc = r.get("ar_all_grid_points_accepted", "")
            if r.get("ar_summary_lower") == r.get("ar_summary_upper"):
                row["ar_honesty_note"] = ("single-grid-point accepted set: only this one grid value was tested and "
                                           "not rejected -- NOT a zero-width continuous CI; interior grid points do "
                                           "not establish exact continuous endpoints")
            elif all_acc == "True":
                row["ar_honesty_note"] = "all grid points accepted -- unbounded/uninformative on this grid, not a tight interval"
            elif bounded == "True":
                row["ar_honesty_note"] = f"interior-bounded across {n_comp} component(s) on the tested grid (grid-truncated approximation, not an exact analytic boundary)"
            else:
                row["ar_honesty_note"] = "boundary-touching or multi-component/errored set -- see ar_extent_unknown_due_to_errors"
            w.writerow(row)
    return out_path


def build_results_md(rows, counts, sample_gate, identity, review):
    by = {(r["design"], r["outcome"]): r for r in rows}
    lines = []
    lines.append("# ACS Twin1 / SameSex2 national housing IV — final results (2026-09-22)\n\n")
    lines.append(f"**Job {review['national_job']}, {review['state'].split(',')[0]}, "
                 f"{review['state'].split(',')[1].strip()}, {review['state'].split(',')[2].strip()}. "
                 f"0 warnings in this national run's log** (the 168 fixest VCOV-not-PSD warnings seen in the "
                 f"prior VT recovery smoke, job 18271160, were a separate, smaller-sample observation and are "
                 f"not claimed resolved, reproduced, or harmless here).\n\n")
    lines.append(f"Code: `{review['code_root'].split('/')[-1]}` (checksummed, read-only). "
                 f"All 18 design×outcome fits reached `full_fit` and pass internal consistency checks "
                 f"(nobs rf==fs==iv, AR error-free, full V symmetric and dimension-matched to named coefficients, "
                 f"V diagonal matches reported SE, primary-receipt coefficients equal final-receipt coefficients). "
                 f"Full per-fit check log: [`national_128g_results/review_summary.json`](national_128g_results/review_summary.json).\n\n")

    lines.append("## Sample construction (reproduces prior construction, not new biology/social-link proof)\n\n")
    lines.append(f"51 state+DC source partitions (2005–2019 ACS 1-year product, `SAMPLE==YEAR*100+1`): "
                 f"{sample_gate['n_input']:,} raw rows read, {sample_gate['n_excluded_year_out_of_range']:,} "
                 f"excluded as out-of-year, {sample_gate['n_excluded_non_acs1yr_product']} excluded as non-1yr "
                 f"product, {sample_gate['n_kept']:,} kept. Mothers age {counts['age_lo']}–{counts['age_hi']} with "
                 f"oldest linked child <18: {counts['n_unique_mothers_sample']:,} unique mothers, "
                 f"{counts['n_households_sample']:,} households. `MOMLOC` links are coresident/social, not "
                 f"confirmed biological histories. These are **contemporaneous cross-sections**, not longitudinal "
                 f"or pre-birth histories.\n\n")
    lines.append(f"- **Twin1**: age-only twin-like proxy (extract27 has no `BIRTHQTR` — not a confirmed-twin test), "
                 f"event age = oldest linked child's age (pooled 0:5), treatment = ≥2 linked children. "
                 f"{counts['Twin1_eligible_N']:,} eligible, {counts['Twin1_proxy_positive_N']:,} instrument-positive.\n")
    lines.append(f"- **SameSex2**: oldest-two linked children same sex, event age = second-oldest child's age "
                 f"(pooled 0:5), treatment = ≥3 linked children — a **third-child margin, a different population "
                 f"from Twin1's second-child margin**. {counts['SameSex2_primary_eligible_N']:,} primary-eligible "
                 f"(of {counts['SameSex2_eligible_pool_N']:,} pool; excludes {counts['SameSex2_primary_age_tie_excluded_N']:,} "
                 f"primary-age ties), {counts['SameSex2_positive_samesex_N']:,} same-sex positive "
                 f"({counts['SameSex2_both_boys_N']:,} both-boys, {counts['SameSex2_both_girls_N']:,} both-girls).\n\n")

    lines.append("## Methods\n\n")
    lines.append("Mother rows are unique by `(YEAR, SAMPLE, SERIAL, PERNUM)`, weighted by mother `PERWT`. "
                 "Standard errors are clustered by household, `(YEAR, SAMPLE, SERIAL)`. Control formula: "
                 "`i(mat_age_at_event) + i(RACE) + i(survey_year) + i(event_age)` — indicator (dummy) sets for "
                 "mother's inferred age at the event (`mat_age_at_event = mother's ACS-interview age − event age`), "
                 "race, ACS survey year, and event age itself (0:5). Confidence intervals are the normal-quantile "
                 "approximation (±1.96×clustered SE) from the reported clustered covariance. The reported F is a "
                 "single-instrument cluster-robust Wald statistic, `F = (FS coefficient / clustered SE)^2` — this "
                 "is **not** a Kleibergen–Paap statistic. Outcomes: `ROOMS` valid codes 1:27,30, capped at 9 for "
                 "the primary outcome; `OWNERSHP` 1=owner/2=renter, recoded to 0/1; `BEDROOMS` (diagnostic only) "
                 "valid codes 1:22, recoded to (code−1) capped at 5. Sample effective size (ESS) was **not saved** "
                 "by the collected receipts and is not reported below — not fabricated, not recomputed from a "
                 "heavy read.\n\n")

    def row_line(name, r, is_iv=False):
        oc_is_own = r["outcome"] == "OWNERSHP_out"
        unit = "pp owner" if oc_is_own else "rooms"
        mult = 100 if oc_is_own else 1
        coef_key = "iv_coef" if is_iv else "rf_coef"
        se_key = "iv_se" if is_iv else "rf_se"
        lo_key = "iv_ci_lower" if is_iv else "rf_ci_lower"
        hi_key = "iv_ci_upper" if is_iv else "rf_ci_upper"
        return (f"| {name} | {f(r[coef_key]) * mult:.4f} | {f(r[se_key]) * mult:.4f} | "
                f"[{f(r[lo_key]) * mult:.4f}, {f(r[hi_key]) * mult:.4f}] | {unit} |\n")

    lines.append("## Table 1 — Reduced form + first stage (PRIMARY objects, pooled event age 0:5)\n\n")
    lines.append("RF = effect of the instrument (Twin1 age-only proxy / SameSex2 oldest-two-same-sex) on the "
                 "housing outcome; interpret as **instrument–outcome association**, not an established causal "
                 "additional-child effect. FS = effect of the instrument on the additional-child treatment "
                 "probability, in **percentage points**.\n\n")
    lines.append("| Design | Outcome | N | HH | Z-pos | RF coef | RF SE | RF 95% CI | Unit | FS coef (pp) | FS SE (pp) | First-stage F |\n")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|\n")
    for d, dname in [("Twin1_pooled0_5", "Twin1"), ("SameSex2_pooled0_5", "SameSex2")]:
        for oc, ocname in [("ROOMS_out", "Rooms"), ("OWNERSHP_out", "Ownership")]:
            r = by[(d, oc)]
            is_own = oc == "OWNERSHP_out"
            mult = 100 if is_own else 1
            unit = "pp owner" if is_own else "rooms"
            lines.append(f"| {dname} | {ocname} | {r['n_usable_prefit']} | {r['n_households_prefit']} | "
                         f"{r['n_instrument_positive']} | {f(r['rf_coef']) * mult:.4f} | {f(r['rf_se']) * mult:.4f} | "
                         f"[{f(r['rf_ci_lower']) * mult:.4f}, {f(r['rf_ci_upper']) * mult:.4f}] | {unit} | "
                         f"{pct(r['fs_coef']):.4f} | {pct(r['fs_se']):.4f} | {f(r['first_stage_F']):.1f} |\n")

    lines.append("\n## Table 2 — 2SLS and Anderson–Rubin (ASSUMPTION-DEPENDENT diagnostics, not causally certified)\n\n")
    lines.append("2SLS is the assumed instrumented **effect of the additional-child treatment** on the outcome. "
                 "Reported only because the first stage is non-degenerate; the exclusion restriction is "
                 "**unresolved** for both designs (see caveats below) so this is a diagnostic, not a causal estimate.\n\n")
    lines.append("| Design | Outcome | 2SLS coef | 2SLS SE | 2SLS 95% CI | Unit | AR accepted set (tested grid) | AR honesty flag |\n")
    lines.append("|---|---|---|---|---|---|---|---|\n")
    for d, dname in [("Twin1_pooled0_5", "Twin1"), ("SameSex2_pooled0_5", "SameSex2")]:
        for oc, ocname in [("ROOMS_out", "Rooms"), ("OWNERSHP_out", "Ownership")]:
            r = by[(d, oc)]
            is_own = oc == "OWNERSHP_out"
            mult = 100 if is_own else 1
            unit = "pp owner" if is_own else "rooms"
            ar_lo, ar_hi = f(r["ar_summary_lower"]) * mult, f(r["ar_summary_upper"]) * mult
            if ar_lo == ar_hi:
                flag = f"SINGLE TESTED GRID POINT ({ar_lo:.4f}) accepted -- not a zero-width CI, interior points untested between grid steps"
            elif r["ar_all_grid_points_accepted"] == "True":
                flag = "all grid points accepted -- unbounded on this grid"
            elif r["ar_fully_interior_bounded"] == "True":
                flag = f"interior-bounded, {r['ar_n_components']} component(s), grid-truncated approximation"
            else:
                flag = "boundary-touching / not fully interior"
            lines.append(f"| {dname} | {ocname} | {f(r['iv_coef']) * mult:.4f} | {f(r['iv_se']) * mult:.4f} | "
                         f"[{f(r['iv_ci_lower']) * mult:.4f}, {f(r['iv_ci_upper']) * mult:.4f}] | {unit} | "
                         f"[{ar_lo:.4f}, {ar_hi:.4f}] | {flag} |\n")

    lines.append("\n## Descriptive interpretation\n\n")
    lines.append("Twin1 and SameSex2 show **opposite-signed** room-count RF (Twin1 positive, SameSex2 negative). "
                 "This is a descriptive contrast across **different samples and margins** (second-child vs. "
                 "third-child; different instrument construction and different populations), not evidence that "
                 "an additional child causally shrinks the home — the SameSex2 result is also consistent with a "
                 "direct room-sharing/allocation response to same-sex composition unrelated to any third-child "
                 "effect (the exclusion restriction concern above). Neither RF is a like-for-like counterfactual "
                 "of the other.\n\n")

    lines.append("## Identification caveats\n\n")
    lines.append("- **RF and FS are the primary, most defensible objects.** 2SLS is an assumption-dependent "
                 "diagnostic, not a causally certified estimate.\n")
    lines.append("- **Exclusion restriction is unresolved for both designs.** Same-sex composition of the oldest "
                 "two children may directly change room-sharing/housing demand independent of any third-child "
                 "effect. Twin-like status associates with maternal health conditions and birth spacing that "
                 "independently affect housing. Neither is addressed by this run.\n")
    lines.append("- **Twin1 is an age-only proxy**, not a confirmed-twin indicator — extract27 has no `BIRTHQTR` "
                 "to bridge to child birth quarter.\n")
    lines.append("- **Twin1 and SameSex2 identify different populations/margins** (second-child vs. third-child) "
                 "and must not be pooled or compared as a single estimate.\n")
    lines.append("- **AR intervals are finite-grid approximations.** Some cells' accepted set is a single tested "
                 "grid point (see the `ar_honesty_note` column in the enhanced CSV) — this is not a zero-width "
                 "continuous confidence interval; untested points between grid steps are not certified rejected "
                 "or accepted.\n")
    lines.append("- No causal certification and no calibration adoption is implied by this run.\n\n")

    lines.append("## Full results table (all 18: 2 designs × 3 outcomes × 3 windows, no cherry-picking)\n\n")
    lines.append("[`national_128g_results/national_18row_table.csv`](national_128g_results/national_18row_table.csv) "
                 "-- original file, **native probability units**, unmodified. "
                 "[`national_128g_results/national_18row_table_with_pp.csv`](national_128g_results/national_18row_table_with_pp.csv) "
                 "-- same data plus added percentage-point and AR-honesty columns (original columns untouched).\n\n")
    lines.append("Brief event-age-3/5 read (from the full CSV, no selection by significance): Twin1 rooms RF is "
                 "positive at both event age 3 and 5, consistent in sign with the pooled estimate; SameSex2 rooms "
                 "RF is negative at both, also consistent with the pooled estimate. Magnitudes and precision vary "
                 "across the smaller event-age-specific subsamples; see the CSV for exact values rather than a "
                 "restated table here.\n\n")
    lines.append("Figure: [`national_128g_results/national_main_results.png`](national_128g_results/national_main_results.png) "
                 "/ [`.pdf`](national_128g_results/national_main_results.pdf)\n\n")

    lines.append("## Prior history (preserved, not overwritten)\n\n")
    lines.append("- VT recovery-smoke review (job 18271160, computational gate only): "
                 "[`recovery_smoke_review/`](recovery_smoke_review/)\n")
    lines.append("- Original national OOM failure and diagnosis (job 18247804, superseded by this run): preserved "
                 "in git history.\n")

    open(os.path.join(RES_DIR, "RESULTS.md"), "w").write("".join(lines))


def build_figure(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    by = {(r["design"], r["outcome"]): r for r in rows}
    designs = [("Twin1_pooled0_5", "Twin1\n(age-only twin-like proxy)"),
               ("SameSex2_pooled0_5", "SameSex2\n(3rd-child margin)")]

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6))
    for ax, (oc, title, ylabel, mult) in zip(axes, [
        ("ROOMS_out", "Rooms (capped at 9)", "Effect on rooms", 1),
        ("OWNERSHP_out", "Ownership", "Effect on P(owner), pp", 100),
    ]):
        xs = [0, 1]
        labels = [d[1] for d in designs]
        rf_y = [f(by[(d[0], oc)]["rf_coef"]) * mult for d in designs]
        rf_lo = [f(by[(d[0], oc)]["rf_ci_lower"]) * mult for d in designs]
        rf_hi = [f(by[(d[0], oc)]["rf_ci_upper"]) * mult for d in designs]
        ax.errorbar([x - 0.08 for x in xs], rf_y,
                    yerr=[[y - l for y, l in zip(rf_y, rf_lo)], [h - y for y, h in zip(rf_y, rf_hi)]],
                    fmt="o", color="#1f77b4", capsize=4, label="RF: instrument–outcome association (95% CI)")
        iv_y = [f(by[(d[0], oc)]["iv_coef"]) * mult for d in designs]
        iv_lo = [f(by[(d[0], oc)]["iv_ci_lower"]) * mult for d in designs]
        iv_hi = [f(by[(d[0], oc)]["iv_ci_upper"]) * mult for d in designs]
        ax.errorbar([x + 0.08 for x in xs], iv_y,
                    yerr=[[y - l for y, l in zip(iv_y, iv_lo)], [h - y for y, h in zip(iv_y, iv_hi)]],
                    fmt="s", color="#d62728", capsize=4, label="2SLS: assumed fertility effect (95% CI, assumption-dependent)")
        ax.axhline(0, color="gray", linewidth=0.8, linestyle="--")
        ax.set_xticks(xs)
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_title(title, fontsize=10)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_xlim(-0.6, 1.6)

    handles, labels_ = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels_, loc="lower center", bbox_to_anchor=(0.5, -0.02), ncol=2, fontsize=8, frameon=False)
    fig.suptitle("National ACS Twin1/SameSex2 housing IV — pooled event age 0:5 (job 18274624)", fontsize=10)
    fig.text(0.5, 0.035,
              "Age-only Twin1 proxy (not confirmed twins); Twin1 and SameSex2 are different populations/margins. "
              "RF/FS primary; 2SLS assumption-dependent, exclusion unresolved.",
              ha="center", fontsize=7, style="italic")
    fig.tight_layout(rect=[0.02, 0.14, 0.98, 0.94])
    fig.savefig(os.path.join(NAT_DIR, "national_main_results.png"), dpi=150, bbox_inches="tight")
    fig.savefig(os.path.join(NAT_DIR, "national_main_results.pdf"), bbox_inches="tight")


if __name__ == "__main__":
    rows, counts, sample_gate, identity, review = load()
    enhance_csv(rows)
    build_results_md(rows, counts, sample_gate, identity, review)
    build_figure(rows)
    print("Rebuilt RESULTS.md, national_18row_table_with_pp.csv, national_main_results.{png,pdf}")
