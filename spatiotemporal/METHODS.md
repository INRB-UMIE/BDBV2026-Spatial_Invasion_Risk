# BDBV 2026 — Spatiotemporal Invasion Forecasting: Methods

*An intuitive but complete description of every method in the pipeline, written so
a first-time reader can follow it end to end. Every model is described with the
**same seven headings** — **Description · Equation · Key details · Assumptions ·
Strengths · Weaknesses · Possible extensions** — and every equation is tied to the
**function that implements it** (file → function), so the maths and the code never
drift apart.*

**Reading guide.** §0 sets up the task and the date-window convention used on every
figure. §1–§4 are the ingredients (data, nowcasting, epidemiological parameters,
mobility). §5 is the primary renewal invasion model — its formulation *and* its
code. §6 is the Bayesian re-implementation of that same model. §7 is the comparator
models. §8–§9 are risk products and evaluation. §10–§12 are trust, outputs, and a
one-line model catalogue.

---

## 0. The forecasting task (what everything is built around)

**Question.** For each of the 519 DRC health zones that has had **no confirmed
BDBV case up to the forecast date**, what is the probability that a person there
experiences their **first symptom onset** in the next **1 week** and within the
next **2 weeks**? This is *spatial invasion* — the arrival of the virus in a
previously-uninfected zone.

**Run mode (important).** The pipeline defaults to a **Bayesian-only run**
(`RUN_FREQUENTIST_MODELS = FALSE` in `run_all.R`): the frequentist penalised-MLE renewal suite and
all comparator/baseline models are **off**, and the Bayesian suite (§6) supplies every featured,
operational product. Where the pipeline still prints a "best frequentist" slot it is filled by the
best-cross-validated **Bayesian** model as a fallback, so both featured slots are Bayesian in the
default run. Sections describing the frequentist suite (§5), the comparators (§7), and the
freq-vs-Bayesian comparison (§6.2) therefore apply only when `RUN_FREQUENTIST_MODELS = TRUE`.

**The target is ONSET-dated, and ascertainment-agnostic.** Two deliberate choices:
1. *Onset, not report.* The event is dated by **symptom onset**, the epidemiologically
   correct clock for transmission — not the (later, delay-distorted) sample or
   report date. The ~15 % of records lacking an onset date are imputed as
   `onset = sample_date − Δ`, where Δ is **drawn per record** from the onset→sample delay
   distribution — the EMPIRICAL delay (bootstrap of the observed complete pairs, **windowed**
   to onset ∈ `[outbreak floor, max(sample) − 5 d]` so the right-truncated recent tail and
   pre-outbreak typos are excluded) when ≥30 exist, else the parametric fallback: the rigorous
   **DHIS2-specific interval-censored MLE** fit (`04c_dhis2_delay_windows.R`; §1) rather than
   the fixed lab rate, and never the delay's mean
   (`IMPUTE_ONSET_FROM_SAMPLE`, seeded, `01_data_prep.R`). Drawing spreads the imputed onsets
   by the delay's true shape instead of piling them on one week. (A single stochastic draw does
   not fully propagate imputation uncertainty; full multiple imputation is the extension.)
   Because recent onsets are not yet fully
   reported, the target is *right-truncated at forecast time* and is nowcast (§2);
   evaluation scores against the eventually-observed first-onset week.
2. *Confirmed-case caveat.* We can only reliably observe **confirmed** cases, so the
   operational target is "first **confirmed-case** onset". The primary probability
   `p_case` is therefore **ascertainment-agnostic** — it is calibrated directly to
   observed first-confirmed-onset frequencies and does *not* divide by any assumed
   reporting rate ρ. The ascertainment-adjusted *infection*-scale score is reported
   separately, and as a **band** across ρ rather than a fixed value (§3.3, §8).

**Why this framing matters.** A zone that is already affected cannot be "newly
invaded", so it carries **no invasion probability** — it is masked (`NA`, not `0`)
everywhere: forecasts, maps, rankings, and evaluation denominators. The
decision-relevant quantity is *where the outbreak jumps next*, which is dominated by
human mobility out of the affected cluster.

**The date-window convention (read every figure this way).** Every forecast is
defined by two windows, and both are now printed in the caption of every
current-forecast figure (via `.window_caption()` in `20_forecast_detail.R`):
- **Fit window** — the span of training data the model saw: from the first
  outbreak week (`min(all_weeks)`) to the **cutoff** (the forecast date,
  `training_cutoff` for the live forecast, or the fold cutoff in cross-validation).
- **Prediction window** — the weeks being forecast: `(cutoff, cutoff + 7·h]` days
  for horizon *h* ∈ {1, 2} weeks. A 1-week forecast issued on a Saturday cutoff
  predicts the following Sun–Sat week.
In cross-validation the cutoff *rolls*: each fold re-fits on `week ≤ cutoff` and
predicts forward, so the fit window grows left-to-right across folds and the
prediction window always sits strictly to its right (this is what makes the
evaluation leakage-free, §9).

**Analysis-date-anchored daily re-issue (`22_daily_reissue.R`).** The operational
product is the probability that each at-risk zone records its first case **within
7 days** and **within 14 days of the analysis date** — rolling cumulative windows
`(d, d+7]` and `(d, d+14]` measured from `forecast_date = ANALYSIS_DATE`, not from
the ISO-week boundary. The renewal engine still runs on the weekly ISO grid
internally; `anchor_windows_from_analysis_date()` converts its weekly expected
introductions to a **uniform daily hazard** (weekly μ / 7) and integrates over the
day-windows. With offset `a = forecast_date − training_cutoff` (0–6, how far into
the current ISO week the analysis date sits), the 7-day window covers the `(6−a)`
remaining days of the current week plus the first `(a+1)` days of the next, using
the workhorse's current-week introduction rate `mu_wk0` (comparators fall back to
the next week's rate). So `p = 1 − exp(−μ_window)`. Re-running daily refits
epinowcast against the new `ANALYSIS_DATE` (§2), so the partial current week's
counts — and every probability — sharpen day by day; each issue **accumulates**
into `daily_invasion_series.rds` (one row per issue-date × method × zone × window)
plus an immutable dated snapshot. Approximation to state honestly: the
uniform-within-week hazard spreads introductions evenly across each ISO week and
scales the partial current week by its remaining fraction.

**The hard truth about sample size.** Across the leave-future-out folds there are
only on the order of **~8–18 realised invasion events from a handful of distinct
zones**, against a per-zone-week base rate under ~1 % at horizon 1 and ~1.5 % at
horizon 2 (data-dependent; the current run's exact base rate and event count are in
`invasion_report.md`), computed per horizon as the invasion prevalence among at-risk
zone-weeks). Every ranking below is therefore **provisional** — differences between
good models are frequently within Monte-Carlo noise, and the pipeline prints the
event count on every evaluation figure.

---

## 1. Data & preprocessing  (`01_data_prep.R`)

**Intuition.** Turn a messy line list of individual cases into a clean matrix of
*confirmed cases per health zone per epidemiological week*, indexed by the date the
illness most plausibly began.

**Key details.**
- Source: the DHIS2 line list (`data/processed/dhis2_linelist_processed/`), read by
  `load_and_prepare_data()` and aggregated by `aggregate_to_zone_week()`.
- Time index = **date of symptom onset** where present, else an **onset imputed from the
  sample date** (`sample_date − Δ`, Δ drawn per record from the onset→sample delay; §0), not
  the raw sample date (onset is the epidemiologically correct clock for transmission;
  ~15 % of records lack it).
- Aggregated to **7-day epidemiological weeks anchored so the current (final) week ends
  exactly on the analysis date** (`WEEK_ANCHOR`; weeks stay a uniform 7 days — only the
  anchor day shifts, so the renewal/R(t)/GT machinery is untouched — and the training
  window closes on the as-of date rather than the preceding ISO Monday, §0). Zero-filled
  across all 519 zones so the space–time grid is always complete.
- **DHIS2-specific onset→sample delay estimation** (`04c_dhis2_delay_windows.R`, the DHIS2
  analog of the lab-linelist `04c_delay_four_windows.R`). The DHIS2 line list reports more
  slowly than the Ituri lab linelist (onset→sample mean ≈ 5.9 d / rate ≈ 0.17 d⁻¹ vs the lab
  4.39 d / 0.228 d⁻¹), so its delay is fit **on the DHIS2 data itself**, windowed to onset ∈
  `[outbreak floor, max(sample) − 5 d]` (dropping the right-truncated tail), by **naive +
  interval-censored MLE** (`fitdistcens`; d=0→[0,0.5], d>0→[d−0.5,d+0.5]; AIC-selected across
  gamma/log-normal/Weibull/exponential — gamma wins) with an optional Bayesian **EpiDist**
  marginal model (`RUN_EPIDIST=TRUE`) that additionally corrects right-truncation. This
  windowed DHIS2 fit drives the onset imputation (§0). *Scope note:* the **nowcast** (§2) and
  the R(t) right-truncation model (§3.2) still use the fixed lab rate `DELAY_ONSET_SAMPLE_RATE`
  (0.228); wiring the DHIS2 rate into those is the documented next step.
- Zone names are harmonised to a canonical 519-zone spine (WorldPop) via an alias
  table; province membership — used to **label every figure by health area**
  (Ituri, Nord-Kivu, Haut-Uele, …) and for Ituri-relative risk — comes from the
  health-zone shapefile's `PROVINCE` field.
- **Line-list-derived weekly covariates** (per zone-week): `confirmed`, `suspected`,
  `total_alerts` (row count = alert volume), `tests_analyzed`, `positivity`
  (= confirmed / tests).
- **Static ancillary covariates** (per zone; WorldPop/CCVI/GRID3/GDP layers):
  `pop_count` / `log_pop`, `ccvi` (socioeconomic deprivation), `gdp_pc`,
  `healthsite_count`, `healthsite_density`, plus OSRM zone-to-zone travel times.
- **Contact-tracing line list** is loaded and used **only** to estimate the serial
  interval (an input to the generation time, §3) — never as a forecast covariate.

**Sitrep reconciliation of the confirmed-case set (`APPEND_SITREP_CONFIRMED`, default on).**
The DHIS2 line list lags the official INSP situation reports for some zones: a zone the sitrep
has already **confirmed** can still appear in the line list as suspect-only, so — because
invasion and every case count in this suite are defined purely from line-list *confirmed* cases
— that zone is silently scored as **never invaded**. In this baseline folder that is especially
consequential: the **retrospective leave-future-out CV that selects the featured model (§9)**
would score every method against a ground truth that omits those real invasions. Before modelling,
`load_linelist()` therefore reconciles the line list up to the sitrep's official cumulative
confirmed counts (`.build_sitrep_confirmed_appends()` in `01_data_prep.R`). Precisely:
  - **The sitrep is a floor, never a ceiling.** For each canonical zone we append the
    **shortfall** = `max(0, sitrep_cumulative_confirmed − linelist_confirmed)` as extra confirmed
    line-list rows, so every zone reaches *at least* its officially confirmed count. Zones where
    the line list already meets or exceeds the sitrep get nothing; **no case is ever removed**.
  - **Source = the CUMULATIVE confirmed file** (`insp_sitrep__cumulative_confirmed_cases__daily.csv`),
    which is authoritative for magnitude; the daily `new_confirmed_cases` stream grossly
    undercounts (e.g. Bunia 66 vs 507) and is **not** used here. Per canonical zone we take the
    **max across spelling variants per date**, the **running maximum over dates** (absorbing the
    sitrep's occasional cumulative dips/revisions), then **diff into dated positive integer
    increments** — one confirmed case per unit increment, dated at the sitrep report date on which
    it appeared. Sitrep `nom`s are canonicalised through the **same alias table** as the line list.
  - **Placement & dating.** The shortfall is filled with the **most recent** sitrep increment dates
    (reporting lag ⇒ the line list already holds the earlier cases). Each appended row carries
    `date_of_sample_collection` = that sitrep date, a **blank onset** (imputed downstream exactly
    like any real onset-less record, §0), `final_mve_case_classification = confirmed_case`, and a
    traceable `alert_id = SITREP-CONF-<zone>-NN`.
  - **Spine & as-of guards.** A sitrep zone not on the canonical WorldPop 519-zone spine is
    **skipped with a warning** (it could not enter the zone-week grid anyway); increments dated
    after `ANALYSIS_DATE` are dropped, so a back-dated / leave-future-out re-run never leaks future
    confirmations — the appended cases enter each CV fold only from the week they were actually
    reported.
  - **Order & determinism.** The append is injected **after** the outbreak filter, unmatched-zone
    exclusion, and aliasing (so appended zones are canonical, in-spine, in-outbreak) and **before**
    onset imputation; the rows are bound at the **end** of the line list, which preserves the
    positional RNG draw sequence of the existing records (bit-reproducible vs the pre-reconciliation
    run for every non-appended row).
  - **Scope & relation to `build_conditional_linelist.R`.** This is a correction to the *official
    record* and is applied in **both** this baseline pipeline and the conditional
    (`spatiotemporal_conditional/`) pipeline — distinct from `build_conditional_linelist.R`, whose
    hardcoded appends are *speculative, prospective* what-ifs for the conditional arm only. The two
    are **idempotent**: those pre-appends raise the line-list confirmed count, so the sitrep
    shortfall for those zones is 0 and they are not re-added.
  - **Effect (current snapshot; recomputed and logged every run).** Against the current line list the
    reconciliation flips **3 zones from uninvaded → invaded** — **Oicha** (Nord-Kivu), **Isiro** and
    **Pawa** (Haut-Uele), each sitrep-confirmed but line-list-zero — and tops up confirmed counts on
    ~15 already-invaded zones, for ~90 appended cases across ~18 zones. Only the 3 flips change
    invasion *status* (and therefore the set of scoreable CV events); the top-ups change case
    *magnitude* (which the count-driven import force `Λ = Σⱼ W[j,i]·Σ_k g(k)·Y_nc[j,t−k]`, §5,
    consumes) but not which zones are invaded. Set `APPEND_SITREP_CONFIRMED = FALSE` to reproduce a
    model selection made on the raw, unreconciled DHIS2 line list (a sitrep-free sensitivity run).

**Data deliberately NOT used by the fits.** Beyond the confirmed-count reconciliation above,
the INSP situation-report streams (`data/insp_sitrep/…`) are **not** otherwise wired into any
spatiotemporal fit — the alert / total-case covariates come from the DHIS2 line list, not the
sitrep. Those weekly alert streams cover only ~6–33 largely already-affected zones over ~2–3
weeks, so they are of limited value for forecasting *unaffected* zones; §13 lists where they
could still help.

**Assumptions.** Onset is a valid transmission clock; the ~15 % of cases without a usable
onset are **imputed as `sample_date − Δ`**, where Δ is **drawn per record** from the fitted
onset→sample delay (the empirical bootstrap of the observed complete pairs when ≥30 exist,
else the parametric Exp(rate) fit) — a single stochastic imputation, not a mean substitution.
An implausible onset — a pre-outbreak year-typo, an onset after its own
sample, or one > 60 d before it (the shared `DELAY_MAX_PLAUSIBLE_DAYS` ceiling) — is treated as
missing and likewise imputed, then clamped
to the outbreak-week floor (the `onset_usable` guard, which stops a typo leaking a spurious
pre-outbreak week into the zone-week grid). The alias table resolves zone identity correctly.
The onset→sample delay is the fixed lab-linelist Exp fit by default, or re-estimated from the
data being processed when `ONSET_SAMPLE_DELAY_SOURCE = "data"` (e.g. for DHIS2 input).
**Strengths.** Onset-based, complete grid, reproducible spine.
**Weaknesses.** Onset-missing cases carry a slight temporal misalignment; geocoding
errors propagate silently; duplicate zone names are resolved only by a parenthetical
province tag.
**Possible extensions.** Onset for the missing ~15 % is already imputed by a per-record draw
from the (empirical) onset→sample delay (§0), not a sample-date substitution. The remaining
extension is **multiple imputation** — repeat the draw M times, run the pipeline per completion,
and pool — to carry the imputation uncertainty into the weekly aggregation and the reported
intervals (a single draw leaves those conditionally slightly narrow).

---

## 2. Nowcasting — reporting-delay / right-truncation correction  (`04_nowcasting.R`, `04b_epinowcast.R`)

**Intuition.** The most recent weeks always look artificially low because cases with
recent onset have not yet been sampled/reported. Nowcasting estimates *how many will
eventually be reported* and scales the recent weeks up — with uncertainty — so
models are not fooled by the dip.

**How the reporting delay is estimated.** The onset→sample delay is fit **once** to
the BDBV-2026 Ituri lab line list by interval-censored MLE
(`pipelines/cfr_analyses/`, exported to `DELAY_ONSET_SAMPLE_*` in `00_config.R`):
an **Exponential** with rate ≈ 0.228 d⁻¹ (mean ≈ 4.39 d), reloaded from the CSV at runtime
(`DELAY_ONSET_SAMPLE_RATE`, the source of truth). This delay drives
*both* nowcast variants below.

### 2a. epinowcast (primary, current forecast)
- **Description.** A Bayesian generative nowcast (`epinowcast`, cmdstanr backend)
  fit on the national daily **onset → sample** reporting triangle (max delay 21 d).
- **Equation.** Expected final onset-week count `λ_w` follows a weekly random walk;
  the reporting-delay CDF `F(d)` (log-normal, chosen because it won a retrospective
  WIS comparison over exponential/gamma) splits `λ_w` across report delays; observed
  cells are Negative-Binomial. The **correction factor** for onset-week *w* is the
  posterior ratio `E[λ_w] / Y^obs_w`, with a credible interval.
- **Key details.** One national factor per week, applied to every zone's observed
  count. The current partial week gets a sensible multiplicative lift
  (e.g. ×2.4, 90 % CrI 1.1–4.4 — note q5/q95 = **90 %**) instead of being discarded.
- **Assumptions.** The onset→sample delay is zone-invariant (a lab/reporting
  property); the weekly random walk captures the true-incidence trend.
- **Strengths.** Proper uncertainty; principled delay model.
- **Weaknesses.** A single national factor cannot capture zone-specific reporting;
  MCMC is slow (~100 s), so it is **not** re-fit inside cross-validation.
- **Possible extensions.** A partial-pooling hierarchical nowcast with province
  random effects, if per-zone triangles become dense enough.

### 2b. Deterministic Exp-CDF correction (fallback + used inside LFO)
- **Description.** A fast closed-form right-truncation correction used inside every
  cross-validation fold (where re-fitting epinowcast is infeasible), via
  `apply_nowcast_correction()`.
- **Equation.** `Y^nc = Y^obs / P(reported by now)`, with
  `P = 1 − exp(−rate · lag)` (rate ≈ 0.228 d⁻¹) from the fitted onset→sample delay. The
  multiplicative correction **factor** `1/P` is capped at **5×** (so `Y^nc ∈ [Y^obs, 5·Y^obs]`) —
  capping the factor, rather than the absolute count, guarantees `Y^nc ≥ Y^obs` for both the
  confirmed and (larger) suspected series.
- **Key details (a subtlety that matters).** The reference date must be the **end**
  of the training week (**cutoff + 7 d**), not its midpoint — the midpoint puts the
  cutoff week at lag 0, which the delay CDF treats as unobserved and zeroes, silently
  collapsing the very week every model forecasts from. This is enforced in the LFO
  loop.
- **Assumptions.** Same fitted onset→sample delay DATA as 2a, but applied through a closed-form
  **Exponential** CDF `1−exp(−rate·lag)` (epinowcast in 2a internally fits a **log-normal** delay,
  WIS-selected); a multiplicative inverse-probability-of-reporting
  correction. **Exception (documented):** when a current week's *national* observed total is
  zero the multiplicative factor cannot lift it, so the epinowcast path adds the national
  posterior-mean deficit back and distributes it across zones by each zone's recent (last ~4
  weeks) case share — an explicit additive step, resting on the assumption that the current
  week's still-unobserved onsets are spatially distributed like the last month's cases. It
  never fabricates cases into never-recently-active zones (their recent share ≈ 0), so the
  at-risk mask is unchanged.
- **Strengths.** Fast, closed-form, **causal** (uses only the delay CDF and counts
  up to the cutoff — no look-ahead), adequate for retrospective scoring.
- **Weaknesses.** No uncertainty; a point correction only.
- **Possible extensions.** Propagate the epinowcast factor's credible interval into
  the fold as a light-weight multiplier prior.

---

## 3. Epidemiological parameters  (`02_epi_params.R`, `00_config.R`)

### 3.1 Generation time (GT) — literature-anchored, with sources

The GT weights how strongly each *past* week of a source zone's incidence
contributes to *today's* import force (§5). It is a **discretised Gamma PMF**
(`make_gt_pmf()`), built daily then binned to weeks (`daily_to_weekly_gt()`).

**There is no Bundibugyo-ebolavirus-specific generation time or serial interval in
the literature.** Towner et al. 2008 characterised the *virus*; MacNeil et al. 2010
and Wamala et al. 2010 reported clinical features and an **incubation** period of
~5.7–7.4 d for the 2007 Uganda outbreak — but neither a generation nor a serial
interval. We therefore proxy the BDBV generation time with the well-characterised
**Zaire ebolavirus serial interval** (standard practice for EVD renewal models; the
serial interval ≈ generation time when latent and infectious profiles are similar,
as for filoviruses). All three profiles trace to peer-reviewed sources:

| Profile | Mean (SD) | Role | Source |
|---|---|---|---|
| **Short** | 12.0 d (6.5) | low-end sensitivity | just below the pooled 95 % CI lower bound (13.2 d) of Nash et al. 2024 |
| **Medium ★** | **15.3 d (9.3)** | **primary** | **WHO Ebola Response Team 2014, NEJM** — serial interval 15.3 d (SD 9.3) |
| **Long** | 18.0 d (10.5) | high-end sensitivity | just above the pooled 95 % CI upper bound (17.5 d) of Nash et al. 2024 |

Sources (URLs verified July 2026):
- **WHO Ebola Response Team, 2014.** *Ebola Virus Disease in West Africa — The First
  9 Months of the Epidemic and Forward Projections.* N Engl J Med 371:1481–1495.
  Serial interval **15.3 d (SD 9.3)**. doi:10.1056/NEJMoa1411100 —
  <https://www.nejm.org/doi/full/10.1056/NEJMoa1411100>
- **Nash et al., 2024.** *Ebola virus disease mathematical models and epidemiological
  parameters: a systematic review.* Lancet Infect Dis 24(10):e647–e657. Pooled
  random-effects **serial interval 15.4 d [95 % CI 13.2–17.5]**.
  doi:10.1016/S1473-3099(24)00374-8 —
  <https://www.thelancet.com/article/S1473-3099(24)00374-8/abstract>
- **Maganga et al., 2014.** *Ebola Virus Disease in the Democratic Republic of
  Congo.* N Engl J Med 371:2083–2091. DRC-specific mean serial interval
  **16.1 ± 4.4 d**. doi:10.1056/NEJMoa1411099 —
  <https://www.nejm.org/doi/full/10.1056/NEJMoa1411099>
- **Van Kerkhove et al., 2015.** *A review of epidemiological parameters from Ebola
  outbreaks to inform early public-health decision-making.* Sci Data 2:150019 —
  the range of published EVD serial-interval estimates the Short/Long scenarios
  bracket. doi:10.1038/sdata.2015.19 —
  <https://www.nature.com/articles/sdata201519>

*If ≥10 contact-tracing transmission pairs are available,
`estimate_si_from_contacts()` fits a Gamma to them by MLE and can replace the proxy —
a data-driven override that is used automatically when it exists.*

- **Description.** A discretised Gamma over daily lags τ = 1…max_tau, normalised.
- **Equation.** `g(τ) = [F_Γ(τ; k, θ) − F_Γ(τ−1; k, θ)] / Σ_{τ'} (·)`, with
  shape/rate from the (mean, SD) above; then weekly `G_w = Σ_{d=7(w−1)+1}^{7w} g(d)`,
  renormalised.
- **Key details.** The medium profile spreads weekly weight ≈ (0.18, 0.35, 0.25,
  0.13, 0.06, …) over the first ~5 weeks — epidemiologically sensible for a ~15-day
  GT.
- **Assumptions.** Serial interval ≈ generation time; Zaire-EBOV proxy transfers to
  BDBV.
- **Strengths.** Sensitivity is quantified across three literature-bracketed
  profiles rather than assumed away; every value traces to a citation.
- **Weaknesses.** No BDBV-2026-specific GT; the proxy could mis-state BDBV's true
  interval.
- **Possible extensions.** Replace with the contact-tracing MLE as pairs accrue;
  carry GT uncertainty as a full distribution rather than three point profiles.

### 3.2 Reproduction number R(t) — how it is computed, and its (limited) role

R(t) is **non-trivial and easy to misuse**, so it is worth being precise about
*which* R is computed *where*, and *what it is (and is not) used for*.

**Two distinct estimators, one narrow purpose.**
- **Display / national R(t)** — `EpiNow2` (a Bayesian renewal model, **not**
  EpiEstim) on the national **confirmed** onset series, with the discretised GT of §3.1 as the
  generation-time distribution. The series is onset-dated, so no reporting-delay convolution is
  applied; instead a **right-truncation model** (`trunc_opts`, a Gamma with mean ≈ 4.4 d from the
  onset→sample delay) nowcasts the incomplete recent weeks — without it the recent tail is
  under-counted and R(t) drops spuriously below 1. Reported with
  **60 % and 90 % credible bands** (q20/q80 and q5/q95 — never 50 %/95 %) and the
  three GT profiles for sensitivity (`rt_national.pdf`). This is the *communicated*
  R(t) and is not fed into the invasion hazard.
- **Projection R** — a fast **Cori-style instantaneous** estimate
  `R̂(t) = I(t) / Σ_k g(k) I(t−k)` on the national weekly incidence, computed inside
  the workhorse/each fold by `estimate_R_local()` (`15_workhorse.R`): it averages
  the last **two reporting-complete** weeks (the most recent onset week is
  truncation-low and is skipped), and is clamped to [0.5, 5]. The *same* estimator
  is used by the Bayesian suite's offset projection (`bayes_forecast_offsets()`),
  so the two paradigms share an identical incidence projection and differ only in
  the *destination* invasion hazard.

**What R is used for — and the key thing it is NOT.** R enters **only** the one-week
projection of *source-zone* incidence needed to build the horizon-2 import force:
`Λ_i(h=2)` needs each source zone's next-week cases, projected as
`Ŷ_j(t+1) = R̂ · [GT-weighted own incidence]_j + imported pressure`. The **invasion
probability itself never multiplies by R** — the import→first-case conversion is the
*separately calibrated* coefficient β (§5.1). This is deliberate: R is a
*self-sustaining transmission* number, whereas invasion is an *importation* event; a
naïve `1 − exp(−R·Λ)` was substantially over-confident (roughly 1–2 orders of magnitude in early testing). For the 1-week horizon R plays no
role at all (Λ(h=1) uses observed data).

**Assumptions / limitations (documented, not hidden).** The projection R is a single
*national* number applied to every source zone's own incidence — it assumes spatially
homogeneous short-term growth, which is a simplification for a sparse, spatially
heterogeneous outbreak, but is defensible because (a) it acts only as a one-week
nuisance projection of *source* incidence, (b) per-zone R with ~a handful of cases per
zone is hopelessly unstable, and (c) the dominant driver of Λ is the mobility routing,
not the exact growth of the source. *Extension:* carry EpiNow2's posterior R (with
uncertainty) into the projection, or a partially-pooled per-province R, once zone-level
incidence is dense enough to estimate them.

### 3.3 Ascertainment / CFR — treated as UNCERTAIN, not fixed

Ascertainment ρ (test positivity) has a nominal value 0.45 and a sensitivity grid
{0.30, 0.45, 0.60} (`ASCERTAINMENT_*`); the delay-adjusted cCFR (0.264) comes from the
CFR pipeline. **Ascertainment is deliberately kept out of the primary forecast** and
carried as a *band* where it is used:
- The primary invasion probability `p_case` is **ascertainment-agnostic** (§0): it is
  P(first confirmed-case onset), calibrated to observed frequencies, and never divides
  by ρ.
- The secondary **infection-scale** score inflates the case hazard to an implied
  *infection* hazard `p_infection = 1 − exp(−μ/ρ)`, and — because ρ is genuinely
  uncertain — is reported as a **band** `[p_infection_lo, p_infection_hi]` spanning the
  grid (ρ = 0.60 gives the lower bound, ρ = 0.30 the upper; `forecast_workhorse()`
  emits `p_infection_lo`/`p_infection_hi`). `p_infection − p_case` flags
  *silent-introduction* risk (zones likely seeded but slow to confirm), now with an
  honest uncertainty range rather than a single assumed ρ.

---

## 4. Mobility matrices (how the virus travels)  (`03_mobility_matrices.R`, `06_simple_models.R`)

**Intuition.** The virus moves with people. An **outward** mobility matrix `W` says
what fraction of travel *out of* zone *j* arrives in zone *i* (`W[j,i]`); it is
row-stochastic (each origin's outflow sums to 1) with a **zero diagonal** (movement
*between* zones only). Inflow to *i* uses column *i*, i.e. `t(W)`.

Variants are built so mobility is an independent axis of sensitivity:

- **M1 / M2** — Flowminder **short-trip** proportions from the epicentre cohort
  (Bunia/Mongbalu/Rwampara), static (M1) or time-evolving (M2). Outbreak-specific
  but pooled over three epicentre rows / ~40 destinations.
- **M3** — full Flowminder national origin–destination matrix. Broadest empirical
  structure; March-2026 baseline, some zones absent.
- **M4 / M4b** — **gravity** GLM (`log pop_i, log pop_j, log dist`) fit to M3 and
  predicted for all 519×519 pairs; M4b swaps the power-law deterrence for an
  exponential kernel (a functional-form sensitivity).
- **M5** — **radiation** model (parameter-free; Simini et al. 2012, Nature
  484:96–100).
- **M6** — OSRM travel-time decay kernels (exponential / power-law).
- **M7** — M3 augmented with IDP displacement flows.
- **M8 ★ (primary outward kernel)** — **composite**: short-trip proportions for the
  epicentre rows, gravity (M4) elsewhere — the most empirically-grounded epicentre
  outflow plus complete national coverage.
- **M9 / M10** — multi-kernel ensemble (mean of M4/M5/M6a elsewhere) / radiation
  composite; robustness checks on the non-epicentre kernel.
- **M13 / M14 — Flowminder cohort composites.** The empirical-CDR analogue of M8/M10:
  the ten Ituri/Nord-Kivu/Tshopo **cohort subscriber-day presence** rows (`build_M_cohort`)
  with gravity (M13) or radiation (M14) elsewhere, in travel-time and road-km (`-dist`)
  variants. Gated by `INCLUDE_COHORT_MODELS` (default on).
- **M15 — Flowminder combined inflow+outflow static.** A base kernel from the *full*
  Flowminder OD data using **both** perspectives: because the inflow and outflow exports are the
  same directed table, the combination that adds information is the **symmetrised total two-way
  flow** `S = O + Oᵀ` (trips out + trips in). Distinct from the directed M3 and materially less
  sparse (origin coverage ≈ 80% → ≈ 84%).
- **M16 — Flowminder cohort + static.** The direct analogue of M13/M14 that pairs the cohort
  presence source rows with the M15 static flows elsewhere — cohort where available, static
  otherwise (the empirical-flow counterpart of cohort+gravity / cohort+radiation).
- **M17 (+ M17-dist) — all-kernel consensus ensemble.** The element-wise mean of every distinct
  base hypothesis {M3 OD, gravity M4, radiation M5, travel-time decay M6a, static M15} (a convex,
  row-stochastic combination) with the empirical cohort/epicentre source rows overlaid — the single
  "combine all mobility matrices" option, extending M9 to the full data set. M15/M16/M17 are gated by
  `INCLUDE_FLOWSTATIC_MODELS` (default on) and carried in the **Bayesian grid** (§6), swept
  exhaustively over GT {short, medium, long} × covariates {none, geo, full}.
- **`-dist` variants (M4-dist / M8-dist / M9-dist / M10-dist / M11-dist; M13/M14-dist; M17-dist)** — each of the
  distinct kernel families rebuilt with the non-epicentre deterrence keyed on OSRM
  **road distance (km)** instead of travel **time (minutes)** — the gravity GLM re-fits on
  km and radiation re-orders intervening opportunities by km (not a rescale). A
  deterrence sensitivity axis (does the ranking survive time→distance?); travel time is
  usually the better human-mobility proxy, so these are robustness checks. Built on demand
  when the OSRM road-distance matrix is present; carried in the **Bayesian grid** (§6), both
  covariate-free and with the geo / full covariate sets. Epicentre short-trip rows are
  distance-agnostic and reused.

### 4.1 M11 — inward / meeting-location force of infection (manuscript-motivated) ★new
- **Description.** A **two-sided, frequency-dependent** contact kernel motivated by
  the meeting-location force-of-infection formulation (Mills 2026, eq. 6–7): rather
  than "cases flow one-way from *j* to *i*", susceptibles from *i* and infectious
  travellers from *k* **co-locate** at shared meeting zones *l*, and transmission
  scales with the *density* of infectious presence there. Built by
  `build_inward_contact_matrix()` (`06_simple_models.R`).
- **Equation.** Let `P` be a **presence matrix** — the fraction of a resident's time
  spent in each zone: diagonal = home retention `home_fraction` (0.70; see below),
  off-diagonal `= (1 − home_fraction)·W`, row-normalised so each resident's presence
  sums to 1. With effective meeting-zone populations `N_eff = Pᵀ N`, the inward FOI
  kernel is
  `M = P · diag(1/N_eff) · Pᵀ`, symmetrised `M ← (M + Mᵀ)/2`, mean-scaled.
  The import force is then `Λ_i = (M · Ỹ)_i` — **exactly the same convolution** as
  the outward model, so `M` drops into `compute_foi()` in place of `t(W)` with **no
  change to the workhorse** (extend, not replace).
- **Key details.** `home_fraction = 0.70` (`MOBILITY_HOME_FRACTION`) is a documented
  modelling assumption (commuting fractions 18–40 % in Mills 2026 ⇒ home fractions
  0.60–0.82; the central 0.70 is used and is sensitivity-analysable). `M` is
  symmetric by construction, so importation pressure between two zones is shared
  rather than purely directional. It reranks risk toward zones that *share meeting
  locations* with affected zones, not merely those *downstream of outward flow*.
- **Assumptions.** Frequency-dependent transmission (per-capita contact at meeting
  zones); presence splits between home and mobility-weighted destinations;
  susceptible and infectious travel with the same presence matrix.
- **Strengths.** Encodes two-sided mixing and crowding (`1/N_eff`) that the one-way
  gravity flow cannot; symmetric and mechanistic; zero-cost to the existing engine.
- **Weaknesses.** `home_fraction` is assumed, not estimated; treats all meeting
  zones as equally mixing; still built on the static outward `W`.
- **Possible extensions.** Estimate `home_fraction` from the short-trip cohort; make
  meeting-zone mixing intensity covariate-dependent (e.g. market/health-facility
  density); time-vary `P`.

*Audit note.* M8's construction (row-stochastic, zero diagonal, epicentre outflow to
the correct observed destinations) and M11's symmetry / row-normalisation / `N_eff`
were verified numerically (M symmetric to ~1e-16; P row-sums = 1).

---

## 5. The mobility-informed renewal invasion model (primary)  (`06_simple_models.R`, `15_workhorse.R`)

This is the centrepiece. §5 gives the **formulation as equations** and the **exact
code path** (Task 1a), then the 7-part description.

### 5.1 Formulation → code, step by step

1. **GT-weighted recent incidence** of each source zone *j* (nowcast-corrected):
   `Ỹ_j(t) = Σ_{k≥1} G_w(k) · Y^nc_{j, t−k}`.
   → the `Σ_k G_weekly[k]·Y(t−k)` loop inside `compute_foi()`.
2. **Import force of infection** into an at-risk zone *i*:
   `Λ_i(t) = Σ_{j≠i} W[j,i] · Ỹ_j(t) = (t(W) · Ỹ)_i` (outward kernel), or
   `Λ_i(t) = (M · Ỹ)_i` (inward kernel M11).
   → `Lambda <- t(W) %*% Y_weighted` in `compute_foi()`. The zero diagonal makes Λ
   **purely external** — exactly the invasion signal for a zone with no local cases.
3. **Calibrated arrival probability — the crucial point.** The import→first-case
   coefficient **β is *learned*, not assumed**, by a complementary-log-log
   regression of realised invasions on the import force with **log Λ as an offset**:
   `cloglog[ P(invade_i) ] = log[−log(1 − p_i)] = β₀ + log Λ_i`
   ⇒ `p_i = 1 − exp(−e^{β₀} · Λ_i) = 1 − exp(−β · Λ_i)`, `β = e^{β₀}`.
   → `fit_import_beta()` / `fit_import_model()` in `15_workhorse.R`
   (`glm(family = binomial("cloglog"), offset = log Λ)`, Firth-penalised via
   `brglm2`). Because β is fit to the *observed* (rare) invasion frequency, the
   probabilities are **calibrated by construction** — no post-hoc Platt scaling.
   *(A naive `1 − exp(−R·Λ)` was substantially over-confident (roughly 1–2 orders of magnitude in early testing), because R is a
   self-sustaining transmission number, not an import scale — hence R does not enter
   the invasion hazard.)*
4. **Two horizons, cumulative hazard.** The 1-week probability uses Λ from observed
   data; for 2 weeks the model projects one week of local + imported incidence
   forward (using `estimate_R_local()`), re-computes Λ, and reports the **cumulative
   hazard** `p(≤h) = 1 − exp(−Σ_{h'≤h} μ_{h'})`, so a first invasion is counted once.
   → `forecast_workhorse()`.
5. **Observation model.** Poisson (default) or **Negative Binomial** (dispersion from
   the *renewal residuals* via `estimate_negbin_dispersion()` /
   `.workhorse_temporal_dispersion()`), for the count quantiles used by WIS. NB is
   the selected primary.
6. **Masking.** Affected zones → `p = NA`, enforced inside the model and every
   downstream product.

### 5.2 The model in seven parts
- **Description.** A zone becomes invaded when infectious pressure arrives from
  affected zones via travel; that pressure is GT-weighted recent incidence routed
  through the mobility matrix, converted to a first-case probability by a single
  calibrated coefficient β.
- **Equation.** `p_i = 1 − exp(−β · Λ_i)`, `Λ_i = Σ_{j≠i} W[j,i] Σ_k G_w(k) Y^nc_{j,t−k}`,
  `β = e^{β₀}` fit by cloglog with offset `log Λ`.
- **Key details.** One free parameter (β) in the base model; NB observation; rolling
  nowcast; cumulative hazard across horizons.
- **Assumptions.** Invasion hazard is proportional to external import force; a single
  national β (spatially homogeneous import→invasion conversion); mobility static over
  the forecast window; no within-zone growth term (correct for zero-case zones).
- **Strengths.** Simple, transparent, mechanistically interpretable; calibrated by
  construction; best discrimination *and* calibration in cross-validation; robust
  because it has essentially one parameter; directly answers the decision question.
- **Weaknesses.** A single national β assumes spatial homogeneity; β is estimated
  from very few events, so its value is uncertain; mobility is static.
- **Possible extensions.** Covariate-modulated β (§5.3); inward FOI (M11, §4.1);
  Bayesian β with a regularising prior (§6); time-varying mobility.

### 5.3 Covariate-modulated β (spatial heterogeneity) + the `d_min` sign question ★
- **Description.** The homogeneous β becomes zone-specific,
  `β_i = exp(β₀ + Σ_m γ_m x_{m,i})`, keeping `log Λ` as a fixed offset so mobility
  drives the *level* while penalised covariates modulate *susceptibility*
  (`fit_import_model()`).
- **Covariates.** An **alert import/local force** (the same GT×mobility construction
  on unconfirmed alerts — a leading indicator), **log-population**, **health-site
  density**, **CCVI** deprivation, and **`d_min`** (travel time to the nearest
  already-affected zone — frontier proximity). Coefficients are **Firth-penalised**
  (`brglm2`) and guarded: any standardised |coef| > 15 is rejected and the fit falls
  back to intercept-only.
- **Why `d_min` carries a *negative* risk effect (the question raised).** `d_min` is
  the travel time from an at-risk zone to the **nearest affected zone**. A **larger**
  `d_min` means the zone is **further from the outbreak frontier**, so it is **less**
  exposed to importation and **less** likely to be invaded — hence the fitted hazard
  ratio is **below 1** (a protective/negative coefficient). It is not that distance
  is protective *per se*; it is that frontier proximity (small `d_min`) is the risk
  factor, and `d_min` measures its opposite. This is the single most interpretable
  covariate and is exactly what an epidemiologist expects: risk decays with distance
  from the active frontier.
- **Why `log_pop` carries a *negative* effect (the second sign question).** In the
  Bayesian fit the population coefficient is credibly **below 1** for M8
  (HR ≈ 0.53, P(HR>1) ≈ 0.002) yet **null** for M11 (HR ≈ 0.80, CrI 0.57–1.17 spans 1). This is
  **not** "big cities are safer" — it is a *conditional-on-the-offset residual*, and
  it is informative. The mobility import force Λ enters as a **fixed offset**, and the
  gravity kernels (M4/M8) build flows from population (`∝ Nᵢ·Nⱼ / d`), so **Λ already
  assigns more absolute import pressure to large zones**. `log_pop` therefore estimates
  what population does *beyond* that offset. A negative residual means the **destination
  hazard is closer to frequency-dependent (per-capita) than the density-dependent
  gravity offset assumes**: conditional on the same import pressure, a first case is
  *less* likely per extra resident. The decisive check is **M11**, whose inward FOI is
  frequency-dependent *by construction* (it divides by the effective meeting-population
  `N_eff`): there `log_pop` collapses to ≈ 1 (no residual), exactly as expected if the
  M8 negative sign is the density-vs-frequency correction. So the two signs are
  mutually consistent and validate the M11 formulation, rather than indicating a bug.
  (A residual selection effect also pushes the same way: the largest zones — Bunia,
  regional capitals — were invaded *early* and leave the at-risk set, so the surviving
  large zones are atypically low-risk.)
- **Important caveat (honest reporting).** On the *full* series the joint covariate
  fit **diverges**: with the mobility offset conditioned on, `alert_local` and test
  `positivity` **quasi-separate** the rare invasions (coefficients → ±∞ ~10¹⁵), trip
  the guard, and the fit collapses to intercept-only — so the live `M8-cov` forecast
  is numerically identical to `M8-med`. The honest way to read the covariates is the
  **association screen** of §5.5 (marginal HRs with adjustment annotations), **not** a
  single unstable joint fit. This is *the* reason the Bayesian suite (§6) exists: a
  weakly-informative prior regularises exactly this separation.
- **Equation.** `p_i = 1 − exp(−e^{β₀ + Σ_m γ_m x_{m,i}} · Λ_i)`.
- **Assumptions / Strengths / Weaknesses.** As §5.2, plus: covariates add
  identifiable susceptibility signal *only in the smaller folds*; on the full series
  they are not identifiable (separation).
- **Possible extensions.** The Bayesian regularised fit (§6); drop endogenous
  covariates (positivity is `confirmed/tests`, circular with the outcome — see §5.4).

### 5.4 Additional independent options (each an *extra* variant, cross-validated)
Each toggles independently; with all off every variant collapses **exactly** to the
base workhorse (verified: identical β and rankings), so these are strict additions.
1. **Reporting-rate structure.** Under-reporting source zones exert more true
   pressure than their counts show. Because `Λ_i = Σ_j W[j,i]·Ỹ_j`, dividing source
   row *j* of `W` by a relative reporting rate `r_j` (health-site density proxy,
   geometric-mean-normalised, bounded [0.25, 4]) is algebraically identical to using
   true incidence `Ỹ_j / r_j`, while observed counts still define the mask/outcome.
2. **Nowcast structure.** `corrected` (delay-adjusted, default) vs `raw`
   (right-censored) — measures what the truncation correction buys.
3. **Completeness-weighted calibration.** Each training week's contribution to the β
   regression is weighted by its nowcast completeness, so recent right-truncated
   weeks count for less.
4. **Endogeneity guard.** Test `positivity = confirmed/tests` is **circular** with the
   invasion outcome and is therefore treated as a *diagnostic* covariate in the
   association screen, **not** a driver of the featured forecast — only exogenous
   covariates (log_pop, CCVI, healthsite_density, `d_min`) legitimately enter β.

### 5.5 Parameter estimates & invasion drivers  (`covariate_associations()`, `plot_model_parameters()`)
- **(a) R(t)** — the renewal reproduction number (EpiNow2, on the onset series)
  plotted with **60 % and 90 % credible bands** (q20/q80 and q5/q95 — *not* 50 %/95 %)
  and an R = 1 line, plus the three GT-profile means (`rt_national.pdf`).
- **(b) Import coefficient + association screen** (`best_model_parameters.pdf`). β₀ =
  `e^{intercept}` is the dominant driver. Because the joint fit diverges (§5.3), each
  candidate covariate is shown as a **marginal** invasion HR per +1 SD with a 95 % CI,
  annotated by what happens **after adjusting for the mobility offset** (significant /
  attenuated / separates). Headline reading: imported alert pressure (HR ≈ 1.9) and
  frontier proximity (`d_min` HR ≈ 0.3) are the strong signals; population has ≈ no
  marginal effect; after conditioning on mobility most covariates attenuate — *where
  the virus can travel from the affected cluster dominates static zone characteristics*.

### 5.6 Ensembles (combining the best members)  (`18_ensemble.R`)
- **Description.** A **pre-specified**, diverse set of strong renewal formulations
  (composite-gravity M8 Poisson & NB, plain gravity M4, ensemble-kernel M9, radiation
  M10, and now inward-FOI M11) combined as the forecast hubs do (Reich 2019; Cramer
  2022; Ray 2023).
- **Equation.** mean = linear opinion pool `p̄ = mean_m p_m` (Vincent averaging of
  count quantiles); median = median probability / median-of-quantiles.
- **Uncertainty (Task 5).** For a rare event with ~17 signals there is no reliable
  single "best" formulation, so the honest, leakage-free band on each zone's
  probability is the **member disagreement** (min / median / max across members),
  drawn by `plot_forecast_uncertainty()` and `plot_spacetime_forecast()`. It is wide
  where the mobility/GT/observation choice matters and narrow where members concur.
- **Strengths / Weaknesses / Extensions.** Robust to any single errant member; no
  selection on test folds (scored in the same LFO). Linear pools are slightly
  under-confident for rare events. Could weight members by CV skill (stacking) — which
  the Bayesian suite does explicitly (§6).

---

## 6. Bayesian renewal suite (brms)  (`21_bayesian_renewal.R`) ★new

- **Description.** A Bayesian re-implementation of the **same** renewal invasion model
  (§5), fit with `brms` / Stan (cmdstanr). It is *not* a different model — it is the
  identical cloglog-with-offset structure, so it re-uses the mobility import force Λ
  verbatim and simply replaces penalised MLE with a posterior. Its purpose is to
  **regularise the quasi-separation** that diverges the frequentist joint covariate
  fit (§5.3) and to produce **proper posterior uncertainty** on both parameters and
  invasion probabilities.
- **Equation.**
  `invaded_i ~ Bernoulli(p_i)`, `cloglog(p_i) = β₀ + Σ_m γ_m x_{m,i} + log Λ_i`
  (`offset(logLam)`), with **weakly-informative priors**
  `β₀ ~ Normal(−3, 2)` (intercept; −3 on the cloglog scale ≈ a low base invasion
  hazard) and `γ_m ~ Normal(0, 1)` on **standardised** covariates (a coefficient of
  ±1 = a hazard-ratio of e^{±1} ≈ 2.7 per SD — regularising but not committal).
  → `fit_bayes_renewal()`.
- **Prediction (offset handled explicitly).** The model is fit with `offset(logLam)`;
  prediction zeroes the offset in the design matrix and **adds the horizon log-Λ back
  manually** (`mu <- exp(sweep(eta0, 2, nd$.off, "+"))`, `p = 1 − exp(−mu)`), then
  accumulates the cumulative hazard across horizons and masks affected zones to `NA`
  — mirroring the frequentist projection exactly (`predict_bayes_invasion()`).
- **Comprehensive grid — the Bayesian analogue of the frequentist suite, ALL four
  axes, `bayes_default_grid()` (single source of truth for the current forecast AND
  cross-validation).** **Several dozen models** — a **curated (fractional), not fully-crossed**
  design (it collapses to a small core when the OSRM road-distance / cohort / flowstatic matrices
  are absent, since those entries auto-drop). It sweeps: **mobility** — the base kernel families
  {gravity M4, composite-gravity M8, multi-kernel ensemble M9, radiation-composite M10,
  inward-FOI M11, **cohort composites M13/M14**, **combined-Flowminder static M15**, **cohort+static
  M16**, **all-kernel consensus ensemble M17**} **plus their OSRM road-distance (`-dist`) twins**
  {M4-dist … M11-dist, M13/M14-dist, M17-dist — the gravity GLM re-fit on kilometres and radiation
  re-ordered by km, §4} (M1/M2/M3/M5/M6/M7/M4b are components or minor variants folded into the
  composites); **generation time**
  {short, medium, long}; **covariates** {none, geo = log_pop/CCVI/`d_min`, **full = + health-site
  density**} — the complete EXOGENOUS set (test `positivity` and the alert covariates
  are excluded as endogenous/leakage, §5.5); and the **observation process** {**cloglog**
  (default, principled) + a **logit-link** sensitivity}. cloglog is not an arbitrary
  choice: it makes the import force a proper **log-cumulative-hazard offset** and yields
  exactly `p = 1 − exp(−β·Λ)`; the logit variant is a robustness check, and prediction
  reconciles both on the **per-week hazard scale** (`μ = exp(η)` for cloglog,
  `μ = −log(1−logit⁻¹(η))` for logit) so the cumulative-survival and ascertainment
  logic are identical. (Poisson vs NB is the *count* observation for the source-incidence
  projection, §3.2 — shared with and varied in the frequentist suite — not an axis of the
  binary invasion likelihood.) The core comprises: 5 base kernels + 5 `-dist` base + 2 GT (on M8)
  + 4 base-covariate (M8/M11 × geo/full) + 10 `-dist`-covariate + 1 logit (M8-geo); plus the
  cohort composites (M13/M14 + `-dist`, none/full [+geo]) and the **Flowminder combined-static
  family** — M15, M16, M17 and M17-dist swept *exhaustively* over GT {short, medium, long} ×
  covariates {none, [geo], full}. **Every
  model in the grid is also cross-validated** on the same leave-future-out folds (lighter
  sampling, `iter = 600`), so the rankings, over-time and detection figures compare the full grid.
- **How to read a model code** (e.g. `Bayes-M10-dist-geo`, shown spelled-out in every figure
  via `.bayes_model_label()`): `Bayes-<mobility>[-dist]-<gt/variant>[-<covariates>][-<link>]` —
  * **mobility**: `M4` = gravity, `M8` = composite-gravity, `M9` = multi-kernel ensemble,
    `M10` = radiation-composite, `M11` = inward meeting-location FOI; a `-dist` token means the
    OSRM **road-distance** twin of that kernel;
  * **gt/variant**: `med`/`short`/`long` = generation time (the M11 base is labelled `-inward`);
  * **covariates**: (none) / `geo` = + log_pop, CCVI, `d_min` / `full` = + health-site density;
  * **link**: (none) = cloglog / `logit` = logit sensitivity.
  So `Bayes-M10-dist-geo` = *radiation-composite mobility on **road distance**, medium GT, geo
  covariates, cloglog*.
- **The FEATURED single Bayesian model** (the one whose maps carry the `bayes_` prefix and
  populate `bayes_risk_scores_all_zones.csv`) is the best cross-validated single model by the
  **calibration-aware CV composite** — the summed within-horizon ranks of AUC-PR skill + mean
  rank-of-truth + log-score, **pooled across both forecast horizons** (lower = better) on the
  leave-future-out folds, the same criterion as
  the frequentist featured model — NOT the loo predictive-stacking weight, and not hardcoded.
  Ensembles are excluded so the featured single model always has a current forecast. In the current
  run this is **`Bayes-M10-med`** (radiation-composite renewal), which also tops the AUC-PR-skill
  leaderboard. The loo predictive-stacking weights instead select the loo-stacked posterior
  *ensemble* (`bayes_ensemble_` products, whose weights concentrate on the `-dist` covariate models
  here — a different structure from the single featured model).
- **Model averaging — robust loo stacking.** The fits are combined by **loo
  predictive stacking** (Yao et al. 2018, Bayesian Analysis 13:917–1007). Because the
  at-risk set (Λ > 0) differs by mobility kernel, all models are first aligned to a
  **common complete-case (zone, week) row-set**; each model's pointwise log-likelihood
  is then evaluated on those aligned rows (`brms::log_lik`, each with *its own* mobility
  offset), collapsed to per-observation ELPD, and stacked with `loo::stacking_weights()`
  on the [obs × model] ELPD matrix. This tolerates a transient cmdstanr chain failure
  (which would otherwise leave models with differing draw counts and break
  `loo_model_weights`); if stacking still cannot run it falls back to an equal-weight
  average, **honestly relabelled** "Bayes-stack (equal-weight)".
- **Reported quantities.** `bayes_param_table()` gives posterior **hazard ratios**
  `e^{γ}` with median and **90 % credible interval** (q05/q95 — *not* 95 %) and the
  posterior probability of direction `P(HR > 1)`; predictions carry posterior mean +
  **90 % CrI**. Convergence is checked by R̂ (fits converged, R̂ ≈ 1.00).
- **Key details.** cmdstanr backend, 2000 iterations / 2 chains, `adapt_delta = 0.9`;
  the LFO closure `make_bayes_lfo_model()` re-fits per fold (`iter = 600` for speed)
  so the Bayesian models are scored in the **same** leakage-free CV as the frequentist
  ones.
- **Assumptions.** As §5, plus the priors above; standardised covariates so the
  `Normal(0,1)` prior is scale-appropriate.
- **Strengths.** Regularises the separation (finite, interpretable coefficients where
  the MLE ran to ±∞); full posterior uncertainty on parameters *and* forecasts;
  principled stacking across structural assumptions; the same mechanistic core so it
  is directly comparable to the frequentist model.
- **Weaknesses.** Slower (MCMC); priors are a modelling choice (deliberately weak, but
  they do shape the tails with so few events); stacking weights themselves are
  uncertain with ~17 events.
- **Possible extensions.** Hierarchical β by province (partial pooling); a spatial CAR
  term on the linear predictor; horseshoe priors for covariate selection; joint
  delay+invasion estimation to propagate GT/nowcast uncertainty into β.

### 6.1 Full map/decision PARITY — one renderer, best-freq + best-bayes, every figure labelled

Every map and decision product is produced by a **single renderer**,
`render_map_suite(rs, model_label, file_prefix)`, called once per **featured** model so
the frequentist and Bayesian paradigms yield an **identical, model-labelled set** of
figures — no map exists for one paradigm without the other, and **no model is chosen
arbitrarily for plotting**:
- **Frequentist** = the skill-selected best renewal model (`primary_method`), no prefix.
- **Bayesian** = the featured single Bayesian model (`best_bayes_method`), chosen by the
  **calibration-aware CV composite** (summed within-horizon ranks of AUC-PR skill + mean
  rank-of-truth + log-score, pooled across both forecast horizons) over the cross-validated single
  models — the same forecast-skill
  criterion as the frequentist featured model, NOT the loo predictive-stacking weight; not
  hardcoded; prefix `bayes_`.
- **Bayesian ensemble** = the loo predictive-stacking posterior mixture (the analogue of the
  frequentist ensemble; the stacking weights select the ensemble, not the featured single
  model), prefix `bayes_ensemble_`.

The suite, for each model, is: invasion choropleths (national + Ituri zoom, h1 & h2),
**per-province** invasion maps (Ituri / Nord-Kivu / Haut-Uele), **prob × vulnerability
bivariate choropleths** (national + each province), ranked at-risk bars, and the
preparedness-priority scatter / component bars / map. **Every figure's title names the
exact model that produced it** (e.g. "…— Frequentist: Renewal-M8-med" /
"…— Bayesian: Bayes-M9-med · multi-kernel ensemble · med GT"). Risk scores come from
`compute_risk_scores()` + `add_risk_indices()` on that model's forecast, saved to
`bayes_risk_scores_current.rds`. The single-model spatial/space-time diagnostics
(`lfo_forecast_vs_outcome`, `spacetime_risk`, `spatial_error_map`) are likewise rendered
for **both** the best frequentist and best Bayesian model (bayes_ prefix), and the
Bayesian LFO models are scored in the *same* `run_invasion_lfo`, so the rankings,
skill-over-time and predicted-vs-observed figures compare both paradigms.

### 6.2 Frequentist vs Bayesian — separation and robustness comparison

The two paradigms are reported **separately** and then **compared**, so agreement is a
robustness check rather than an assumption:
- **Ranking on the same folds** (`plot_freq_bayes_ranking()`,
  `freq_vs_bayes_ranking_h1.pdf`): AUC-PR skill per method coloured by inferential
  family, on the identical leave-future-out folds. The two families fit the *same*
  cloglog renewal model; only the estimation (penalised MLE vs posterior) differs.
- **Per-zone agreement** (`plot_freq_bayes_agreement()`,
  `freq_vs_bayes_agreement_h1.pdf`): the frequentist featured `p_case` vs the Bayesian
  stacked `p_invasion` for every at-risk zone, with the Spearman correlation and the
  y = x line. On-diagonal points = the paradigms concur; systematic departures localise
  where the prior/regularisation matters.

**Bayesian figures.** `bayes_parameters.pdf` — a posterior figure restricted to the
**best-fitting subset** of the suite (top few models by loo predictive-stacking weight),
since the full suite's spelled-out model names are far too long to serve as y-axis ticks.
Each shown model gets a short code (B1, B2, …) on the axes and is spelled out in a custom
**model-key** legend panel (assumptions + stacking weight). **(a)** the import coefficient
β₀ = exp(intercept) for each shown model and **(b)** covariate hazard ratios (90 % CrI)
for the shown models that include covariates (`.bayes_model_label()`);
`bayes_invasion_uncertainty_h1.pdf` /
`bayes_stacked_invasion_uncertainty_h1.pdf` (top zones' posterior invasion probability
with 90 % CrI, coloured by province), `bayes_stacking_weights.pdf`, the decision maps
of §6.1, and `bayes_parameters.csv`.

**Bayesian posterior-distribution & sensitivity products (★new).** Three further figures make the
Bayesian uncertainty and its assumptions explicit:
- `bayes_posterior_densities.pdf` — the **full posterior densities** (ridgelines, 5/50/95 % lines)
  of β₀ per model and of the covariate effects, i.e. the whole distribution rather than only the
  median + interval.
- `bayes_gt_posterior.pdf` (+ `.csv`) — a **generation-time preference**: the GT is otherwise a
  fixed assumption, so the featured model is refit over a grid of GT means and each is weighted by
  its leave-one-out predictive density combined with the SI literature prior (Normal(15.3, 1.5)).
  This is a **loo-predictive pseudo-BMA+ preference weight, NOT a jointly-estimated GT posterior**
  (the identical-observations precondition is enforced; Pareto-k̂ is checked), and the prior is
  informative so the result is prior-dominated — read as "which GT the invasion data prefer".
- `bayes_nowcast_sensitivity.pdf` (+ `.csv`) — a **sensitivity of β₀ to the two-stage nowcast
  INPUT**: the featured model is refit on raw / epinowcast-corrected / fast-delay-CDF training
  counts. This probes robustness to the *choice* of point nowcast; it does **not** propagate the
  nowcast's own uncertainty (a joint delay+invasion model, or pooling nowcast posterior draws,
  would — the honest caveat on the two-stage plug-in), and β₀ shifts partly reflect at-risk-set
  and outcome recomposition across variants.

---

## 7. Comparator models (kept, cross-validated on the same folds)

Each is a genuine independent check on the renewal model, scored in the identical LFO.

### 7.1 hhh4 endemic–epidemic  (`07_hhh4_model.R`)
- **Description.** The Held/Meyer `surveillance::hhh4` count model with
  autoregressive, neighbourhood (mobility-weighted) and endemic components.
- **Equation.** `λ_{it} = e_{it}·ν_{it} + φ_{it}·λ_{i,t−1} + τ_{it}·Σ_{j≠i} w_{ji} λ_{j,t−1}`
  (endemic + autoregressive + neighbour-driven), with log-linear predictors on each.
- **Key details / Assumptions.** Population offsets; neighbourhood weights from the
  mobility structure; NB counts. A rich model for a sparse signal.
- **Strengths.** Mature, well-tested framework; near-calibrated; independent check.
- **Weaknesses.** Heavier, more parameters, occasional convergence issues.
- **Extensions.** Random effects by zone; power-law neighbourhood weights.

### 7.2 Gravity-decay baseline (B4)  (`05_baseline_models.R`)
- **Description / Equation.** Invasion hazard `p_i = 1 − exp(−k · Σ_j N_j^β d_{ij}^{−γ} Y_j)` —
  mobility as pure gravity, no renewal/GT. The gravity **shape** exponents (β, γ) are **fixed**
  (population 0.5, distance −1), while a single overall **scale k** is fit by MLE to the observed
  invasions (a grid search maximising the Bernoulli likelihood).
- **Strengths / Weaknesses.** Simple mechanistic baseline; only its overall scale is calibrated
  (the gravity shape is fixed, not fit) and it ignores infectious timing (no GT weighting).
  **Extension:** add the GT weighting and let the mobility routing (not raw gravity) carry the
  force → recovers the renewal model.

### 7.3 Distance-weighted baseline (B1)  (`05_baseline_models.R`)
- **Description / Equation.** `μ_i = Σ_j d_{ij}^{−α} Y_j / Σ_j d_{ij}^{−α}` — "risk =
  nearby recent cases".
- **Strengths / Weaknesses.** The honest spatial null every real model must beat; full
  519-zone support; no mobility structure, no calibration (poor log-score).

### 7.4 Stochastic spatial SEIR (C6)  (`08_stochastic_seir.R`)
- **Description.** A full metapopulation compartmental simulation; invasion read from
  *new infections* (not ascertained cases, which lag).
- **Equation.** Per-zone S→E→I→R with `σ = 1/incubation` (5.3 d) and `γ = 1/infectious`
  (mean infectious period ≈ 8 d — a fixed modelling choice, *not* derived from the cCFR;
  the cCFR instead **partitions** the infectious outflow into death `γ_D = cCFR·γ` and
  recovery `γ_R = (1−cCFR)·γ`), between-zone force of infection via the mobility
  matrix; stochastic (Gillespie/τ-leap) simulation.
- **Strengths / Weaknesses.** Mechanistic, propagates stochastic uncertainty; many
  assumed parameters for a ~15-event task — over-engineered, noisiest, least
  calibrated here. **Extension:** fit rather than assume the mixing/infectious
  parameters if data allow.

### 7.5 Models removed after honest attempts (kept out of all figures/rankings)
zone-specific-R renewal S2 (identical to national-R by construction),
free-coefficient covariate GLM S3 (coefficients explode — the problem the Bayesian
suite solves), rare-event logistic S4 (a double-logistic bug pinned it near 0.5),
INLA spatio-temporal C4 (predictive explosion), ZINB C5, and the
KNN/population/persistence baselines B2/B5/B6 (no invasion signal, or never
cross-validated). Guiding principle: **only show models that demonstrably work.**

---

## 8. Risk scores & preparedness priority (at-risk zones only)  (`compute_risk_scores()`, `20_forecast_detail.R`)

All are computed **only** for zones with zero confirmed cases to date (`NA`
otherwise), from the featured model's expected first-case count `μ_i` and calibrated
probability `p_i`:
1. **Absolute invasion probability** `p_case = 1 − exp(−μ_i)` (ascertainment-agnostic,
   §0/§3.3).
2. **Relative risk within each province of interest** `RR_p = μ_i / mean(μ over at-risk
   zones in province p)` (+ share, rank) for **Ituri, Nord-Kivu and Haut-Uele**
   (`PROVINCES_OF_INTEREST`; Task 7). A zone can top its province yet be modest
   nationally, and vice-versa. `compute_risk_scores()` generalises the old Ituri-only
   score with a per-province loop (`.prov_suffix`: `rr_ituri`, `rr_nordkivu`,
   `rr_hautuele`).
3. **Relative risk nationwide** `RR_nat = μ_i / mean(μ over all at-risk)` (+ share, rank).
4. **Infection vs case as a BAND** `p_inf = 1 − exp(−μ_i/ρ)` reported as
   `[p_infection_lo, p_infection_hi]` across ρ ∈ {0.30, 0.45, 0.60} rather than one
   fixed ρ (Task 8); `p_inf − p_case` flags **silent-introduction** risk with honest
   ascertainment uncertainty.
5. **0–1 relative-risk index** `rr01 = p_i / max(p_i)` within a geography — national
   *and* within each province of interest (`rr01_nat`, `rr01_ituri`, `rr01_nordkivu`,
   `rr01_hautuele`) — a bounded, comparable score for maps and tables.

**Preparedness priority (vulnerability- & capacity-adjusted).**
- **Vulnerability & capacity index** `V ∈ [0,1]` = equal-weight mean of **four**
  **percentile-ranked** pillars (higher = worse-off; Task 2 adds the fourth):
  * **surveillance gap** `1 − rank(healthsite_density)` — detection reach;
  * **healthcare-load gap** `1 − rank(healthsite_count/pop)` — facilities per capita;
  * **healthcare-ACCESS gap** `rank(travel_time to nearest health-facility zone)` — road
    travel time (OSRM) to the nearest zone with any facility (0 for zones with their
    own), a distinct axis from density/per-capita (`compute_vulnerability_index(..,
    osrm_mat)`; the axis is dropped, and V averages the remaining three, if no OSRM
    matrix is supplied);
  * **social vulnerability** `rank(CCVI)`.
- **Priority** `= rr01_nat × V`, rescaled to [0,1]. **Multiplicative**: a zone must be
  *both* at material invasion risk *and* vulnerable/under-resourced to rank high.
- **Bivariate prob × vulnerability choropleth (Task 1).** Alongside the single priority
  number, `plot_prob_vuln_choropleth()` maps the two dimensions the priority multiplies
  on **one** 4×4 bivariate colour grid (invasion-prob quartile × vulnerability quartile),
  nationally and per province — so the operational "hot corner" (likely-invaded *and*
  under-resourced) is legible without collapsing to a scalar.
- **Per-province masked invasion maps (Task 7).** `plot_province_risk_maps()` renders a
  zoomed invasion choropleth for each province of interest (Ituri, Nord-Kivu, Haut-Uele)
  plus a national panel.
- **Strengths / Weaknesses.** Relative and 0–1 scores are robust to residual absolute
  miscalibration and directly operational; the ascertainment band is now explicit; V's
  four pillars are equal-weighted by choice (re-weightable); the surveillance pillar is
  a facility-density proxy (dedicated testing data cover only ~19/519 zones).

---

## 9. Evaluation (leave-future-out, invasion-focused)  (`16_invasion_eval.R`, `19_spacetime_eval.R`, `20_forecast_detail.R`)

**Intuition.** Stand at each past week, forecast the next 1–2 weeks with only the
data available then, and score against what actually happened — with metrics that
make sense for a **rare binary** event.

### 9.1 What an LFO fold *is* (Task 1c)
A **fold** = one historical **cutoff date**. Leave-**future**-out means: train on
`week ≤ cutoff`, forecast `(cutoff, cutoff + 7h]`, score against realised first
cases. The cutoff **rolls forward** week by week, so the training (fit) window
**expands** and the evaluation (prediction) window always sits strictly in its
future — this is the temporal analogue of cross-validation folds, and it is what
guarantees no look-ahead. A cutoff qualifies once its shortest-horizon evaluation
week is reporting-complete; each fold scores **only** horizons that are themselves
complete (a per-horizon completeness guard recovers recent cutoffs where h=1 is
evaluable but h=2 is not). Each fold's training data is nowcast **as of its own
cutoff + 7 d** (rolling data-vintage). ~5 folds survive for h=1. `run_invasion_lfo()`
drives the loop; `make_bayes_lfo_model()` slots the Bayesian models into the same
loop.

### 9.1a Analysis-date-anchored daily backtest (`run_invasion_lfo_daily`)
The weekly LFO in §9.1 issues at the **clean week boundary** against ISO-week
targets. The **live** forecast instead trains on the *partial* current week, issues
**mid-week**, and scores rolling windows measured from the analysis date (§daily
re-issue). `run_invasion_lfo_daily()` validates exactly that. For each historical
cutoff week *W* it re-issues on `W + {0, 3, 6}` days (Mon/Thu/Sun by default), with:
- **partial issue-week in training** — the model's last training week is exactly
  *W*, nowcast **as of the issue date `d`**, matching the live state;
- **analysis-date-anchored windows** — the weekly forecast is re-anchored by
  `anchor_windows_from_analysis_date()` to `(d, d+7]` and `(d, d+14]`, identical to
  the live product;
- **onset-defined outcome at daily resolution** — a zone "invades within *N* days"
  iff its first-ever confirmed **onset** (`date_index`) falls in `(d, d+N]`; the
  at-risk set is zones with no confirmed onset on/before `d`. A window is scored only
  once it is old enough (`analysis_date − (d+N) ≥ min_eval_age_days`);
- **leakage-free training counts** — TRAINING is re-aggregated from the line list
  censored to each issue date's sample-observation date (`reaggregate_asof()`), never
  the final weekly counts; a `week_spine` zero-fills the cutoff week so the weekly
  horizon alignment never slips. The OUTCOME uses final onset dates (what truly
  happened), which is correct.

Cost note: epinowcast is **not** refit per issue-fold (intractable across
issues × folds × models); the deterministic per-fold correction (§2b) is used for
scoring, while the *live* forecast refits epinowcast once per day. Results pool to
`lfo_daily_results.rds` with `issue_date`, `cutoff`, and `window_days`/`lead_days`
retained so skill reads against intra-week lead time. Toggle with `RUN_DAILY_LFO`.
Pooling caveat: a fixed event recurs across issue offsets, so an aggregate metric
should be computed per offset or cluster-bootstrapped by `cutoff`.

### 9.2 The outcome, the at-risk set, and how β is pooled
- **Outcome.** A zone is an "invasion" at horizon *h* iff its **first-ever** confirmed
  case falls in `(cutoff, cutoff + 7h]` — cumulative, matched to the model's
  cumulative probability, counted once (`is_new_invasion`).
- **At-risk only, pooled across folds.** Because events are so few, metrics are
  computed on the *pooled* at-risk zone-weeks (not per-fold-then-averaged, which is
  degenerate), using the models' **raw** probabilities.
- **β pooling.** β is pooled over week-transitions (each week's arrivals are Poisson
  with rate β·Λ, so the cumulative over *h* weeks is exact), optionally
  completeness-weighted.

### 9.3 How uncertainty is generated (Task 1c)
Three distinct, honestly-labelled sources:
1. **Forecast uncertainty** = pre-specified **ensemble member spread** (min/median/max
   across mobility/GT/observation members) — the frequentist band (§5.6).
2. **Posterior uncertainty** = the **Bayesian 90 % credible interval** on β and on each
   invasion probability (§6).
3. **Sampling uncertainty on the metrics** = a **cluster bootstrap resampling zones**
   (not zone-weeks), so the CI respects the fact that a zone's weeks are correlated.
All three uncertainty **bands** — the ensemble spread, the Bayesian credible interval, and
the metric cluster-bootstrap — use the q5/q95 = **90 %** / q20/q80 = **60 %** convention
(matching the forecast-count quantiles and R(t) bands); these bands are never labelled 50 %.
The **one** deliberate exception is the frequentist **covariate-association screen** (§5.5)
and the params-/β-over-folds traces, where each regression hazard ratio is reported with a
**standard 95 % Wald confidence interval** (the conventional coefficient interval), clearly
labelled as such.

### 9.4 Metrics (and why AUC-PR is reported as a *skill ratio*) (Task 1b)
- **Discrimination.** AUC-ROC and **Average Precision** (the correct step-wise AP, not
  the interpolated AUC-PR that flatters constant predictors at low prevalence),
  reported as **AUC-PR skill = AP / base_rate**. *Why divide by the base rate?*
  Average Precision equals the base rate for a random ranker, so raw AP is
  incomparable across horizons with different prevalence. Dividing by the per-horizon
  base rate (invasion prevalence among at-risk zone-weeks, under ~1 % at h1 and ~1.5 % at h2; data-dependent)
  makes **1 = no skill** and, e.g., 5 = "five times better than chance at
  concentrating invasions at the top of the list", comparable across horizons.
- **Probabilistic accuracy.** log-score, Brier skill, calibration-in-the-large, and a
  reliability diagram (`plot_reliability()`, `17_invasion_viz.R`).
- **Operational.** Precision/Recall/hit-rate@K and the **mean rank of invaded zones**
  (all at-risk zones ranked by p descending, **ties share the average rank** so
  hundreds of zeros do not get arbitrary positions; the rank of *every* invaded zone
  is averaged within a fold, then across folds).
- **An intuitive headline metric for the imbalanced task (Task 6).** For a
  ~0.7 %-prevalence event, raw accuracy and even AUC-ROC are misleading (both are high
  for a do-nothing model). The pipeline therefore reports two honest, plain-language
  summaries:
  * **Detection-vs-budget curve** (`compute_detection_curve()`,
    `detection_vs_budget_h{1,2}.pdf`) — *"if we actively monitor the top-K highest-risk
    zones each week, what fraction of the zones that actually get invaded do we catch?"*
    i.e. **sensitivity (recall) at a fixed weekly alert budget**, pooled over folds,
    with the matching precision. A good model reaches a high catch-rate at small K; this
    is the single most decision-legible picture of skill.
  * **Balanced skill at the Youden-optimal threshold** (`invasion_balance_metrics()`,
    `invasion_balanced_skill.csv`) — **balanced accuracy** = ½(sensitivity+specificity),
    **F1**, and **Matthews correlation (MCC)**, which — unlike accuracy/AUC-ROC — are not
    inflated by the ~99.3 % of zone-weeks with no invasion, giving one honest number that
    explicitly balances sensitivity and specificity.
- **WIS (count score) is deliberately NOT a selector.** Over 519 mostly-zero zones it
  rewards "predict zero everywhere", crowning a no-signal persistence baseline; it is
  computed **exactly per Bracher et al. 2021** (interval score with weights
  `w_k = α_k/2`, `w_0 = 1/2`; levels q5/q95 = 90 %, q20/q80 = 60 %) and kept only as a
  secondary count diagnostic.

### 9.5 Evaluation over time (Task 4)
Beyond the pooled score, the pipeline now traces performance **across folds**:
- **AUC-PR skill & hit@K over time** — `plot_skill_over_time()` plots per-fold skill
  vs forecast date (does discrimination strengthen as the outbreak matures?).
- **Predicted vs observed over folds** — `plot_predobs_over_folds()` overlays, per
  fold, the **mean predicted probability**, the **observed invasion fraction**, and
  the **mean predicted probability at the zones that actually invaded**: a direct,
  fold-by-fold calibration + discrimination read (well-calibrated ⇒ mean predicted ≈
  observed fraction; discriminating ⇒ the invaded-zone line sits above the mean).
- **Parameter estimates over time** — `plot_params_over_time()` re-fits the
  identifiable exogenous covariate cloglog GLM **at each fold cutoff** (causally,
  training only on `week ≤ cutoff`) and plots each hazard ratio with its 95 % CI over
  time, showing whether/when each driver stabilises.
- **Import coefficient β₀ and reporting completeness over folds (Task 4)** —
  `compute_beta_over_folds()`/`plot_beta_over_folds()` refit the intercept-only cloglog
  (with the log Λ offset) at each cutoff to trace the calibrated import→invasion
  coefficient **β₀ = exp(intercept)** (with 95 % CI), alongside the mean **reporting
  completeness** (observed ÷ nowcast-corrected counts) applied that fold — so one can
  see whether the calibration and the truncation correction are stable as the outbreak
  accrues.
- **Reporting rate across space (Task 4)** — `plot_reporting_rate_map()` maps the
  per-zone relative reporting-rate proxy the model uses to up-weight under-ascertained
  source zones (health-site density, geometric-mean-normalised, bounded [0.25, 4]);
  red = likely under-reporting (true import pressure exceeds observed counts).
- **Lead time (anticipation)** — for every newly affected zone, how many weeks *before*
  its first case it entered the top-K risk list; **spatial error map** —
  hit/miss/false-alarm/correct-negative per zone; **space-time risk tile** — top zones
  × cutoff with realised first cases marked.

**Strengths.** Honest, task-matched, leakage-free, with explicit uncertainty and event
counts on every figure. **Weaknesses.** ~8–17 events is a genuinely tiny sample — CIs
are wide and rankings indicative; ground truth is the current snapshot, so the newest
weeks are excluded rather than corrected; per-fold R(t) is a fast renewal ratio.

### 9.6 How the featured models are selected  (`08c`)
Every model is scored in the same LFO on at-risk zones, pooled across folds, and
ranked by a **calibration-aware composite** (sum of within-horizon ranks of AUC-PR
skill + mean-rank-of-invaded + log-score, pooled across both forecast horizons) — never a single
metric, never count-WIS.
Two models are featured: the **best over all families** and the **best renewal model**
(so the featured mechanistic product is always renewal-based). With ~8–17 events the
leaders are within CV noise, which is why the ensemble spread and Bayesian posterior
are reported alongside the single best model.

---

## 10. What to trust, and what not to

**Trust:** the *ranking* of at-risk zones, the *relative-risk* and 0–1 scores, the
vulnerability-adjusted **priority ordering** (a resource-targeting lens, not a
probability), that mobility-informed renewal models clearly beat distance/persistence
baselines on discrimination, the **direction** of the identifiable effects (frontier
proximity `d_min` lowers risk, deprivation raises it), and the masked risk maps.

**Treat as provisional:** the exact *absolute* probabilities (calibrated to a handful
of events; NB is best but still ~2–3× over the realised rate), the *magnitude* of the
priority score and its pillar weighting, the alert/positivity covariate effects (not
identifiable on the full series — the frequentist featured model runs intercept-only
live; the Bayesian suite regularises them but with wide posteriors), and fine-grained
differences between good models. As more of the outbreak accrues, the same
cross-validation sharpens automatically.

---

## 11. Outputs: figures, tables & reports  (`spatiotemporal/outputs/`)

Affected zones are masked (grey / NA) throughout; figures are colour-blind-safe,
label zones by **province / health area**, and current-forecast figures print their
**fit + prediction date windows** (§0).

**Forecast & risk (featured model).** `invasion_risk_map_h{1,2}.pdf`,
`risk_scores_bars_h1.pdf`, `forecast_uncertainty_h1.pdf`,
`spacetime_forecast_uncertainty_h1.pdf`, `forecast_ensemble_panel_h1.pdf`,
`risk_table_{ituri,national}.csv`.
**Priority (§8).** `priority_scatter_h1.pdf`, `priority_bars_h1.pdf`,
`priority_map_h1.pdf`, `priority_table.csv`.
**Bayesian (§6).** `bayes_parameters.pdf`, `bayes_posterior_densities.pdf`,
`bayes_gt_posterior.pdf`/`.csv`, `bayes_nowcast_sensitivity.pdf`/`.csv`,
`bayes_invasion_uncertainty_h1.pdf`,
`bayes_stacked_invasion_uncertainty_h1.pdf`, `bayes_stacking_weights.pdf`,
`bayes_discrimination_summary_h{1,2}.pdf`, `bayes_parameters.csv`.
**Parameters & selection (§5.5, §9.6).** `rt_national.pdf`,
`best_model_parameters.pdf`, `covariate_associations.rds`, `model_specification.md`,
`invasion_report.md`.
**Evaluation over time (§9.5).** `skill_over_time_*.pdf`,
`predicted_vs_observed_over_folds_h1.pdf`, `params_over_time.pdf`, `lead_time.pdf`,
`spatial_error_map_h{1,2}.pdf`, `spacetime_risk_h1.pdf`, plus `invasion_evaluation.csv`,
`reliability_h1.pdf`, `discrimination_summary_h{1,2}.pdf`, and the
`skill_over_time.csv` / `lead_time.csv` / `zone_spatial_error.csv` tables.

---

## 12. Model catalogue (one-line reference)

The **featured** model is **selected by the cross-validation composite each run** — AUC-PR skill +
mean rank-of-truth + log-score, summed within-horizon ranks pooled across both forecast horizons,
over the cross-validated single
Bayesian models (and, when `RUN_FREQUENTIST_MODELS = TRUE`, the frequentist renewal by the same
composite) — NOT the loo-stacking weight, and NOT hardcoded. In the current (Bayesian-only) run the
featured model is the radiation-composite renewal **`Bayes-M10-med`** (which also tops the AUC-PR-
skill leaderboard); the loo-stacking weights instead define the separate loo-stacked ensemble (its
weights concentrate on the `-dist` covariate models, e.g. `Bayes-M10-dist-geo` 0.27). The exact
choices, with all parameters, are in `invasion_report.md` / `model_specification.md`. M8 below is the
recommended **default** kernel, not necessarily the headline model.

| Model | File → function | Core equation | Role |
|---|---|---|---|
| Renewal M8 (base) | `15_workhorse.R` → `forecast_workhorse` | `p = 1 − e^{−β Λ}`, `Λ = t(W)·Ỹ` | default kernel |
| Renewal M8-cov | `fit_import_model` | `β_i = e^{β₀+Σγx}` | covariate-heterogeneity variant (exogenous covariates only) |
| Renewal M11 (inward) | `build_inward_contact_matrix` | `Λ = M·Ỹ`, `M = P diag(1/N_eff) Pᵀ` | manuscript-motivated |
| Bayesian renewal | `21_bayesian_renewal.R` → `fit_bayes_suite` | cloglog + offset, `γ~N(0,1)` | regularised + posterior |
| Ensemble (mean/median) | `18_ensemble.R` | `p̄ = mean/median p_m` | robust combination |
| hhh4 | `07_hhh4_model.R` | endemic+AR+neighbour | comparator |
| Gravity B4 / Distance B1 | `05_baseline_models.R` | gravity / `d^{−α}` | baselines |
| Stochastic SEIR C6 | `08_stochastic_seir.R` | metapop S→E→I→R | mechanistic comparator |

---

## 13. Where sitrep (and other unused) data could still help

As of the sitrep reconciliation (§1), the sitrep's **cumulative confirmed count** *is* now used —
as a floor that tops the line-list confirmed set up to the official count. Its **other** streams
remain outside the fits. Highest-value additions, in order:
1. **Suspected cases as a leading indicator** (`new_suspected_cases`, reported before
   lab confirmation — a timelier invasion early-warning than the line-list alert).
2. **Contacts traced / isolated → surveillance pillar** (a *direct* capacity measure
   that would replace the facility-density proxy for the zones they cover).
3. **Hospitalisations → healthcare pillar** (burden/utilisation signal).
4. **Independent case series for nowcasting** (`new_confirmed_cases` cross-checking the
   line-list truncation correction).
Caveat: current coverage is ~6–33 mostly already-affected zones over ~2–3 weeks, so
these help most for *response capacity of the source cluster* and *near-frontier early
warning*, less as covariates across the 519-zone at-risk grid.
