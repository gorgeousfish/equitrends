# equitrends

Equivalence Tests for Pre-Trends in DiD Estimation

[![Stata](https://img.shields.io/badge/Stata-16%2B-blue)](https://www.stata.com/)
![Version](https://img.shields.io/badge/version-0.1.0-informational)
![License](https://img.shields.io/badge/license-AGPL--3.0-blue)

![Reversing the Burden of Proof](images/image.png)

EQUITRENDS is a Stata package implementing equivalence tests for pre-trends in Difference-in-Differences (DiD) designs, based on [Dette & Schumann (2024)](https://doi.org/10.1080/07350015.2024.2308121), published in the *Journal of Business & Economic Statistics*.

## Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Getting help](#getting-help)
- [Data requirements](#data-requirements)
- [Quick start](#quick-start)
- [Command reference](#command-reference)
- [Test selection guide](#test-selection-guide)
- [Methodology](#methodology)
- [RMS test alpha restriction](#rms-test-alpha-restriction)
- [Stored results](#stored-results)
- [Visualization](#visualization)
- [Examples](#examples)
- [Citation](#citation)
- [Authors](#authors)
- [License](#license)
- [Future roadmap](#future-roadmap)
- [Related packages](#related-packages)
- [Recommended resources](#recommended-resources)
- [Support](#support)

## Overview

Standard pre-trend tests in DiD designs test the null hypothesis of *exact* parallel trends (H₀: β = 0). This approach suffers from fundamental limitations:

1. **Failure to reject ≠ Evidence in favor**: Low statistical power may prevent detection of actual violations. Failing to reject does not support the parallel trends assumption.
2. **Conditional bias amplification**: Roth (2022) demonstrates that conditioning on passing traditional pre-tests can *amplify* DiD bias when violations exist.
3. **No explicit threshold**: Traditional tests provide no framework for determining what constitutes a "negligible" deviation from parallel trends.

EQUITRENDS implements *equivalence* tests that **reverse the burden of proof**: the null hypothesis is that deviations are *large* (H₀: ‖β‖ ≥ threshold), and rejection provides statistical evidence that deviations are *small*. This approach:

- Requires explicit justification of the equivalence threshold
- Controls Type I error (falsely concluding equivalence)
- Increases statistical power with sample size
- Allows researchers to quantify the smallest threshold at which equivalence holds

### Key features

- **Three equivalence hypotheses**: maximum, mean, and RMS (Dette & Schumann, 2024, Section 3.1)
- **Minimum equivalence threshold**: compute the smallest threshold at which equivalence can be concluded
- **Multiple inference methods for the maximum test**: IU (analytical), spherical bootstrap, and wild bootstrap
- **Visualization**: coefficient plots with equivalence bounds (`equivtest_plot`)

## Requirements

- **Stata 16.0** or higher
- **Python (Stata Python integration)** with **numpy** and **scipy** for bootstrap-based maximum tests:
  - `equivtest ..., type(max) method(boot)`
  - `equivtest ..., type(max) method(wild)`

**Note**: The IU method for `type(max)` and the `type(mean)` / `type(rms)` tests do not require Python. Bootstrap methods (`method(boot)` and `method(wild)`) require Python with **numpy** and **scipy**.

## Installation

### From GitHub

```stata
* Install the package
net install equitrends, from("https://raw.githubusercontent.com/gorgeousfish/equitrends/main") replace
```

### Local installation

```stata
* Install the package
net install equitrends, from("/path/to/equitrends-main") replace
```

### Loading example data

```stata
* Load the bundled example dataset
equitrends_data, clear
```

The dataset is bundled with the package and loaded from the local installation directory.

## Getting help

After installation, access the built-in Stata documentation:

```stata
help equitrends       // Package overview and command list
help equivtest        // Main unified testing interface
help equivtest_plot   // Visualization options
help maxequivtest     // Maximum test details
help meanequivtest    // Mean test details
help rmsequivtest     // RMS test details
help equivsim         // Monte Carlo simulation
help equitrends_data  // Load example datasets
```

## Data requirements

**Before using this package, ensure your data meets the following requirements:**

| Requirement                             | Description                                                                                                                                              |
| :-------------------------------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Panel structure**               | Panel data with individual (`id`) and time (`time`) identifiers. **Both balanced and unbalanced panels are supported.**                        |
| **Minimum pre-treatment periods** | At least one pre-treatment period (*T* ≥ 1). More periods increase statistical power.                                                                 |
| **Treatment group indicator**     | Binary variable coded as 0 (control) or 1 (treated).                                                                                                     |
| **Block adoption design**         | All treated units must receive treatment at the same time. Staggered adoption requires cohort-specific analysis (see Dette & Schumann, 2024, Section 5). |
| **Complete time-group cells**     | Each time period must contain observations in both treatment and control groups.                                                                         |

### Unbalanced panels

The package automatically detects and handles unbalanced panels (where individuals have different numbers of observed time periods). When an unbalanced panel is detected:

- `e(is_balanced)` returns 0
- `e(T_min)` and `e(T_max)` store the minimum and maximum number of periods per individual
- Output displays "Panel type: Unbalanced" with the period range

No special syntax is required—simply run the command as usual:

```stata
xtset id time  // Shows "unbalanced" if panel is unbalanced
equivtest y, type(max) id(id) group(treat) time(time) pretreatment(1 2 3 4 5) baseperiod(5)
```

### Data completeness

When some time periods lack observations for either group, the placebo regression cannot be estimated correctly. In such cases:

- The command will issue an error or produce missing values
- Check that each time-group cell contains at least one observation

**Common causes:**

- Missing values in treatment, outcome, or identifier variables
- Sample restrictions that eliminate entire time periods for one group

**Solution:** Ensure at least one observation per time-treatment cell, or restrict analysis to time periods with complete coverage. Use `xtset` to verify panel structure before running equivalence tests.

## Quick start

### Empirical example: Di Tella & Schargrodsky (2004)

This example replicates the empirical application in Dette & Schumann (2024, Section 7) using the Buenos Aires crime data from Di Tella & Schargrodsky (2004). The dataset contains monthly car theft counts for 876 Buenos Aires city blocks (April–December 1994), of which 37 blocks received police protection after a July terrorist attack.

```stata
* Load and prepare data
equitrends_data, clear
drop if mes == 72 | mes == 73                        // Remove post-treatment months
drop if mes > 7                                      // Keep pre-treatment only
rename (observ totrob mes) (ID Y period)
gen G = (distanci == 0)                              // Treatment indicator
xtset ID period

* Maximum test (IU method, cluster-robust SE)
equivtest Y, type(max) method(iu) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) vce(cluster) cluster(ID)

* Mean test
equivtest Y, type(mean) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) vce(cluster) cluster(ID)

* RMS test
equivtest Y, type(rms) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) seed(2024)

* Visualize results
equivtest_plot, ci
```

### Bootstrap methods (maximum test only)

```stata
* Spherical bootstrap (assumes spherical errors; Theorem 1)
equivtest Y, type(max) method(boot) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) nboot(1000) seed(12345)

* Wild bootstrap (recommended for non-spherical errors; Remark 1(c))
equivtest Y, type(max) method(wild) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) nboot(1000) seed(12345)
```

### Control variables (conditional parallel trends)

The `x()` option includes additional control variables in the TWFE placebo regression (Dette & Schumann, 2024, Section 5). Time-invariant covariates are absorbed by the individual fixed effects during double demeaning. For the full conditional PTA specification (Eq. 5.5), construct covariate-by-time interactions in your data and pass them via `x()`.

```stata
equivtest Y, type(max) method(iu) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) x(edpub estserv banco) ///
    vce(cluster) cluster(ID)
```

## Command reference

| Command             | Description                                         |
| :------------------ | :-------------------------------------------------- |
| `equivtest`       | Unified interface for all three tests (recommended) |
| `maxequivtest`    | Maximum absolute coefficient test (IU/Boot/Wild)    |
| `meanequivtest`   | Mean coefficient test                               |
| `rmsequivtest`    | Root mean square (RMS) test                         |
| `equivtest_plot`  | Coefficient plot with equivalence bounds            |
| `equivsim`        | Monte Carlo simulation for power analysis           |
| `equitrends_data` | Load bundled example datasets                       |

### Unified syntax

```stata
equivtest depvar, type(max|mean|rms) id(varname) group(varname) time(varname) [options]
```

### Core options

| Option                    | Description                                                    |
| :------------------------ | :------------------------------------------------------------- |
| `type(max/mean/rms)`    | Test type (required)                                           |
| `id(varname)`           | Panel identifier (required)                                    |
| `group(varname)`        | Treatment group indicator 0/1 (required); also accepts `g()` |
| `time(varname)`         | Time period variable (required); also accepts `period()`     |
| `threshold(#)`          | Equivalence threshold; omit to compute minimum threshold       |
| `alpha(#)`              | Significance level; default 0.05                               |
| `pretreatment(numlist)` | Pre-treatment periods to include                               |
| `baseperiod(#)`         | Base period for placebo construction                           |
| `x(varlist)`            | Control variables                                              |

Options specific to `type(max)`:

| Option                   | Description                                  |
| :----------------------- | :------------------------------------------- |
| `method(iu/boot/wild)` | Inference method; default `iu`             |
| `nboot(#)`             | Bootstrap replications; default 1000         |
| `seed(#)`              | Random seed for bootstrap                    |
| `nodots`               | Suppress bootstrap progress display          |

Options specific to `type(rms)`:

| Option          | Description                                                 |
| :-------------- | :---------------------------------------------------------- |
| `nolambda(#)` | Number of subsamples used for self-normalization; default 5 |
| `seed(#)`     | Random seed for subsampling                                 |

Robust SE options (available for IU/mean; not for bootstrap or RMS):

| Option               | Description                                              |
| :------------------- | :------------------------------------------------------- |
| `vce(vcetype)`     | Variance estimator; see table below                      |
| `cluster(varname)` | Cluster variable (required for cluster-robust VCE types) |

**Variance estimator types (`vcetype`):**

| Type                | Description                                                                        |
| :------------------ | :--------------------------------------------------------------------------------- |
| `ols`             | Homoskedastic OLS variance (default)                                               |
| `robust`/`hc1`  | HC1 heteroskedasticity-robust (White, 1980)                                        |
| `hc2`             | HC2 leverage-adjusted (MacKinnon & White, 1985); better finite-sample properties   |
| `hc3`             | HC3 more conservative leverage adjustment (Davidson & MacKinnon, 1993)             |
| `hac`             | Arellano (1987) HAC estimator for panel data                                   |
| `cluster`/`cr0` | CR0 cluster-robust without small-sample adjustment; requires `cluster()`         |
| `cr1`             | CR1 cluster-robust with G/(G-1) adjustment (Stata default); requires `cluster()` |
| `hc1_cluster`     | HC1 cluster-robust with small-sample adjustment; requires `cluster()`            |

## Test selection guide

EQUITRENDS offers three equivalence tests with different properties. Choose based on your research context:

| Feature                    | Maximum test               | Mean test                   | RMS test                   |
| :------------------------- | :------------------------- | :-------------------------- | :------------------------- |
| **Hypothesis**       | max\|βₜ\| < δ           | \|β̄\| < τ               | β_RMS < ζ                |
| **Measures**         | Largest single violation   | Average violation           | Root mean square           |
| **Cancellation**     | No                         | Yes (opposing signs cancel) | No                         |
| **Sensitivity**      | Any single large deviation | Systematic directional bias | Balanced across deviations |
| **Conservativeness** | Most conservative          | Least conservative          | Moderate                   |

### Recommendations

1. **Maximum test (`type(max)`)**: Start here as the default, conservative choice.

   - Detects any single large violation
   - Use `method(iu)` for analytical inference or `method(wild)` for non-spherical errors
   - Recommended when you want to rule out *any* substantial pre-trend violation
2. **Mean test (`type(mean)`)**: Use when violations are expected to be monotone (same sign).

   - More powerful when deviations are directionally consistent
   - **Caution**: Opposing violations may cancel out, leading to false equivalence
3. **RMS test (`type(rms)`)**: Use as a general-purpose alternative.

   - Balances sensitivity across all placebo coefficients
   - No cancellation problem
   - Self-normalized (no variance estimation required)

### Interpreting minimum thresholds

When `threshold()` is omitted, EQUITRENDS reports the smallest equivalence threshold (δ\*, τ\*, or ζ\*) at which equivalence can be concluded at the specified significance level. Compare this to your estimated treatment effect:

- **δ\* << estimated ATT**: Strong evidence for negligible pre-trends
- **δ\* ≈ estimated ATT**: Pre-trend violations may explain the treatment effect
- **δ\* >> estimated ATT**: Insufficient evidence for parallel trends; consider alternative designs

## Methodology

Let $\beta = (\beta_1,\ldots,\beta_T)'$ denote the vector of placebo (pre-treatment) coefficients from the TWFE placebo regression (Dette & Schumann, 2024, Eq. (2.5)). EQUITRENDS implements three equivalence hypotheses (Section 3.1):

1. **Maximum deviation (Eq. (3.1))**:

$$
H_0: \|\beta\|_{\infty} \ge \delta \quad \text{vs.} \quad H_1: \|\beta\|_{\infty} < \delta, \qquad \|\beta\|_{\infty}=\max_{l\in\{1,\ldots,T\}}|\beta_l|
$$

2. **Mean deviation (Eq. (3.2))**:

$$
\bar{\beta}=\frac{1}{T}\sum_{l=1}^{T}\beta_l, \qquad H_0: |\bar{\beta}| \ge \tau \quad \text{vs.} \quad H_1: |\bar{\beta}| < \tau
$$

3. **RMS deviation (Eq. (3.3))**:

$$
\beta_{\mathrm{RMS}}=\sqrt{\frac{1}{T}\sum_{l=1}^{T}\beta_l^2}, \qquad H_0: \beta_{\mathrm{RMS}} \ge \zeta \quad \text{vs.} \quad H_1: \beta_{\mathrm{RMS}} < \zeta
$$


### Inference methods for the maximum test

- **IU (Intersection-Union, analytical)**: For each placebo coefficient *t* = 1, ..., *T*, the test rejects H0 iff all |beta_t| < Q(alpha), where Q denotes the alpha-quantile of the folded normal distribution with mean delta and variance sigma_tt/n (Dette & Schumann, 2024, Eq. (4.4)). Computationally attractive but conservative for large *T*. Folded normal CDF/quantiles are implemented in Mata.
- **Bootstrap** (`method(boot)`): Generates bootstrap samples under the constraint on beta using constrained OLS, then computes the empirical alpha-quantile as the critical value (Dette & Schumann, 2024, Theorem 1). Assumes spherical errors. More powerful than IU for *T* > 1. Requires Python with numpy and scipy.
- **Wild bootstrap** (`method(wild)`): Replaces i.i.d. bootstrap errors with Rademacher-weighted residuals, making the test robust to heteroskedasticity and serial correlation (Dette & Schumann, 2024, Remark 1(c)). Recommended for non-spherical errors. Requires Python with numpy and scipy.

### Inference for the mean test

The mean test rejects H0 whenever the absolute sample mean of placebo coefficients falls below the alpha-quantile of the folded normal with mean tau and variance 1'Sigma1/(nT^2) (Dette & Schumann, 2024, Eq. (4.12)).

### Inference for the RMS test

The RMS test uses a self-normalized statistic based on subsampling (Dette & Schumann, 2024, Theorem 2). It rejects H0 whenever beta_RMS^2 < zeta^2 + Q_W(alpha) * V_n, where Q_W(alpha) is the alpha-quantile of the limiting distribution (a functional of Brownian motion) and V_n is computed from subsample estimates (Eq. (4.18)). This test is pivotal and does not require variance estimation.

## RMS test alpha restriction

The RMS test supports only:

$$
\alpha \in \{0.01, 0.025, 0.05, 0.1, 0.2\}
$$

This reflects the implementation based on critical values for the limiting distribution in Dette & Schumann (2024, Theorems 2-3).

## Stored results

`equivtest` stores results in `e()` (see also `src/sthlp/equivtest.sthlp`).

### Scalars

| Result                     | Description                                                                        |
| :------------------------- | :--------------------------------------------------------------------------------- |
| `e(N)`                   | Number of observations                                                             |
| `e(N_g)`                 | Number of individuals (panels)                                                     |
| `e(no_placebos)`         | Number of placebo coefficients (*T*)                                             |
| `e(alpha)`               | Significance level used                                                            |
| `e(base_period)`         | Base period for placebo construction                                               |
| `e(is_balanced)`         | 1 if balanced panel, 0 otherwise                                                   |
| `e(threshold_specified)` | 1 if `threshold()` was specified, 0 otherwise                                    |
| `e(threshold)`           | Equivalence threshold value (if specified)                                         |
| `e(reject)`              | 1 if H₀ rejected, 0 otherwise (if threshold specified)                            |
| `e(min_threshold)`       | Minimum threshold δ\*/τ\*/ζ\* at which equivalence holds (if threshold omitted) |

Type-specific scalars for `type(max)`:

| Result                | Description                                                    |
| :-------------------- | :------------------------------------------------------------- |
| `e(max_abs_coef)`   | Maximum absolute placebo coefficient                           |
| `e(nboot)`          | Number of bootstrap replications (`method(boot)` or `wild`)  |
| `e(boot_critical)`  | Bootstrap critical value (`method(boot)` or `wild`)          |

Type-specific scalars for `type(mean)`:

| Result                     | Description                                        |
| :------------------------- | :------------------------------------------------- |
| `e(abs_mean_placebo)`    | Absolute value of mean placebo coefficient          |
| `e(var_mean_placebo)`    | Variance of mean placebo coefficient                |
| `e(se_mean_placebo)`     | Standard error of mean placebo coefficient          |
| `e(p_value)`             | p-value (if threshold specified)                    |
| `e(mean_critical_value)` | Critical value for mean test (if threshold specified) |

Type-specific scalars for `type(rms)`:

| Result                    | Description                                    |
| :------------------------ | :--------------------------------------------- |
| `e(rms_placebo_coefs)`  | Root mean square of placebo coefficients        |
| `e(nolambda)`           | Number of lambda subsamples                     |
| `e(rms_critical_value)` | RMS critical value (if threshold specified)     |

### Macros

| Result                    | Description                                                             |
| :------------------------ | :---------------------------------------------------------------------- |
| `e(cmd)`                | Command name (`equivtest`)                                            |
| `e(cmdline)`            | Full command as typed                                                   |
| `e(type)`               | Test type: `max`, `mean`, or `rms`                                  |
| `e(method)`             | Inference method: `iu`, `boot`, or `wild` (only for `type(max)`)  |
| `e(depvar)`             | Dependent variable name                                                 |
| `e(idvar)`              | Panel identifier variable                                               |
| `e(groupvar)`           | Treatment group variable                                                |
| `e(timevar)`            | Time variable                                                           |
| `e(vce)`                | Variance estimator type (when applicable)                               |
| `e(clustvar)`           | Cluster variable (if cluster-robust VCE)                                |
| `e(preperiods)`         | Specified pre-treatment periods (if `pretreatment()` used)            |
| `e(placebo_coef_names)` | Names of placebo coefficient periods                                    |

### Matrices

| Result                      | Description                                                                          |
| :-------------------------- | :----------------------------------------------------------------------------------- |
| `e(b_placebo)`            | Placebo coefficient vector (*T* × 1); available for IU max and mean tests |
| `e(V_placebo)`            | Placebo variance-covariance matrix (*T* × *T*); available for IU max/mean tests |
| `e(se_placebo)`           | Placebo standard errors (*T* × 1); available for IU max test                      |
| `e(IU_critical_values)`   | Critical values per coefficient (type=max, method=iu, threshold specified)           |
| `e(min_equiv_thresholds)` | Minimum thresholds per coefficient (type=max, method=iu, threshold not specified)    |

## Visualization

`equivtest_plot` creates coefficient plots of placebo coefficients with equivalence bounds. Run it after `equivtest` to visualize results.

### Basic syntax

```stata
equivtest_plot [, options]
```

### Main options

| Option           | Description                                                                                      |
| :--------------- | :----------------------------------------------------------------------------------------------- |
| `threshold(#)` | Equivalence threshold for horizontal lines; defaults to `e(threshold)` or `e(min_threshold)` |
| `ci`           | Display confidence intervals (IU method only)                                                    |
| `level(#)`     | Confidence level; default 95                                                                     |
| `connect`      | Connect points with lines                                                                        |
| `noline`       | Suppress zero reference line                                                                     |
| `nobase`       | Suppress base period reference point                                                             |
| `nothreshold`  | Suppress equivalence threshold lines                                                             |

### Style options

| Option                               | Description                              |
| :----------------------------------- | :--------------------------------------- |
| `msymbol(symbolstyle)`             | Marker symbol; default `O` (circle)    |
| `msize(markersizestyle)`           | Marker size; default `medium`          |
| `mcolor(colorstyle)`               | Marker color; default `navy`           |
| `basemsymbol(symbolstyle)`         | Base period marker symbol; default `S` |
| `basemcolor(colorstyle)`           | Base period marker color                 |
| `threshlcolor(colorstyle)`         | Threshold line color; default `red`    |
| `threshlwidth(linewidthstyle)`     | Threshold line width; default `medium` |
| `threshlpattern(linepatternstyle)` | Threshold line pattern; default `dash` |
| `cilcolor(colorstyle)`             | CI line color; default `navy`          |
| `cilwidth(linewidthstyle)`         | CI line width; default `medium`        |

### Output options

| Option                     | Description                       |
| :------------------------- | :-------------------------------- |
| `title(string)`          | Graph title                       |
| `subtitle(string)`       | Graph subtitle                    |
| `xtitle(string)`         | X-axis title                      |
| `ytitle(string)`         | Y-axis title                      |
| `xlabel(rule_or_values)` | X-axis labels                     |
| `ylabel(rule_or_values)` | Y-axis labels                     |
| `note(string)`           | Graph note                        |
| `scheme(schemename)`     | Graph scheme                      |
| `saving(filename)`       | Save graph to file                |
| `replace`                | Replace existing file when saving |
| `name(windowname)`       | Graph window name                 |

### Examples

```stata
* Basic plot after running equivtest
equivtest_plot

* Plot with 95% confidence intervals
equivtest_plot, ci

* Publication-quality plot
equivtest_plot, ci level(95) title("Pre-trend Analysis") ///
    msymbol(O) mcolor(navy) threshlpattern(dash) scheme(s2color)

* Save graph to file
equivtest_plot, ci saving(pretrend_plot) replace
```

## Examples

### Simulated panel data

This example demonstrates the package using simulated data with parallel trends satisfied.

```stata
* Generate simulated panel data
clear
set seed 12345
set obs 1000
gen id = ceil(_n/10)
gen time = mod(_n-1, 10) + 1
gen treat = (id <= 50)

* Generate outcome: parallel trends hold in pre-treatment (time <= 5)
* Treatment effect of 0.5 in post-treatment (time > 5)
gen y = rnormal() + 0.1*time + treat*(time > 5)*0.5

* Run maximum test with IU method
equivtest y, type(max) method(iu) id(id) group(treat) time(time) ///
    pretreatment(1 2 3 4) baseperiod(5)

* Display minimum threshold
display "Minimum equivalence threshold (delta*): " %6.4f e(min_threshold)
```

### Interpreting results

When running `equivtest` without specifying a `threshold()`, the command reports the minimum equivalence threshold at which the null hypothesis of non-equivalence can be rejected:

```
Equivalence Test for Pre-Trends (Maximum)
──────────────────────────────────────────────────────────────────────────────
Type:             max                    Observations:       1000
Method:           iu                     Individuals:         100
VCE:              ols                    Placebo coefs:         4
──────────────────────────────────────────────────────────────────────────────

Hypothesis Test:
  H0: max|placebo effect| >= delta  (non-equivalence)
  H1: max|placebo effect| <  delta  (equivalence)

Minimum Equivalence Threshold Search
──────────────────────────────────────────────────────────────────────────────
        Period |  Abs. Estimate   Std. Error   Min. Threshold
───────────────+──────────────────────────────────────────────────────────────
    placebo_1  |       0.034521     0.100234       0.199543
    placebo_2  |       0.012345     0.098765       0.174832
    placebo_3  |       0.056789     0.101234       0.223456
    placebo_4  |       0.023456     0.099876       0.187654
──────────────────────────────────────────────────────────────────────────────

Minimum equivalence threshold delta* =   0.2235
```

**Interpretation**: At alpha = 0.05, we can conclude that the maximum absolute placebo coefficient is less than delta\*. If this threshold is small relative to your estimated treatment effect, you have evidence supporting negligible pre-trend violations.

### Testing with a pre-specified threshold

```stata
* Test whether max violation < 0.2
equivtest y, type(max) method(iu) id(id) group(treat) time(time) ///
    pretreatment(1 2 3 4) baseperiod(5) threshold(0.2)

* Check rejection decision
if e(reject) == 1 {
    display "Equivalence concluded: max|beta_t| < 0.2"
}
else {
    display "Cannot conclude equivalence at threshold 0.2"
}
```

## Citation

If you use this package in your research, please cite both the methodology paper and the Stata implementation:

**APA Format:**

> Dette, H., & Schumann, M. (2024). Testing for Equivalence of Pre-Trends in Difference-in-Differences Estimation. *Journal of Business & Economic Statistics*, 42(4), 1289–1301. https://doi.org/10.1080/07350015.2024.2308121
>
> Cai, X., & Xu, W. (2025). *Equitrends: Stata module for equivalence tests for pre-trends in DiD* (Version 0.1.0) [Computer software]. GitHub. https://github.com/gorgeousfish/equitrends

**BibTeX:**

```bibtex
@article{dette2024testing,
  title={Testing for Equivalence of Pre-Trends in Difference-in-Differences Estimation},
  author={Dette, Holger and Schumann, Martin},
  journal={Journal of Business \& Economic Statistics},
  volume={42},
  number={4},
  pages={1289--1301},
  year={2024},
  publisher={Taylor \& Francis},
  doi={10.1080/07350015.2024.2308121}
}

@software{equitrends2025stata,
  title={Equitrends: Stata module for equivalence tests for pre-trends in DiD},
  author={Cai, Xuanyu and Xu, Wenli},
  year={2025},
  version={0.1.0},
  url={https://github.com/gorgeousfish/equitrends}
}
```

## Authors

**Stata Implementation:**

- **Xuanyu Cai**, City University of Macau
  Email: [xuanyuCAI@outlook.com](mailto:xuanyuCAI@outlook.com)
- **Wenli Xu**, City University of Macau
  Email: [wlxu@cityu.edu.mo](mailto:wlxu@cityu.edu.mo)

**Methodology:**

- **Holger Dette**, Department of Mathematics, Ruhr University Bochum
- **Martin Schumann**, School of Business and Economics, Maastricht University

## License

AGPL-3.0 License. See [LICENSE](LICENSE).

## Future roadmap

The development team is evaluating the following extensions for future versions:

- **Staggered Adoption Support**: Extending the equivalence testing framework to staggered treatment adoption designs via cohort stacking (Dette & Schumann, 2024, Section 5).
- **Threshold Sensitivity Analysis**: Visualizing how the test decision and minimum threshold vary across a continuous range of equivalence thresholds.
- **Additional Bootstrap Methods**: Implementing multiplier bootstrap and subsampling alternatives.
- **Bootstrap Parallelization**: Multi-core bootstrap computation via the `parallel` package for reduced runtime on large datasets.

## Related packages

- **EquiTrends (R)**: [TiesBos/EquiTrends](https://github.com/TiesBos/EquiTrends) — Original R implementation of equivalence tests for pre-trends (Dette & Schumann, 2024)
- **pretest (Stata)**: [gorgeousfish/pretest](https://github.com/gorgeousfish/pretest) — Conditional extrapolation pre-test with bias-adjusted confidence intervals for DiD (Mikhaeil & Harshaw, 2025)

## Recommended resources

For beginners in causal inference and econometrics:

- [Causal Inference for the Brave and True](https://matheusfacure.github.io/python-causality-handbook/landing-page.html) — An excellent introductory tutorial on causal inference by Matheus Facure
- [Causal Inference for the Brave and True (Chinese Edition)](https://ci-book.huangwz.com/intro) — Chinese translation by Wenzhe Huang and Wenli Xu
- [What&#39;s Trending in Difference-in-Differences? A Synthesis of the Recent Econometrics Literature](https://doi.org/10.1016/j.jeconom.2023.03.008) — Comprehensive review by Roth, Sant'Anna, Bilinski & Poe (2023)

## Support

For questions and bug reports, please open an issue on [GitHub](https://github.com/gorgeousfish/equitrends/issues).
