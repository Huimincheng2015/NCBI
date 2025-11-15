# Timing Modification Documentation

## Overview

This modification separates **estimation time** from **cross-validation overhead time** in the hyperparameter tuning process for graphon estimation methods (GGCV and ECV).

## Key Changes

### 1. New Timed Functions (`cv_functions_timed.R`)

#### `BM_timed()` - Timed version of GGCV cross-validation
**Timing breakdown:**
- **Estimation time**: Time spent fitting the graphon model to each of K folds
- **CV overhead time**:
  - Fold construction (masking/perturbation)
  - Rescaling predictions
  - Computing validation loss

**Returns:**
```r
list(
  cv_score = <MSE value>,
  time_estimation = <total estimation time across all folds>,
  time_cv_overhead = <total CV overhead time>,
  time_total = <sum of above>
)
```

#### `CV_timed()` - Timed version of ECV cross-validation
**Timing breakdown:**
- **Estimation time**: Time spent fitting the graphon model to each of K imputed networks
- **CV overhead time**:
  - Fold construction (masking entries)
  - SVD imputation for network completion
  - Computing validation loss

**Returns:**
```r
list(
  cv_score = <MSE value>,
  time_estimation = <total estimation time across all folds>,
  time_cv_overhead = <total CV overhead time>,
  time_total = <sum of above>
)
```

### 2. Modified Main Script (`model_selection_accuracy_timed.R`)

#### New Data Columns
The result dataframe now includes:
- `time_BM_total` - Total time for GGCV (estimation + CV overhead)
- `time_BM_estimation` - Estimation time only for GGCV
- `time_BM_cv_overhead` - CV overhead time only for GGCV
- `time_CV_total` - Total time for ECV (estimation + CV overhead)
- `time_CV_estimation` - Estimation time only for ECV
- `time_CV_cv_overhead` - CV overhead time only for ECV

These times are summed across all hyperparameter values tested.

#### New Output Files

1. **`par_sel_ns_timed.csv`** - Main results file with detailed timing columns

2. **`../time_results/time_NS_total.pdf`** - Total runtime comparison (original metric)

3. **`../time_results/time_NS_cv_overhead_only.pdf`** - **NEW: CV overhead only comparison**
   - Shows only the cross-validation overhead time
   - Excludes model estimation time
   - Direct comparison of GGCV vs ECV CV overhead

4. **`../time_results/time_breakdown_stacked.pdf`** - **NEW: Stacked bar chart breakdown**
   - Shows both estimation and CV overhead times
   - Faceted by method (GGCV vs ECV)
   - Grouped by sample size

5. **`../time_results/timing_breakdown_detailed.csv`** - **NEW: Detailed timing data**
   - Tabular format of all timing breakdowns
   - Useful for further analysis

## Interpretation

### CV Overhead Time Components

**GGCV (BM method):**
- Fold construction: Randomly perturbing K validation sets with Bernoulli(θ) samples
- Validation: Rescaling predictions and computing MSE on validation sets

**ECV (CV-imputation method):**
- Fold construction: Masking validation sets
- **SVD imputation**: Completing the masked network (dominant cost)
- Validation: Computing MSE on validation sets

### Expected Results

**ECV should have higher CV overhead than GGCV because:**
1. SVD imputation is computationally expensive (O(n³) per fold)
2. Rank estimation via `rankMatrix()` adds overhead
3. GGCV only does simple perturbation (O(m) where m = #edges)

### Formula Verification

For each method:
```
CV_overhead_time = Total_time - Estimation_time
```

Where:
- **Total_time** = Complete hyperparameter tuning runtime
- **Estimation_time** = Sum of model fitting times across all folds and hyperparameters
- **CV_overhead_time** = Fold construction + validation overhead

## Usage

### Running the Modified Code

```r
source("model_selection_accuracy_timed.R")
```

### Accessing Timing Results

```r
# Load results
dat = read.csv("par_sel_ns_timed.csv")

# Extract CV overhead for a specific run
sample_run = dat[dat$graphon_id == 5 & dat$sample_size == 200 & dat$repeat_i == 1, ]

# GGCV timing
ggcv_total = sample_run$time_BM_total[1]
ggcv_estimation = sample_run$time_BM_estimation[1]
ggcv_cv_overhead = sample_run$time_BM_cv_overhead[1]

cat("GGCV Total:", ggcv_total, "secs\n")
cat("GGCV Estimation:", ggcv_estimation, "secs\n")
cat("GGCV CV Overhead:", ggcv_cv_overhead, "secs\n")
cat("Verification:", ggcv_estimation + ggcv_cv_overhead == ggcv_total, "\n")
```

## Implementation Notes

1. **Timing granularity**: We use `Sys.time()` and `difftime()` for millisecond-level precision

2. **Aggregation**: Times are summed across:
   - All K folds in cross-validation
   - All hyperparameter values tested (h_lis)

3. **Parallel execution**: Each `repeat_i` runs independently, so timing is per-repeat

4. **Validation**: The sum `time_estimation + time_cv_overhead` should equal `time_total` (within rounding error)

## Verification Checklist

- [ ] `time_BM_total ≈ time_BM_estimation + time_BM_cv_overhead`
- [ ] `time_CV_total ≈ time_CV_estimation + time_CV_cv_overhead`
- [ ] CV overhead time for ECV > CV overhead time for GGCV (expected due to SVD)
- [ ] Timing values are non-negative
- [ ] PDF plots generated successfully

## Troubleshooting

**Issue**: Total time doesn't match sum of components
- **Cause**: Minor rounding errors or timing measurement overhead
- **Solution**: Check if difference < 0.01 seconds (acceptable)

**Issue**: CV overhead time is 0
- **Cause**: `run_ecv = FALSE` so ECV is not executed
- **Solution**: Set `run_ecv = TRUE` in the main script

**Issue**: Missing plots
- **Cause**: Directory `../time_results/` doesn't exist
- **Solution**: Create directory: `dir.create("../time_results", recursive=TRUE)`
