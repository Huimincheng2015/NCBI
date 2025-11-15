# Summary of Timing Modifications

## Files Created

### 1. `cv_functions_timed.R`
**Purpose:** New implementations of BM and CV functions with detailed timing tracking

**Functions:**
- `BM_timed()` - GGCV with separate estimation/overhead timing
- `CV_timed()` - ECV with separate estimation/overhead timing

**Key differences from original:**
- Returns a list with timing breakdown instead of just CV score
- Tracks time for each fold's estimation and validation separately
- Aggregates timing across all folds

### 2. `model_selection_accuracy_timed.R`
**Purpose:** Modified main script using the timed CV functions

**Key modifications:**
```r
# OLD: Called BM() and CV(), only tracked total time
BM_lis[i] = BM(...)
time_our[i] = total_time

# NEW: Call BM_timed() and CV_timed(), capture breakdown
BM_result = BM_timed(...)
BM_lis[i] = BM_result$cv_score
time_BM_total[i] = BM_result$time_total
time_BM_estimation[i] = BM_result$time_estimation
time_BM_cv_overhead[i] = BM_result$time_cv_overhead
```

**New data columns:**
- `time_BM_total`, `time_BM_estimation`, `time_BM_cv_overhead`
- `time_CV_total`, `time_CV_estimation`, `time_CV_cv_overhead`

**New visualizations:**
1. Total time plot (same as before)
2. **CV overhead only plot** (MAIN NEW CONTRIBUTION)
3. Stacked breakdown plot (estimation vs overhead)

### 3. `TIMING_MODIFICATION_README.md`
**Purpose:** Complete documentation of the modifications

**Contents:**
- Detailed explanation of timing breakdown
- Interpretation guide
- Usage examples
- Verification checklist

### 4. `test_timing_functions.R`
**Purpose:** Minimal test to verify timing functions work correctly

**Tests:**
- BM_timed returns valid results
- CV_timed returns valid results
- Total = estimation + overhead (verification)
- Timing values are sensible

### 5. `CHANGES_SUMMARY.md` (this file)
**Purpose:** High-level overview of all changes

## What You Asked For

✅ **Separately record:**
   - (a) Time spent on graphon estimation/model fitting
   - (b) Time spent on cross-validation overhead

✅ **Compute CV-only runtime:**
   - Formula: `CV_overhead = Total_time - Estimation_time`
   - Implemented in both `BM_timed()` and `CV_timed()`

✅ **Produce new figure:**
   - File: `../time_results/time_NS_cv_overhead_only.pdf`
   - Compares GGCV vs ECV CV overhead across sample sizes

## Additional Contributions

Beyond what was requested, I also provided:

1. **Stacked breakdown plot** showing both components visually
2. **Detailed CSV export** of all timing data for further analysis
3. **Test script** to verify implementation correctness
4. **Comprehensive documentation** with interpretation guide

## How to Use

### Quick Start
```r
# Run the modified version
source("model_selection_accuracy_timed.R")

# The script will generate:
# - par_sel_ns_timed.csv (results with timing)
# - ../time_results/time_NS_cv_overhead_only.pdf (your requested plot)
# - ../time_results/time_breakdown_stacked.pdf (bonus visualization)
```

### Verify Installation
```r
# Test that timing functions work
source("test_timing_functions.R")
# Should print timing breakdowns and verification
```

### Analyze Results
```r
# Load results
dat = read.csv("par_sel_ns_timed.csv")

# Get CV overhead for GGCV at n=200
dat_200 = dat[dat$sample_size == 200,]
mean(dat_200$time_BM_cv_overhead)  # GGCV overhead
mean(dat_200$time_CV_cv_overhead)  # ECV overhead
```

## Key Insights

### Expected Results

**GGCV CV Overhead includes:**
- Fold construction: O(m) where m = number of edges
- Perturbation: O(m/K) per fold
- Validation: O(n²) per fold

**ECV CV Overhead includes:**
- Fold construction: O(m) masking
- **SVD imputation: O(n³) per fold** ← Dominant cost!
- Validation: O(n²) per fold

**Prediction:** ECV overhead >> GGCV overhead (especially for large n)

### Verification

The code includes automatic verification:
```r
time_total ≈ time_estimation + time_cv_overhead
```

If this doesn't hold (within ~0.01 sec), there's a timing bug.

## File Structure

```
/home/user/NCBI/
├── cv_functions_timed.R              # New timed CV functions
├── model_selection_accuracy_timed.R  # Modified main script
├── test_timing_functions.R           # Test script
├── TIMING_MODIFICATION_README.md     # Full documentation
├── CHANGES_SUMMARY.md                # This file
└── par_sel_ns_timed.csv              # Output (generated when run)

../time_results/                       # (created by script)
├── time_NS_cv_overhead_only.pdf      # ⭐ YOUR REQUESTED PLOT
├── time_NS_total.pdf                 # Total time comparison
├── time_breakdown_stacked.pdf        # Stacked breakdown
└── timing_breakdown_detailed.csv     # Detailed timing data

../box_results_ns/                     # (existing output directory)
├── h_paper.pdf
└── ...
```

## Comparison with Original Code

| Aspect | Original | Modified |
|--------|----------|----------|
| BM function | Returns CV score only | Returns {score, times} |
| CV function | Returns CV score only | Returns {score, times} |
| Timing data | Total time only | Total + estimation + overhead |
| Output CSV | `par_sel_ns.csv` | `par_sel_ns_timed.csv` |
| Plots | Total time plot | Total + CV-only + breakdown |
| Documentation | In-code comments | Full README |

## Next Steps

1. **Create required directories:**
   ```r
   dir.create("../time_results", recursive=TRUE)
   dir.create("../box_results_ns", recursive=TRUE)
   ```

2. **Ensure dependencies are available:**
   - All source files (ICE.R, function.R, network_generate.R, graphon-master/)
   - All required R packages

3. **Run the test:**
   ```r
   source("test_timing_functions.R")
   ```

4. **Run the full analysis:**
   ```r
   source("model_selection_accuracy_timed.R")
   ```

5. **Check the CV-overhead plot:**
   ```r
   system("open ../time_results/time_NS_cv_overhead_only.pdf")
   ```

## Questions?

Refer to `TIMING_MODIFICATION_README.md` for detailed documentation.
