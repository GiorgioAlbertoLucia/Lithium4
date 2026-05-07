# Femtoscopy Analysis Framework

A unified framework for femtoscopy correlation function analysis, including significance calculation and yield estimation.

## Overview

This framework provides a modular, organized codebase for two related analyses:
1. **Significance Calculation** - Computes chi2 distributions and statistical significance of correlation function features
2. **Yield and Upper Limit** - Estimates particle yields and their uncertainties through toy Monte Carlo

## Framework Structure

### Header Files (Shared Utilities)

#### `FemtoConfig.h`
Central configuration file containing:
- Common constants (k* ranges, chi2 histogram parameters)
- Enumeration for sampling methods (Poissonian, Gaussian)
- Data structures (InputDataPaths, IterationResult)

#### `HistogramUtils.h`
Histogram manipulation utilities:
- `loadHistograms()` - Load same/mixed event histograms from files
- `loadSameMixedWithCorrection()` - Load and apply correction factors
- `computeCorrelationFunction()` - Calculate C(k*) = Same/Mixed with error propagation
- `sampleHistogramPoisson()` - Poisson sampling for toy MC
- `sampleHistogramGaussian()` - Gaussian sampling for toy MC

#### `Chi2Utils.h`
Chi2 and significance calculation utilities:
- `computeChi2()` - Calculate chi2 between observed and expected distributions
- `fillRunningChi2()` - Compute cumulative chi2 vs k*
- `fillWindowChi2()` - Compute chi2 in sliding windows
- `createSignificanceGraphs()` - Generate p-value and significance graphs

#### `RooFitUtils.h`
RooFit model preparation utilities (for yield calculation):
- `prepareSignalModel()` - Load and setup signal PDF from histogram
- `prepareBackgroundModel()` - Load and setup background PDF from histogram
- `performFit()` - Perform fits with diagnostic plotting

### Analysis Programs

#### `computeSignificance.cpp`
Computes statistical significance of correlation function features.

**Purpose**: Determine if observed correlation function deviates significantly from null hypothesis

**Method**:
1. Load observed correlation functions (Matter and Antimatter)
2. Generate toy MC samples using Poisson sampling
3. Compute chi2 distributions for null hypothesis
4. Calculate running chi2 and windowed chi2
5. Determine p-values and significance vs k*

**Output**: `output/significance.root` containing:
- Chi2 distributions for signal region and background region
- Running chi2 distributions for each k* bin
- Graphs of p-value and significance vs k*
- Separate results for Matter and Antimatter

#### `compute_yield_and_upper_limit.cpp`
Estimates particle yields through template fitting.

**Purpose**: Extract signal yield and uncertainties from correlation function

**Method**:
1. Load same/mixed event histograms
2. Prepare signal and background templates (RooFit PDFs)
3. For each toy MC iteration:
   - Sample correlation function (Poisson or Gaussian)
   - Pre-fit background in flat region
   - Fit signal + background in full range
   - Extract yield
4. Fit yield distribution to estimate mean and uncertainty

**Output**: `yield_upper_limit_[Matter/Antimatter].root` containing:
- Raw yield distribution
- Chi2 distribution from fits
- Diagnostic plots of fits (every N iterations)
- High-yield fit examples

## Building the Code

### Prerequisites
- ROOT (tested with ROOT 6.x)
- RooFit (included with ROOT)
- C++17 compatible compiler

### Compilation

```bash
# Build all programs
make

# Build only significance calculation
make computeSignificance

# Build only yield calculation
make compute_yield_and_upper_limit

# Clean build artifacts
make clean
```

## Running the Analysis

### Significance Calculation

```bash
./computeSignificance
```

Configuration is in the `SignificanceConfig` namespace within the source file. Key parameters:
- `N_ITERATIONS` - Number of toy MC iterations (default: 1,000,000)
- `N_BINS` - Number of k* bins (must match input histogram)
- `N_BINS_WINDOW` - Size of sliding window for windowed chi2

### Yield and Upper Limit

```bash
./compute_yield_and_upper_limit
```

Configuration is in the `YieldConfig` namespace. Key parameters:
- `N_ITERATIONS` - Number of toy MC iterations (default: 1,000)
- `SAMPLING_METHOD` - Poissonian or Gaussian
- `KSTAR_MIN/MAX` - Full fit range
- `PREFIT_MIN/MAX` - Background-only fit range (flat region)

## Customization

### Modifying File Paths

Edit the configuration namespaces in each `.cpp` file:

**For Significance**:
```cpp
namespace SignificanceConfig {
    InputDataPaths paths = {
        .inputMixedFile = "your/path/to/mixed.root",
        // ... other paths
    };
}
```

**For Yield**:
```cpp
namespace YieldConfig {
    const char* CORRELATION_FILE = "your/path/to/correlation.root";
    const char* SIGNAL_FILE = "your/path/to/signal.root";
    // ... other paths
}
```

### Adjusting Analysis Parameters

Common parameters can be modified in `FemtoConfig.h`:
```cpp
namespace Config {
    const double KSTAR_MIN = 0.02;  // Lower k* bound
    const double KSTAR_MAX = 0.4;   // Upper k* bound
    // ... other parameters
}
```

### Changing Sampling Method

In `compute_yield_and_upper_limit.cpp`:
```cpp
// Choose one:
SamplingMethod SAMPLING_METHOD = SamplingMethod::POISSONIAN;
// or
SamplingMethod SAMPLING_METHOD = SamplingMethod::GAUSSIAN;
```

## Code Organization Benefits

### Before (Original Code)
- Duplicated functions (e.g., `computeCorrelationFunction` in both files)
- Inconsistent implementations
- Hard to maintain consistency across analyses
- No clear separation of concerns

### After (Unified Framework)
- **Shared utilities** in header files eliminate duplication
- **Consistent implementations** across analyses
- **Modular design** allows easy extension
- **Clear separation** between configuration, utilities, and analysis logic
- **Easier maintenance** - fix bugs in one place
- **Better documentation** through organized structure

## Key Features

### Shared Functionality
- Common histogram operations
- Consistent correlation function computation
- Unified chi2 calculation methods
- Standard error propagation

### Flexible Configuration
- Easy to switch between Matter/Antimatter
- Adjustable iteration counts
- Configurable k* ranges
- Multiple sampling methods

### Diagnostic Output
- Comprehensive ROOT files with all intermediate results
- Optional diagnostic plotting (controlled by iteration interval)
- Automatic detection and saving of high-yield events

## Future Extensions

The modular structure makes it easy to add:
- New sampling methods (add to `SamplingMethod` enum and `HistogramUtils.h`)
- Additional chi2 calculations (add to `Chi2Utils.h`)
- Different fit models (extend `RooFitUtils.h`)
- Alternative significance metrics (modify `Chi2Utils.h`)

## Troubleshooting

### Common Issues

1. **ROOT not found**: Ensure ROOT is properly installed and `root-config` is in PATH
2. **RooFit missing**: RooFit should be included with ROOT installation
3. **File not found errors**: Check file paths in configuration namespaces
4. **Histogram binning mismatch**: Ensure `N_BINS` matches input histogram binning

### Debug Mode

To enable verbose output, add print statements or use ROOT's logging:
```cpp
gErrorIgnoreLevel = kInfo;  // Show all messages
```

## Contact

For questions about the framework structure or implementation, please refer to the inline documentation in the header files.
