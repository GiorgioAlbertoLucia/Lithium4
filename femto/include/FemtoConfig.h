/**
 * FemtoConfig.h
 * 
 * Configuration and data structures shared across femtoscopy analysis tools
 */

#ifndef FEMTO_CONFIG_H
#define FEMTO_CONFIG_H

#include <string>
#include <unordered_map>

namespace FemtoAnalysis {

// Sampling method for toy MC generation
enum class SamplingMethod {
    POISSONIAN,
    GAUSSIAN
};

// Common configuration parameters
namespace Config {
    // File paths - these should be adjusted per analysis
    struct FilePaths {
        const char* correlationFile;
        const char* signalFile;
        const char* backgroundFile;
        const char* chi2File;
        const char* outputFile;
    };
    
    // Analysis parameters
    const double KSTAR_MIN = 0.02;
    const double KSTAR_MAX = 0.4;
    const double PREFIT_MIN = 0.2;
    const double PREFIT_MAX = 0.4;
    
    // Reference k* values for normalization
    const double REFERENCE_KSTAR_FOR_BKG_NORMALIZATION = 0.3;
    const double REFERENCE_KSTAR_FOR_SIG_NORMALIZATION = 0.07;
    
    // Chi2 histogram parameters
    const int NBINS_CHI2 = 1600;
    const float CHI2_MAX_VALUE = 160;
    
    // Iteration parameters (can be overridden per analysis)
    const int DEFAULT_N_ITERATIONS = 1000000;
    const int DEFAULT_PRINT_INTERVAL = 10000;
    const double HIGH_YIELD_THRESHOLD = 1000.0;
}

// Structure for input data paths (used in significance computation)
struct InputDataPaths {
    const char* inputMixedFile;
    const char* inputMixedNameMatter;
    const char* inputMixedNameAntimatter;
    const char* inputCorrectionFile;
    const char* inputCorrectionNameAntimatter;
    const char* inputCorrectionNameMatter;
    const char* inputChi2File;
    const char* inputChi2NameAntimatter;
    const char* inputChi2NameMatter;
};

// Structure for iteration results
struct IterationResult {
    double yield;
    double chi2;
};

} // namespace FemtoAnalysis

#endif // FEMTO_CONFIG_H
