# Makefile for Femtoscopy Analysis Tools
# Requires ROOT to be installed and configured

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2 -I.
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs) -lRooFit -lRooFitCore

# Targets
TARGETS = computeSignificance compute_yield_and_upper_limit

# Header files
HEADERS = include/FemtoConfig.h include/HistogramUtils.h include/Chi2Utils.h include/RooFitUtils.h

# Default target
all: $(TARGETS)

# Build computeSignificance
computeSignificance: computeSignificance.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $< $(ROOTLIBS)

# Build compute_yield_and_upper_limit
compute_yield_and_upper_limit: compute_yield_and_upper_limit.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $< $(ROOTLIBS)

# Clean build artifacts
clean:
	rm -f $(TARGETS)

# Phony targets
.PHONY: all clean

# Help target
help:
	@echo "Femtoscopy Analysis Tools - Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  all                            - Build all programs (default)"
	@echo "  computeSignificance            - Build significance calculation program"
	@echo "  compute_yield_and_upper_limit  - Build yield and upper limit program"
	@echo "  clean                          - Remove built executables"
	@echo "  help                           - Show this help message"
	@echo ""
	@echo "Usage:"
	@echo "  make                           - Build all programs"
	@echo "  make computeSignificance       - Build only significance program"
	@echo "  make clean                     - Clean build artifacts"
