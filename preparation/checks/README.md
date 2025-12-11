# Processed Datasets Overview

This directory contains several processed datasets produced with different selection strategies for comparison.

## **Standard Selection**
File: **`same_event_all_pass1_pass4.root`**

Uses:
- **Central TPC calibration**
- **Hard TOF cut**
- **PID in tracking corrections**

---

## **Alternative Datasets**
The datasets below each differ from the standard selection by one specific modification.

### **Summary Table**

| File | Difference From Standard | Notes |
|------|--------------------------|-------|
| **`same_event_all_pass1_pass4_no_pid_tracking.root`** | No PID in tracking correction | He3 tracked as He3, H3, He4; no pt correction applied; ≥2.5 GeV/c still requires He3 hypothesis |
| **`same_event_all_pass1_pass4_no_triton.root`** | No triton mass hypothesis included | He3 tracked as He3; no pt correction needed |
| **`same_event_all_pass1_pass4_tof_veto.root`** | Conditional TOF veto | Proton TOF cut applied only if TOF is reached |
| **`same_event_pass1_pass4_offline_tpc.root`** | Offline TPC calibration | Uses `Common.h` parametrisation instead of the TTree `fNSigmaTPCHe3` column |
| **`same_event_pass1_pass4_nclstpc.root`** | Cluster-count requirement | Applies **fNClsTPCHe3 > 110** |
| **`same_event_pass1_pass4_nclstpc_tofvetohe3.root`** | Cluster-count + TOF mass veto | fNClsTPCHe3 > 110 **and** m_TOF > 2 GeV/c^2 for tracks reaching TOF |
| **`mixed_event_pass1_pass4_nclstpc_5me_batches.root`** | Mixed-event extension | Same cluster cuts as above, includes **256** and **3112** mixed-event batches |
| **`same_event_pass1_pass4_nocuts.root`** | No cuts during event mixing | Should not affect same-event results |

---

## **Notes**
- These variations allow systematic comparison of PID strategies, TPC calibrations, TOF logic, and track-quality cuts.
- Mixed-event datasets include optional batch extensions for improved statistics.
