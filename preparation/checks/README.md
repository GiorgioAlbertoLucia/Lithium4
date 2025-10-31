This folder contains several processed dataset, using  different strategies for comparison.
The "standard" selection, adopted in same_event_all_pass1_pass4.root, uses
* central TPC calibration
* hard cut on TOF
* PID in tracking corrections

Other dataset differ by one of these selections
- same_event_all_pass1_pass4_no_pid_tracking.root
    * contains He3 tracked as He3, H3, He4. No correction is applied to the reconstructed pt based on the pid in tracking hypothesis (after pt > .5, He3 hypotehesis is still required)
- same_event_all_pass1_pass4_no_triton.root
    * contains He3 tracked as He3. No correction is applied (but none is really required)
- same_event_all_pass1_pass4_tof_veto.root
    * proton TOF cut is only applied if the hadron reaches the TOF, otherwise it is skipped
- same_event_pass1_pass4_offline_tpc.root
    * He3 TPC calibration using the parameterisation in Common.h (instead of using the fNSigmaTPCHe3 column as already defined in the TTree)
