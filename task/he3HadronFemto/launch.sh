# inputs
# data: /data/galucia/lithium4/LHC24_pass1_skimmed/LHC24ag_pass1_skimmed_slice.root
# data: /Users/glucia/Projects/ALICE/data/lithium/same/LHC24ag_pass1_skimmed_test.root
# data: /data/galucia/lithium_local/local_task/LHC23zzh_pass4_slice0.root (slice1, slice2)
# mc: /data/galucia/lithium4/MC/AO2D_injectedli4.root
# mc: /data/galucia/lithium_local/MC/raw/LHC25a4_001_raw.root
# mc: /data/galucia/lithium_local/MC/raw/LHC25a4_002_raw.root
# mc: /data/galucia/lithium_local/MC/raw/LHC25a4_003_raw.root
# mc: /data/galucia/lithium_local/MC/raw/LHC25a4_004_raw.root

LOGFILE="output.log"
CONF="-b --configuration json://configuration.json"
#CONF="-b --configuration json://configuration_mc.json"
OUTPUT_DIR="OutputDirector.json"
#OUTPUT_DIR="OutputDirector_mc.json"



o2-analysis-lf-nucleiqc $CONF |
    
    #o2-analysis-mc-converter $CONF|
    o2-analysis-mccollision-converter $CONF --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000|

    o2-analysis-pid-tof-merge $CONF |
    #o2-analysis-tracks-extra-v002-converter $CONF |
    o2-analysis-multcenttable $CONF |
    o2-analysis-event-selection-service $CONF |
    o2-analysis-pid-tpc-service $CONF |
    o2-analysis-ft0-corrected-table $CONF --aod-writer-json $OUTPUT_DIR --aod-file @input_data.txt > $LOGFILE

# report the status of the workflow
rc=$?
if [ $rc -eq 0 ]; then
    echo "Workflow finished successfully"
else
    echo "Error: Workflow failed with status $rc"
    echo "Check the log file for more details: $LOGFILE"
    #exit $rc
fi

    