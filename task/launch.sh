# inputs
# data: /data/galucia/lithium4/LHC24_pass1_skimmed/LHC24ag_pass1_skimmed_slice.root
# data: /Users/glucia/Projects/ALICE/data/lithium/same/LHC24ag_pass1_skimmed_test.root
# data: /data/galucia/lithium_local/local_task/LHC23zzh_pass4_slice0.root (slice1, slice2)
# mc: /data/galucia/lithium4/MC/AO2D_injectedli4.root

LOGFILE="output.log"
#CONF="-b --configuration json://configuration.json"
CONF="-b --configuration json://configuration_mc.json"
#OUTPUT_DIR="OutputDirector.json"
OUTPUT_DIR="OutputDirector_mc.json"


o2-analysis-lf-he3hadronfemto $CONF --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000|

    # converters
    #o2-analysis-tracks-extra-converter $CONF --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000|
    o2-analysis-tracks-extra-v002-converter $CONF|
    #o2-analysis-mc-converter $CONF|
    #o2-analysis-bc-converter $CONF --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000|
    #o2-analysis-mccollision-converter $CONF --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000|

    # standard wagons
    o2-analysis-timestamp $CONF --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000|
    o2-analysis-multiplicity-table $CONF --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000|
    o2-analysis-centrality-table $CONF --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000|
    o2-analysis-event-selection $CONF --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000|
    o2-analysis-track-propagation $CONF --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000|
    o2-analysis-trackselection $CONF --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000|
    o2-analysis-pid-tof-full $CONF --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000 |
    o2-analysis-pid-tof-base $CONF --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000|
    o2-analysis-pid-tpc-base $CONF --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000|
    o2-analysis-ft0-corrected-table $CONF --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000|
    o2-analysis-pid-tpc $CONF --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000 --aod-writer-json $OUTPUT_DIR --aod-file @input_data.txt > $LOGFILE


# report the status of the workflow
rc=$?
if [ $rc -eq 0 ]; then
    echo "Workflow finished successfully"
else
    echo "Error: Workflow failed with status $rc"
    echo "Check the log file for more details: $LOGFILE"
    exit $rc
fi

    