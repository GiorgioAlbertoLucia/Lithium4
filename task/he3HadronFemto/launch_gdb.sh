LOGFILE="output.log"
CONF=(-b --configuration json://configuration_mc.json --shm-segment-size 750000000000)
OUTPUT_DIR="OutputDirector_mc.json"

GDB_WRAP="gdb -batch -ex run -ex bt -ex quit --args"
CONF_DEBUG=("${CONF[@]}" --child-driver "$GDB_WRAP")

o2-analysis-lf-he3hadronfemto "${CONF_DEBUG[@]}" |
    o2-analysis-pid-tof-merge "${CONF_DEBUG[@]}" |
    o2-analysis-trackselection "${CONF_DEBUG[@]}" |
    o2-analysis-propagationservice "${CONF_DEBUG[@]}" |
    o2-analysis-multcenttable "${CONF_DEBUG[@]}" |
    o2-analysis-event-selection-service "${CONF_DEBUG[@]}" |
    o2-analysis-pid-tpc-service "${CONF_DEBUG[@]}" |
    o2-analysis-ft0-corrected-table "${CONF_DEBUG[@]}" --aod-writer-json $OUTPUT_DIR --aod-file @input_data.txt > $LOGFILE 2>&1

rc=$?
if [ $rc -eq 0 ]; then
    echo "Workflow finished successfully"
else
    echo "Error: Workflow failed with status $rc"
    echo "Check the log file for more details: $LOGFILE"
fi