#!/bin/bash

SHINE_DIR=/afs/cern.ch/work/v/viklochk/public/Shine_Legacy/shine_install

function initShine {
	if ! source "$SHINE_DIR/scripts/env/lxplus_32bit_slc6.sh"; then exit 1; fi
	eval $($SHINE_DIR/bin/shine-offline-config --env-sh)	
}


timestamp_start=$(date +%Y%m%d_%H%M%S)

run=${1}
shift 1
drop_dir=${1}
shift 1
in_dir=${1}

####################
# Init Shine
####################

echo
echo "Configuring Shine... "
initShine

echo
echo Merging input...
hadd run_${run}.root ${in_dir}/*${run}*x*.root

log_file=log_$run.txt

echo
echo Starting calibration...
./vD -i run_${run}.root --root-output &> $log_file

cp *output.root $drop_dir/run_${run}_calib.root
cp run_${run}.root $drop_dir/



EXITCODE=$?

echo `ls -lrt`

timestamp_end=$(date +%Y%m%d_%H%M%S)
echo "TIME: ""$timestamp_start"" -- ""$timestamp_end"

exit $EXITCODE
