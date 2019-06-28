#!/bin/bash

#############################################################
echo ""
echo "========================= START ========================="

# AMOUNT OF CLUSTERS THIS JOB WILL RUN ON
NUM_CLUSTERS=1

# ROOT DIR IN HIERARCHY
# PLACE FOR SMALL LOGS 
RECO_WORK_DIR=/afs/cern.ch/work/v/viklochk/public/test/

# LOCATION OF OUTPUT *.root files
# PLACE FOR DATATREE
RECO_DROP_DIR=/eos/user/v/viklochk/public/vD/

calib_dir=/eos/user/v/viklochk/public/vD/calibrated/

EXECUTABLE=/afs/cern.ch/work/v/viklochk/public/vD/TpcCalibration/NewAnalyzer/vD

CONDOR_CMD=condor_submit
CMD=condor_run.sh
QUEUE=microcentury

# espresso     = 20 minutes
# microcentury = 1 hour
# longlunch    = 2 hours
# workday      = 8 hours
# tomorrow     = 1 day
# testmatch    = 3 days
# nextweek     = 1 week

runlist=$1
shift 1
commit=$1
if [[ -z $commit ]]; then commit=no_commit; fi

#############################################################

timestamp=$(date +%Y%m%d_%H%M%S)

work_dir=$RECO_WORK_DIR/$commit/$timestamp
mkdir -p "$work_dir"
echo "[INFO] Output directory: $work_dir"

batch_dir=$work_dir/batch
mkdir -p $batch_dir
echo "[INFO] Batch directory: $batch_dir"

drop_dir=$RECO_DROP_DIR/$commit/$timestamp
mkdir -p $drop_dir

echo "[INFO] Drop directory: $drop_dir"

rsync -r ./condor_run.sh   $batch_dir
rsync -r ./condor.sub.in.head $batch_dir
rsync -r $EXECUTABLE  $batch_dir
rsync -r $runlist $batch_dir

chmod -R 0755 $batch_dir

##############################
#	INSIDE BATCH
#############################
cd $batch_dir

transfer_input_files=$EXECUTABLE

#
# CLEAN-UP BATCH
#
rm -f ./*.sub

#
# POPULATE HEAD OF SUB
#

execution_prefix=`head /dev/urandom | tr -dc A-Za-z0-9 | head -c8`

sub_head_name=.${execution_prefix}_condor.sub.head
cat ./condor.sub.in.head | sed \
	-e "s~@EXECUTABLE@~${batch_dir}/${CMD}~g" \
	-e "s~@QUEUE@~${QUEUE}~g" \
	-e "s~@TRANSFER_INPUT_FILES@~${transfer_input_files}~g" \
	-e "s~@TRANSFER_OUTPUT_FILES@~${transfer_output_files}~g" \
> $sub_head_name

runlist=`basename $runlist`

pwd

run=""
i=0
while read run; do
  batch_ind=`expr $i % $NUM_CLUSTERS`
  sub_batch_name=.${execution_prefix}_condor_${batch_ind}.sub

  if ! [[ -e $sub_batch_name ]]; then cat $sub_head_name > $sub_batch_name; fi 

  echo "$i:$batch_ind >--------------------"$run"--------------------"

  initial_dir=$work_dir/$run
  mkdir -p $initial_dir
  echo "[INFO] InitialDir $initial_dir"

  condor=$initial_dir/condor_${run}.sub

  cat $sub_head_name > $condor

  echo | tee -a $condor $sub_batch_name
  echo "initialdir  = $initial_dir"				              | tee -a $condor $sub_batch_name > /dev/null
  echo "arguments    = $run $drop_dir $calib_dir"             | tee -a $condor $sub_batch_name > /dev/null
  echo "queue"                                                | tee -a $condor $sub_batch_name > /dev/null

  let i+=1 
done <$runlist

#
# SUBMIT TASK
#
echo Submitting jobs...
echo

for sub in `find ./ -type f -name *.sub`
do
	echo Submitting ${sub}...
	until $CONDOR_CMD $sub; do echo "[WARN] Retry"; done
done

echo Cleaning up...
echo rm -f .${execution_prefix}*
rm -r .${execution_prefix}*

echo "------------------------------------------------"
cd -

#############################################################

echo "========================= END ========================="
echo ""
