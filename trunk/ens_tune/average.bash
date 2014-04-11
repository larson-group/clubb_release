#!/bin/bash

# Script to average the netCDF outputs of a tuning run

for CASE in 404 409
do
 for NAME in bomex combined dycoms2_rf01
 do

  if [ ${NAME} = 'bomex' ]; then
    EXP=(bomex)
    t1=(180)
    t2=(360)
  elif [ ${NAME} = 'combined' ]; then
    EXP=(bomex dycoms2_rf01)
    t1=(180 180)
    t2=(360 240)
  elif [ ${NAME} = 'dycoms2_rf01' ]; then
    EXP=(dycoms2_rf01)
    t1=(180)
    t2=(240)
  else
    exit 0
  fi

  NEXP=${#EXP[@]}
  RUN_CASE=${NAME}_${CASE}
  destdir=/home/cjg/calibration/revised/paper/figures_${CASE}/data

  cd ${RUN_CASE}

# loop over best NEXPS experiments and perform time averaging

  NEXPS=40
  list=`head -${NEXPS} ${RUN_CASE}_summary.dat | cut -c20-23`
  for x in $list
  do

    echo Averaging results for ${RUN_CASE}_${x}
    mkdir -p ${destdir}/${RUN_CASE}/${RUN_CASE}_${x}

    cd ${RUN_CASE}_${x}
    dmget *.nc
    for (( n=0; n < ${NEXP}; n++ ))
    do
      for file in zt zm sfc
      do
       ncra -d time,${t1[$n]},${t2[$n]} -F \
            ${EXP[$n]}_${file}.nc ${EXP[$n]}_${file}_avg.nc
       mv ${EXP[$n]}_${file}_avg.nc ${destdir}/${RUN_CASE}/${RUN_CASE}_${x}/
      done
    done
    cd ..

  done

# Copy summary file over

  cp ${RUN_CASE}_summary.dat ${destdir}/${RUN_CASE}/

  cd ..

 done
done

