#!/bin/bash
#Script to download Stokes Drift Data for Jan 2000 - Dec 2014 and split it into single time files
n=0

for y in {2002..2012}
do
  export year=`printf "%04d\n" $y`
  echo "getting files at year $y"
  
  for m in {1..12}
  do
    export month=`printf "%02d\n" $m`
    echo "getting files at month $m"
    wget ftp://ftp.ifremer.fr/ifremer/ww3/HINDCAST/GLOBAL/${year}_CFSR/wnd/ww3.${year}${month}_wnd.nc
	
    for k in {0..300}
    do 
      export name=`printf "%06d\n" $n`
      ncks -d time,$k,$k ww3.${year}${month}_wnd.nc WaveWatchWind${name}.nc
	if [ -e "WaveWatchWind${name}.nc" ]
	then
          n=`expr $n + 1`
        fi
    done 
  done
done

for y in {2013..2014}
do
  export year=`printf "%04d\n" $y`
  echo "getting files at year $y"

  for m in {1..12}
  do
    export month=`printf "%02d\n" $m`
    echo "getting files at month $m"
    wget ftp://ftp.ifremer.fr/ifremer/ww3/HINDCAST/GLOBAL/${year}_CFSR/wnd/glob_30m.${year}${month}_wnd.nc

    for k in {0..300}
    do
      export name=`printf "%06d\n" $n`
      ncks -d time,$k,$k glob_30m.${year}${month}_wnd.nc WaveWatchWind${name}.nc
	if [ -e "WaveWatchWind${name}.nc" ]
	then
          n=`expr $n + 1`
        fi
    done
  done
done

