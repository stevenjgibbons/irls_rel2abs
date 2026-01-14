#!/bin/sh
set -x
xfile=relative_x.txt
yfile=relative_y.txt
irls_rel2abs 6  2 < $xfile | grep ABS > absolute_x.txt
irls_rel2abs 6  2 < $yfile | grep ABS > absolute_y.txt
outfile=DPRK_xy.txt
if test -r $outfile
then
  rm $outfile
fi
touch $outfile
for ev in 1 2 3 4 5 6
do
  xval=`awk '$2 == ev {print $3}' ev=$ev absolute_x.txt`
  yval=`awk '$2 == ev {print $3}' ev=$ev absolute_y.txt`
  echo "DPRK$ev  $xval  $yval" 
  echo "DPRK$ev  $xval  $yval"  >> $outfile
done
