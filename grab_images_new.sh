#!/bin/bash

# Copy dicom of a scan series from the scanner and convert to nifti

exam=$1
series=$2
if (( $# < 3))
then
  outname="anat"
else
  outname=$3
fi

rm $outname.nii.gz 2> /dev/null # remove existing nifti
dcmpath=$(ssh sdc@cnimr "pathExtract $exam $series | sed -n '2 p' ") # print the second line of the returned text
echo $dcmpath
seriesdir=$(dirname $dcmpath)
echo $seriesdir
mkdir tmp
scp sdc@cnimr:$seriesdir/* tmp/
dcm2niix -z y -b n -o $(pwd) -f $outname tmp/
rm -r tmp
