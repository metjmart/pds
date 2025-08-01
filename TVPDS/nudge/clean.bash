#! /bin/bash

###--- Remove sym links to wrfinput files and wrffdda_d01 file

###--- Determine how many wrfinput sym links exist
nsl=`ls t[0-9]* | wc -l`
syms=`ls t[0-9]*`

###--- Remove the sym links
for i in $(seq 1 1 ${nsl})
do
	cur=`echo ${syms} | cut -f ${i} -d " "`
	echo "Removing ${cur}"
	rm ${cur}
done

###--- Remove the wrffdda_d01 file
echo "Removing wrffdda_d01"
rm wrffdda_d01
