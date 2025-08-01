#! /bin/bash

##############################################
# Run the WRF NetCDF nudging file generator. #
##############################################
#
#
# Sample Usage:
#               ./nudge.bash wrfinput_d01_t1 wrfinput_d01_t2 wrfinput_d01_t3 wrfinput_d01_t4 180
#
# where: wrfinput_d01_tX are the wrfinput_d01 nudging envioronmental condition files and 180 is the
#        desired time interval between the wrfinput_d01 times
#
#

start=`date`

#############
# Greetings #
#############
echo "---------------"
echo "NetNudge driver"
echo "---------------"

######################
# Define environment #
######################
export FC=ifort
export APP=NetNudge.exe
# export INPUT="wrfinput_d01_t1 wrfinput_d01_t2"
###--- Allow for 8 days at a 3 hr nudging interval (requires 82 inputs: 81 times and then the N hr interval)
export INPUT="${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10} \
               ${11} ${12} ${13} ${14} ${15} ${16} ${17} ${18} ${19} ${20} \
               ${21} ${22} ${23} ${24} ${25} ${26} ${27} ${28} ${29} ${30} \
               ${31} ${32} ${33} ${34} ${35} ${36} ${37} ${38} ${39} ${40} \
               ${41} ${42} ${43} ${44} ${45} ${46} ${47} ${48} ${49} ${50} \
               ${51} ${52} ${53} ${54} ${55} ${56} ${57} ${58} ${59} ${60} \
               ${61} ${62} ${63} ${64} ${65} ${66} ${67} ${68} ${69} ${70} \
               ${71} ${72} ${73} ${74} ${75} ${76} ${77} ${78} ${79} ${80} \
               ${81} ${82}"
export OUTPUT=wrffdda_d01

##################
# Run unit tests #
##################
./${APP} ${INPUT}
if [ -e ${OUTPUT} ];then
        echo "Output file: wrffdda_d01"
else
    	echo "Problems! No output file."
fi

#########
# Ciao! #
#########
end=`date`
echo "done!"
echo "Start time was ${start} and End time was ${end}"
exit

