#! /bin/bash

################################################
# Build the WRF NetCDF nudging file generator. #
################################################

######################
# Define environment #
######################
export FC=ifort
export APP=NetNudge.exe

############
# Clean up #
############
echo "Cleaning out old build"
rm -f ${APP} *.mod *.o ${OUTPUT}
sleep 2

########################
# Build new executable #
########################
echo "Building new executable"
${FC} `nc-config --fflags` netNudge.f90 -o ${APP} `nc-config --flibs`

if [ ! -e ${APP} ];then
        echo "Problems encountered during build! Bailing out."
        exit
fi

# Done
exit

