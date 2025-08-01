#!/bin/make
#################################################
# Makefile to build a Fortran 90/95 executable. #
#################################################

##############
# File names #
##############
DIR=.
NAME=sheargeneric
COMMAND=shg
OUTPUT=${COMMAND}.out

####################
# Fortran compiler #
####################
# FC=pgf90
FC=ifort

##########################
# Fortran compiler flags #
##########################
#FFLAGS=-O3 -fastsse -Mvect=noaltcode -Msmartalloc -Mprefetch=distance:8 -Mfprelaxed -Mfixed
FFLAGS=-O3 
# -byteswapio -Mfree -Minfo=all =Mneginfo=al

# FFLAGS=-O0 -g -check all -warn all -traceback -fp-stack-check -fpe0
# FFLAGS=-O0 -g
# FFLAGS=-O3 -r8 -xHost -ip -fp-model precise -fp-speculation=strict -fp-stack-check -fpconstant -ftz -align all -m64
# FFLAGS=-O3 -xHost -ip -fp-model fast=2 -no-heap-arrays -no-prec-div -no-prec-sqrt -ftz -align all -m64
# FFLAGS=

#######################
# Fortran source code #
#######################
SOURCES=sheargeneric.f makevtbar.f makevtofrz.f smoothvort.f hermite3.f \
	makephyd.f makepofrz.f inittofrz.f makevortex.f \
	updatetofrz.f evalgraderr.f getdataxyz.f getdatazlev.f \
	makepofyz.f makeuofyz.f makeshear.f inittofyz.f initqofyz.f \
	updatetofyz.f updateqofyz.f evalgraderryz.f getrh.f getq.f	\
	coriolis.f beta.f t2th.f t2tv.f linterp.f lextrap.f \
	interpsounding.f

#############################################
# Object code (SOURCES but with .o suffix): #
#############################################
OBJECTS=${SOURCES:.f=.o}

################
# Dependencies #
################
${NAME}:	${OBJECTS}

#############################################
# General make target for user's program -  #
# "make" with no option specified compiles, #
# links, and builds the executable program. #
#############################################
${NAME}:
	${FC} -c ${FFLAGS} ${SOURCES}
	$(FC) ${FFLAGS} -o ${COMMAND} ${OBJECTS}

######################################################
# "make run" executes the program in the background. #
######################################################
.PHONY:	run
run:	${COMMAND}
	${DIR}/${COMMAND}
#	${DIR}/${COMMAND} > ${DIR}/${OUTPUT} 2>&1 &

###########################################################
# "make mem" checks the memory while running the program. #
###########################################################
.PHONY:	mem
mem:	${COMMAND}
	valgrind --leak-check=full --dsymutil=yes ${DIR}/${COMMAND}

################################################
# "make clean" removes objects and executables #
################################################
.PHONY:	clean
clean:
	\rm -f *.o ${COMMAND}

#############################################################
# "make real clean" removes objects, executables and output #
#############################################################
.PHONY:	realclean
realclean:
	\rm -f *.o *genmod* ${COMMAND} *.out input_sounding*; \rm -rf *.dSYM
