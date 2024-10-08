#!/bin/bash
################
# WHERE WE ARE #
################
QSLEHOME=$(pwd)
#############
## VERSION ##
#############
echo "Please, select what to do:"
OPTIONS="build clean exit"
select opt in $OPTIONS; do
	if [ "$opt" = "build" ]; then
		CC="gcc"
		F77="gfortran"
		CPP="g++"
		COMPOPT="-O2"
		break
	elif [ "$opt" = "clean" ]; then
		echo "This operation cleans all libraries and executables leaving only the source code. Continue? [yes | no] "
		read answer
		if [ "$answer" == "yes" ]; then
			cd ${QSLEHOME}/lib/zmatlib
			make clean
			
			cd ${QSLEHOME}/lib/dite2lib
			make clean
			
			cd ${QSLEHOME}/src
			make clean

			cd ${QSLEHOME}
			rm -rf ./run_* ./qsle ./Makefile.in

			echo ""
			echo "QSLE build dir completely clean."
			echo ""
			exit
		else
			echo ""
			echo "Nothing to be done. Build script stopped."
			echo ""
			exit
		fi
	elif [ "$opt" = "exit" ]; then
		echo "Nothing to be done. Build script stopped."
		exit
	else
		echo "Bad option - please answer 1 to build, 2 to clean, or 3 to exit"
	fi
done
###########################################
# PREPARATION OF THE make.COMPILERS FILE ##
###########################################
cat > ./Makefile.in << EOF
CC      = ${CC}
CPP     = ${CPP}
F77     = ${F77}
COMPOPT = ${COMPOPT}
LIBS=-std=c++11 -lzmat -ldite2 -llapack -lblas -larpack -lsuperlu -fopenmp
LL=-L$(pwd)/lib/zmatlib/ -L$(pwd)/lib/dite2lib/
II=-I$(pwd)/lib/zmatlib/include/ -I$(pwd)/lib/armadillo-12.8.0/include/ -I$(pwd)/src/include/ -I$(pwd)/lib/dite2lib/include/

EOF
#############################
# BUILD THE ZMATLIB LIBRARY #
#############################
cd ${QSLEHOME}/lib/zmatlib/
make -j2
ZMATPATH=$(pwd)
# Check if ZMAT has been built
log=`find -name "libzmat.a"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build the ZMAT library (lib/zmatlib/). Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
#############################
# BUILD THE DITE2 LIBRARY #
#############################
cd ${QSLEHOME}/lib/dite2lib/
make -j2
# Check if ZMAT has been built
log=`find -name "libdite2.a"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build the ZMAT library (lib/dite2lib/). Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
#######################
# BUILD QSLE PARALLEL #
#######################
cd ${QSLEHOME}/src/
make
# Check if QSLE has been built
cd ${QSLEHOME}
log=`find -name "qsle"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build the QSLE stand alone (src). Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
############################
# GENERATE EXAMPLE SCRIPTS #
############################
cat > ./run_qsle.sh << EOF
#!/bin/bash

###########################################################
# This bash script can be used to run QSLE:               #
# - Step 1: copy it in your working directory;            #
# - Step 2: change NUMBER with the number of cpus to be   #
#           used;                                         #
# - Step 3: change QSLE_INPUT_FILE with the actual name   #
#           of the input file;                            #
# - Step 4: run it as a bash script.                      #
#                                                         #
# N.B.: if the QSLE directory is moved to a different to  #
# path after the program has been compiled, the paths in  #
# the following commands must be changed accordingly to   #
# the new location of the files                           #
###########################################################

export OMP_NUM_THREADS=NUMBER
export QSLEINFO=${QSLEHOME}/src/qsle_info
${QSLEHOME}/qsle QSLE_INPUT_FILE
EOF
chmod u+x ./run_qsle.sh

#######
# END #
#######
echo ""
echo "==============================================================="
echo "QSLE succesfully compiled! "
echo ""
echo "A calculation should be done in 3 steps:"
echo "1) copy the 'run_qsle.sh' bash script in your working directory"
echo "2) change the parameters in 'run_qsle.sh' following the"
echo "   instructions in the bash script"
echo "3) run the QSLE simulation with the command:"
echo "      ./run_qsle input.inp"
echo ""
echo "For further information see 'QSLE_manual.pdf' in the doc folder"
echo ""
echo "Please note that the bash scripts point to the actual position"
echo "of this directory:"
echo ${QSLEHOME}
echo "If this directory is moved to another location, the scripts "
echo "should be modified accordingly"
echo "==============================================================="
echo ""

