#!/bin/bash

START=0
END=4

for (( i=$START; i<=$END; i++ ))
do
	mkdir case$i
	mkdir case$i/CPS_data
	mkdir case$i/CPS_pics
	sed "/case0/s/case0/case$i/" para_template.txt > ./case$i/para.txt
done

# change the parameters for each case

STEP=0.2   # this is the step for thermodynamic coupling strength

for (( i=$START; i<=$END; i++ ))
do
	# information about each case
	Lambda=$(bc <<< $(($i))*$STEP-0.4)
	echo "Heterogeneous membrane and solvent with Lambda=$Lambda" >> ./case$i/info.txt

	# change the parameters specifically in each case
	awk -v Lambda="${Lambda}" '$1=="Lambda" {$2 = Lambda}
	{print}' ./case$i/para.txt > ./case$i/temp.txt
	mv ./case$i/temp.txt ./case$i/para.txt
done







####################below this line, do not touch!!!#############################



# change the directory to be current working directory in job scheduler script
awk -v awkvar=`pwd` '$1=="cd" {$2 = awkvar} {print}' run_coupling_template > temp
mv temp run_coupling_template

# create job scheduler script for each case and store into each directory
for (( i=$START; i<=$END; i++ ))
do
	awk -v awkvar="./case$i/para.txt" '$1=="./cpl_PS" {$2 = awkvar} {print}' run_coupling_template > ./case$i/run_coupling
done

# make files
make -f Makefile.txt cpl_ps
make -f Makefile.txt clean

# submit jobs
for (( i=$START; i<=$END; i++ ))
do
	sbatch ./case$i/run_coupling
done




