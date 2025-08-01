#!/bin/bash


workdir=`pwd`
indir=$workdir/inp
parafile=$indir/LSQ_1D.txt
datadir=`cat $parafile | awk 'NR==1{print $1}'`
datafile=$workdir/$datadir
weightfile=`cat $parafile | awk 'NR==3{print $1}'`


echo  "1D inversion on the fault plane"
if [ ! -f ${parafile} ]
then
echo "lack of ${parafile}"
exit
fi

# Generating output files and coping src files into it
outputfiles=`awk '{if(NR == 2)print $1}' ${parafile}` 
if [ ! -d $outputfiles ]
then
mkdir -p $outputfiles
fi
curdir=`pwd`
sed -i '15c '"${curdir}"'' $parafile

cp ${parafile} $outputfiles/LSQ_1D.txt
cp ./src/*.m  $outputfiles/
cp $weightfile  $outputfiles/
echo "$workdir" >> $outputfiles/LSQ_1D.txt

# Run inversion and plot
curdir=`pwd`
cd $outputfiles
########################################## 1D inversion #############################################
matlab -nodesktop -nosplash -nodisplay -r "main_LSQ_1D ;quit"
######################################### Plot results ##############################################
cp ../bin/plot_LSQ_1D.sh .
bash plot_LSQ_1D.sh ${parafile}
rm *.m *.sh
cd ..
