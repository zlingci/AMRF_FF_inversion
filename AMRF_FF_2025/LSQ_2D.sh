#!/bin/bash

workdir=`pwd`
indir=$workdir/inp
parafile=$indir/LSQ_2D.txt
datadir=`cat $parafile | awk 'NR==1{print $1}'`
datafile=$workdir/$datadir
weightfile=`cat $parafile | awk 'NR==3{print $1}'`
modelfile=`cat $parafile | awk 'NR==4{print $1}'`

echo  "2D inversion on the fault plane"
if [ ! -f ${parafile} ]
then
echo "lack of ${parafile} "
exit
fi


# Generating output files and coping src files into it
outputfiles=`awk '{if(NR == 2)print $1}' ${parafile}`
if [ ! -d $outputfiles ]
then
mkdir -p $outputfiles
fi
curdir=`pwd`
sed -i '17c '"${curdir}"'' $parafile
cp ${parafile} $outputfiles/LSQ_2D.txt
cp ./src/*.m  $outputfiles/
cp $weightfile  $outputfiles/
cp ./inp/*${modelfile}* $datafile/
echo "$workdir" >> $outputfiles/LSQ_2D.txt

cd $outputfiles
########################################## 2D inversion #############################################
matlab -nodesktop -nosplash -nodisplay -r "main_LSQ_2D ;quit"
######################################### Plot results ##############################################
cp ../bin/plot_LSQ_2D.sh .
bash plot_LSQ_2D.sh ${parafile}
rm *.m *.sh
cd ..
