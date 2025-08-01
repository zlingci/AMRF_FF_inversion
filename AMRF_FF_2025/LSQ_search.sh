#!bin/bash

workdir=`pwd`
indir=$workdir/inp
parafile=$indir/LSQ_search.txt
datadir=`cat $parafile | awk 'NR==1{print $1}'`
datafile=$workdir/$datadir
weightfile=`cat $parafile | awk 'NR==3{print $1}'`


echo  "1D inversion searching on the fault plane"
if [ ! -f ${parafile} ]
then
echo "lack of ${parafile}"
exit
fi

# Generating output files and coping src files into it
outputfiles=`awk '{if(NR == 2)print $1}' ${parafile}`
curdir=`pwd`
sed -i '15c '"${curdir}"'' $parafile

if [ ! -d $outputfiles ]
then
mkdir -p $outputfiles
fi
cp ${parafile} $outputfiles/LSQ_search.txt
cp ./src/*.m  $outputfiles/
cp $weightfile  $outputfiles/
echo "$workdir" >> $outputfiles/LSQ_search.txt

# Run inversion and plot
cd $outputfiles
########################################## 1D inversion #############################################
matlab -nodesktop -nosplash -nodisplay -r "main_LSQ_search ;quit"
#################################### calculate normal distribution ####################################
matlab -nodesktop -nosplash -nodisplay -r "norm_fitting ;quit"
######################################### Plot results ##############################################
bash ../bin/plot_LSQ_search.sh $parafile
rm *.m
cd $curdir
