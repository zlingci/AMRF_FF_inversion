#!/bin/bash
datafiles=`pwd` #`awk '{if(NR == 2)print $1}' LSQ_2D.txt`
depfile=${datafiles}/dep.txt
distfile=${datafiles}/dist.txt
ptakeofffile=${datafiles}/takeoff_p.txt
stakeofffile=${datafiles}/takeoff_s.txt
#if [ -f $ptakeofffile -a -f $stakeofffile ]
#then
#exit
#fi
tempsfile=${datafiles}/temps
temppfile=${datafiles}/tempp
modelfile=`awk '{if(NR == 1)print $1}' $depfile`
depnum=`wc $depfile | awk '{print $1}'`
distnum=`wc $distfile | awk '{print $1}'`

if [ -f $ptakeofffile ] 
then
rm $ptakeofffile
fi
if [ -f $stakeofffile ]
then
rm $stakeofffile
fi

i=1
while [ $i -le $distnum ]
do
dist=`awk -v s=$i '{if(NR==s)print $1}' $distfile`
if [ -f $tempsfile ]
then
rm $tempsfile
fi
if [ -f  $temppfile ]
then
rm $temppfile
fi
j=2
while [ $j -le $depnum ]
do
dep=`awk -v s=$j '{if(NR==s)print $1}' $depfile` 
as=`taup_time -ph s,Sg,S -h $dep -km $dist -mod $modelfile | awk 'NR==6 {print $6}'`
ap=`taup_time -ph p,Pg,P -h $dep -km $dist -mod $modelfile | awk 'NR==6 {print $6}'`
echo $as >> $tempsfile
echo $ap >> $temppfile
j=`expr $j + 1`
done
echo `cat $tempsfile` >> $stakeofffile
echo `cat $temppfile` >> $ptakeofffile
i=`expr $i + 1`
done
rm $tempsfile $temppfile
