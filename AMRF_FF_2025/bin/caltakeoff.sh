#!/bin/bash

workdir=`pwd`
indir=$workdir/inp
infile=$indir/LSQ_1D.txt

datadir=`cat $infile | awk 'NR==1{print $1}'`
datafiles=$workdir/$datadir

if [ "$1"x = x ]
then
echo "Please give th row number in $model, such as: sh run_dec_ye.sh CRUST1.0"
exit
fi
modelname=$1
model=${indir}/${modelname}.nd

echo egfname station cha dist azimuth takeoff weight > weightall.txt
# comp=("zp" "rp" "zs" "rs" "ts")
comp[0]='zp'
comp[1]='rp'
comp[2]='zs'
comp[3]='rs'
comp[4]='ts'
# channel=("Z" "R" "Z" "R" "T")
channel[0]='Z'
channel[1]='R'
channel[2]='Z'
channel[3]='R'
channel[4]='T'
# wave=("p,Pg,Pn,P" "p,Pg,Pn,P" "s,Sg,S" "s,Sg,S" "s,Sg,S")
wave[0]='p,Pg,Pn,P'
wave[1]='p,Pg,Pn,P'
wave[2]='s,Sg,S'
wave[3]='s,Sg,S'
wave[4]='s,Sg,S'
for ((i_comp=0;i_comp<5;i_comp+=1))
do
############################################ cal takeoff ##############################################
stationf=${comp[$i_comp]}.txt
cat $datafiles/*.${comp[$i_comp]}.txt > $stationf
snum=`wc $stationf | awk '{print $1}'`
i=1
while [ $i -le $snum ]
do
egfname=`awk -v s=$i '{if(NR == s)print $1}' $stationf`
sname=`awk -v s=$i '{if(NR == s)print $2}' $stationf`
cha=`awk -v s=$i '{if(NR == s)print $3}' $stationf`
datafile=$datafiles/$egfname.${sname}.${cha}.${comp[$i_comp]}.out
if [ ! -f $datafile ]
then
datafile=$datafiles/$egfname.${cha}.${comp[$i_comp]}.out
if [ ! -f $datafile ]
then
i=` expr $i + 1`
continue
fi
fi
sta=`saclst kstnm f $datafile | awk '{print $2}'`
evdp=`saclst evdp f $datafile | awk '{print $2}'`
dist=`saclst dist f $datafile | awk '{print $2}'`
az=`saclst az f $datafile | awk '{print $2}'`
ap=`taup_time -ph ${wave[$i_comp]} -h $evdp -km $dist -mod $model | awk 'NR==6 {print $6}'`
echo $egfname $sta ${comp[$i_comp]} $dist $az $ap "1.0" >> weightall.txt
i=` expr $i + 1`
done
rm $stationf
done

