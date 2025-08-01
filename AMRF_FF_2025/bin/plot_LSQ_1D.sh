#!/bin/bash

parfile=$1
datadir=`pwd`

AMRFfile=${datadir}/AMRF.txt
AMRF_estfile=${datadir}/AMRF_est.txt
moment_AMRFfile=${datadir}/moment_AMRF.txt
moment_AMRF_estfile=${datadir}/moment_AMRF_est.txt
MRFd_estfile=${datadir}/MRFd_est.txt
D0_estfile=${datadir}/D0_est.txt
VRfile=${datadir}/VR.txt
parameterfile=${datadir}/parameter.txt

gmt gmtset ANNOT_FONT_SIZE_PRIMARY 10p
gmt gmtset LABEL_FONT_SIZE 12p
gmt gmtset BASEMAP_TYPE plain
gmt gmtset FRAME_PEN 0.75p,black
gmt gmtset TICK_LENGTH -5p
gmt gmtset ANNOT_OFFSET_PRIMARY 3p
gmt gmtset LABEL_OFFSET 3p

PS=${datadir}/temp.ps
PS1=${datadir}/temp_zp.ps
PS2=${datadir}/temp_rp.ps
PS3=${datadir}/temp_zs.ps
PS4=${datadir}/temp_rs.ps
PS5=${datadir}/temp_ts.ps
######################################################################################
############################### plot MRFd_est #######################################
#####################################################################################
Vr=`awk '{if(NR == 9)print $1}' $parfile`
ts=`awk '{print $1}' $MRFd_estfile|sort -n|tail -1`
L_start=`awk '{if(NR == 1)print $1}' $D0_estfile`
L_end=`awk '{print $1}' $D0_estfile|tail -1`
L2=`awk '{if(NR == 2)print $1}' $D0_estfile`
dL=`bc << EOF
scale=6
(${L2} - ${L_start})
EOF
`
y_start=`bc << EOF
sacle=6
(${L_start} - ${dL})
EOF
`
y_end=`bc << EOF
sacle=6
(${L_end} + 1.2*${dL})
EOF
`
x_start=0
x_end=`bc << EOF
scale=6
(${L_end}/${Vr} + $ts)
EOF
`
J1=X7.0c/18c
R1=$x_start/$x_end/$y_start/$y_end
Xa1=-1c
Ya1=-1.3c
gmt psxy -R$R1 -J$J1 -T -K > $PS
max_MRFd=`awk '{$1=""; print $0}' $MRFd_estfile|xargs -n1|sort -n|tail -1`
snum=`awk '{print NF}' $MRFd_estfile | awk '{if(NR == 1)print}'`
i=$snum
while [ $i -ge 2 ]
do
k=`expr $i - 1`
L=`awk -v s=$k '{if(NR == s)print $1}' $D0_estfile`
L1=`awk -v s=$k '{if(NR == s && $1 < 0)print $1*-1; else if(NR == s && $1 >= 0)print $1}' $D0_estfile`
time_start=`bc << EOF
scale=6
($L1/${Vr})
EOF
`
time_end=`bc << EOF
scale=6
(${time_start} + ${ts})
EOF
`
gmt psxy -R$R1 -J$J1 -Xa$Xa1 -Ya$Ya1 -W0.5p,black,solid -O -K >> $PS << EOF
0 $L
$time_start $L
EOF
gmt psxy -R$R1 -J$J1 -Xa$Xa1 -Ya$Ya1 -W0.5p,black,solid -O -K >> $PS << EOF
$time_end $L
$x_end $L
EOF
awk -v s=$i -v n=$max_MRFd -v b=$dL -v m=$Vr -v o=$L -v o1=$L1 '{print $1+o1/m,$s*2*b/n+o}' $MRFd_estfile | gmt psxy -R$R1 -J$J1 -Xa$Xa1 -Ya$Ya1 -W0.5p,red,solid -G150/150/150 -O -K >> $PS
i=`expr $i - 1`
done
gmt psbasemap -R$R1 -J$J1 -Xa$Xa1 -Ya$Ya1 -Ba1f0.5:"Time (s)":/a0.4f0.2:"Distance along rupture direction (km)":SWne -K -O >> $PS
######################################################################################
################ plot AMRF, moment_AMRF and AMRF_est, moment_AMRF_est ################
######################################################################################
cp $PS $PS1
cp $PS $PS2
cp $PS $PS3
cp $PS $PS4
cp $PS $PS5
duration=`awk '{print $1}' $AMRF_estfile|sort -n|tail -1`
x_start=`bc << EOF
scale=3
($duration*0.1*-1)
EOF
`
x_end=`bc << EOF
scale=3
($duration*1.1)
EOF
`
J2=X7.0c/18c
R2=$x_start/$x_end/-5/420
Xa2=7.7c
Ya2=-1.3c
J3=X7.0c/18c
R3=0/$duration/-5/460
Xa3=15.2c
Ya3=-1.3c
max_AMRF_moment=`awk '{$1=""; print $0}' $moment_AMRFfile|xargs -n1|sort -n|tail -1`
max_AMRF=`awk '{$1=""; print $0}' $AMRFfile|xargs -n1|sort -n|tail -1`
snum=`awk '{print NF}' $AMRFfile | awk '{if(NR == 1)print}'`
i=2
while [ $i -le $snum ]
do
k=`expr $i - 1`
cha=`awk -v s=$k '{if(NR == s)print $3}' $VRfile`
vr=`awk -v s=$k '{if(NR == s)print $7}' $VRfile`
az=`awk -v s=$k '{if(NR == s)print $4}' $VRfile`
az=`printf '%-.2f' $az`
vr=`printf '%-.1f' $vr`
if [ "$cha" = "zp" ]
then
PStemp=$PS1 
state1=a
elif [ "$cha" = "rp" ]
then
PStemp=$PS2
state2=a
elif [ "$cha" = "zs" ]
then
PStemp=$PS3
state3=a
elif [ "$cha" = "rs" ]
then
PStemp=$PS4
state4=a
else
PStemp=$PS5
state5=a
fi
gmt psxy -R$R2 -J$J2 -Xa$Xa2 -Ya$Ya2 -W0.5p,red,solid -O -K >> $PStemp << EOF
$x_start $az
0 $az
EOF
gmt psxy -R$R2 -J$J2 -Xa$Xa2 -Ya$Ya2 -W0.5p,red,solid -O -K >> $PStemp << EOF
$duration $az
$x_end $az
EOF
awk -v s=$i -v n=$max_AMRF -v m=$az '{print $1,$s*80/n+m}' $AMRFfile | gmt psxy -R$R2 -J$J2 -Xa$Xa2 -Ya$Ya2 -W0.5p,black,solid -O -K >> $PStemp
awk -v s=$i -v n=$max_AMRF -v m=$az '{print $1,$s*80/n+m}' $AMRF_estfile | gmt psxy -R$R2 -J$J2 -Xa$Xa2 -Ya$Ya2 -W0.5p,red,solid -O -K >> $PStemp
gmt pstext -R0/1/-5/420 -J$J2 -Xa$Xa2 -Ya$Ya2 -N -K -O >> $PStemp << EOF
0.85 $az 10 0 5 ML @;0/0/0;${vr}%@;;
EOF
awk -v s=$i -v n=$max_AMRF_moment -v m=$az '{print $1,$s*80/n+m}' $moment_AMRFfile | gmt psxy -R$R3 -J$J3 -Xa$Xa3 -Ya$Ya3 -W0.5p,black,solid -O -K >> $PStemp
awk -v s=$i -v n=$max_AMRF_moment -v m=$az '{print $1,$s*80/n+m}' $moment_AMRF_estfile | gmt psxy -R$R3 -J$J3 -Xa$Xa3 -Ya$Ya3 -W0.5p,red,solid -O -K >> $PStemp
i=` expr $i + 1`
done
VR=`awk '{if(NR == 1)print $2}' $parameterfile`
if [ "$state1"x = ax ]
then
gmt pstext -R0/1/0/1 -J$J2 -Xa$Xa2 -Ya$Ya2 -N -K -O >> $PS1 << EOF
0.40 0.96 10 0 5 ML @;255/0/0;VR ${VR}%@;;
EOF
gmt psbasemap -R$R2 -J$J2 -Xa$Xa2 -Ya$Ya2 -Ba1f0.5:"Time (s)":/a60f30:"Azimuth (degree)":SWne -K -O >> $PS1
gmt psbasemap -R$R3 -J$J3 -Xa$Xa3 -Ya$Ya3 -Ba1f0.5:"Time (s)":/a60f30:"Azimuth (degree)":Swne -K -O >> $PS1
gmt psxy -R$R1 -J$J1 -O -T >> $PS1
cat $PS1 >> ${datadir}/all.ps
fi
if [ "$state2"x = ax ]
then
gmt pstext -R0/1/0/1 -J$J2 -Xa$Xa2 -Ya$Ya2 -N -K -O >> $PS2 << EOF
0.40 0.96 10 0 5 ML @;255/0/0;VR ${VR}%@;;
EOF
gmt psbasemap -R$R2 -J$J2 -Xa$Xa2 -Ya$Ya2 -Ba1f0.5:"Time (s)":/a60f30:"Azimuth (degree)":SWne -K -O >> $PS2
gmt psbasemap -R$R3 -J$J3 -Xa$Xa3 -Ya$Ya3 -Ba1f0.5:"Time (s)":/a60f30:"Azimuth (degree)":Swne -K -O >> $PS2
gmt psxy -R$R1 -J$J1 -O -T >> $PS2
cat $PS2 >> ${datadir}/all.ps
fi
if [ "$state3"x = ax ]
then
gmt pstext -R0/1/0/1 -J$J2 -Xa$Xa2 -Ya$Ya2 -N -K -O >> $PS3 << EOF
0.40 0.96 10 0 5 ML @;255/0/0;VR ${VR}%@;;
EOF
gmt psbasemap -R$R2 -J$J2 -Xa$Xa2 -Ya$Ya2 -Ba1f0.5:"Time (s)":/a60f30:"Azimuth (degree)":SWne -K -O >> $PS3
gmt psbasemap -R$R3 -J$J3 -Xa$Xa3 -Ya$Ya3 -Ba1f0.5:"Time (s)":/a60f30:"Azimuth (degree)":Swne -K -O >> $PS3
gmt psxy -R$R1 -J$J1 -O -T >> $PS3
cat $PS3 >> ${datadir}/all.ps
fi
if [ "$state4"x = ax ]
then
gmt pstext -R0/1/0/1 -J$J2 -Xa$Xa2 -Ya$Ya2 -N -K -O >> $PS4 << EOF
0.40 0.96 10 0 5 ML @;255/0/0;VR ${VR}%@;;
EOF
gmt psbasemap -R$R2 -J$J2 -Xa$Xa2 -Ya$Ya2 -Ba1f0.5:"Time (s)":/a60f30:"Azimuth (degree)":SWne -K -O >> $PS4
gmt psbasemap -R$R3 -J$J3 -Xa$Xa3 -Ya$Ya3 -Ba1f0.5:"Time (s)":/a60f30:"Azimuth (degree)":Swne -K -O >> $PS4
gmt psxy -R$R1 -J$J1 -O -T >> $PS4
cat $PS4 >> ${datadir}/all.ps
fi
if [ "$state5"x = ax ]
then
gmt pstext -R0/1/0/1 -J$J2 -Xa$Xa2 -Ya$Ya2 -N -K -O >> $PS5 << EOF
0.40 0.96 10 0 5 ML @;255/0/0;VR ${VR}%@;;
EOF
gmt psbasemap -R$R2 -J$J2 -Xa$Xa2 -Ya$Ya2 -Ba1f0.5:"Time (s)":/a60f30:"Azimuth (degree)":SWne -K -O >> $PS5
gmt psbasemap -R$R3 -J$J3 -Xa$Xa3 -Ya$Ya3 -Ba1f0.5:"Time (s)":/a60f30:"Azimuth (degree)":Swne -K -O >> $PS5
gmt psxy -R$R1 -J$J1 -O -T >> $PS5
cat $PS5 >> ${datadir}/all.ps
fi
#####################################################################################################
figuredir=../figure/
echo ${datadir}
if [ ! -d $figuredir ]
then
mkdir $figuredir
fi
ps2pdf ${datadir}/all.ps $figuredir/plot_1D.pdf
evince $figuredir/plot_1D.pdf &
rm gmt.* $PS $PS1 $PS2 $PS3 $PS4 $PS5 ${datadir}/all.ps
