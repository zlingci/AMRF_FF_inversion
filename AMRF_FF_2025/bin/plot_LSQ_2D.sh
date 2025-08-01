#!/bin/bash

parfile=$1
datadir=`pwd`

MRFfile=${datadir}/AMRF.txt
MRF_estfile=${datadir}/AMRF_est.txt
moment_MRFfile=${datadir}/moment_AMRF.txt
moment_MRF_estfile=${datadir}/moment_AMRF_est.txt
MRF_est_normalfile=${datadir}/AMRF_est_normal.txt
D0_estfile=${datadir}/D0_est.txt
slip_estfile=${datadir}/time_slip.txt
VRfile=${datadir}/VR.txt
parameterfile=${datadir}/parameter.txt
takeofffile1=${datadir}/takeoff_s.txt
takeofffile2=${datadir}/takeoff_p.txt

gmt gmtset ANNOT_FONT_SIZE_PRIMARY 10p
gmt gmtset LABEL_FONT_SIZE 12p
gmt gmtset BASEMAP_TYPE plain
gmt gmtset FRAME_PEN 0.75p,black
gmt gmtset TICK_LENGTH -5p
gmt gmtset ANNOT_OFFSET_PRIMARY 3p
gmt gmtset LABEL_OFFSET 3p
gmt gmtset COLOR_FOREGROUND 115/0/0
gmt gmtset COLOR_BACKGROUND blue 
PS=${datadir}/temp.ps
PS1=${datadir}/temp_p.ps
PS2=${datadir}/temp_s.ps
######################################################################################
############################ plot D0 and MRFd_est ###################################
#####################################################################################
echo "plot coseimic slip"
VR=`awk '{if(NR == 1)print $2}' $parameterfile`
VR=`printf '%-.1f' $VR`
dep_centroid=`awk '{if(NR == 5)print $2}' $parameterfile`
Vr=`awk -v s=9 '{if(NR == s)print $1}' $parfile`
source_dep=`awk -v s=11 '{if(NR == s)print $1}' $parfile`
dip=`awk -v s=8 '{if(NR == s)print $2}' $parfile`
strike=`awk -v s=8 '{if(NR == s)print $1}' $parfile`
rake=`awk -v s=8 '{if(NR == s)print $3}' $parfile`
dx=`awk -v s=12 '{if(NR == s)print $3}' $parfile`
dy=`awk -v s=13 '{if(NR == s)print $3}' $parfile`
dt=`awk -v s=6 '{if(NR == s)print $2}' $parfile`
tri_num=`awk -v s=14 '{if(NR == s)print $1}' $parfile`
tri_point=`awk -v s=14 '{if(NR == s)print $2}' $parfile`
x_start=`awk -v s=12 -v m=$dx '{if(NR == s)print $1-m/2}' $parfile`
x_end=`awk -v s=12 -v m=$dx '{if(NR == s)print $2+m/2}' $parfile`
x_start_1=`awk -v s=12 '{if(NR == s)print $1}' $parfile`
x_end_1=`awk -v s=12 '{if(NR == s)print $2}' $parfile`
y_start1=`awk -v s=13 -v m=$dy '{if(NR == s)print $1-m/2}' $parfile`
y_end1=`awk -v s=13 -v m=$dy '{if(NR == s)print $2+m/2}' $parfile`
y_start1_1=`awk -v s=13 '{if(NR == s)print $1}' $parfile`
y_end1_1=`awk -v s=13 '{if(NR == s)print $2}' $parfile`
y_start2=`awk -v s=13 '{if(NR == s)print $1}' $parfile`
y_start2=`bc -l << EOF
scale=3
(($y_start2-$dy/2)*s($dip*3.1415/180)+${source_dep})
EOF
`
y_end2=`awk -v s=13 '{if(NR == s)print $2}' $parfile`
y_end2=`bc -l << EOF
scale=3
(($y_end2+$dy/2)*s($dip*3.1415/180)+${source_dep})
EOF
`
Jy=9 #5.8
ratio=`bc -l << EOF
scale=3
(($y_end1 - $y_start1)*s($dip*3.1415/180)/($x_end - $x_start))
EOF
`
Jx=`bc <<EOF
scale=3
($Jy/$ratio)
EOF
`
if [ $(echo "$Jx > 11" | bc) = 1 ]
then
Jx=11
Jy=`bc <<EOF
scale=3
($Jx*$ratio)
EOF
`
fi
J1=X11c/-9c
Xa1=-1.0c
Ya1=8.5c
R1=$x_start_1/$x_end_1/$y_start1_1/$y_end1_1
R1_1=$x_start/$x_end/$y_start1/$y_end1
R1_2=$x_start/$x_end/$y_start2/$y_end2
#max_D0=2.2
max_D0=`awk '{if(NR >= 4)print $5}' $D0_estfile|sort -g|tail -1|awk '{print $1*1.001}'`
echo ${max_D0}
dD0=`bc << EOF
scale=14
($max_D0/200)
EOF
`
min_D0=`bc << EOF
scale=14
(-1*$max_D0/10)
EOF
`
min_D0=0.0
cpts=../bin/sss.cpt
gmt makecpt -C$cpts -T$min_D0/$max_D0/$dD0 > ${datadir}/temp.cpt

awk '{if(NR >= 4)print $1,$2,$5}' $D0_estfile | gmt xyz2grd -R$R1 -G${datadir}/temp.grd -I$dx/$dy
gmt psxy -R$R1_1 -J$J1 -T -K > $PS
gmt grdimage -R$R1_1 -J$J1 ${datadir}/temp.grd -Xa$Xa1 -Ya$Ya1 -C${datadir}/temp.cpt -K -O >> $PS
w=`bc << EOF
scale=3
($Jx*0.6)
EOF
`
y=`bc << EOF
scale=3
($Jy + 1.0)
EOF
`
x=`bc << EOF
scale=3
($Jy*2)
EOF
`
cptlocx=$(echo "$Jx/2" | bc -l)
gmt psscale -Xa$Xa1 -Ya$Ya1 -D${cptlocx}c/-1.1c/${w}c/0.2ch -E -C${datadir}/temp.cpt -Ba0.4f0.2/:"Coseismic slip (m)": -K -O >> $PS
J11=X5.3c/${Jy}c
Xa11=11c
gmt pstext -R0/1/0/1 -J${J11} -Xa$Xa11 -Ya$Ya1 -K -O >> $PS << EOF
0.15 0.95 10 0 5 BL Vr = ${Vr} km/s; VR = ${VR}%
0.15 0.85 10 0 5 BL Strike = ${strike}\260
0.15 0.75 10 0 5 BL Dip = ${dip}\260; Rake = ${rake}\260
0.15 0.65 10 0 5 BL H@-0@- = $source_dep km; H@-c@- = $dep_centroid km
EOF
######################################## plot MRFd_est ############################################
snum=`awk '{if(NR == 3)print $1}' ${slip_estfile}`
i=1
k=4
while [ $i -le $snum ]
do
subnum=`awk -v s=$k '{if(NR == s)print $9}' ${slip_estfile}`
ss=`expr $k + 1`
ee=`expr $ss + $subnum - 1`
if [ $i -eq 1 ]
then
awk -v s=$ss -v e=$ee '{if(NR>=s && NR<=e)print $3}' ${slip_estfile} > temp.txt
else
awk -v s=$ss -v e=$ee '{if(NR>=s && NR<=e)print $3}' ${slip_estfile} >> temp.txt
fi
k=`expr $k + $subnum + 1`
i=`expr $i + 1`
done
max_subMRF=`awk '{printf "%f\n", $1}' temp.txt| sort -n|tail -1`
ts=`bc << EOF
scale=6
(($tri_num+1)*($tri_point-1)*$dt/2)
EOF
`
i=1
k=4
while [ $i -le $snum ]
do
subnum=`awk -v s=$k '{if(NR == s)print $9}' ${slip_estfile}`
x_shift=`awk -v s=$k -v n=$dx '{if(NR == s)print $2-n/2}' ${slip_estfile}`
y_shift=`awk -v s=$k -v n=$dy '{if(NR == s)print $3+n/2}' ${slip_estfile}`
ss=`expr $k + 1`
ee=`expr $ss + $subnum - 1`
awk -v e=$ee -v s=$ss -v n11=$x_shift -v n12=$dx -v n21=$y_shift -v n22=$dy -v m2=$max_subMRF -v m1=$ts '{if(NR>=s && NR<=e)print $1*n12/m1+n11,$3*n22*-1/m2+n21}' ${slip_estfile} | gmt psxy -R$R1_1 -J$J1 -Xa$Xa1 -Ya$Ya1 -W0.01p,black,solid+s -G150/150/150 -O -K >> $PS
k=`expr $k + $subnum + 1`
i=`expr $i + 1`
done
gmt psxy -R$R1_1 -J$J1 -Xa$Xa1 -Ya$Ya1 -Sa0.4c -Gred -W0.5p,black -O -K >> $PS << EOF
0 0
EOF
gmt psbasemap -R$R1_2 -J$J1 -Xa$Xa1 -Ya$Ya1 -B/a2f1:"Depth (km)":E -K -O >> $PS
gmt psbasemap -R$R1_1 -J$J1 -Xa$Xa1 -Ya$Ya1 -Ba2f1:"Distance along strike (km)":/a2f1:"Distance along dip (km)":SWn -K -O >> $PS

######################################################################################
################################## plot MRF_normal_est ##############################
######################################################################################
echo "plot MRF"
moment1=`awk '{if(NR == 3)print $2}' $parameterfile | awk -F 'e' '{print $1}'`
moment2=`awk '{if(NR == 3)print $2}' $parameterfile | awk -F 'e' '{print $2+0}'`
moment_1=`awk '{if(NR == 4)print $2}' $parameterfile | awk -F 'e' '{print $2-0}'`
Centriod_time=`awk '{if(NR == 8)print $2}' $parameterfile`
Duration=`awk '{if(NR == 14)print $2}' $parameterfile`
Centriod_time=`printf "%-.2f" $Centriod_time`
Duration=`printf "%-.2f" $Duration`
Mw=`bc -l << EOF
scale=4
((${moment2}+7+l($moment1)/l(10))*2/3 - 10.73)
EOF
`
Mw=`printf "%-.2f" $Mw`
x_start=`awk '{print $1}' $MRF_est_normalfile|sort -n|awk '{if(NR == 1)print}'`
x_end=`awk '{print $1}' $MRF_est_normalfile|sort -n|tail -1`
y_start=0
y_end=`awk '{printf "%f\n", $2}' $MRF_est_normalfile|sort -n|tail -1 | awk -v s=${moment_1}  '{print $1*1.15*10^(-s)}'`
J3=X11c/7c
R3=$x_start/$x_end/$y_start/$y_end
Xa3=-1.2c
Ya3=-1.6c
awk -v s=${moment_1} '{print $1,$2*10^(-s)}' $MRF_est_normalfile | gmt psxy -R$R3 -J$J3 -Xa$Xa3 -Ya$Ya3 -W1.0p,black,solid+s -G220/220/220 -O -K >> $PS
y1=$y_start
y2=`bc << EOF
scale=6
($y_end*0.2)
EOF
`
gmt psxy -R$R3 -J$J3 -Xa$Xa3 -Ya$Ya3 -W1.5p,red,solid -O -K >> $PS << EOF 
$Centriod_time $y1 
$Centriod_time $y2
EOF
gmt psxy -R$R3 -J$J3 -Xa$Xa3 -Ya$Ya3 -W1.5p,blue,solid -O -K >> $PS << EOF 
$Duration $y1 
$Duration $y2
EOF
gmt pstext -R0/1/0/1 -J$J3 -Xa$Xa3 -Ya$Ya3 -K -O >> $PS << EOF
0.52 0.90 10 0 5 BL M@-0@- = $moment1 x 10@+$moment2@+
0.52 0.80 10 0 5 BL M@-w@- = $Mw
0.52 0.70 10 0 5 BL Duration = $Duration s
0.52 0.60 10 0 5 BL T@-c@- = ${Centriod_time} s
EOF
gmt psbasemap -R$R3 -J$J3 -Xa$Xa3 -Ya$Ya3 -Ba4f2:"Time (s)":/a20f10:"MRF (x10@+${moment_1}@+ Nm/s)":SWne -K -O >> $PS

######################################################################################
################ plot AMRF and AMRF_est ################
######################################################################################
echo "plot AMRF"
cp $PS $PS1
cp $PS $PS2
plnum=29
plnum1=`bc << EOF
scale=2
($plnum + 0.5)
EOF
`
duration=`awk -v s=7 '{if(NR == s)print $1}' $parfile`
x_start=`bc << EOF
scale=3
($duration*0.15*-1)
EOF
`
x_end=`bc << EOF
scale=3
($duration*1.00)
EOF
`
J5=X4.7c/15.c
R5=$x_start/$x_end/-0.30/0.38
Xa5=12.8c
Ya5=-1.6c
J6=X4.9c/19.1c
R6=$x_start/$x_end/-0.5/$plnum1
Xa61=19.1c
Ya6=-1.6c
R7=0/1/-0.5/$plnum1
J7=X1.6c/19.1c
Xa71=17.5c
Ya7=-1.6c
snum=`awk '{if(NR == 3)print $1}' $MRFfile`
i=1
k=4
while [ $i -le $snum ]
do
subnum=`awk -v s=$k '{if(NR == s)print $5}' $MRFfile`
ss=`expr $k + 1`
ee=`expr $ss + $subnum - 1`
if [ $i -eq 1 ]
then
awk -v s=$ss -v e=$ee '{if(NR>=s && NR<=e)print $2}' $MRFfile > temp1.txt
else
awk -v s=$ss -v e=$ee '{if(NR>=s && NR<=e)print $2}' $MRFfile >> temp1.txt
fi
k=`expr $k + $subnum + 1`
i=`expr $i + 1`
done
max_AMRF=`awk '{printf "%f\n", $1}' temp1.txt|sort -n |tail -1`
i=1
ns=0
np=0
while [ $i -le $snum ]
do
k=$i
cha=`sort -k4 -n $VRfile | awk -v s=$k '{if(NR == s)print $3}'`
egfname=`sort -k4 -n $VRfile | awk -v s=$k '{if(NR == s)print $1}'`
sta=`sort -k4 -n $VRfile | awk -v s=$k '{if(NR == s)print $2}'`
vr=`sort -k4 -n $VRfile | awk -v s=$k '{if(NR == s)print $6}'`
az=`sort -k4 -n $VRfile | awk -v s=$k '{if(NR == s)print $4}'`
dpr=`sort -k4 -n $VRfile | awk -v s=$k '{if(NR == s)print $7}'`
index=`grep -n "$egfname $sta $cha" $VRfile | awk -F ':' '{print $1}'`
echo $dpr
az=`printf '%-.2f' $az`
vr=`printf '%-.1f' $vr`
if [ "$cha" = "zs" -o "$cha" = "rs" -o "$cha" = "ts" ]
then
PStemp=$PS1
ns=`expr $ns + 1`
state=1
else
PStemp=$PS2
np=`expr $np + 1`
state=0
fi
gmt psxy -R$R5 -J$J5 -Xa$Xa5 -Ya$Ya5 -W0.5p,red,solid -O -K >> $PStemp << EOF
$x_start $dpr
0 $dpr
EOF
gmt psxy -R$R5 -J$J5 -Xa$Xa5 -Ya$Ya5 -W0.5p,red,solid -O -K >> $PStemp << EOF
$duration $dpr
$x_end $dpr
EOF
j=1
m=4
while [ $j -le $index ]
do
subnum=`awk -v s=$m '{if(NR == s)print $5}' $MRFfile`
m=`expr $m + $subnum + 1`
j=`expr $j + 1`
done
ee=`expr $m - 1`
ss=`expr $m - $subnum`
awk -v e=$ee -v s=$ss -v n=$max_AMRF -v m=$dpr '{if(NR>=s && NR<=e)print $1,$2*0.10/n+m}' $MRFfile | gmt psxy -R$R5 -J$J5 -Xa$Xa5 -Ya$Ya5 -W0.5p,black,solid -O -K >> $PStemp
awk -v e=$ee -v s=$ss -v n=$max_AMRF -v m=$dpr '{if(NR>=s && NR<=e)print $1,$2*0.10/n+m}' $MRF_estfile | gmt psxy -R$R5 -J$J5 -Xa$Xa5 -Ya$Ya5 -W0.5p,red,solid -O -K >> $PStemp
####################################################################################################
if [ $state -eq 0 ]
then
if [ $np -le $plnum ]
then
k=`expr $np - 1`
Xa6=$Xa61
Xa7=$Xa71
elif [ $(echo "$np <= 2*$plnum" | bc) = 1 ]
then
k=`expr $np - $plnum - 1`
Xa6=$Xa62
Xa7=$Xa72
else
k=`expr $np - $plnum - $plnum  - 1`
Xa6=$Xa63
Xa7=$Xa73
fi
else
if [ $ns -le $plnum ]
then
k=`expr $ns - 1`
Xa6=$Xa61
Xa7=$Xa71
elif [ $(echo "$ns <= 2*$plnum" | bc) = 1 ]
then
k=`expr $ns - $plnum - 1`
Xa6=$Xa62
Xa7=$Xa72
else
k=`expr $ns - $plnum - $plnum - 1`
Xa6=$Xa63
Xa7=$Xa73
fi
fi
gmt psxy -R$R6 -J$J6 -Xa$Xa6 -Ya$Ya6 -W0.5p,red,solid -O -K >> $PStemp << EOF
$x_start $k
0 $k
EOF
gmt psxy -R$R6 -J$J6 -Xa$Xa6 -Ya$Ya6 -W0.5p,red,solid -O -K >> $PStemp << EOF
$duration $k
$x_end $k
EOF
max_AMRF1=`awk -v s=$ss -v e=$ee '{if(NR>=s && NR<=e)printf "%f\n",$2}' $MRFfile|sort -n|tail -1`
awk -v s=$ss -v e=$ee -v n=$max_AMRF1 -v m=$k '{if(NR>=s && NR<=e)print $1,$2/n+m}' $MRFfile | gmt psxy -R$R6 -J$J6 -Xa$Xa6 -Ya$Ya6 -W0.5p,black,solid -O -K >> $PStemp
awk -v s=$ss -v e=$ee -v n=$max_AMRF1 -v m=$k '{if(NR>=s && NR<=e)print $1,$2/n+m}' $MRF_estfile | gmt psxy -R$R6 -J$J6 -Xa$Xa6 -Ya$Ya6 -W0.5p,red,solid -O -K >> $PStemp
k1=`bc << EOF
scale=3
($k + 0.4)
EOF
`
k2=`bc << EOF
scale=3
($k + 0.45)
EOF
`
gmt pstext -R0/1/-0.5/$plnum1 -J$J6 -Xa$Xa6 -Ya$Ya6 -N -K -O >> $PStemp << EOF
0.81 $k1 8 0 5 ML  ${vr}%
EOF
az=`printf '%-.1f' $az`
gmt pstext -R0/1/-0.5/$plnum1 -J$J7 -Xa$Xa7 -Ya$Ya7 -N -K -O >> $PStemp << EOF
0.95 $k 8 0 5 MR ${sta}
EOF
gmt pstext -R0/1/-0.5/$plnum1 -J$J7 -Xa$Xa7 -Ya$Ya7 -N -K -O >> $PStemp << EOF
0.95 $k2 8 0 5 MR ${az}\260
EOF
i=`expr $i + 1`
done
if [ $ns -ge 1 ]
then
gmt pstext -R0/1/0/1 -J$J5 -Xa$Xa5 -Ya$Ya5 -N -K -O >> $PS1 << EOF
0.55 0.96 10 0 5 BL VR ${VR}% (S)
EOF
gmt psbasemap -R$R5 -J$J5 -Xa$Xa5 -Ya$Ya5 -Ba4f2:"Time (s)":/a0.1f0.02:"Directivity parameter (${strike}\260)":SWne -K -O >> $PS1
gmt psbasemap -R$R6 -J$J6 -Xa$Xa61 -Ya$Ya6 -Ba4f2:"Time (s)":/a100Swen -K -O >> $PS1
if [ $ns -gt $plnum ]
then
gmt psbasemap -R$R6 -J$J6 -Xa$Xa62 -Ya$Ya6 -Ba4f2:"Time (s)":/a100Swen -K -O >> $PS1
fi
if [ $(echo "$ns > 2*$plnum" | bc) = 1 ]
then
gmt psbasemap -R$R6 -J$J6 -Xa$Xa63 -Ya$Ya6 -Ba4f2:"Time (s)":/a100Swen -K -O >> $PS1
fi
gmt psxy -R$R1_1 -J$J1 -O -T >> $PS1
cat $PS1 >> ${datadir}/all.ps
fi
if [ $np -ge 1 ]
then
gmt pstext -R0/1/0/1 -J$J5 -Xa$Xa5 -Ya$Ya5 -N -K -O >> $PS2 << EOF
0.55 0.96 10 0 5 BL @;255/0/0;VR ${VR}% (P)@;;
EOF
gmt psbasemap -R$R5 -J$J5 -Xa$Xa5 -Ya$Ya5 -Ba4f2:"Time (s)":/a0.1f0.02:"Directivity parameter (${az}\260)":SWne -K -O >> $PS2
gmt psbasemap -R$R6-J$J6 -Xa$Xa61 -Ya$Ya6 -Ba4f2:"Time (s)":/a100Swen -K -O >> $PS2
if [ $np -gt $plnum ]
then
gmt psbasemap -R$R6 -J$J6 -Xa$Xa62 -Ya$Ya6 -Ba4f2:"Time (s)":/a100Swne -K -O >> $PS2
fi
if [ $(echo "$np > 2*$plnum" | bc) = 1 ]
then
gmt psbasemap -R$R6 -J$J6 -Xa$Xa63 -Ya$Ya6 -Ba4f2:"Time (s)":/a100Swne -K -O >> $PS2
fi
gmt psxy -R$R1_1 -J$J1 -O -T >> $PS2
#cat $PS2 >> ${datadir}/all.ps
gs -q -dNOPAUSE -dBATCH -sDEVICE=ps2write -sOutputFile=${datadir}/all.ps -f $PS1 $PS2
fi
##################################################################################################
figuredir=../figure
if [ ! -d $figuredir ]
then
mkdir $figuredir
fi

ps2pdf ${datadir}/all.ps $figuredir/plot_2D.pdf
evince $figuredir/plot_2D.pdf &
rm gmt.* ${datadir}/temp.grd ${datadir}/temp.cpt $PS $PS1 $PS2 ${datadir}/all.ps
