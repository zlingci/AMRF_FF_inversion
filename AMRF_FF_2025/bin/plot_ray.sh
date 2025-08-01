#!/bin/bash

cha1=$1
cha=$2
Eevent=egf
Tevent=mainshock
duration=12
dir_data=dec_results
Vr_direction=137

if [ "$2"x = x ]
then
echo "Please give the wave (P,S) and component (rp,zp,rs,zs,ts) such as: sh plot_stf.sh S ts"
exit
fi

sfile=${dir_data}/$Eevent.${cha}.txt
if [ ! -f $sfile ]
then
echo "Parameter file $sfile does not exist , Please check it"
exit
fi

#####################################################################################################
#gmtset ANOT_FONT_SIZE  10 LABEL_FONT_SIZE 12 HEADER_FONT_SIZE 14p
#gmtset ANNOT_OFFSET_PRIMARY 0.075c ANNOT_OFFSET_SECONDARY 0.1c LABEL_OFFSET 0.075c HEADER_OFFSET -0.15c
#gmtset TICK_LENGTH -0.10c LABEL_FONT Helvetica ANNOT_FONT Helvetica 
gmt gmtset MAP_FRAME_TYPE plain
gmt gmtset MAP_FRAME_PEN 0.75p,black
gmt gmtset MAP_TICK_LENGTH_PRIMARY -5p
gmt gmtset MAP_ANNOT_OFFSET 8p
gmt gmtset MAP_TICK_PEN 0.75p,black

stf_start=`bc << EOF
scale=2
(-0.1*$duration)
EOF
`
stf_end=$duration
if [ "$cha1" = "S" ]
then
wave_start=0.0
wave_end=110
if [ "$cha" != "ts" -a "$cha" != "zs" -a "$cha" != "rs" ]
then
echo "Please give the corrected compomnet (ts,zs,rs) for S wave at \$2"
exit
fi
elif [ "$cha1" = "P" ]
then
wave_start=0.0
wave_end=35
if [ "$cha" != "zp" -a "$cha" != "rp" ]
then
echo "Please give the corrected compomnet (rp,zp) for P wave at \$2"
exit
fi
else
echo "Uncorreted wave: $cha1, please give the P or S at \$1"
exit
fi
snum=`wc $sfile | awk '{print $1}'`
PS=temp.ps
################################### plot az-stf and wave ##############################################
rm max_stf.txt
if [ "$cha1" = "S" ]
then
Rup=0.31
Rdown=-0.30
else
Rup=0.24
Rdown=-0.21
fi
R=${stf_start}/${stf_end}/$Rdown/$Rup
R1=${wave_start}/${wave_end}/-1/$snum
J1=X8/19c
Xa1=-0.5c
Ya1=-1c
J2=X2.5c/19c
Xa2=-3.3c
Ya2=-1c
J3=X8c/19c
Xa3=10.3c
Ya3=-1c
J4=X2.7c/19c
Xa4=7.8c
Ya4=-1c
J5=X8c/19c
Xa5=18.7c
Ya5=-1c
gmt psxy -R$R -J$J1 -T -K > $PS
i=1
while [ $i -le $snum ]
do
net=`sort -k 6 -n $sfile | awk -v s=$i '{if(NR == s)print $2}'`
sname=`sort -k 6 -n $sfile | awk -v s=$i '{if(NR == s)print $3}'`
dist=`sort -k 6 -n  $sfile | awk -v s=$i '{if(NR == s)print $5}'`
az=`sort -k 6 -n $sfile | awk -v s=$i '{if(NR == s)print $6}'`
fitting=`sort -k 6 -n $sfile | awk -v s=$i '{if(NR == s)print $9}'`
dist=`printf '%-.1f' $dist`
az=`printf '%-.1f' $az`
fitting=`printf '%-.1f' $fitting`
dpar=`sort -k 6 -n $sfile | awk -v s=$i '{if(NR == s)print $7}'`
dpar=`bc -l << EOF
scale=4
(-1*$dpar*c(($az-${Vr_direction})*3.1415/180))
EOF
`
dpar=`printf %.4f $dpar`
if [ $(echo "$dpar <= 0"| bc) = 1 ]
then
dpar=`bc << EOF
sacle=4
($dpar*-1)
EOF
`
k=-1
else	
k=1
fi
######################################################### plot stf ###############################
Decon=$dir_data/${Eevent}.${net}.$sname.$cha.out
ptime=`saclst a f $Decon | awk '{print $2}'`
sac << EOD
cut a ${stf_start} ${stf_end}
r $Decon
ch allt (0- &1,a&)  iztype io
taper
w temp.sac
q
EOD
ind=$dpar
inpsac=temp.sac
outfil=tmp.txt
./bin/sac2wiggle1 -I$inpsac -O$outfil $ind
max_stf=`sort -rk 3 -n $outfil | awk -v s=1 '{if(NR == s)print $3}'`
echo $sname $max_stf >> max_stf.txt
#max_stf=1000
if [ "$cha1" = "S" ]
then
awk -v n=${max_stf} -v s=$k  '{print $1,$3*0.12/400+$2*s}' $outfil | gmt psxy -R$R -J$J1 -Xa$Xa1 -Ya$Ya1 -W1p,black,solid -O -K >> $PS
else
awk -v n=${max_stf} -v s=$k  '{print $1,$3*0.08/400+$2*s}' $outfil | gmt psxy -R$R -J$J1 -Xa$Xa1 -Ya$Ya1 -W1p,black,solid -O -K >> $PS
fi
########################################### plot wave ############################################ 
ind1=` expr $i - 1`
inpsac=$dir_data/${Eevent}.${net}.$sname.$cha.target1
outfil=tmp.txt
./bin/sac2wiggle -I$inpsac -O$outfil $ind1
scale=1.4
awk '{print $1,$2,$4}' $outfil | gmt pswiggle -R$R1 -J$J3 -Xa$Xa3 -Ya$Ya3 -Z$scale -W0.5p,black,solid -O -K >> $PS
inpsac=$dir_data/${Eevent}.${net}.$sname.$cha.predicted
outfil=tmp.txt
./bin/sac2wiggle -I$inpsac -O$outfil $ind1
scale=1.4
awk '{print $1,$2,$4}' $outfil | gmt pswiggle -R$R1 -J$J3 -Xa$Xa3 -Ya$Ya3 -Z$scale -W0.5p,red,solid -O -K >> $PS
inpsac=$dir_data/${Eevent}.${net}.$sname.$cha.egf1
outfil=tmp.txt
./bin/sac2wiggle -I$inpsac -O$outfil $ind1
scale=1.4
awk '{print $1,$2,$4}' $outfil | gmt pswiggle -R$R1 -J$J5 -Xa$Xa5 -Ya$Ya5 -Z$scale -W0.5p,100/100/100,solid -O -K >> $PS
########################################## text ##################################################3
x=` bc << EOF
sacle=2
($ind1 + 0)
EOF
`
x1=` bc << EOF
sacle=2
($ind1 + 0.35)
EOF
`
gmt pstext -R0/1/-1/$snum -J$J3 -Xa$Xa3 -Ya$Ya3 -F+f10,5+jLM -K -O >> $PS << EOF
0.02 $x1 @;0/0/0;${fitting}%@;;
EOF
gmt pstext -R0/1/-1/$snum -J$J4 -Xa$Xa4 -Ya$Ya4 -F+f10,5+jLM -K -O >> $PS << EOF
0.01 $x @;0/0/0;${sname}_${cha}@;;
EOF
dpar=`bc << EOF
scale=4
($dpar*$k)
EOF
`
dpar=`printf "%.3f" $dpar`
echo $dpar
gmt pstext -R0/1/-1/$snum -J$J4 -Xa$Xa4 -Ya$Ya4 -F+f10,5+jLM -K -O >> $PS << EOF
0.5 $x @;0/0/0;${dpar}@;;
EOF
i=` expr $i + 1`
done

x=` bc << EOF
sacle=2
($snum - 0.5)
EOF
`
gmt pstext -R0/1/-1/$snum -J$J4 -Xa$Xa4 -Ya$Ya4 -F+f10,5+jLM -K -O >> $PS << EOF
0.01 $x @;0/0/0;Sta@;;
EOF
gmt pstext -R0/1/-1/$snum -J$J4 -Xa$Xa4 -Ya$Ya4 -F+f10,5+jLM -K -O >> $PS << EOF
0.5 $x @;0/0/0;Az@;;
EOF
gmt psxy -R$R -J$J1 -Xa$Xa1 -Ya$Ya1 -W1p,red,dashed -O -K >> $PS << EOF
0 $Rdown
0 $Rup
EOF
#gmt psxy -R$R -J$J1 -Xa$Xa1 -Ya$Ya1 -W2p,red,dashed -O -K >> $PS << EOF
#2.015 -0.30
#7.276 0.31
#EOF
#gmt psxy -R$R -J$J1 -Xa$Xa1 -Ya$Ya1 -Sc0.25c -W1p,black,solid -Gred -O -K >> $PS << EOF
#2.3 -0.267
#7.0 0.278
#4.603 0.000
#EOF
#gmt psxy -R$R -J$J1 -Xa$Xa1 -Ya$Ya1 -W2p,blue,dashed -O -K >> $PS << EOF
#2.982 -0.30
#7.011 0.31
#EOF
#gmt psxy -R$R -J$J1 -Xa$Xa1 -Ya$Ya1 -Sc0.25c -W1p,black,solid -Gblue -O -K >> $PS << EOF
#3.2 -0.267
#6.8 0.278
#4.964 0.000
#EOF
gmt psbasemap -R$R -J$J1 -Xa$Xa1 -Ya$Ya1 -Bxaf+l"Time (s)" -Bya0.1f0.02+l"Directivity parameter (@~\146@~@-r@- = ${Vr_direction})" -BWSen -K -O >> $PS
gmt psbasemap -R$R1 -J$J3 -Xa$Xa3 -Ya$Ya3 -Bxaf+l"Time (s)" -BwSen -K -O >> $PS
gmt psbasemap -R$R1 -J$J5 -Xa$Xa5 -Ya$Ya5 -Bxaf+l"Time (s)" -BwSen -K -O >> $PS
gmt psxy -R$R -J$J1 -O -T >> $PS 
if [ ! -d figure ]
then
mkdir -p figure
fi
ps2pdf ${PS} figure/${Tevent}_${Eevent}_${cha}.ray.pdf
evince figure/${Tevent}_${Eevent}_${cha}.ray.pdf &
rm $PS gmt.* tmp.* tempall.ps temp.sac max_stf.txt
