#!/bin/bash


cha1=$1
cha=$2
Eevent=egf
Tevent=mainshock
duration=12
dir_data=dec_results

if [ "$2"x = x ]
then
echo "Please give the wave (P,S) and component (rp,zp,rs,zs,ts) such as: sh plot_stf.sh S ts "
exit
fi

sfile=${dir_data}/$Eevent.${cha}.txt
if [ ! -f $sfile ]
then
echo "Parameter file $sfile does not exist , Please check it"
exit
fi
###################################################################################################
#gmtset ANOT_FONT_SIZE  10 LABEL_FONT_SIZE 12 HEADER_FONT_SIZE 14p
#gmtset ANNOT_OFFSET_PRIMARY 0.075c ANNOT_OFFSET_SECONDARY 0.1c LABEL_OFFSET 0.075c HEADER_OFFSET -0.15c
#gmtset TICK_LENGTH -0.10c LABEL_FONT  Helvetica ANNOT_FONT Helvetica

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
plnum=`wc $sfile | awk '{print $1}'`
PS=temp.ps
##############################################################################################
####################################### plot single wave and stf #############################
#############################################################################################
J1=X10c/15c
Xa1=7c
Ya1=0c
J2=X15c/3c
Xa2=6c
Ya2=15c
J3=X5c/15c
Xa3=2c
Ya3=0c
R=0/1/0/1
if [ -f ${sfile} ]
then
snum=`wc $sfile | awk '{print $1}'`
else
snum=0
fi
i=1
while [ $i -le $snum ]
do
echo $i
net=`awk -v s=$i '{if(NR == s)print $2}' $sfile`
sname=`awk -v s=$i '{if(NR == s)print $3}' $sfile`
dist=`awk -v s=$i '{if(NR == s)print $5}' $sfile`
az=`awk -v s=$i '{if(NR == s)print $6}' $sfile`
fitting=`awk -v s=$i '{if(NR == s)print $9}' $sfile`
dist=`printf '%-.1f' $dist`
az=`printf '%-.1f' $az`
fitting=`printf '%-.1f' $fitting`
gmt psxy -R$R -J$J1 -T -K > $PS
Tfile=$dir_data/${Eevent}.${net}.${sname}.${cha}.target
Tfile1=$dir_data/${Eevent}.${net}.${sname}.${cha}.target1
Tfile2=$dir_data/${Eevent}.${net}.${sname}.${cha}.predicted
Efile=$dir_data/${Eevent}.${net}.${sname}.${cha}.egf
Efile1=$dir_data/${Eevent}.${net}.${sname}.${cha}.egf1
Decon=$dir_data/${Eevent}.${net}.${sname}.${cha}.out
s_start=`saclst b f ${Tfile} | awk '{print $2}'`
s_end=`saclst e f ${Tfile} | awk '{print $2}'`
R=$s_start/$s_end/0/20
gmt psbasemap -R$R -J$J1 -Xa$Xa1 -Ya$Ya1 -Bxaf+l"Time (s)" -BSwen -K -O >> $PS
scale=0.8
ind=18
inpsac=$Tfile
outfil=tmp.txt
./bin/sac2wiggle -I$inpsac -O$outfil $ind
awk '{print $1,$2,$4}' $outfil | gmt pswiggle -R$R -J$J1 -Xa$Xa1 -Ya$Ya1 -Z$scale -W0.5p,black,solid -BwSne -O -K >> $PS
ind=14
inpsac=$Tfile1
outfil=tmp.txt
./bin/sac2wiggle -I$inpsac -O$outfil $ind
awk '{print $1,$2,$4}' $outfil | gmt pswiggle -R$R -J$J1 -Xa$Xa1 -Ya$Ya1 -Z$scale -W0.5p,black,solid -BwSne -O -K >> $PS
ind=14
inpsac=$Tfile2
outfil=tmp.txt
./bin/sac2wiggle -I$inpsac -O$outfil $ind
awk '{print $1,$2,$4}' $outfil | gmt pswiggle -R$R -J$J1 -Xa$Xa1 -Ya$Ya1 -Z$scale -W0.5p,red,solid -BwSne -O -K >> $PS
ind=10
inpsac=$Efile
outfil=tmp.txt
./bin/sac2wiggle -I$inpsac -O$outfil $ind
awk '{print $1,$2,$4}' $outfil | gmt pswiggle -R$R -J$J1 -Xa$Xa1 -Ya$Ya1 -Z$scale -W0.5p,black,solid -BwSne -O -K >> $PS
ind=6
inpsac=$Efile1
outfil=tmp.txt
./bin/sac2wiggle -I$inpsac -O$outfil $ind
awk '{print $1,$2,$4}' $outfil | gmt pswiggle -R$R -J$J1 -Xa$Xa1 -Ya$Ya1 -Z$scale -W0.5p,black,solid -BwSne -O -K >> $PS
ind=2
inpsac=$Decon
outfil=tmp.txt
./bin/sac2wiggle -I$inpsac -O$outfil $ind
awk '{print $1,$2,$4}' $outfil | gmt pswiggle -R$R -J$J1 -Xa$Xa1 -Ya$Ya1 -Z$scale -W0.5p,black,solid -BwSne -O -K >> $PS
############################################### text ###################################
color=0/0/0
gmt pstext -R0/1/0/1 -J$J1 -Xa$Xa1 -Ya$Ya1 -F+f15,5 -K -O >> $PS << EOF
0.05 0.95 @;$color;@S@;;
EOF
gmt pstext -R0/1/0/1 -J$J2 -Xa$Xa2 -Ya$Ya2 -F+f15,5+jLM -K -O >> $PS << EOF
0.05 0.5 @;$color;@${sname}-${cha} Dist=${dist} Az=${az} fitting=${fitting}% @;;
EOF
gmt pstext -R0/1/0/20 -J$J3 -Xa$Xa3 -Ya$Ya3 -F+f12,5+jLM -K -O >> $PS << EOF
0.5 2 @;$color;STF@;;
0.5 6 @;$color;EGF(filter)@;;
0.5 10 @;$color;EGF(raw)@;;
0.5 14 @;$color;Event(filter)@;;
0.5 13 @;255/0/0;Event(syn)@;;
0.5 18 @;$color;Event(raw)@;;
EOF
gmt psxy -R$R -J$J1 -O -T >> $PS 
cat $PS >> tempall.ps
rm $PS
i=` expr $i + 1`
done
##########################################################################################
################################# plot az-stf and wave ###################################
#########################################################################################
snum=`wc $sfile | awk '{print $1}'`
R=${stf_start}/${stf_end}/-1/$plnum
R1=${wave_start}/${wave_end}/-1/$plnum
J1=X5.5c/19c
Xa1=1c
Ya1=-1c
J2=X2.8c/19c
Xa2=-1.9c
Ya2=-1c
J3=X8c/19c
Xa3=7c
Ya3=-1c
J4=X8c/19c
Xa4=15.5c
Ya4=-1c
PS=temp.ps
gmt psxy -R$R -J$J1 -T -K > $PS
gmt psbasemap -R$R -J$J1 -Xa$Xa1 -Ya$Ya1 -Bxaf+l"Time (s)" -BwSen -K -O >> $PS
gmt psbasemap -R$R1 -J$J3 -Xa$Xa3 -Ya$Ya3 -Bxaf+l"Time (s)" -BwSen -K -O >> $PS
i=1
k=0
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
color=blue
######################################## plot stf ################################################ 
Decon=$dir_data/${Eevent}.${net}.${sname}.$cha.out
ptime=`saclst a f $Decon | awk '{print $2}'`
sac << EOD
cut a ${stf_start} ${stf_end}
r $Decon
ch allt (0- &1,a&)  iztype io
taper
w temp.sac
q
EOD
ind=`bc << EOF
scale=0
($i - $plnum*$k - 1)
EOF
`
inpsac=temp.sac
outfil=tmp.txt
./bin/sac2wiggle -I$inpsac -O$outfil $ind
max_stf=`sort -rk 3 -n $outfil | awk -v s=1 '{if(NR == s)print $3}'`
awk -v n=${max_stf} '{print $1,$3/n+$2}' $outfil | gmt psxy -R$R -J$J1 -Xa$Xa1 -Ya$Ya1 -W0.1p,black,solid -G$color -O -K >> $PS
awk -v n=${max_stf} '{print $1,$3/n+$2}' $outfil | gmt psxy -R$R -J$J1 -Xa$Xa1 -Ya$Ya1 -W1p,black,solid -O -K >> $PS
gmt psxy -R0/1/-1/$plnum -J$J1 -Xa$Xa1 -Ya$Ya1 -W0.5p,black,solid -O -K >> $PS << EOF
0 $ind
1 $ind
EOF
############################################### plot wave ########################################
target1=$dir_data/${Eevent}.${net}.$sname.$cha.target1
target2=$dir_data/${Eevent}.${net}.$sname.$cha.predicted
egf=$dir_data/${Eevent}.${net}.$sname.$cha.egf1
ind=`bc << EOF
scale=0
($i - $plnum*$k - 1)
EOF
`
inpsac=$target1
outfil=tmp.txt
./bin/sac2wiggle -I$inpsac -O$outfil $ind
scale=1.4
awk '{print $1,$2,$4}' $outfil | gmt pswiggle -R$R1 -J$J3 -Xa$Xa3 -Ya$Ya3 -Z$scale -W0.5p,black,solid -O -K >> $PS
inpsac=$target2
outfil=tmp.txt
./bin/sac2wiggle -I$inpsac -O$outfil $ind
scale=1.4
awk '{print $1,$2,$4}' $outfil | gmt pswiggle -R$R1 -J$J3 -Xa$Xa3 -Ya$Ya3 -Z$scale -W0.5p,red,solid -O -K >> $PS
inpsac=$egf
outfil=tmp.txt
./bin/sac2wiggle -I$inpsac -O$outfil $ind
scale=1.4
awk '{print $1,$2,$4}' $outfil | gmt pswiggle -R$R1 -J$J4 -Xa$Xa4 -Ya$Ya4 -Z$scale -W0.5p,100/100/100,solid -O -K >> $PS
################################################### text ###########################################
x=` bc << EOF
sacle=2
($ind + 0.0)
EOF
`
x1=` bc << EOF
sacle=2
($ind + 0.35)
EOF
`
gmt pstext -R0/1/-1/$plnum -J$J2 -Xa$Xa2 -Ya$Ya2 -F+f10,5+jLM -K -O >> $PS << EOF
0.01 $x @;0/0/0;${sname}_${cha}@;;
EOF
gmt pstext -R0/1/-1/$plnum -J$J2 -Xa$Xa2 -Ya$Ya2 -F+f10,5+jLM -K -O >> $PS << EOF
0.5 $x @;0/0/0;${az}@;;
EOF
gmt pstext -R0/1/-1/$plnum -J$J3 -Xa$Xa3 -Ya$Ya3 -F+f10,5+jLM -K -O >> $PS << EOF
0.02 $x1 @;0/0/0;${fitting}%@;;
EOF
if [ $(echo "$i - $k*$plnum == $plnum" |bc) = 1 -o $i -eq $snum ]
then
k=`expr $k + 1`
x=` bc << EOF
sacle=2
($plnum - 0.5)
EOF
`
gmt pstext -R0/1/-1/$plnum -J$J2 -Xa$Xa2 -Ya$Ya2 -F+f10,5+jLM -K -O >> $PS << EOF
0.01 $x @;0/0/0;Sta@;;
EOF
gmt pstext -R0/1/-1/$plnum -J$J2 -Xa$Xa2 -Ya$Ya2 -F+f10,5+jLM -K -O >> $PS << EOF
0.5 $x @;0/0/0;Az@;;
EOF
gmt psxy -R$R -J$J1 -Xa$Xa1 -Ya$Ya1 -W0.5p,red,dashed -O -K >> $PS << EOF
0 -1
0 $plnum
EOF
gmt psbasemap -R$R -J$J1 -Xa$Xa1 -Ya$Ya1 -Bxaf+l"Time (s)" -BwSen -K -O >> $PS
gmt psbasemap -R$R1 -J$J3 -Xa$Xa3 -Ya$Ya3 -Bxaf+l"Time (s)" -BwSen -K -O >> $PS
gmt psbasemap -R$R1 -J$J4 -Xa$Xa4 -Ya$Ya4 -Bxaf+l"Time (s)" -BwSen -K -O >> $PS
gmt psxy -R$R -J$J1 -O -T >> $PS
cat $PS >> tempall.ps
rm $PS
if [ $i -ne $snum ]
then
gmt psxy -R$R -J$J1 -T -K > $PS
fi
fi
i=` expr $i + 1`
done
if [ ! -d figure ]
then
mkdir -p figure
fi
ps2pdf tempall.ps figure/${Tevent}_${Eevent}_${cha}.pdf
evince figure/${Tevent}_${Eevent}_${cha}.pdf & 
rm tempall.ps gmt.* tmp.* temp.sac 
