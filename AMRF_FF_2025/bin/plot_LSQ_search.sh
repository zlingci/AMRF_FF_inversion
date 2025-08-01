#!bin/sh

parfile=$1
datadir=`pwd`

resfile=${datadir}/VR.txt
resfile_use=${datadir}/VR_use.txt
ellipse=${datadir}/ellipse.txt
mean_var=${datadir}/mean_var.txt
normal_speed=${datadir}/normal_speed.txt
normal_direction=${datadir}/normal_direction.txt

gmt gmtset ANNOT_FONT_SIZE_PRIMARY 12p
gmt gmtset LABEL_FONT_SIZE 14p
gmt gmtset BASEMAP_TYPE plain
gmt gmtset FRAME_PEN 0.75p,black
gmt gmtset TICK_LENGTH -5p
gmt gmtset ANNOT_OFFSET_PRIMARY 3p
gmt gmtset LABEL_OFFSET 4p
gmt gmtset COLOR_FOREGROUND 115/0/0
gmt gmtset COLOR_BACKGROUND darkblue
##################################################################################################
#################plot variance redution varing with rupture speed and ruputre direction ########## 
##################################################################################################
dip=`awk -v s=8 '{if(NR == s)print $2}' $parfile`
strike=`awk -v s=8 '{if(NR == s)print $1}' $parfile`
dangle=`awk -v s=10 '{if(NR == s)print $3}' $parfile`
angle_shift=`awk -v s=10 '{if(NR == s)print $4}' $parfile`
echo $angle_shift
#angle_start=`awk -v s=10 -v n=$dangle -v m=$angle_shift '{if(NR == s)print $1+m-n/2}' $parfile`
#angle_start1=`awk -v s=10 -v m=$angle_shift '{if(NR == s)print $1+m}' $parfile`
#angle_end=`awk -v s=10 -v n=$dangle -v m=$angle_shift '{if(NR == s)print $2+m-n/2}' $parfile`
#angle_end1=`awk -v s=10 -v m=$angle_shift '{if(NR == s)print $2+m}' $parfile`
angle_start=`awk -v s=10 -v n=$dangle -v m=$angle_shift '{if(NR == s)print $1-n/2}' $parfile`
angle_start1=`awk -v s=10 -v m=$angle_shift '{if(NR == s)print $1}' $parfile`
angle_end=`awk -v s=10 -v n=$dangle -v m=$angle_shift '{if(NR == s)print $2-n/2}' $parfile`
angle_end1=`awk -v s=10 -v m=$angle_shift '{if(NR == s)print $2}' $parfile`
dvs=`awk -v s=9 '{if(NR == s)print $3}' $parfile`
vs_start=`awk -v s=9 -v n=$dvs '{if(NR == s)print $1-n/2}' $parfile`
vs_start1=`awk -v s=9 '{if(NR == s)print $1}' $parfile`
vs_end=`awk -v s=9 -v n=$dvs '{if(NR == s)print $2+n/2}' $parfile`
vs_end1=`awk -v s=9 '{if(NR == s)print $2}' $parfile`
R1=$angle_start/$angle_end/$vs_start/$vs_end
R11=$angle_start1/$angle_end1/$vs_start1/$vs_end1
J1=X15c/9c
Xa1=-1c
Ya1=5.3c
PS=tamp.ps
gmt psxy -R$R1 -J$J1 -T -K > $PS
max_var=`awk '{print $3}' $resfile|sort -n|tail -1|awk '{print $1*1.001}'`
min_var=`bc << EOF
scale=6
($max_var - 0.3*$max_var)
EOF
`
dvar=`bc << EOF
scale=6
(12/200)
EOF
`
gmt makecpt -Cjet -T$min_var/$max_var/0.01 > temp.cpt
awk -v m=$angle_shift '{if($1 < m)print $1+360,$2,$3;else print $1,$2,$3}' $resfile | gmt xyz2grd -R$R11 -Gtemp.grd -I$dangle/$dvs
gmt grdimage -R$R1 -J$J1 temp.grd -Xa$Xa1 -Ya$Ya1 -Ctemp.cpt -K -O >> $PS
awk '{print $1,$2}' $ellipse | gmt psxy -R$R1 -J$J1 -Xa$Xa1 -Ya$Ya1 -W2p,255/255/255,dashed+s -O -K >> $PS
awk '{print $1}' $resfile_use | gmt pshistogram -R$angle_start/$angle_end/0/130 -J$J1 -Xa$Xa1 -Ya$Ya1 -W5 -L2p,white -O -K >> $PS
awk '{print $2-0.001}' $resfile_use | gmt pshistogram -R$vs_start/$vs_end/0/120 -J$J1 -Xa$Xa1 -Ya$Ya1 -L2p,white -A -W0.2 -O -K >> $PS
max1=`awk '{print $2}' $normal_speed|sort -n|tail -1|awk '{print $1*1.001}'`
awk -v s=$max1 '{print $2/s,$1}' $normal_speed | gmt psxy -R0/4.5/$vs_start/$vs_end -J$J1 -Xa$Xa1 -Ya$Ya1 -W2p,255/255/255+s -O -K >> $PS
max2=`awk '{print $2}' $normal_direction|sort -n|tail -1|awk '{print $1*1.001}'`
awk -v s=$max2 '{print $1,$2/s}' $normal_direction | gmt psxy -R$angle_start/$angle_end/0/4.5 -J$J1 -Xa$Xa1 -Ya$Ya1 -W2p,255/255/255+s -O -K >> $PS
awk '{print $1,$2}' $mean_var | gmt psxy -R$R1 -J$J1 -Xa$Xa1 -Ya$Ya1 -Sa0.5c -G255/255/255 -O -K >> $PS
gmt psscale -Xa$Xa1 -Ya$Ya1 -D15.5c/4.5c/8c/0.5c -Ctemp.cpt -Ba4f2:"Variance reduction (%)": -K -O >> $PS
gmt psbasemap -R$R1 -J$J1 -Xa$Xa1 -Ya$Ya1 -Ba60f30:"Rupture direction (\260)":/a1f0.2:"Rupture speed (km/s)":SWne -K -O >> $PS
direction=`awk '{printf "%.1f",$1}' $mean_var`
speed=`awk '{printf "%.1f",$2}' $mean_var`
direction_var=`awk '{printf "%.1f",$3*2}' $mean_var`
speed_var=`awk '{printf "%.1f",$4*2}' $mean_var`
gmt pstext -R0/1/0/1 -J$J1 -Xa$Xa1 -Ya$Ya1 -N -K -O >> $PS << EOF
0.75 0.95 15 0 5 ML @;255/0/0; $direction\260 \261 $direction_var\260 @;;
0.75 0.88 15 0 5 ML @;255/0/0; $speed \261 $speed_var km/s @;;
EOF
#####################################################################################################
figuredir=../figure
if [ ! -d $figuredir ]
then
mkdir $figuredir
fi
ps2pdf $PS $figuredir/plot_search.pdf
evince $figuredir/plot_search.pdf &
rm $PS temp.grd temp.cpt
rm gmt.*
