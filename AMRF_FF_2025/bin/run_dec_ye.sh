#!/bin/bash

workdir=`pwd`
indir=$workdir/inp/
Tevent=mainshock
Eevent=egf
parfile=${indir}/parameter_${Eevent}.txt

if [ ! -d "dec_results/" ]
then
mkdir dec_results/
fi

if [ "$1"x = x ]
then
echo "Please give th row number in $parafile, such as: sh run_dec_ye.sh 2"
exit
fi

k=$1
snum=$1
while [ $k -le $snum ]
do
## reading parameter file
net=`awk -v s=$k '{if(NR == s)print $1}' $parfile`
sta=`awk -v s=$k '{if(NR == s)print $2}' $parfile`
f1=`awk -v s=$k '{if(NR == s)print $3}' $parfile`
f2=`awk -v s=$k '{if(NR == s)print $4}' $parfile`
dt=`awk -v s=$k '{if(NR == s)print $5}' $parfile`
tshift=`awk -v s=$k '{if(NR == s)print $6}' $parfile`
err=`awk -v s=$k '{if(NR == s)print $7}' $parfile`
state_zp=`awk -v s=$k '{if(NR == s)print $8}' $parfile`
state_rp=`awk -v s=$k '{if(NR == s)print $9}' $parfile`
state_zs=`awk -v s=$k '{if(NR == s)print $10}' $parfile`
state_rs=`awk -v s=$k '{if(NR == s)print $11}' $parfile`
state_ts=`awk -v s=$k '{if(NR == s)print $12}' $parfile`
pstart=`awk -v s=$k '{if(NR==s)print $13}' $parfile`
pend=`awk -v s=$k '{if(NR==s)print $14}' $parfile`
pshift=`awk -v s=$k '{if(NR==s)print $15}' $parfile`
sstart=`awk -v s=$k '{if(NR==s)print $16}' $parfile`
send=`awk -v s=$k '{if(NR==s)print $17}' $parfile`
sshift=`awk -v s=$k '{if(NR==s)print $18}' $parfile`

# deconvolution for S wave in .[ZRT] components
wave="s"
for component in "z" "r" "t"
do
comp=`echo $component | tr 'a-z' 'A-Z'`
Twavefile=mainshock_data/$net.${sta}.BH${comp}.sac
Ewavefile=egf_data/$net.${sta}.BH${comp}.sac  
dist=`saclst dist f $Twavefile | awk '{print $2}'`
az=`saclst az f $Twavefile | awk '{print $2}'`
state_flag=`eval echo '$state_'"$component""$wave"`
time_start=`eval echo '$'"$wave"'start'`
time_end=`eval echo '$'"$wave"'end'`
time_shift=`eval echo '$'"$wave"'shift'`
if [ $state_flag -eq 1 ]
then
Tstart=$time_start
Tend=$time_end
Estart=`bc << EOF
scale=3
($time_start + $time_shift)
EOF
`
Eend=`bc << EOF
scale=3
($time_end + $time_shift)
EOF
`
echo $Estart $Eend $Tstart $Tend

# cut wave
obs0=inp0.sac
gf0=gf0.sac
obs=inp.sac
gf=gf.sac
sac << EOD
cuterr fillz
cut a $Tstart $Tend 
r $Twavefile
rtr
rmean
rtr
taper
ch allt (0- &1,b&) iztype io
interpolate delta $dt
hp n 4 c $f1
w $obs0
cut a $Estart $Eend
r $Ewavefile
rtr
rmean
rtr
taper
ch allt (0- &1,b&)  iztype io
interpolate delta $dt
taper
hp n 4 c $f1
w $gf0
quit
EOD
$workdir/bin/taper -I$gf0 -O$gf
$workdir/bin/taper -I$obs0 -O$obs

# deconvolution
$workdir/bin/iterdeconfd.ye << end
$obs
$gf
10000        * nbumps
$tshift      * phase delay for result
$err         * min error improvement to accept
0.000001 $f2 4 1
1            * 1 allows negative bumps
0            * 0 form minimal output (1) will output lots of files
end
sac << EOD
r decon.out
ch allt (0- &1,b&) iztype io
w decon.out
quit
EOD
cp $obs0 dec_results/${component}${wave}.target
cp $gf0 dec_results/${component}${wave}.egf
cp numerator dec_results/${component}${wave}.target1
cp denominator dec_results/${component}${wave}.egf1
cp predicted dec_results/${component}${wave}.predicted
cp decon.out dec_results/${component}${wave}.out
fi
done

# deconvolution for P wave in .[ZR] components
wave="p"
for component in "z" "r"
do
comp=`echo $component | tr 'a-z' 'A-Z'`
Twavefile=mainshock_data/$net.${sta}.BH${comp}.sac
Ewavefile=egf_data/$net.${sta}.BH${comp}.sac  
dist=`saclst dist f $Twavefile | awk '{print $2}'`
az=`saclst az f $Twavefile | awk '{print $2}'`
state_flag=`eval echo '$state_'"$component""$wave"`
time_start=`eval echo '$'"$wave"'start'`
time_end=`eval echo '$'"$wave"'end'`
time_shift=`eval echo '$'"$wave"'shift'`
if [ $state_flag -eq 1 ]
then
Tstart=$time_start
Tend=$time_end
Estart=`bc << EOF
scale=3
($time_start + $time_shift)
EOF
`
Eend=`bc << EOF
scale=3
($time_end + $time_shift)
EOF
`
echo $Estart $Eend $Tstart $Tend

# cut wave
obs0=inp0.sac
gf0=gf0.sac
obs=inp.sac
gf=gf.sac
sac << EOD
cuterr fillz
cut a $Tstart $Tend 
r $Twavefile
rtr
rmean
rtr
taper
ch allt (0- &1,b&) iztype io
interpolate delta $dt
hp n 4 c $f1
w $obs0
cut a $Estart $Eend
r $Ewavefile
rtr
rmean
rtr
taper
ch allt (0- &1,b&)  iztype io
interpolate delta $dt
taper
hp n 4 c $f1
w $gf0
quit
EOD
$workdir/bin/taper -I$gf0 -O$gf
$workdir/bin/taper -I$obs0 -O$obs

# deconvolution
$workdir/bin/iterdeconfd.ye << end
$obs
$gf
10000        * nbumps
$tshift      * phase delay for result
$err         * min error improvement to accept
0.000001 $f2 4 1
1            * 1 allows negative bumps
0            * 0 form minimal output (1) will output lots of files
end
sac << EOD
r decon.out
ch allt (0- &1,b&) iztype io
w decon.out
quit
EOD
cp $obs0 dec_results/${component}${wave}.target
cp $gf0 dec_results/${component}${wave}.egf
cp numerator dec_results/${component}${wave}.target1
cp denominator dec_results/${component}${wave}.egf1
cp predicted dec_results/${component}${wave}.predicted
cp decon.out dec_results/${component}${wave}.out
fi
done

echo "sta:" $sta
echo "fh:" $f2 "Hz"
k=` expr $k + 1`
done
rm c00000 gf0.sac inp0.sac decon.out  denominator denominator1 gf.sac inp.sac  numerator1 numerator predicted temp.out
