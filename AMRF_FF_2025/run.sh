
modelname=yunnanEYA
workdir=`pwd`


### 0. deconvolution for each station (each line in the parameter file):
# sh bin/run_dec_ye.sh $NR   # $NR>=2

### 1. deconvolution for all data
sh bin/runs_dec_ye.sh $modelname

### 2. cd dec_results/ and manually mark the onsets of AMRFs(*.out) as timemark A
###    Here the manually picked data are provided in 'dec_results_picked/'
rm -r dec_results/
cp -r dec_results_picked/ dec_results/

### 3. plot figures for deconvolution
sh bin/plot_stf.sh S ts
sh bin/plot_az.sh S ts
sh bin/plot_ray.sh S ts

### 4. prepare data for finite fault inversion
if [ ! -d data_inv ]; then
mkdir data_inv
fi
cp dec_results/*.out data_inv/
cp dec_results/*.txt data_inv/
cd data_inv/
ls -l *out | awk '{print $9}' | awk -F. '{print "mv "$0" "$1"."$3"."$4"."$5}' | sh
cd ..

### 5.calculate directivity for the mainshock
### uses grid search, takes long time
# sh LSQ_search.sh

### 6. calculate takeoff for finite fault inversion
sh bin/caltakeoff.sh $modelname
### run 1D inversion with parameter file 'inp/LSQ_1D.txt'
sh LSQ_1D.sh

### 7. run 2D inversion with parameter file 'inp/LSQ_2D.txt'
sh LSQ_2D.sh


