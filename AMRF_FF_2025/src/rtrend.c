#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include"sac.h"

typedef   int   BOOLEAN;
#define   TRUE   1
#define   FALSE   0 

float Squre_sum(float *, int );
float Sum_Average(float *, int );
float X_Y_By(float *, float *, int );
float Squre_sum(float *, int );

int main(int argc, char **argv) { 
 int i;
 int nrtrend = 0;
 float *data,*time;
 char sacf[128],outf[128];
 SACHEAD hd;
 BOOLEAN   error = FALSE;
 float tlength = 0;
 float K = 0;               //拟合直线的斜率
 float R = 0;               //拟合直线的截距
 float time_sum_average = 0;   //数组 X[N] 个元素求和 并求平均值
 float data_sum_average = 0;   //数组 Y[N] 个元素求和 并求平均值
 float time_square_sum = 0;    //数组 X[N] 个个元素的平均值
 float time_multiply_data = 0;    //数组 X[N]和Y[N]对应元素的乘机

 /* input parameters */
 for (i=1; !error && i < argc; i++) {
        if (argv[i][0] == '-') {
                switch(argv[i][1]) {
                case 'I':
                  sscanf(&argv[i][2],"%s",sacf);
                  break;
                case 'O':
                  sscanf(&argv[i][2],"%s",outf);
                  break;
                default:
                  error = TRUE;
                  break;
                }
        }
 }
 
 if (argc == 1 || error ) {
     fprintf(stderr,"convert sac file to x-var-y ascii file. \n\
usage: %s -Isacfile -Ooutfile var\n",argv[0]);
                return -1;
 }

 if ( (data=read_sac(sacf,&hd)) == NULL) {
    fprintf(stderr,"error occur when open sac file1\n");
 }
 
 tlength = hd.a - 0.5;
 nrtrend = tlength/hd.delta; 
 time = (float *)malloc(nrtrend*sizeof(float));
 for (i=0;i<nrtrend;i++) {
     time[i] = i*hd.delta; 
 }
 
 time_sum_average= Sum_Average(time,nrtrend);
 data_sum_average = Sum_Average(data,nrtrend);
 time_square_sum = Squre_sum(time,nrtrend);
 time_multiply_data = X_Y_By(time,data,nrtrend);
 K = (time_multiply_data - nrtrend * time_sum_average * data_sum_average)/(time_square_sum - nrtrend * time_sum_average*time_sum_average);
 R = data_sum_average - K * time_sum_average;
 printf("K = %f\n",K);
 printf("R = %f\n",R);
 
 for (i=0;i<hd.npts;i++){
    data[i] = data[i] - (K*i*hd.delta + R);       
 }

 if ( (error = write_sac(outf,hd,data)) == -1) {
    fprintf(stderr,"error occur when write sac file\n");
 }

 return 0;
}

float Sum_Average(float d[],int N)
{
   unsigned int i=0;
   float z=0;
   for(i=0;i<N;i++) {
      z = z + d[i];
   }
   z = z/N;
   return z;
}

float X_Y_By(float m[],float n[],int N)
{
   unsigned int i=0;
   float z=0;
   for(i=0;i<N;i++) {
      z = z + m[i]*n[i];
   }
   return z;
 }

float Squre_sum(float c[], int N)
{
   unsigned int i=0;
   float z=0;
   for(i=0;i<N;i++) {
    z = z + c[i]*c[i];
   }
   return z;
}

