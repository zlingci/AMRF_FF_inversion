#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include"sac.h"

typedef   int   BOOLEAN;
#define   TRUE   1
#define   FALSE   0 
int main(int argc, char **argv) {
 int i;
 int nmean = 0;
 float *data;
 float mean = 0;
 float tlength = 0;
 char sacf[128],outf[128];
 SACHEAD hd;
 BOOLEAN   error = FALSE;

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

 tlength = hd.a - 0.2;
 nmean = tlength/hd.delta;
 for (i=0;i<nmean;i++){
     mean += data[i];
 }
 mean = mean/nmean;
 for (i=0;i<hd.npts;i++){
    data[i] = data[i] - mean;
 } 

if ( (error = write_sac(outf,hd,data)) == -1) {
  fprintf(stderr,"error occur when write sac file\n");
}

 return 0;
}


