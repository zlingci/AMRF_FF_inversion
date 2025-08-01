/*************************************************************
 * *     Name:       sac2wiggle.c:
 * *
 * *     Description:      convert sac file to 3 columns data [x y var]
 * *
 * *     Author:           Can Ge, July, 2012 at WhIGG
 * *
 * * modification history:
 * **************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "sac.h"

typedef   int   BOOLEAN;
#define   TRUE   1
#define   FALSE   0 
int main(int argc, char **argv) {
 int i,k,npts;
 float tlength,*r,*data,t,dm,dm1,ave;
 char sacf[128],outf[128],var[128];
 FILE *fp;
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
	}else {
		strcpy(var, argv[i]);
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
	fp = fopen(outf,"w");
        dm=0.0;
	for (i=0;i<hd.npts;i++){
		if( fabs(data[i]) >=dm) dm=fabs(data[i]);
	}
        dm1=hd.dist;
        ave=0.0;
	for (i=0;i<40;i++){
		ave+=data[i];
	}
        ave=ave/40.0;
	for (i=0;i<hd.npts;i++){
		t=hd.delta*i+hd.b;
		fprintf(fp,"%f   %s   %f   %f   %f   %f \n",t,var,data[i],data[i]/dm,fabs(data[i]),data[i]*dm1);
	}
/*        fprintf(fp,"> max =  %f \n",dm);*/
      fclose(fp);


 return 0;
}


