/*
 *   This file is part of TISEAN
 *
 *   Copyright (c) 1998-2007 Rainer Hegger, Holger Kantz, Thomas Schreiber
 *
 *   TISEAN is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   TISEAN is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with TISEAN; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
/*Author: Rainer Hegger. Last modified, Jul 13, 2010 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include "routines/tsa.h"

#define WID_STR "Estimates the time delayed mutual information\n\t\
of one or two data columns"


char *file_out=NULL,stout=1,dimset=0;
char *infile=NULL;
unsigned long length=ULONG_MAX,exclude=0;
unsigned int dim=1;
unsigned int verbosity=0xff;
long partitions=16,corrlength=20;
long **array,*h1,*h11,**h2;
char *column=NULL;

void show_options(char *progname)
{
  what_i_do(progname,WID_STR);
  fprintf(stderr," Usage: %s [Options]\n\n",progname);
  fprintf(stderr," Options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted"
          " as a possible"
          " datafile.\nIf no datafile is given stdin is read. Just - also"
          " means stdin\n");
  fprintf(stderr,"\t-l # of points to be used [Default is all]\n");
  fprintf(stderr,"\t-x # of lines to be ignored [Default is 0]\n");
  fprintf(stderr,"\t-c columns to read  [Default is 1]\n");
  fprintf(stderr,"\t-b # of boxes [Default is 16]\n");
  fprintf(stderr,"\t-D max. time delay [Default is 20]\n");
  fprintf(stderr,"\t-o output file [-o without name means 'datafile'.mut;"
	  "\n\t\tNo -o means write to stdout]\n");
  fprintf(stderr,"\t-V verbosity level [Default is 1]\n\t\t"
          "0='only panic messages'\n\t\t"
          "1='+ input/output messages'\n");
  fprintf(stderr,"\t-h  show these options\n");
  fprintf(stderr,"\n");
  exit(0);
}

void scan_options(int n,char** in)
{
  char *out;

  if ((out=check_option(in,n,'l','u')) != NULL)
    sscanf(out,"%lu",&length);
  if ((out=check_option(in,n,'x','u')) != NULL)
    sscanf(out,"%lu",&exclude);
  if ((out=check_option(in,n,'c','s')) != NULL)
    column=out;
  if ((out=check_option(in,n,'b','u')) != NULL)
    sscanf(out,"%lu",&partitions);
  if ((out=check_option(in,n,'D','u')) != NULL)
    sscanf(out,"%lu",&corrlength);
  if ((out=check_option(in,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(in,n,'o','o')) != NULL) {
    stout=0;
    if (strlen(out) > 0)
      file_out=out;
  }
}

double make_cond_entropy(long t)
{
  long i,j,hi,hii,count=0,start,stop;
  double hpi,hpj,pij,cond_ent=0.0,norm;

  for (i=0;i<partitions;i++) {
    h1[i]=h11[i]=0;
    for (j=0;j<partitions;j++)
      h2[i][j]=0;
  }
  if (t < 0) {
    start=0;
    stop=length+t;
    t=-t;
  }
  else {
    start=t;
    stop=length-t;
  }

  for (i=start;i<stop;i++) {
    hii=array[0][i];
    hi=array[1][i+t];
    h1[hi]++;
    h11[hii]++;
    h2[hi][hii]++;
    count++;
  }

  norm=1.0/(double)count;
  cond_ent=0.0;

  for (i=0;i<partitions;i++) {
    hpi=(double)(h1[i])*norm;
    if (hpi > 0.0) {
      for (j=0;j<partitions;j++) {
	hpj=(double)(h11[j])*norm;
	if (hpj > 0.0) {
	  pij=(double)h2[i][j]*norm;
	  if (pij > 0.0)
	    cond_ent += pij*log(pij/hpj/hpi);
	}
      }
    }
  }

  return cond_ent;
}

int main(int argc,char** argv)
{
  char stdi=0;
  long tau,i,j;
  double **series,min1,min2,interval1,interval2,shannon,condent;
  FILE *file;
  
  if (scan_help(argc,argv))
    show_options(argv[0]);
  
  scan_options(argc,argv);
#ifndef OMIT_WHAT_I_DO
  if (verbosity&VER_INPUT)
    what_i_do(argv[0],WID_STR);
#endif

  infile=search_datafile(argc,argv,NULL,verbosity);
  if (infile == NULL)
    stdi=1;

  if (file_out == NULL) {
    if (!stdi) {
      check_alloc(file_out=(char*)calloc(strlen(infile)+5,(size_t)1));
      strcpy(file_out,infile);
      strcat(file_out,".mut");
    }
    else {
      check_alloc(file_out=(char*)calloc((size_t)10,(size_t)1));
      strcpy(file_out,"stdin.mut");
    }
  }
  if (!stout)
    test_outfile(file_out);

  if (column == NULL) {
    series=(double**)get_multi_series(infile,&length,exclude,&dim,"",
				      dimset,verbosity);
  }
  else {
    series=(double**)get_multi_series(infile,&length,exclude,&dim,column,
                                      dimset,verbosity);
  }

  check_alloc(h1=(long *)malloc(sizeof(long)*partitions));
  check_alloc(h11=(long *)malloc(sizeof(long)*partitions));
  check_alloc(h2=(long **)malloc(sizeof(long *)*partitions));
  for (i=0;i<partitions;i++) 
    check_alloc(h2[i]=(long *)malloc(sizeof(long)*partitions));
  check_alloc(array=(long **)malloc(sizeof(long*)*2));
  check_alloc(array[0]=(long *)malloc(sizeof(long)*length));
  check_alloc(array[1]=(long *)malloc(sizeof(long)*length));
  if (dim == 1) {
    rescale_data(series[0],length,&min1,&interval1);
    for (i=0;i<length;i++)
      if (series[0][i] < 1.0)
	array[0][i]=array[1][i]=(long)(series[0][i]*(double)partitions);
      else
	array[0][i]=array[1][i]=partitions-1;
    free(series[0]);
  }
  else {
    rescale_data(series[0],length,&min1,&interval1);
    rescale_data(series[1],length,&min2,&interval2);
    for (j=0;j<2;j++) {
      for (i=0;i<length;i++)
	if (series[j][i] < 1.0)
	  array[j][i]=(long)(series[j][i]*(double)partitions);
	else
	  array[j][i]=partitions-1;
      free(series[j]);
    }
  }
  free(series);

  shannon=make_cond_entropy(0);
  if (corrlength >= length)
    corrlength=length-1;

  if (!stout) {
    file=fopen(file_out,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",file_out);
    for (tau=-corrlength;tau<=corrlength;tau++) {
      condent=make_cond_entropy(tau);
      fprintf(file,"%ld %e %e\n",tau,condent,condent/log((double)partitions));
      fflush(file);
    }
    fclose(file);
  }
  else {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Writing to stdout\n");
    for (tau=-corrlength;tau<=corrlength;tau++) {
      condent=make_cond_entropy(tau);
      fprintf(stdout,"%ld %e %e\n",tau,condent,condent/log((double)partitions));
      fflush(stdout);
    }
  }

  return 0;
}

