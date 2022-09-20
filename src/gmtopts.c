/* @(#)gmtopts.c
 */

#ifndef BUFSIZ
#define BUFSIZ 1024
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gmtopts.h"

/*
 *
 */
int parse_R_option(char *item, double *w, double *e, double *s, double *n)
{

	/*  char text[BUFSIZ]; */
	char string[BUFSIZ];
	int ll = 0, nn = 0;

	/* Parse the -R option.  Full syntax:  -Rg or -Rd or -R[g|d]w/e/s/n[r] */

	/* global range format ? */
	if (item[0] == 'g' || item[0] == 'd') {
		if (item[0] == 'g') {	/* -Rg is shorthand for -R0/360/-90/90 */
			/*g_flag = 1; */
			fprintf(stderr, "Error: -Rg shorthand is not supported");
			abort();

		} else {	/* -Rd is horthand for -R-180/+180/-90/90 */
			/*      Rd_flag = 1;  *//* global flag */
			fprintf(stderr, "Error: -Rd shorthand is not supported");
			abort();
		}
		*w = -180.0, *e = 180.0;
		*s = -90.0;
		*n = +90.0;
		if (!item[1])
			return (0);
		strcpy(string, &item[1]);
	} else {
		strcpy(string, &item[0]);
	}

	/* Remove the trailing r */
	if (string[strlen(string) - 1] == 'r') {
		ll = 1;
		string[strlen(string) - 1] = '\0';
	}

	if (ll) {
		nn = sscanf(string, "%lf/%lf/%lf/%lf", w, s, e, n);
	} else {
		nn = sscanf(string, "%lf/%lf/%lf/%lf", w, e, s, n);
	}

	if (4 != nn) {
		fprintf(stderr, "Parse error option R: <%s>\n", string);
		abort();
	}

	/* If region is given then we must have w < e and s < n */
	if (*w > *e || *s > *n) {
		fprintf(stderr, "Region error w=%f e=%f s=%f n=%f \n", *w, *e,
			*s, *n);
		abort();
	}

	return EXIT_SUCCESS;
}



int parse_I_option (char *item, double *dlat, double *dlon) {

  int nn=0;
  
  /* Parse the -I option.  Full syntax:  -Idx/[dy] */
  nn = sscanf(item,"%lf/%lf",dlat,dlon);

  if (  2 != nn ){
    nn = sscanf(item,"%lf",dlat);
    if (  1 != nn ){
      fprintf(stderr,"Parse error option I: <%s>\n",item);    
      abort ();
    } else {
      *dlon = *dlat;
    }
  } 
  
  return EXIT_SUCCESS;
}

