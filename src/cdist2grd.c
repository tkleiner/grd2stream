/*
 * cdist2grd
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

// #include <getopt.h>
#include <math.h>
#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>

#include "getopt.h"
#include "debug_printf.h"
#include "gmtopts.h"
#include "grdio.h"

/**
 * print usage information and exit
 */
static void usage(void);

/**
 * print program version and exit
 */
static void version(void);

/* for logging */
const char *program_name = "cdist2grd";

/**
 * index in the primary 2D data fields
 * ixy = j * nx + i;
 */
#ifndef IJ2IND
#define IJ2IND(i, j, nx) (((size_t)(j)) * ((size_t)(nx)) + (size_t)(i))
#endif

#ifndef MIN
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#endif

#ifndef MAX
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#endif

/*****************************************************************************/
int main(int argc, char **argv) {

  double *p_x = NULL;
  double *p_y = NULL;
  double *p_grd = NULL;

  char *p_file_name = NULL;
  char *p_grd_name = NULL;

  FILE *fp = NULL;
  char line[BUFSIZ];

  int oc, err = 0;
  int R_opt = 0;
  int I_opt = 0;
  int verbose = 0;

  double xmin = 0.0, xmax = 0.0, xinc = 0.0, xlen = 0.0, xpos = 0.0,
         xdist = 0.0;
  double ymin = 0.0, ymax = 0.0, yinc = 0.0, ylen = 0.0, ypos = 0.0,
         ydist = 0.0;
  double dist = 0.0;
  size_t nx = 0, ny = 0, nn = 0;
  size_t i = 0, j = 0, k = 0;
  unsigned int nlines = 0;

  /* parse commandline args */
  while ((oc = getopt(argc, argv, "VhR:I:v")) != -1)
    switch (oc) {
    case 'R':
      /* Region */
      R_opt = 1;
      parse_R_option((char *)optarg, &xmin, &xmax, &ymin, &ymax);
      break;
    case 'I':
      /* incerement */
      I_opt = 1;
      parse_I_option((char *)optarg, &xinc, &yinc);
      break;
    case 'h':
      /* help */
      usage();
      break;
    case 'v':
      /* version */
      version();
      break;
    case 'V':
      /* verbose opttion */
      verbose++;
      break;
    case '?':
      fprintf(stderr, "Unknown option `-%c'.\n", optopt);
      return 1;
    default:
      abort();
    }

  if ((argc - optind) == 1) {
    p_grd_name = (char *)argv[optind];
  } else if ((argc - optind) == 2) {
    p_file_name = (char *)argv[optind];
    p_grd_name = (char *)argv[optind + 1];
  } else {
    usage();
  }

  xlen = xmax - xmin;
  ylen = ymax - ymin;
  if (xlen * ylen <= 0.0) {
    fprintf(stderr, "xmin: %.3f xmax: %.3f\n", xmin, xmax);
    fprintf(stderr, "ymin: %.3f ymax: %.3f\n", ymin, ymax);
    return EXIT_FAILURE;
  }

  if (remainder(xlen, xinc) != 0.0) {
    fprintf(stderr, "xlen: %.3f xinc: %.3f\n", xlen, xinc);
    return EXIT_FAILURE;

  } else {
    nx = (size_t)(xlen / xinc) + 1;
  }

  if (remainder(ylen, yinc) != 0.0) {
    fprintf(stderr, "ylen: %.3f yinc: %.3f\n", ylen, yinc);
    return EXIT_FAILURE;
  } else {
    ny = (size_t)(ylen / yinc) + 1;
  }

  if (verbose) {
    fprintf(stderr, "Input:\n");
    fprintf(stderr, "xmin: %.3f xmax: %.3f x_inc: %f nx: %lu\n", xmin, xmax,
            xinc, nx);
    fprintf(stderr, "ymin: %.3f ymax: %.3f y_inc: %f ny: %lu\n", ymin, ymax,
            yinc, ny);
    fprintf(stderr, "verbose: %u\n", verbose);
  }

  p_x = (double *)calloc(nx, sizeof(double));
  for (k = 0; k < nx; k++)
    p_x[k] = k * xinc + xmin;

  p_y = (double *)calloc(ny, sizeof(double));
  for (k = 0; k < ny; k++)
    p_y[k] = k * yinc + ymin;

  nn = ny * nx; /* number of points in grd */
  p_grd = (double *)calloc(nn, sizeof(double));
  for (k = 0; k < nn; k++) {
    p_grd[k] = 9.999e32;
  }

  if (p_file_name == NULL) { /* Just read standard input */
    fp = stdin;
    if (verbose)
      fprintf(stderr, "Reading from standard input\n");
  } else {
    if ((fp = fopen(p_file_name, "r")) == NULL) {
      fprintf(stderr, "Cannot open file %s\n", p_file_name);
      return EXIT_FAILURE;
    } else {
      if (verbose)
        fprintf(stderr, "Working on file %s\n", p_file_name);
    }
  }

  while (!feof(fp) && (fgets(line, BUFSIZ, fp) != 0)) {

    nlines++;
    if (verbose) {
      fprintf(stderr, "Processing input line: %u ...\n", nlines);
    }

    if (line[0] == '#' || line[0] == '>' || line[0] == '\n') {
      continue; /* continue with next point in file */
    }

    if (2 != sscanf(line, "%lf %lf", &xpos, &ypos)) {
      fprintf(
          stderr,
          "ERROR: Mismatch between actual and expected fields near line %u\n",
          nlines);
      return EXIT_FAILURE;
    }

    if (verbose) {
      fprintf(stderr, "       x0=%.3f y0=%.3f\n", xpos, ypos);
    }

    /*****************************************************************
     * the work is starting here
     ******************************************************************/

    for (j = 0; j < ny; j++) {
      ydist = fabs(ypos - p_y[j]);

      for (i = 0; i < nx; i++) {
        xdist = fabs(xpos - p_x[i]);
        dist = sqrt(xdist * xdist + ydist * ydist);
        k = IJ2IND(i, j, nx);
        p_grd[k] = MIN(p_grd[k], dist);
      }
    }
  }

  /* write blank array  */
  grdwrite(p_grd_name, nx, ny, p_x, p_y, (void *)p_grd, NC_DOUBLE);

  if (fp != stdin)
    fclose(fp);
  if (p_x)
    free(p_x);
  if (p_y)
    free(p_y);
  if (p_grd)
    free(p_grd);
  return EXIT_SUCCESS;
} /* main */

/**
 *
 */
void usage(void) {
  fprintf(stderr,
          "NAME:\n"
          "  %s - calculates a Cartesian distance grid to the point coordinates\n",
          program_name);
  fprintf(stderr,"USAGE:\n"
          "  %s -Rw/e/s/n[r] -Idx/[dy] in.xy out.grd \n\n"
          "  in.xy     two column ASCII file containing (x,y) coordinates\n"
          "  out.grd   GMT grid file conatining the Cartesian distances to the points in in.xy\n"
          ,program_name);
  fprintf(stderr,"\nOPTIONS:\n"
          "  -Rw/e/s/n  w, e, s, and n specify the Region of interest (GMT)\n"
          "  -Idx/[dy]  dx [and optionally dy] is the grid spacing (GMT)\n"
          "  -V         verbose output\n"
          "  -v         version\n"
          "  -h         help\n\n");
  fprintf(stderr,"\nEXAMPLE:\n"
          "  inc=-I1000/1000\n"
          "  reg=-R33000/565800/-1800000/-1023600\n"
          "  blockmean -V $reg $inc RES_profiles.xy > _tmp.d\n"
          "  %s -V $reg $inc _tmp.d $grd\n"
          "  \n",program_name);
  exit(0);
}

/**
 *
 */
void version(void) {
  fprintf(stderr, "This is %s version %s (%s).\n", PACKAGE_NAME,
          PACKAGE_VERSION, __DATE__);
  exit(0);
}
