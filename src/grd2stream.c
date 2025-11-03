#ifndef LAST_UPDATE
#define LAST_UPDATE "Time-stamp: <2025-10-28 10:02:27 (tkleiner)>"
#endif

/*
 * grd2stream
 * reads two 2-D gridded files which represents the  x-  and  y-
 * components  of a vector field and produces a streamline polygon with
 * starting at point x0,y0 reading from stdin or file
 */

/*
 * @todo: get a better guess for step size
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "none"
#endif

#ifndef PACKAGE_NAME
#define PACKAGE_NAME "grd2stream"
#endif
#include <ctype.h>
#include <math.h>
#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
//#include <unistd.h>  // for system getopt.h
#include <errno.h>
#include <string.h>
#include <limits.h> // for UINT_MAX


#include "debug_printf.h"
#include "getopt.h"
#include "grdio.h"
#include "log.h"

#if ENABLE_GMT_API
#include <gmt.h>
#define GRDREAD grdread_gmt
#else
#define GRDREAD grdread
#endif

#define MAXSTEPS 10000
#define SUBSTEPS 2

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

/** bi-linear interpolation
 *
 * using grid locations (i, j), (i + 1, j), (i, j + 1) and (i + 1, j + 1)
 *
 */
static int interp2(size_t nx, size_t ny, double *p_x, double *p_y, double *p_vx, double *p_vy, double xi, double yi,
                   double *p_vxi, double *p_vyi, double eps);

/**
 * Given an array xx[1..n], and given a value x, returns a value j such that x is between xx[j] and xx[j+1].
 * xx must be monotonic, either increasing or decreasing. j=0 or j=n is returned to indicate that x is out of range.
 * This includes the situation, when x == xx[j], and thus, xx[j + 1] is just the next neighbour
 *
 * See "Numerical Recipies in C" 2nd Edition section "3.4 How to Search an Ordered Table"
 */
static void locate(double *xx, size_t n, double x, size_t *j);

/**
 * print usage information and exit
 */
static void usage(void);

/**
 * print program version and exit
 */
static void version(void);

/* make global to use in print functions below */
int verbose = 0;
int log_breaks = 0;

// Global format strings with defaults
const char *fmt_f = "%.3f";  // for coordinates
const char *fmt_g = "%.3g";  // for other quantities
char msg[256];               // contains fprintf/printf string with formats

// Accepts formats containing '%' and one of 'f', 'e', or 'g'
static int is_valid_float_format(const char *fmt) {
  if (!fmt || !strchr(fmt, '%'))
    return 0;
  // Check for at least one type specifier
  const char *types = "fegEG";
  for (const char *p = types; *p; ++p) {
    if (strchr(fmt, *p))
      return 1;
  }
  return 0;
}

// Set float format (coordinates, general floats)
void set_float_format(const char *fmt) {
  if (is_valid_float_format(fmt)) {
    fmt_f = fmt;
  } else {
    fprintf(stderr, "Warning: invalid float format '%s'. Using default '%s'.\n", fmt, fmt_f);
  }
}

// Set general format (velocities, general numbers)
void set_general_format(const char *fmt) {
  if (is_valid_float_format(fmt)) {
    fmt_g = fmt;
  } else {
    fprintf(stderr, "Warning: invalid general format '%s'. Using default '%s'.\n", fmt, fmt_g);
  }
}

// Log messages
void log_break_nan(double x, double y, double x0, double y0) {
  if (log_breaks) {
    snprintf(msg, sizeof(msg),
             "Stop: NaN found for velocity vx or vy at position (%s, %s) "
             "for seed (%s, %s).\n",
             fmt_f, fmt_f, fmt_f, fmt_f);
    fprintf(stderr, msg, x, y, x0, y0);
  }
}

void log_break_zero(double x, double y, double x0, double y0) {
  if (log_breaks) {
    snprintf(msg, sizeof(msg),
             "Stop: Zero velocity magnitude found at position (%s, %s) "
             "for seed (%s, %s).\n",
             fmt_f, fmt_f, fmt_f, fmt_f);
    fprintf(stderr, msg, x, y, x0, y0);
  }
}

void log_break_dx(double x, double x0, double y0) {
  if (log_breaks) {
    snprintf(msg, sizeof(msg),
             "Stop: Estimated x + dx (%s) is outside valid region for seed "
             "(%s, %s).\n",
             fmt_f, fmt_f, fmt_f);
    fprintf(stderr, msg, x, x0, y0);
  }
}

void log_break_dy(double y, double x0, double y0) {
  if (log_breaks) {
    snprintf(msg, sizeof(msg),
             "Stop: Estimated y + dy (%s) is outside valid region for seed "
             "(%s, %s).\n",
             fmt_f, fmt_f, fmt_f);
    fprintf(stderr, msg, y, x0, y0);
  }
}

void log_break_stepsize(double stpsz, double xinc, double yinc, double x0, double y0) {
  if (log_breaks) {
    snprintf(msg, sizeof(msg), "Stop: Step size too small %s << (%s, %s) for seed (%s, %s)\n", fmt_f, fmt_f, fmt_f,
             fmt_f, fmt_f);
    fprintf(stderr, msg, stpsz, xinc, yinc, x0, y0);
  }
}

void log_break_maxiter(double x0, double y0) {
  if (log_breaks) {
    snprintf(msg, sizeof(msg), "Stop: Maximum number of iterations reached for seed (%s, %s)\n", fmt_f, fmt_f);
    fprintf(stderr, msg, x0, y0);
  }
}

void log_break_maxtime(double x0, double y0, unsigned long iter) {
  if (log_breaks) {
    snprintf(msg, sizeof(msg),
             "Stop: Maximum integration time reached for seed (%s, %s) "
             "after %%u iterations\n",
             fmt_f, fmt_f);
    fprintf(stderr, msg, x0, y0, (unsigned int)iter);
  }
}

void log_break_delta(double x0, double y0, unsigned long iter, double delta) {
  if (log_breaks) {
    snprintf(msg, sizeof(msg), "Stop: Invalid delta = %s for seed (%s, %s) after %%u iterations\n", fmt_g, fmt_f,
             fmt_f);
    fprintf(stderr, msg, delta, x0, y0, (unsigned int)iter);
  }
}

void test_locate(void) {
  double xiv[] = {3.5, 2.0, 6.0, -1.0};    // test locations xi
  double xv[] = {0., 1., 2., 3., 4., 5.};  // x-vector
  const size_t nx = sizeof(xv) / sizeof(double);
  double xi;
  size_t ixel;

  for (int i = 0; i < sizeof(xiv) / sizeof(double); i++) {
    xi = xiv[i];
    // printf("search xi = %.3f in x[0] = %.3f < xi < x[%lu] = %.3f\n", xi, xv[0], nx - 1, xv[nx - 1]);
    snprintf(msg, sizeof(msg), "search xi = %s in x[0] = %s < xi < x[%%lu] = %s\n", fmt_f, fmt_f, fmt_f);
    printf(msg, xi, xv[0], (unsigned long)(nx - 1), xv[nx - 1]);

    locate(xv, nx, xi, &ixel);
    // printf("found xi = %.3f -> x[%lu] = %.3f\n", xi, ixel, xv[ixel]);
    snprintf(msg, sizeof(msg), "found xi = %s -> x[%%lu] = %s\n", fmt_f, fmt_f);
    printf(msg, xi, (unsigned long)ixel, xv[ixel]);
  }
}

/** Helper: parse a double safely
 * For floating-point options (-d, -e, -T, etc.)
 */
static double parse_double(const char *s, const char *optname) {
  char *end;
  errno = 0;
  double val = strtod(s, &end);
  if (end == s || *end != '\0' || errno == ERANGE) {
    fprintf(stderr, "Invalid numeric argument for %s: '%s'\n", optname, s);
    exit(EXIT_FAILURE);
  }
  return val;
}

/** Helper: parse a long safely
 * for integer options (-k, -n, etc.)
 * Replace atoi() with strtol() or strtoul() (for unsigned)
 */
static long parse_long(const char *s, const char *optname) {
  char *end;
  errno = 0;
  long val = strtol(s, &end, 10);
  if (end == s || *end != '\0' || errno == ERANGE) {
    fprintf(stderr, "Invalid integer argument for %s: '%s'\n", optname, s);
    exit(EXIT_FAILURE);
  }
  return val;
}


/** Helper: parse a long and verify it fits in unsigned int */
static unsigned int parse_uint(const char *s, const char *optname) {
  char *end;
  errno = 0;
  long val = strtol(s, &end, 10);

  if (end == s || *end != '\0') {
    fprintf(stderr, "Invalid integer argument for %s: '%s'\n", optname, s);
    exit(EXIT_FAILURE);
  }
  if (errno == ERANGE) {
    fprintf(stderr, "Numeric value out of range for %s: '%s'\n", optname, s);
    exit(EXIT_FAILURE);
  }
  if (val < 0) {
    fprintf(stderr, "Negative value not allowed for %s: '%s'\n", optname, s);
    exit(EXIT_FAILURE);
  }
  if ((unsigned long)val > UINT_MAX) {
    fprintf(stderr, "Value too large for %s (max %u): '%s'\n",
            optname, UINT_MAX, s);
    exit(EXIT_FAILURE);
  }

  return (unsigned int)val;
}


/* for logging */
const char *program_name = PACKAGE_NAME;

/*****************************************************************************/
int main(int argc, char **argv) {
  size_t nx = 0, ny = 0;
  double xmin, xmax, ymin, ymax, x_inc, y_inc;

  unsigned long i, j, iter, maxiter = MAXSTEPS;
  unsigned int nlines = 0, npoly = 0;

  double y0 = 0.0, x0 = 0.0; /* start positions */
  double yi, xi;
  double yt, xt;
  double vxi = 0.0, vyi = 0.0, uv; /* vx, vy and magnitude of velocity */

  double dx0, dx1, dx2, dx3 = 0.0;
  double dy0, dy1, dy2, dy3 = 0.0;
  double dx, dy;
  double ex = 0.0, ey = 0.0; /* local error */

  double lim;

  double itime;                            /* stream-line integration time in time units given by thge
                                              velocity field */
  double maxtime;                          /* same units as itime, less equal zero indicates an error */
  double dist, dout, delta, delta_initial; /* unit m ???*/
  double dir = 1.0;                        /* direction (1.. forward, -1..backward) */
  unsigned int freq = 2;                   /* default freq = 2 samples per grid cell */

  double *p_x = NULL;
  double *p_y = NULL;
  double *p_vx = NULL;
  double *p_vy = NULL;

  /* x and y coordinates of the coarse grid */
  double *p_xc = NULL;
  double *p_yc = NULL;
  int *blank = NULL;
  size_t nbx = 0, nby = 0;
  size_t ib = 0, jb = 0;
  double xb_inc, yb_inc;
  double density = .1; /* 10% */

  char *p_vx_name = NULL;
  char *p_vy_name = NULL;
  char *p_file_name = NULL;

  size_t nxm = 0, nym = 0;
  size_t im = 0, jm = 0;
  double xm_inc = 0.0, ym_inc = 0.0;
  char *p_mask_name = NULL;
  char *p_blank_name = NULL;
  double *p_xm = NULL;
  double *p_ym = NULL;
  double *p_mask = NULL; /* read as double for now */
  double mask_val = 0.0;

  FILE *fp = NULL;
  char line[BUFSIZ];

  int oc, err = 0;
  int b_opt = 0; /* backward steps */
  int d_opt = 0; /* delta */
  int k_opt = 4; /* Runge Kutta 4 */
  int l_opt = 0; /* print also u,v in columns 4,5 */
  int t_opt = 0; /* print integration time in column 6*/
  int T_opt = 0; /* stop after given time */

  /* experimental options */
  int L_opt = 0;     /* 3col input */
  int D_opt = 0;     /* check distance */
  int B_opt = 0;     /* write blank file: this is experimental */
  int M_opt = 0;     /* read MASK */
  double eps = 1e-3; /* relative distance error for a point to be considered at the grid */

  /* Initialize logging facility */
  /*  log_initialize(opt_nodaemon ? LOG_TO_STDERR : LOG_TO_SYSLOG);*/

#if 0
  (void)log_initialize(LOG_TO_STDERR);
  (void)log_set_debug(DEBUG_TRACE_NONE); /* or DEBUG_TRACE_ALL,... */
#endif

  /* parse commandline args */
  while ((oc = getopt(argc, argv, "tbd:lhvk:n:f:VLDrM:B:T:e:F:G:")) != -1)
    switch (oc) {
      case 'e':
        eps = parse_double(optarg, "-e");
        break;
      case 't':
        /* report time */
        t_opt = 1;
        break;
      case 'T':
        /* maximum integration time */
        T_opt = 1;
        maxtime = parse_double(optarg, "-T");
        break;
      case 'b':
        /* go backward */
        b_opt = 1;
        break;
      case 'l':
        /* long output */
        l_opt = 1;
        break;
      case 'L':
        /* long input */
        L_opt = 1;
        break;
      case 'D':
        /* check spacing */
        D_opt = 1;
        break;
      case 'd':
        /* output step size */
        d_opt = 1;
        // dout = (double)atof(optarg);
        dout = parse_double(optarg, "-d");
        break;
      case 'F':
        /* float format */
        set_float_format(optarg);
        break;
      case 'G':
        /* data format */
        set_general_format(optarg);
        break;
      case 'h':
        /* help */
        usage();
        break;
      case 'v':
        /* version */
        version();
        break;
      case 'k':
        /* stepping */
        k_opt = (int)parse_long(optarg, "-k");
        break;
      case 'n':
        /* maximum number of steps allowed */
        // maxiter = (unsigned int) parse_long(optarg, "-n");
        maxiter = parse_uint(optarg, "-n");
        break;
      case 'f':
        /* read initial points from file */
        p_file_name = optarg;
        break;
      case 'M':
        /* read mask */
        M_opt = 1;
        p_mask_name = optarg;
        break;
      case 'B':
        /* read mask */
        B_opt = 1;
        p_blank_name = optarg;
        break;
      case 'V':
        /* verbose opttion */
        verbose++;
        break;
      case 'r':
        /* verbose opttion */
        log_breaks = 1;
        break;
      case '?':
        fprintf(stderr, "Unknown option `-%c'.\n", optopt);
        return 1;
      default:
        abort();
    }

  if ((argc - optind) != 2)
    usage();

  p_vx_name = argv[optind];
  p_vy_name = argv[optind + 1];

  /* Order of integration */
  if ((k_opt < 1 && k_opt > 4) || (k_opt == 3)) {
#if 0
    (void) log_printf(OM_LOG_CRITICAL, "Wrong selection option k_opt");
    (void) log_finalize();
#endif
    usage();
  }

  /* u,v output and integration time output at the same time is invalid */
  if (t_opt > 0 && l_opt > 0) {
    usage();
  }

  // set direction (upstream or downstream)
  dir = (b_opt > 0) ? -1.0 : 1.0;

  /* read the grid files */
  err = GRDREAD(p_vx_name, &nx, &ny, &p_x, &p_y, &p_vx);
  err += GRDREAD(p_vy_name, &nx, &ny, &p_x, &p_y, &p_vy);
  if (M_opt > 0 && p_mask_name != NULL) {
    err += GRDREAD(p_mask_name, &nxm, &nym, &p_xm, &p_ym, &p_mask);
    xm_inc = p_xm[1] - p_xm[0];
    ym_inc = p_ym[1] - p_ym[0];
    if (verbose) {
      fprintf(stderr, "Input mask:\n");
      snprintf(msg, sizeof(msg), "xmin: %s xmax: %s x_inc: %s nx: %%lu\n", fmt_f, fmt_f, fmt_f);
      fprintf(stderr, msg, p_xm[0], p_xm[nxm - 1], xm_inc, (unsigned long)nxm);
      snprintf(msg, sizeof(msg), "ymin: %s ymax: %s y_inc: %s ny: %%lu\n", fmt_f, fmt_f, fmt_f);
      fprintf(stderr, msg, p_ym[0], p_ym[nym - 1], ym_inc, (unsigned long)nym);
    }
  }
  if (err != 0) {
    fprintf(stderr, "error reading grd files\n");
    exit(EXIT_FAILURE);
  }

  /* grid limits */
  xmin = p_x[0];
  xmax = p_x[nx - 1];
  ymin = p_y[0];
  ymax = p_y[ny - 1];
  x_inc = (xmax - xmin) / (double)(nx - 1);
  y_inc = (ymax - ymin) / (double)(ny - 1);

  /*
   * Blank array: This is the heart of the algorithm. It begins life
   * zeroed, but is set to one when a streamline passes through each
   * box. Then streamlines are only allowed to pass through zeroed
   * boxes. The lower resolution of this grid determines the
   * approximate spacing between trajectories.
   */
  nbx = (int)(nx * density + 1);
  nby = (int)(ny * density + 1);

  blank = (int *)calloc(nby * nbx, sizeof(int));
  /* box spacing */
  xb_inc = (p_x[nx - 1] - p_x[0]) / (double)(nbx - 1);
  yb_inc = (p_y[ny - 1] - p_y[0]) / (double)(nby - 1);

  p_xc = (double *)calloc(nbx, sizeof(double));
  p_yc = (double *)calloc(nby, sizeof(double));

  for (i = 0; i < nbx; i++) {
    p_xc[i] = p_x[0] + i * xb_inc;
  }
  for (j = 0; j < nby; j++) {
    p_yc[j] = p_y[0] + j * yb_inc;
  }
  for (i = 0; i < nbx * nby; i++) {
    blank[i] = IO_NAN;
  }

  /* default freq = 2 samples per grid cell */
  delta_initial = MIN(x_inc, y_inc) / ((double)freq);

  /*
   * TODO: makes no sense at the moment
   */
  if (d_opt) {
    freq = 1;
  } else {
    freq = 1;
    dout = MIN(x_inc, y_inc) / ((double)5);
  }
  delta_initial = dout / ((double)freq);

  if (verbose) {
    fprintf(stderr, "Input:\n");
    snprintf(msg, sizeof(msg), "xmin: %s xmax: %s x_inc: %s nx: %%lu\n", fmt_f, fmt_f, fmt_f);
    fprintf(stderr, msg, xmin, xmax, x_inc, nx);
    snprintf(msg, sizeof(msg), "ymin: %s ymax: %s y_inc: %s ny: %%lu\n", fmt_f, fmt_f, fmt_f);
    fprintf(stderr, msg, ymin, ymax, y_inc, ny);
    snprintf(msg, sizeof(msg), "d_out: %s d_inc: %s RK: %%d freq: %%u\n", fmt_f, fmt_f);
    fprintf(stderr, msg, dout, delta_initial, k_opt, freq);
    fprintf(stderr, "verbose: %u\n", verbose);
  }

  if (verbose > 1) {
    fprintf(stderr, "Coarse grid:\n");
    snprintf(msg, sizeof(msg), "xmin: %s xmax: %s x_inc: %s nx: %%lu\n", fmt_f, fmt_f, fmt_f);
    fprintf(stderr, msg, p_xc[0], p_xc[nbx - 1], xb_inc, nbx);
    snprintf(msg, sizeof(msg), "ymin: %s ymax: %s y_inc: %s ny: %%lu\n", fmt_f, fmt_f, fmt_f);
    fprintf(stderr, msg, p_yc[0], p_yc[nby - 1], yb_inc, nby);
  }

  if ((delta_initial > x_inc) || (delta_initial > y_inc)) {
    snprintf(msg, sizeof(msg), "WARN: Step size too large: %s > (%s, %s)\n", fmt_f, fmt_f, fmt_f);
    fprintf(stderr, msg, delta_initial, x_inc, y_inc);
  }

  if (p_file_name == NULL) { /* Just read standard input */
    fp = stdin;
    if (verbose)
      fprintf(stderr, "Reading from standard input\n");
  } else {
    if ((fp = fopen(p_file_name, "r")) == NULL) {
      fprintf(stderr, "Cannot open file %s\n", p_file_name);
      return EXIT_FAILURE;
    }
    if (verbose)
      fprintf(stderr, "Working on file %s\n", p_file_name);
  }

  while (!feof(fp) && (fgets(line, BUFSIZ, fp) != 0)) {
    nlines++;
    if (verbose) {
      fprintf(stderr, "Processing input line: %u ...\n", nlines);
    }

    if (line[0] == '#' || line[0] == '>' || line[0] == '\n') {
      /* if(verbose) fprintf(stderr,"skip header\n"); */
    } else {
      if (L_opt) {
        if (3 != sscanf(line, "%lf %lf %lf", &x0, &y0, &dir)) {
          fprintf(stderr,
                  "ERROR: Mismatch between actual and expected fields near "
                  "line %u\n",
                  nlines);
          return EXIT_FAILURE;
        }
        dir = (dir < 0.0) ? -1.0 : 1.0;

      } else {
        if (2 != sscanf(line, "%lf %lf", &x0, &y0)) {
          fprintf(stderr,
                  "ERROR: Mismatch between actual and expected fields near "
                  "line %u\n",
                  nlines);
          return EXIT_FAILURE;
        }
        if (verbose > 1) {
          snprintf(msg, sizeof(msg), "       x0=%s y0=%s\n", fmt_f, fmt_f);
          fprintf(stderr, msg, x0, y0);
        }
      }

      /*****************************************************************
       *
       * the work is starting here
       *
       ******************************************************************/

      if ((x0 > p_x[nx - 1]) || (y0 > p_y[ny - 1]) || (x0 < p_x[0]) || (y0 < p_y[0])) {
        if (log_breaks) {
          snprintf(msg, sizeof(msg),
                   "Stop: Seed (%s, %s) is outside valid region "
                   "(xmin=%s, xmax=%s, ymin=%s, ymax=%s)!\n",
                   fmt_f, fmt_f, fmt_f, fmt_f, fmt_f, fmt_f);
          fprintf(stderr, msg, x0, y0, p_x[0], p_x[nx - 1], p_y[0], p_y[ny - 1]);
        }
        continue; /* continue with next point in file */
      }

      /* insert ogr2ogr header: https://docs.generic-mapping-tools.org/6.3/cookbook/ogrgmt-format.html */
      if (npoly == 0) {
        printf(
            "# @VGMT1.0 @GLINESTRING \n"
            "# @Nname|id\n"
            "# @Tstring|integer\n"
            "# FEATURE_DATA\n");
      }

      npoly++;
      printf(">\n# @D\"streamline %u\"|%u\n", npoly, npoly);

      dist = 0.0;
      itime = 0.0;
      delta = delta_initial; /* reset delta for next streamline, important for
                                -T option */

      /* simple step iteration */
      xi = x0;
      yi = y0;

      /*
       * Points at current streamline
       */
      for (iter = 0; iter < maxiter; iter++) {
        // if (iter == 0) {
        //   dx0 = dx1 = dx2 = dx3 = 0.0;
        //   dy0 = dy1 = dy2 = dy3 = 0.0;
        // }

        /*
         * STEP 0 (get initial velocity data at requested point)
         */
        (void)interp2(nx, ny, p_x, p_y, p_vx, p_vy, xi, yi, &vxi, &vyi, eps);
        if (verbose > 1) {
          snprintf(msg, sizeof(msg), "#S: x = %s, y = %s, vx = %s, vy = %s\n", fmt_f, fmt_f, fmt_g, fmt_g);
          fprintf(stderr, msg, xi, yi, vxi, vyi);
        }

        /* check if mask reached */
        if (M_opt > 0) {
          im = (size_t)(((xi - p_xm[0]) / xm_inc));
          jm = (size_t)(((yi - p_ym[0]) / ym_inc));
          mask_val = p_mask[jm * nxm + im];
        }

        /* this is the output */
        if (!(iter % freq)) {
          if (l_opt) {
            snprintf(msg, sizeof(msg), "%s %s %s %s %s\n", fmt_f, fmt_f, fmt_f, fmt_g, fmt_g);
            printf(msg, xi, yi, dist, vxi, vyi);
          } else if (t_opt) {
            snprintf(msg, sizeof(msg), "%s %s %s %s %s %s\n", fmt_f, fmt_f, fmt_f, fmt_g, fmt_g, fmt_f);
            printf(msg, xi, yi, dist, vxi, vyi, itime);
          } else if (M_opt > 0) {
            snprintf(msg, sizeof(msg), "%s %s %s %s\n", fmt_f, fmt_f, fmt_f, fmt_f);
            printf(msg, xi, yi, dist, mask_val);
          } else {
            snprintf(msg, sizeof(msg), "%s %s %s\n", fmt_f, fmt_f, fmt_f);
            printf(msg, xi, yi, dist);
          }
        }

        if (isnan(mask_val) || mask_val == 0.0) {
          /* continue until valid mask value is found */
        } else {
          /* break valid mask point */
          /* print start position and mask of final position */
          printf("#M# %.3f %.3f %.3f\n", x0, y0, mask_val);
          break;
        }

        if (isnan(vxi) || isnan(vyi)) {
          log_break_nan(xi, yi, x0, y0);
          if (M_opt) {
            snprintf(msg, sizeof(msg), "#M# %s %s NaN\n", fmt_f, fmt_f);
            printf(msg, x0, y0);
          }
          break;
        }
        if ((uv = sqrt(vxi * vxi + vyi * vyi)) <= 0.0) {
          log_break_zero(xi, yi, x0, y0);
          if (M_opt) {
            snprintf(msg, sizeof(msg), "#M# %s %s NaN\n", fmt_f, fmt_f);
            printf(msg, x0, y0);
          }
          break;
        }

        /* prepare the advance */
        if (T_opt) {
          /* next time would be t + dt, dt = delta/uv */
          double time_excess = itime + delta / uv - maxtime;
          if (time_excess > 0.0) {
            /* reduce time step length to match target time
             * new_dt = dt - time_excess = delta/uv - time_excess
             * new_delta = new_dt * uv
             */
            if (verbose) {
              snprintf(msg, sizeof(msg), "#T# reducing delta = %s ", fmt_g);
              printf(msg, delta);
            }
            delta = (delta / uv - time_excess) * uv;
            if (verbose) {
              snprintf(msg, sizeof(msg), " -> %s to match maxtime = %s \n", fmt_f, fmt_g);
              printf(msg, delta, maxtime);
            }
          } else if (time_excess < 0.0) {
            /* regular time step required */
          } else {
            log_break_maxtime(x0, y0, iter);
            break;
          }
        }

        if (delta <= 0.0) {
          log_break_delta(x0, y0, iter, delta);
          break;
        }

        /* advance now */
        dist += (delta * dir);
        dx0 = dir * delta * vxi / uv; /* if vxi > 0 && dir <= 0, then dx0 is already negative */
        dy0 = dir * delta * vyi / uv; /* if vyi > 0 && dir <= 0, then dy0 is already negative */
        itime += delta / uv;

        /*
         * check initial guess
         */
        if (x0 + dx0 > p_x[nx - 1] || x0 + dx0 < p_x[0]) {
          log_break_dx(x0 + dx0, x0, y0);
          if (M_opt) {
            snprintf(msg, sizeof(msg), "#M# %s %s NaN\n", fmt_f, fmt_f);
            printf(msg, x0, y0);
          }
          break;
        }
        if (y0 + dy0 > p_y[ny - 1] || y0 + dy0 < p_y[0]) {
          log_break_dy(y0 + dy0, x0, y0);
          if (M_opt) {
            snprintf(msg, sizeof(msg), "#M# %s %s NaN\n", fmt_f, fmt_f);
            printf(msg, x0, y0);
          }
          break;
        }
        if (verbose > 2) {
          snprintf(msg, sizeof(msg), "#\tRK0: x=%s, y=%s, vx = %s, vy=%s, dx=%s, dy=%s\n", fmt_f, fmt_f, fmt_g, fmt_g,
                   fmt_f, fmt_f);
          fprintf(stderr, msg, xi, yi, vxi, vyi, dx0, dy0);
        }

        /*  check stepsize  */
        if (iter > 0 && k_opt == 4) {
          ex = dx3 / 6.0 - vxi * delta / 6.0;
          ey = dy3 / 6.0 - vyi * delta / 6.0;
        }

        /* early break for Euler method. Continue with next i */
        if (k_opt == 1) {
          xi += dx0;
          yi += dy0;
          continue;
        }

        /*
         * RK-STEP 1
         */
        xt = xi + dx0 / 2.0;
        yt = yi + dy0 / 2.0;
        (void)interp2(nx, ny, p_x, p_y, p_vx, p_vy, xt, yt, &vxi, &vyi, eps);
        if (isnan(vxi) || isnan(vyi)) {
          log_break_nan(xt, yt, x0, y0);
          if (M_opt) {
            snprintf(msg, sizeof(msg), "#M# %s %s NaN\n", fmt_f, fmt_f);
            printf(msg, x0, y0);
          }
          break;
        }
        if ((uv = sqrt(vxi * vxi + vyi * vyi)) <= 0.0) {
          log_break_zero(xt, yt, x0, y0);
          if (M_opt) {
            snprintf(msg, sizeof(msg), "#M# %s %s NaN\n", fmt_f, fmt_f);
            printf(msg, x0, y0);
          }
          break;
        }
        dx1 = dir * delta * vxi / uv;
        dy1 = dir * delta * vyi / uv;
        if (verbose > 2) {
          snprintf(msg, sizeof(msg), "#\tRK1: x=%s, y=%s, vx = %s, vy=%s, dx=%s, dy=%s\n", fmt_f, fmt_f, fmt_g, fmt_g,
                   fmt_f, fmt_f);
          fprintf(stderr, msg, xt, yt, vxi, vyi, dx1, dy1);
        }

        /*
         * RK-STEP 2
         */
        xt = xi + dx1 / 2.0;
        yt = yi + dy1 / 2.0;
        (void)interp2(nx, ny, p_x, p_y, p_vx, p_vy, xt, yt, &vxi, &vyi, eps);
        if (isnan(vxi) || isnan(vyi)) {
          log_break_nan(xt, yt, x0, y0);
          if (M_opt) {
            snprintf(msg, sizeof(msg), "#M# %s %s NaN\n", fmt_f, fmt_f);
            printf(msg, x0, y0);
          }
          break;
        }
        if ((uv = sqrt(vxi * vxi + vyi * vyi)) <= 0.0) {
          log_break_zero(xt, yt, x0, y0);
          if (M_opt) {
            snprintf(msg, sizeof(msg), "#M# %s %s NaN\n", fmt_f, fmt_f);
            printf(msg, x0, y0);
          }
          break;
        }
        dx2 = dir * delta * vxi / uv;
        dy2 = dir * delta * vyi / uv;
        if (verbose > 2) {
          snprintf(msg, sizeof(msg), "#\tRK2: x=%s, y=%s, vx = %s, vy=%s, dx=%s, dy=%s\n", fmt_f, fmt_f, fmt_g, fmt_g,
                   fmt_f, fmt_f);
          fprintf(stderr, msg, xt, yt, vxi, vyi, dx2, dy2);
        }

        /*
         * RK-STEP 3
         */
        xt = xi + dx2;
        yt = yi + dy2;
        (void)interp2(nx, ny, p_x, p_y, p_vx, p_vy, xt, yt, &vxi, &vyi, eps);
        if (isnan(vxi) || isnan(vyi)) {
          log_break_nan(xt, yt, x0, y0);
          if (M_opt) {
            snprintf(msg, sizeof(msg), "#M# %s %s NaN\n", fmt_f, fmt_f);
            printf(msg, x0, y0);
          }
          break;
        }
        if ((uv = sqrt(vxi * vxi + vyi * vyi)) <= 0.0) {
          log_break_zero(xt, yt, x0, y0);
          if (M_opt) {
            snprintf(msg, sizeof(msg), "#M# %s %s NaN\n", fmt_f, fmt_f);
            printf(msg, x0, y0);
          }
          break;
        }
        dx3 = dir * delta * vxi / uv;
        dy3 = dir * delta * vyi / uv;
        if (verbose > 2) {
          snprintf(msg, sizeof(msg), "#\tRK3: x=%s, y=%s, vx = %s, vy=%s, dx=%s, dy=%s\n", fmt_f, fmt_f, fmt_g, fmt_g,
                   fmt_f, fmt_f);
          fprintf(stderr, msg, xt, yt, vxi, vyi, dx3, dy3);
        }

        /*
         * final RK-STEP update
         * x_{n+1} = x_{n} + h/6 (k_1 + 2*k_2 + 2*k_3 + k_4)
         *         = x_{n} + dx
         */
        dx = (dx0 / 6.0 + dx1 / 3.0 + dx2 / 3.0 + dx3 / 6.0);
        dy = (dy0 / 6.0 + dy1 / 3.0 + dy2 / 3.0 + dy3 / 6.0);

        /*
         * check final step
         */
        if ((xi + dx > p_x[nx - 1]) || (xi + dx < p_x[0])) {
          log_break_dx(xi + dx, x0, y0);
          if (M_opt) {
            snprintf(msg, sizeof(msg), "#M# %s %s NaN\n", fmt_f, fmt_f);
            printf(msg, x0, y0);
          }
          break;
        }
        if ((yi + dy > p_y[ny - 1]) || (yi + dy < p_y[0])) {
          log_break_dy(yi + dy, x0, y0);
          if (M_opt) {
            snprintf(msg, sizeof(msg), "#M# %s %s NaN\n", fmt_f, fmt_f);
            printf(msg, x0, y0);
          }
          break;
        }

        /* final update of position */
        xi += dx;
        yi += dy;

        /*
         * check if already visited
         *
         * is not working now for backward streamlines! Why?
         */

        /* compute index in blank array */
        ib = (size_t)(((xi - p_x[0]) / xb_inc));
        jb = (size_t)(((yi - p_y[0]) / yb_inc));
        blank[jb * nbx + ib] = npoly;

        /*
         * START TESTING HERE
         */

#if 0
        if (verbose) {
          fprintf(stderr, " xi=%15.5f yi=%15.5f\n", xi, yi);
          fprintf(stderr, " ib=%d jb=%d\n", ib, jb);
          fprintf(stderr, " xb=%15.5f yb=%15.5f\n", p_xc[ib], p_yc[jb]);
        }

        if (D_opt) {
          if (ib < 0 || ib > nbx - 1) {
            fprintf(stdout, "# error: %lu %lu %f\n", ib, nbx, (xi - p_x[0]));
            break;
          }
          if (jb < 0 || jb > nby - 1) {
            fprintf(stdout, "# error: %lu %lu %f\n", ib, nby, (yi - p_y[0]));
            break;
          }

          if ((blank[jb * nbx + ib] > 0) && (blank[jb * nbx + ib] != npoly)) {
            if (verbose) {
              fprintf(stdout, "# blank: %lu %lu %d\n", ib, jb, blank[jb * nbx + ib]);
              snprintf(msg, sizeof(msg), "# Stop: x = %s, y = %s already visited by %%d\n", fmt_f, fmt_f);
              fprintf(stdout, msg, xi, yi, id);
            }
            break;

          } else {
            blank[jb * nbx + ib] = npoly;
          }
        }
#endif
        /*
         * END TESTING HERE
         */
        if (verbose > 1) {
          snprintf(msg, sizeof(msg), "#E: x = %s, y = %s, dx = %s, dy = %s\n", fmt_f, fmt_f, fmt_f, fmt_f);
          fprintf(fp, msg, xi, yi, dx, dy);
        }

        lim = sqrt(dx * dx + dy * dy);
        if (lim * 1000.0 < (MIN(x_inc, y_inc))) {
          log_break_stepsize(lim, x_inc, y_inc, x0, y0);
          if (M_opt) {
            snprintf(msg, sizeof(msg), "#M# %s %s NaN\n", fmt_f, fmt_f);
            printf(msg, x0, y0);
          }
          break;
        }
      }
      /* end of work for one streamline */
      if (iter == maxiter) {
        log_break_maxiter(x0, y0);
        if (M_opt) {
          snprintf(msg, sizeof(msg), "#M# %s %s NaN\n", fmt_f, fmt_f);
          printf(msg, x0, y0);
        }
      }
    }
  }

  /* write blank array for later usage */
  if (B_opt && NULL != p_blank_name) {
    grdwrite(p_blank_name, nbx, nby, p_xc, p_yc, (void *)blank, NC_INT);
  }

  if (fp != stdin)
    fclose(fp);
  if (p_x)
    free(p_x);
  if (p_y)
    free(p_y);
  if (p_xc)
    free(p_xc);
  if (p_yc)
    free(p_yc);
  if (p_vx)
    free(p_vx);
  if (p_vy)
    free(p_vy);

  if (p_mask)
    free(p_mask);
  if (p_xm)
    free(p_xm);
  if (p_ym)
    free(p_ym);

  if (blank)
    free(blank);

#if 0
  log_finalize();
#endif

  return EXIT_SUCCESS;
} /* main */

/*****************************************************************************
 * interp2_fwd results in nan if any of the four pixel is nan
 * See also https://github.com/JuliaMath/Interpolations.jl/issues/192
 *          https://github.com/scipy/scipy/issues/11381
 */
int interp2(size_t nx, size_t ny, double *p_x, double *p_y, double *p_vx, double *p_vy, double xi, double yi,
            double *p_vxi, double *p_vyi, double eps) {
  size_t ixel = 0, iyel = 0;   /* where we are */
  size_t ixel1 = 0, iyel1 = 0; /* the next */
  size_t i00, i01, i10, i11;
  double p1, p2, q1, q2;
  double xinc, yinc;
  double xerr, yerr;

  debug_printf(DEBUG_INFO, "search xi = %.3f in x[%lu] = %.3f < xi < x[%lu] = %.3f\n", xi, 0, p_x[0], nx - 1,
               p_x[nx - 1]);
  debug_printf(DEBUG_INFO, "search yi = %.3f in y[%lu] = %.3f < yi < y[%lu] = %.3f\n", yi, 0, p_y[0], ny - 1,
               p_y[ny - 1]);

  if (xi >= p_x[nx - 1]) {
    /* right margin */
    ixel = nx - 1;
  } else if (xi <= p_x[0]) {
    /* left margin */
    ixel = 0;
  } else {
    locate(p_x, nx, xi, &ixel);
  }
  debug_printf(DEBUG_INFO, "found xi = %.3f -> x[%lu] = %.3f\n", xi, ixel, p_x[ixel]);

  if (yi >= p_y[ny - 1]) {
    /* top margin */
    iyel = ny - 1;
  } else if (yi <= p_y[0]) {
    /* bottom margin */
    iyel = 0;
  } else {
    locate(p_y, ny, yi, &iyel);
  }
  debug_printf(DEBUG_INFO, "found yi = %.3f -> y[%lu] = %.3f\n", yi, iyel, p_y[iyel]);

  ixel1 = MIN(ixel + 1, nx - 1);
  iyel1 = MIN(iyel + 1, ny - 1);

  xinc = fabs(p_x[ixel1] - p_x[ixel]);
  yinc = fabs(p_y[iyel1] - p_y[iyel]);
  debug_printf(DEBUG_INFO, "found xinc = %.3f, yinc = %.3f\n", xinc, yinc);

  /* index in two-dimensional data array */
  i00 = (iyel) * (nx) + ixel;
  i01 = (iyel) * (nx) + ixel1;
  i10 = (iyel1) * (nx) + ixel;
  i11 = (iyel1) * (nx) + ixel1;

  /* check if we are on the grid */
  xerr = fabs(xi - p_x[ixel]);
  yerr = fabs(yi - p_y[iyel]);
  debug_printf(DEBUG_INFO, "Found xerr = %.3f, yerr = %.3f\n", xerr, yerr);
  if (xerr < eps * xinc && yerr < eps * yinc) {
    debug_printf(DEBUG_INFO, "Point is on grid -> early exit interp2(...)\n");
    (*p_vxi) = p_vx[i00];
    (*p_vyi) = p_vy[i00];
    return EXIT_SUCCESS;
  }

  debug_printf(DEBUG_INFO, "vx[i00] =  %.3g\n", p_vx[i00]);
  debug_printf(DEBUG_INFO, "vx[i10] =  %.3g\n", p_vx[i10]);
  debug_printf(DEBUG_INFO, "vx[i01] =  %.3g\n", p_vx[i01]);
  debug_printf(DEBUG_INFO, "vx[i11] =  %.3g\n", p_vx[i11]);

  debug_printf(DEBUG_INFO, "vy[i00] =  %.3g\n", p_vy[i00]);
  debug_printf(DEBUG_INFO, "vy[i10] =  %.3g\n", p_vy[i10]);
  debug_printf(DEBUG_INFO, "vy[i01] =  %.3g\n", p_vy[i01]);
  debug_printf(DEBUG_INFO, "vy[i11] =  %.3g\n", p_vy[i11]);

  // linear interpolation between the lower two quadrants along the X axis
  p1 = p_vx[i00] + (xi - p_x[ixel]) * (p_vx[i10] - p_vx[i00]) / (p_x[ixel + 1] - p_x[ixel]);
  // linear interpolation between the upper two quadrants along the X axis
  q1 = p_vx[i01] + (xi - p_x[ixel]) * (p_vx[i11] - p_vx[i01]) / (p_x[ixel + 1] - p_x[ixel]);
  // linear interpolation between the two virtual points p1 and q1 along the Y axis
  (*p_vxi) = p1 + (yi - p_y[iyel]) * (q1 - p1) / (p_y[iyel + 1] - p_y[iyel]);

  p2 = p_vy[i00] + (xi - p_x[ixel]) * (p_vy[i10] - p_vy[i00]) / (p_x[ixel + 1] - p_x[ixel]);

  q2 = p_vy[i01] + (xi - p_x[ixel]) * (p_vy[i11] - p_vy[i01]) / (p_x[ixel + 1] - p_x[ixel]);

  (*p_vyi) = p2 + (yi - p_y[iyel]) * (q2 - p2) / (p_y[iyel + 1] - p_y[iyel]);

  debug_printf(DEBUG_INFO, "p1 =  %.3g -> q1 = %.3g\n", p1, q1);
  debug_printf(DEBUG_INFO, "p2 =  %.3g -> q2 = %.3g\n", p2, q2);

  return EXIT_SUCCESS;
} /* interp2 */

void locate(double *xx, size_t n, double x, size_t *j) {
  /* fast search
   * table lookup */
  size_t iu, im, il;
  int ascnd;

  il = 0;
  iu = n - 1;
  ascnd = (xx[n - 1] >= xx[0]);
  while (iu - il > 1) {
    im = (iu + il) >> 1;
    if (x >= xx[im] == ascnd)
      il = im;
    else
      iu = im;
  }

  if (x == xx[0])
    *j = 0;
  else if (x == xx[n - 1])
    *j = n - 1;
  else
    *j = il;

} /* locate */

/**
 *
 */
void usage(void) {
  fprintf(stderr,
          "NAME:\n"
          "  %s - generate stream lines based on GMT grid files\n",
          program_name);
  fprintf(stderr,
          "USAGE:\n"
          "  %s v_x.grd v_y.grd -f xyfile \n\n"
          "  v_x.grd & v_y.grd   grid files with the 2 vector components\n"
          "  xyfile              two column ASCII file containing seed point coordinates (x0,y0)\n",
          program_name);
  fprintf(stderr,
          "\nOPTIONS:\n"
          "  -b                  backward steps\n"
          "  -d inc              stepsize\n"
          "  -f xyfile           file with seed points\n"
          /* "  -k                  select stepping method (default: RK4)\n" */
          "  -l                  output format: 'x y dist v_x v_y' (5 columns)\n"
          "  -t                  output format: 'x y dist v_x v_y time' (6 columns)\n"
          "  -T maxtime          maximum integration time (default: none)\n"
          "  -n maxsteps         maximum number of steps (default: %d)\n"
          "  -V                  verbose output\n"
          "  -r                  report why a streamline stopped to stderr "
          "(default: off)\n"
          "  -F fmt_f            format for coordinates (default: %s)\n"
          "  -G fmt_g            format for data (default: %s)\n"
          "  -v                  version\n"
          "  -h                  help\n\n",
          MAXSTEPS, fmt_f, fmt_g);
  fprintf(stderr,
          "\nDESCRIPTION:\n"
          "  %s - reads (x0,y0) pairs from standard input or xyfile (-f option)\n"
          "  and generates polylines in multiple segment mode each starting at x0,y0.\n"
          "  Output: 'x y dist' (3 columns) to stdout.\n",
          program_name);
  fprintf(stderr,
          "\nNOTES:\n"
          "  Units of the x- and y-direction must match the spatial unit of the "
          "velocity\n"
          "  to avoid unexpected results, e.g. m and m/a. \n"
          "  If the velocity is given in e.g. m/a, m/d or m/s \n"
          "  '-T 1' stops after one year, one day or 1 second respectively.\n"
          "  Units are not converted in %s, thus keep the units consistent!\n"
          "\n"
          "  If you have more than one seed point, '-f xyfile' option has better performance,\n"
          "  because the grid files need to be read only once.\n",
          program_name);
#if ENABLE_GMT_API
  fprintf(stderr,
          "\n"
          "  Use 'gmt grdconvert' without further arguments to "
          "retrieve a list of supported grid formats\n");
#endif

  fprintf(stderr,
          "\nEXAMPLE:\n"
          "  echo \"0 0\" | %s vx.grd vy.grd | psxy -m -R -J ... \n"
          "  \n",
          program_name);

  exit(0);
}

/**
 *
 */
void version(void) {
  unsigned int gmt_major, gmt_minor, gmt_patch;

  fprintf(stderr, "This is %s version %s (%s)", PACKAGE_NAME, PACKAGE_VERSION, __DATE__);

#if ENABLE_GMT_API
  (void)GMT_Get_Version(NULL, &gmt_major, &gmt_minor, &gmt_patch);
  fprintf(stderr, ", GMT API version %u.%u.%u", gmt_major, gmt_minor, gmt_patch);
#endif
  fprintf(stderr, ".\n");

  exit(0);
}
