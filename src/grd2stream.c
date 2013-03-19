#ifndef LAST_UPDATE
#define LAST_UPDATE "Time-stamp: <2013-03-19 11:45:56 (tkleiner)>"
#endif

/*
 * grd2stream
 * reads two 2-D gridded files which represents the  x-  and  y-
 * components  of a vector field and produces a stream line polygon with
 * starting at point x0,y0 reading from stdin or file
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.2.2"
#endif

#ifndef PACKAGE_NAME
#define PACKAGE_NAME "grd2stream"
#endif


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "grdio.h"
#include "debug_printf.h"
#include "log.h"


#define MAXSTEPS 10000
#define SUBSTEPS 2

#define MIN(X,Y) ( ((X) < (Y)) ? (X) : (Y) )
#define MAX(X,Y) ( ((X) > (Y)) ? (X) : (Y) )

#define SQRT sqrtf


/**
 *
 */
static int interp2(size_t nx, size_t ny, float *p_x, float *p_y,
                   float *p_vx, float *p_vy, float xi, float yi,
                   float *p_vxi,float *p_vyi);

/**
 * returns a value j such that x is between xx[j] and xx[j+1]
 */
static void locate(float* xx, size_t n, float x, size_t *j);

/**
 * print usage information and exit
 */
static void usage(void);

/**
 * print program version and exit
 */
static void version(void);


/* for logging */
const char *program_name = PACKAGE_NAME;

/*****************************************************************************/
int main( int argc, char** argv )
{
  size_t nx=0,ny=0;
  float xmin,xmax,ymin,ymax,x_inc,y_inc;

  unsigned long i,cnt, n=MAXSTEPS;
  unsigned int nlines=0,npoly=0;

  float y0=0.0f,x0=0.0f;   /* start positions */
  float yi=0.0f,xi=0.0f;
  float yt=0.0f,xt=0.0f;
  float vxi=0.0f,vyi=0.0f,uv=0.0f;

  float dx0=0.0f,dx1=0.0f,dx2=0.0f,dx3=0.0f;
  float dy0=0.0f,dy1=0.0f,dy2=0.0f,dy3=0.0f;
  float ex=0.0f,ey=0.0f; /* local error */
  float dx=0.0f,dy=0.0f; /* local error */
  float lim=0.0f;

  float dist,dout,delta = 1000.0; /* unit m ???*/
  float dir=1.0;              /* direction (1.. forward, -1..backward) */
  unsigned int freq = 1;

  float *p_x=NULL;
  float *p_y=NULL;
  float *p_vx=NULL;
  float *p_vy=NULL;

  char* p_vx_name = NULL;
  char* p_vy_name = NULL;
  char* p_file_name = NULL;

  FILE* fp = NULL;
  char line[BUFSIZ];

  int oc,err=0;
  int verbose = 0;
  int b_opt = 0; /* backward steps */
  int d_opt = 0; /* delta */
  int k_opt = 4; /* Runge Kutta 4 */
  int l_opt = 0; /* print also u,v in columns 4,5 */

  /* Initialize logging facility */
  /*  log_initialize(opt_nodaemon ? LOG_TO_STDERR : LOG_TO_SYSLOG);*/

#if 0
  (void)log_initialize(LOG_TO_STDERR);
  (void)log_set_debug(DEBUG_TRACE_NONE); /* or DEBUG_TRACE_ALL,... */
#endif

  /* parse commandline args */
  while ((oc = getopt (argc, argv, "bd:lhvk:n:f:V")) != -1)
    switch (oc) {
    case 'b': 
      /* go backward */
      b_opt=1; break;
    case 'l': 
      /* long output */
      l_opt=1; break;
    case 'd': 
      /* output step size */
      d_opt=1; dout= (float) atof(optarg); break;
    case 'h': 
      /* help */
      usage(); break;
    case 'v':
      /* version */
      version(); break;
    case 'k':
      /* stepping */
      k_opt = atoi(optarg); break;
    case 'n': 
      /* maximum number of steps allowed */
      n = (unsigned int)atoi(optarg); break; 
    case 'f': 
      /* read initial points from file */
      p_file_name = (optarg); break; 
    case 'V': 
      /* verbose opttion */
      verbose++; break;
    case '?':
      fprintf (stderr, "Unknown option `-%c'.\n", optopt);
      return 1;
    default:
      abort ();
    }
  
  if ((argc - optind) != 2) usage();
  p_vx_name = (char*)argv[optind];
  p_vy_name = (char*)argv[optind+1];
  
  /* Order of integration */
  if ( (k_opt < 1 && k_opt > 4) || (k_opt == 3) ) {
#if 0
    (void) log_printf(OM_LOG_CRITICAL, "Wrong selection option k_opt");
    (void) log_finalize();
#endif
    usage();
  }

  
  if(b_opt) {
    dir = -1.0;
  }

  /* read the grid files */
  err  = grdread (p_vx_name, &nx, &ny, &p_x, &p_y, &p_vx);
  err += grdread (p_vy_name, &nx, &ny, &p_x, &p_y, &p_vy);
  if (err != 0) {
    fprintf(stderr,"error reading grd files\n");
    exit(EXIT_FAILURE);
  }

  /* grid limits */
  xmin = p_x[0];
  xmax = p_x[nx-1];
  ymin = p_y[0];
  ymax = p_y[ny-1];
  x_inc= (xmax - xmin)/ ( (float) (nx - 1) );
  y_inc= (ymax - ymin)/ ( (float) (ny - 1) );

  /* default freq = 2 samples per grid cell */
  freq = 2;
  delta = MIN( x_inc, y_inc ) / ( (float) freq );

  /*
   * TODO: makes no sense at the moment
   */
  if (d_opt) {
    freq = 1;
  } else {
    freq = 1;
    dout = MIN( x_inc, y_inc ) / ( (float) 5 );
  }
  delta = dout / ( (float) freq );

  
  if (verbose) {
    fprintf(stderr,"xmin: %.3f xmax: %.3f x_inc: %.3f nx: %lu\n",
            xmin,xmax,x_inc,nx);
    fprintf(stderr,"ymin: %.3f ymax: %.3f y_inc: %.3f ny: %lu\n",
            ymin,ymax,y_inc,ny);
    fprintf(stderr,"d_out: %.3f d_inc: %.3f RK: %d freq: %u\n",
            dout,delta,k_opt,freq);
  }
  

  if ( (delta > x_inc) || (delta > y_inc) ) {
    fprintf(stderr,"WARN: Stepsize to large: %.3f > (%.3f, %.3f)\n",
            delta,x_inc,y_inc);
  }


  if (p_file_name == NULL) {   /* Just read standard input */
    fp = stdin;
    if (verbose) fprintf (stderr, "Reading from standard input\n");
  } else {
    if ( (fp = fopen (p_file_name, "r")) == NULL) {
      fprintf (stderr, "Cannot open file %s\n",p_file_name);
      return EXIT_FAILURE;
    } else {
      if(verbose) fprintf (stderr, "Working on file %s\n", p_file_name);
    }
  }


  while (!feof(fp) && (fgets (line, BUFSIZ, fp) != 0) ) {

    nlines++;
    if (verbose) {
      fprintf(stderr,"Processing input line: %u ...\n",nlines);
    }

    if( line[0] == '#' || line[0] == '>' || line[0]== '\n' ) {
      /* if(verbose) fprintf(stderr,"skip header\n"); */
    } else {

      if(2 != sscanf (line, "%f %f", &x0,&y0)){
        fprintf(stderr,
                "ERROR: Mismatch between actual and expected fields near line %u\n",
                nlines);
        return EXIT_FAILURE;
      }

      /*****************************************************************
       *
       * the work is starting here
       *
       ******************************************************************/
      if ( (x0 > p_x[nx-2]) || (y0 > p_y[ny-2]) ||
           (x0 < p_x[0]) || (y0 < p_y[0]) ) {
        fprintf(stderr,"Error: P(%.3f, %.3f) is outside valid region.\n",x0,y0);
        fprintf(stderr,"... using next\n");
        continue; /* continue with next point in file */
      }

      npoly++;
      printf("> streamline: %u\n",npoly);
      dist = 0.0;
      
      /* simple step iteration */
      xi = x0;
      yi = y0;
      i = 0;
      
      /*
       * Points at current streamline
       */
      for ( i = 0; i <= n; i++ ) {
        
        if ( i == 0 ) {
          dx0 =  dx1 = dx2 = dx3 = 0.0f;
          dy0 =  dy1 = dy2 = dy3 = 0.0f;
        }

        /*
         * STEP 0 (get inital velocity data at requested point)
         */
        (void)interp2(nx,ny,p_x,p_y,p_vx,p_vy,xi,yi,&vxi,&vyi);
        debug_printf(DEBUG_INFO,
                     "STEP0: xi = %.3f, yi = %.3f, vxi = %.3f, vyi = %.3f\n",xi,yi,vxi,vyi);

        if (verbose>1){
          fprintf(stderr,
                  "# x = %.3f, y = %.3f, vx = %.3f, vy = %.3f\n",xi,yi,vxi,vyi);
        }

        if (! (i % freq) ) {
          if(l_opt) {
            printf("%.3f %.3f %.3f %.3f %.3f\n",xi,yi,dist,vxi,vyi);
          } else {
            printf("%.3f %.3f %.3f\n",xi,yi,dist);
          }
        }

        dist += (delta * dir);

        if(isnan(vxi) || isnan(vyi)) break;
        if ( (uv = SQRT(vxi*vxi + vyi*vyi)) <= 0.0f ) break;

        dx0 = dir * delta * vxi/uv;
        dy0 = dir * delta * vyi/uv;

        
        if (verbose>2) {
          fprintf(stderr,"#\t x0=%15.3f, y0=%15.3f,",xi,yi);
          fprintf(stderr," dx0=%15.3f, dy0=%15.3f\n",dx0,dy0);
        }

        /*  check stepsize  */
        if ( i > 0 && k_opt == 4) {
          ex = dx3/6.0f - vxi*delta/6.0f;
          ey = dy3/6.0f - vyi*delta/6.0f;
          /*      fprintf(stderr,"%d %.6f\n",i,SQRT(ex*ex + ey*ey)); */
        }

        /* early break for Euler method. Continue with next i */
        if ( k_opt == 1 ) {
          xi += dx0;
          yi += dy0;
          continue;
        }

        /*
         * RK-STEP 1
         */
        xt = xi + dir*dx0/2.0f;
        yt = yi + dir*dy0/2.0f;
        (void)interp2(nx,ny,p_x,p_y,p_vx,p_vy,xt,yt,&vxi,&vyi);
        if(isnan(vxi) || isnan(vyi)) break;
        if ( (uv = SQRT(vxi*vxi + vyi*vyi)) <= 0.0f ) break;
        dx1 = dir * delta * vxi/uv;
        dy1 = dir * delta * vyi/uv;

        if (verbose > 2) {
          fprintf(stderr,"#\t x1=%15.3f, y1=%15.3f,",xt,yt);
          fprintf(stderr," dx1=%15.3f, dy1=%15.3f\n",dx1,dy1);
        }

        /*
         * RK-STEP 2
         */
        xt = xi + dx1 / 2.0f;
        yt = yi + dy1 / 2.0f;
        (void)interp2(nx,ny,p_x,p_y,p_vx,p_vy,xt,yt,&vxi,&vyi);

        if(isnan(vxi) || isnan(vyi)) break;
        if ( (uv = SQRT(vxi*vxi + vyi*vyi)) <= 0.0f ) break;
        dx2 = dir * delta * vxi/uv;
        dy2 = dir * delta * vyi/uv;

        if (verbose > 2) {
          fprintf(stderr,"#\t x2=%15.3f, y2=%15.3f,",xt,yt);
          fprintf(stderr," dx2=%15.3f, dy2=%15.3f\n",dx2,dy2);
        }


        /*
         * RK-STEP 3
         */
        xt = xi + dx2;
        yt = yi + dy2;
        (void)interp2(nx,ny,p_x,p_y,p_vx,p_vy,xt,yt,&vxi,&vyi);

        if(isnan(vxi) || isnan(vyi)) break;
        if ( (uv = SQRT(vxi*vxi + vyi*vyi)) <= 0.0f ) break;
        dx3 = dir * delta * vxi/uv;
        dy3 = dir * delta * vyi/uv;

        if(verbose>2) {
          fprintf(stderr,"#\t x3=%15.3f, y3=%15.3f,", xt,yt);
          fprintf(stderr," dx3=%15.3f, dy3=%15.3f\n",dx3,dy3);
        }

        /*
         * RK-STEP update
         */
        dx = ( dx0/6.0f + dx1/3.0f + dx2/3.0f + dx3/6.0f );
        dy = ( dy0/6.0f + dy1/3.0f + dy2/3.0f + dy3/6.0f );

        /*  fprintf(stderr,"#->   x=%15.5f,\t y=%15.5f\n", xi,yi); */

        xi += dx;
        yi += dy;

        if(verbose>1) {
          fprintf(stderr,"# x = %.3f, y = %.3f,", xi,yi);
          fprintf(stderr," dx = %.3f, dy = %.3f\n", dx,dy);
        }

        lim=SQRT(dx*dx+dy*dy);
        if (lim*1000.0f < (MIN(x_inc,y_inc))) {
          if (verbose>1) {
            fprintf(stderr,"# error: stepsize to small %.3f < (%.3f,%.3f)\n",lim,x_inc,y_inc);
          }
          break;
        }


      }
      /*
       * end of work
       */
    }


  }

  if (fp != stdin) fclose(fp);

  if(p_x) free(p_x);
  if(p_y) free(p_y);
  if(p_vx) free(p_vx);
  if(p_vy) free(p_vy);

#if 0
  log_finalize();
#endif

  return  EXIT_SUCCESS;
} /* main */



/*****************************************************************************/
int interp2(size_t nx, size_t ny, float* p_x, float* p_y,
            float* p_vx, float* p_vy,
            float xi, float yi, float* p_vxi, float* p_vyi){

  size_t ixel=0,iyel=0; /* where we are */
  size_t i00,i01,i10,i11;
  float p1=0.0f,p2=0.0f,q1=0.0f,q2=0.0f;

  locate(p_x, nx, xi, &ixel);
  debug_printf(DEBUG_INFO,"found xi = %.3f -> x[%lu] = %.3f\n",
               xi,ixel,p_x[ixel]);
  locate(p_y, ny, yi, &iyel);
  debug_printf(DEBUG_INFO,"found yi = %.3f -> y[%lu] = %.3f\n",
               yi,iyel,p_y[iyel]);

  i00 = (iyel)*(nx)+ixel;
  i01 = (iyel)*(nx)+ixel+1;
  i10 = (iyel+1)*(nx)+ixel;
  i11 = (iyel+1)*(nx)+ixel+1;

  debug_printf(DEBUG_INFO,"vx[i00] =  %.3f\n",p_vx[i00]);
  debug_printf(DEBUG_INFO,"vx[i10] =  %.3f\n",p_vx[i10]);
  debug_printf(DEBUG_INFO,"vx[i01] =  %.3f\n",p_vx[i01]);
  debug_printf(DEBUG_INFO,"vx[i11] =  %.3f\n",p_vx[i11]);

  debug_printf(DEBUG_INFO,"vy[i00] =  %.3f\n",p_vy[i00]);
  debug_printf(DEBUG_INFO,"vy[i10] =  %.3f\n",p_vy[i10]);
  debug_printf(DEBUG_INFO,"vy[i01] =  %.3f\n",p_vy[i01]);
  debug_printf(DEBUG_INFO,"vy[i11] =  %.3f\n",p_vy[i11]);

  p1 = p_vx[i00]
    + (xi - p_x[ixel]) * (p_vx[i10] - p_vx[i00])/(p_x[ixel+1]-p_x[ixel]);

  q1 = p_vx[i01]
    + (xi - p_x[ixel]) * (p_vx[i11] - p_vx[i01])/(p_x[ixel+1]-p_x[ixel]);

  (*p_vxi) = p1 + (yi-p_y[iyel])*(q1-p1)/(p_y[iyel+1]-p_y[iyel]);


  p2 = p_vy[i00]
    + (xi - p_x[ixel]) * (p_vy[i10] - p_vy[i00])/(p_x[ixel+1]-p_x[ixel]);

  q2 = p_vy[i01]
    + (xi - p_x[ixel]) * (p_vy[i11] - p_vy[i01])/(p_x[ixel+1]-p_x[ixel]);

  (*p_vyi) = p2 + (yi-p_y[iyel])*(q2-p2)/(p_y[iyel+1]-p_y[iyel]);

  debug_printf(DEBUG_INFO,"p1 =  %.3f -> q1 = %.3f\n",p1,q1);
  debug_printf(DEBUG_INFO,"p2 =  %.3f -> q2 = %.3f\n",p2,q2);


  /* fprintf(stderr,"# (xi=%.3f, yi=%.3f) => (ui=%.3f, vi=%.3f) \n", */
  /*     xi,yi,(*p_vxi),(*p_vyi) ); */




  return EXIT_SUCCESS;
} /* interp2 */


/**
 *
 */
void locate(float* xx, size_t n, float x, size_t *j)
{
#if 0
  /* fast search
   * table lookup */
  size_t ju,jm,jl;
  int ascnd;
  jl=0;
  ju=n+1;
  ascnd=(xx[n] >= xx[1]);
  while (ju-jl > 1) {
    jm=(ju+jl) >> 1;
    if (x >= xx[jm] == ascnd)
      jl=jm;
    else
      ju=jm;
  }
  if (x == xx[1]) *j=1;
  else if(x == xx[n]) *j=n-1;
  else *j=jl;
#else
  size_t i;
  i=0;
  while (xx[i] <= x) {
    i++;
  }
  *j = i-1;
#endif
} /* locate */


/**
 *
 */
void usage(void)
{
  fprintf(stderr,"NAME:\n"
          "  %s - generate stream lines based on GMT grid files\n",program_name);
  fprintf(stderr,"USAGE:\n"
          "  %s v_x.grd v_y.grd -f xyfile \n\n"
          "  v_x.grd & v_y.grd   grid files with the 2 vector components\n"
          "  xyfile              two column ASCII file containing (x0,y0)\n",program_name);
  fprintf(stderr,"\nOPTIONS:\n"
          "  -b                  backward steps\n"
          "  -d inc              stepsize\n"
          /* "  -k                  select stepping method (default: RK4)\n" */
          "  -l                  long output format: 'x y dist v_x v_y' (5 cols)\n"
          "  -n maxsteps         maximum number of steps (default: %d)\n"
          "  -V                  verbose output\n"
          "  -v                  version\n"
          "  -h                  help\n\n",MAXSTEPS);
  fprintf(stderr,"\nDESCRIPTION:\n"
          "  %s - reads (x0,y0) pairs from standard input or xyfile (-f option)\n"
          "  and generates polylines in multiple seqment mode each starting at x0,y0.\n"
          "  Output: 'x y dist' (3 cols) to stdout.\n",program_name);
  fprintf(stderr,"\nEXAMPLE:\n"
          "  echo \"0 0\" | %s vx.grd vy.grd | psxy -m -R -J ... \n"
          "  \n",program_name);
  exit(0);
}


/**
 *
 */
void version(void) {
  fprintf(stderr,"This is %s version %s (%s).\n",
          PACKAGE_NAME,PACKAGE_VERSION, __DATE__);
  exit(0);
}
