#ifndef LAST_UPDATE
#define LAST_UPDATE "Time-stamp: <2016-04-12 08:57:56 (tkleiner)>"
#endif

/*
 * grd2stream
 * reads two 2-D gridded files which represents the  x-  and  y-
 * components  of a vector field and produces a stream line polygon with
 * starting at point x0,y0 reading from stdin or file
 */

/*
 * @todo: fix points outside grids in output
 * @todo: get a better gues for step size
 * @todo: 
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.2.3"
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
static int interp2(size_t nx, size_t ny, double *p_x, double *p_y,
                   double *p_vx, double *p_vy, double xi, double yi,
                   double *p_vxi,double *p_vyi);

/**
 * returns a value i such that x is between xx[i] and xx[i+1] 
 * for non-equidistant grids
 */
static void locate(double* xx, size_t n, double x, size_t *j);


/**
 * old version (obsolete)
 */
void locate_slow(double* xx, size_t n, double x, size_t *j);


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

void log_break_nan(double x, double y, double x0, double y0)
{
  if (log_breaks) {
    fprintf(stderr,
            "Stop: NaN found for velocity vx or vy at position (%.3f, %.3f) "
            "for seed (%.3f, %.3f).\n",x,y,x0,y0);
  }
}


void log_break_zero(double x, double y, double x0, double y0)
{
  if (log_breaks) {
    fprintf(stderr,
            "Stop: Zero velocity magnitude found at position (%.3f, %.3f) "
            "for seed (%.3f, %.3f).\n",x,y,x0,y0);
  }
}


void log_break_dx(double x, double x0, double y0)
{
  if (log_breaks){
    fprintf(stderr,
            "Stop: Estimated x + dx (%.3f) is outside valid region for seed (%.3f, %.3f).\n",
            x,x0,y0);
  }
}

void log_break_dy(double y, double x0, double y0)
{
  if (log_breaks){
    fprintf(stderr,
            "Stop: Estimated y + dy (%.3f) is outside valid region for seed (%.3f, %.3f).\n",
            y,x0,y0);
  }
}

void log_break_stepsize(double stpsz, double xinc, double yinc, double x0, double y0)
{
  if (log_breaks){
    fprintf(stderr,
            "Stop: Stepsize to small %.3f << (%.3f,%.3f) for seed (%.3f, %.3f)\n",
            stpsz,xinc,yinc,x0,y0);
  }
}

void log_break_maxiter(double x0, double y0)
{
  if (log_breaks){
    fprintf(stderr,
            "Stop: Maximum number of iterations reached for seed (%.3f, %.3f)\n",
            x0,y0);
  }
}



/* for logging */
const char *program_name = PACKAGE_NAME;


/*****************************************************************************/
int main( int argc, char** argv )
{
  size_t nx=0,ny=0;
  double xmin,xmax,ymin,ymax,x_inc,y_inc;

  unsigned long i, j, cnt, iter, maxiter = MAXSTEPS;
  unsigned int nlines=0,npoly=0;

  double y0=0.0,x0=0.0;   /* start positions */
  double yi=0.0,xi=0.0;
  double yt=0.0,xt=0.0;
  double vxi=0.0,vyi=0.0,uv=0.0;

  double dx0=0.0,dx1=0.0,dx2=0.0,dx3=0.0;
  double dy0=0.0,dy1=0.0,dy2=0.0,dy3=0.0;
  double ex=0.0,ey=0.0; /* local error */
  double dx=0.0,dy=0.0; /* local error */
  double lim=0.0;

  double dist,dout,delta = 1000.0; /* unit m ???*/
  double dir=1.0;              /* direction (1.. forward, -1..backward) */
  unsigned int freq = 1;

  double *p_x=NULL;
  double *p_y=NULL;
  double *p_vx=NULL;
  double *p_vy=NULL;


  
  /* x and y coordinates of the coarse grid */
  double *p_xc = NULL;
  double *p_yc = NULL;
  int *blank = NULL;
  size_t nbx=0,nby=0;
  size_t ib=0,jb=0;
  double xb_inc = 0.0, yb_inc = 0.0;
  double density = .1; /* 10% */

  char* p_vx_name = NULL;
  char* p_vy_name = NULL;
  char* p_file_name = NULL;


  size_t nxm=0,nym=0;
  size_t im=0,jm=0;
  double xm_inc = 0.0, ym_inc = 0.0;
  char* p_mask_name = NULL;
  double *p_xm=NULL;
  double *p_ym=NULL;
  double *p_mask=NULL; /* read as double for now */
  double mask_val = 0.0;


  
  FILE* fp = NULL;
  char line[BUFSIZ];

  int oc,err=0;
  int b_opt = 0; /* backward steps */
  int d_opt = 0; /* delta */
  int k_opt = 4; /* Runge Kutta 4 */
  int l_opt = 0; /* print also u,v in columns 4,5 */

  /* experimental options */
  int L_opt = 0; /* 3col input */
  int D_opt = 0; /* check distance */
  int M_opt = 0; /* read MASK */


  /* Initialize logging facility */
  /*  log_initialize(opt_nodaemon ? LOG_TO_STDERR : LOG_TO_SYSLOG);*/

#if 0
  (void)log_initialize(LOG_TO_STDERR);
  (void)log_set_debug(DEBUG_TRACE_NONE); /* or DEBUG_TRACE_ALL,... */
#endif

  /* parse commandline args */
  while ((oc = getopt (argc, argv, "bd:lhvk:n:f:VLDrM:")) != -1)
    switch (oc) {
    case 'b': 
      /* go backward */
      b_opt=1; break;
    case 'l': 
      /* long output */
      l_opt=1; break;
    case 'L': 
      /* long input */
      L_opt=1; break;
    case 'D': 
      /* check spacing */
      D_opt=1; break;
    case 'd': 
      /* output step size */
      d_opt=1; dout= (double) atof(optarg); break;
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
      maxiter = (unsigned int)atoi(optarg); break; 
    case 'f': 
      /* read initial points from file */
      p_file_name = (optarg); break; 
    case 'M': 
      /* read mask */
      M_opt=1;
      p_mask_name = (optarg); break; 
    case 'V': 
      /* verbose opttion */
      verbose++; break;
    case 'r': 
      /* verbose opttion */
      log_breaks=1; break;
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
  if ( M_opt > 0 && p_mask_name != NULL) {
    err  += grdread (p_mask_name, &nxm, &nym, &p_xm, &p_ym, &p_mask);
    xm_inc = p_xm[1]- p_xm[0];
    ym_inc = p_ym[1]- p_ym[0];
    if (verbose) {
      fprintf(stderr,"Input mask:\n");
      fprintf(stderr,"xmin: %.3f xmax: %.3f x_inc: %f nx: %lu\n",
              p_xm[0],p_xm[nxm-1],xm_inc,nxm);
      fprintf(stderr,"ymin: %.3f ymax: %.3f y_inc: %f ny: %lu\n",
              p_ym[0],p_ym[nym-1],ym_inc,nym);
    }
  }
  if (err != 0) {
    fprintf(stderr,"error reading grd files\n");
    exit(EXIT_FAILURE);
  }
  

  /* grid limits */
  xmin = p_x[0];
  xmax = p_x[nx-1];
  ymin = p_y[0];
  ymax = p_y[ny-1];
  x_inc= (xmax - xmin)/ ( (double) (nx - 1) );
  y_inc= (ymax - ymin)/ ( (double) (ny - 1) );


  /* 
   * Blank array: This is the heart of the algorithm. It begins life
   * zeroed, but is set to one when a streamline passes through each
   * box. Then streamlines are only allowed to pass through zeroed
   * boxes. The lower resolution of this grid determines the
   * approximate spacing between trajectories.
   */  
  nbx = (int)(nx*density+1);
  nby = (int)(ny*density+1);

  blank = (int *)calloc(nby * nbx, sizeof(int));
  /* box spacing */
  xb_inc = (p_x[nx-1] - p_x[0]) / ( (double)(nbx-1) );
  yb_inc = (p_y[ny-1] - p_y[0]) / ( (double)(nby-1) );

  if(p_xc == NULL) p_xc = (double*) calloc(nbx, sizeof(double));
  if(p_yc == NULL) p_yc = (double*) calloc(nby, sizeof(double));

  for (i = 0; i<nbx; i++) {
    p_xc[i] = p_x[0] + i * xb_inc;
  }
  for (j = 0; j<nby; j++) {
    p_yc[j] = p_y[0] + j * yb_inc;
  }
  for (i = 0; i<nbx*nby; i++) {
    blank[i] = IO_NAN; 
  }
  

  /* default freq = 2 samples per grid cell */
  freq = 2;
  delta = MIN( x_inc, y_inc ) / ( (double) freq );

  /*
   * TODO: makes no sense at the moment
   */
  if (d_opt) {
    freq = 1;
  } else {
    freq = 1;
    dout = MIN( x_inc, y_inc ) / ( (double) 5 );
  }
  delta = dout / ( (double) freq );

  
  if (verbose) {
    fprintf(stderr,"Input:\n");
    fprintf(stderr,"xmin: %.3f xmax: %.3f x_inc: %f nx: %lu\n",
            xmin,xmax,x_inc,nx);
    fprintf(stderr,"ymin: %.3f ymax: %.3f y_inc: %f ny: %lu\n",
            ymin,ymax,y_inc,ny);
    fprintf(stderr,"d_out: %.3f d_inc: %.3f RK: %d freq: %u\n",
            dout,delta,k_opt,freq);
    fprintf(stderr,"verbose: %u\n",verbose);
  }
  
  if (verbose>1) {
    fprintf(stderr,"Coarse grid:\n");
    fprintf(stderr,"xmin: %.3f xmax: %.3f x_inc: %f nx: %lu\n",
            p_xc[0],p_xc[nbx-1],xb_inc,nbx);
    fprintf(stderr,"ymin: %.3f ymax: %.3f y_inc: %f ny: %lu\n",
            p_yc[0],p_yc[nby-1],yb_inc,nby);
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

      if (L_opt) {

        if(3 != sscanf (line, "%lf %lf %lf", &x0,&y0,&dir)){
          fprintf(stderr,
                  "ERROR: Mismatch between actual and expected fields near line %u\n",
                  nlines);
          return EXIT_FAILURE;
        } else {
          if ( dir<0.0 ) {
            dir = -1.0;
          } else {
            dir = 1.0;
          }
        }
      
      } else {

        if(2 != sscanf (line, "%lf %lf", &x0,&y0)){
          fprintf(stderr,
                  "ERROR: Mismatch between actual and expected fields near line %u\n",
                  nlines);
          return EXIT_FAILURE;
        } else {
          if (verbose > 1) {
            fprintf(stderr,"       x0=%.3f y0=%.3f\n",x0,y0);
          }
        }

      }

      /*****************************************************************
       *
       * the work is starting here
       *
       ******************************************************************/
      
      if ( (x0 > p_x[nx-1]) || (y0 > p_y[ny-1]) ||
           (x0 < p_x[0]) || (y0 < p_y[0]) ) {
        if (log_breaks) {
          fprintf(stderr,
                  "Stop: Seed (%.3f, %.3f) is outside valid region "
                  "(xmin=%.3f, xmax=%.3f, ymin=%.3f, ymax=%.3f)!\n",
                  x0, y0, p_x[0], p_x[nx-1], p_x[0], p_x[nx-1]);
        }
        continue; /* continue with next point in file */
      }

      npoly++;
      printf("> streamline: %u\n",npoly);
      dist = 0.0;
      
      /* simple step iteration */
      xi = x0;
      yi = y0;
      iter = 0;
      
      /*
       * Points at current streamline
       */
      for ( iter = 0; iter < maxiter; iter++ ) {
        
        if ( iter == 0 ) {
          dx0 =  dx1 = dx2 = dx3 = 0.0;
          dy0 =  dy1 = dy2 = dy3 = 0.0;
        }

        /*
         * STEP 0 (get inital velocity data at requested point)
         */
        (void)interp2(nx,ny,p_x,p_y,p_vx,p_vy,xi,yi,&vxi,&vyi);
        if (verbose>1){
          fprintf(stderr,
                  "# x = %.3f, y = %.3f, vx = %.3f, vy = %.3f\n",xi,yi,vxi,vyi);
        }
        /* check if mask reached */
        if (M_opt > 0) {
          im = (size_t)( ((xi - p_xm[0]) / xm_inc));
          jm = (size_t)( ((yi - p_ym[0]) / ym_inc));
          mask_val = p_mask[jm*nxm + im];
        }

        /* this is the output */
        if (! (iter % freq) ) {
          if(l_opt) {
            printf("%.3f %.3f %.3f %.3f %.3f\n",xi,yi,dist,vxi,vyi);
          } else {
            if (M_opt > 0) {
              printf("%.3f %.3f %.3f %.3f\n",xi,yi,dist,mask_val);
            } else {
              printf("%.3f %.3f %.3f\n",xi,yi,dist);
            }
          }
        }

        if (isnan(mask_val) || mask_val == 0.0) {
          /* continue until valid mask value is found */
        } else {
          /* break valid mask point */
          /* print start position and mask of final position */
          printf("#M# %.3f %.3f %.3f\n",x0,y0,mask_val);
          break;
        }

        
        if(isnan(vxi) || isnan(vyi)) {
          log_break_nan(xi,yi,x0,y0);
          if(M_opt)printf("#M# %.3f %.3f NaN\n",x0,y0);
          break;
        }
        if ( (uv = SQRT(vxi*vxi + vyi*vyi)) <= 0.0 ) {
          log_break_zero(xi,yi,x0,y0);
          if(M_opt)printf("#M# %.3f %.3f NaN\n",x0,y0);
          break;
        }

        /* advance now */
        dist += (delta * dir);
        dx0 = dir * delta * vxi/uv;
        dy0 = dir * delta * vyi/uv;


        /*
         * check initial guess
         */
        if ( (x0+dx0 > p_x[nx-1]) || (x0+dx0 < p_x[0]) ) {
          log_break_dx(x0+dx0, x0, y0);
          if(M_opt)printf("#M# %.3f %.3f NaN\n",x0,y0);
          break;
        }
        if ( (y0+dy0 > p_y[ny-1]) || (y0+dy0 < p_y[0]) ) {
          log_break_dy(y0+dy0, x0, y0);
          if(M_opt)printf("#M# %.3f %.3f NaN\n",x0,y0);
          break;
        }


        if (verbose>2) {
          fprintf(stderr,"#\t x0=%15.3f, y0=%15.3f,",xi,yi);
          fprintf(stderr," dx0=%15.5f, dy0=%15.5f\n",dx0,dy0);
        }

        /*  check stepsize  */
        if ( iter > 0 && k_opt == 4) {
          ex = dx3/6.0f - vxi*delta/6.0f;
          ey = dy3/6.0f - vyi*delta/6.0f;
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
        if(isnan(vxi) || isnan(vyi)) {
          log_break_nan(xt,yt,x0,y0);
          if(M_opt)printf("#M# %.3f %.3f NaN\n",x0,y0);
          break;
        }
        if ( (uv = SQRT(vxi*vxi + vyi*vyi)) <= 0.0 ) {
          log_break_zero(xt,yt,x0,y0);
          if(M_opt)printf("#M# %.3f %.3f NaN\n",x0,y0);
          break;
        }
        dx1 = dir * delta * vxi/uv;
        dy1 = dir * delta * vyi/uv;

        if (verbose > 2) {
          fprintf(stderr,"#\t x1=%15.3f, y1=%15.3f,",xt,yt);
          fprintf(stderr," dx1=%15.5f, dy1=%15.5f\n",dx1,dy1);
        }

        /*
         * RK-STEP 2
         */
        xt = xi + dx1 / 2.0f;
        yt = yi + dy1 / 2.0f;
        (void)interp2(nx,ny,p_x,p_y,p_vx,p_vy,xt,yt,&vxi,&vyi);

        if(isnan(vxi) || isnan(vyi)) {
          log_break_nan(xt,yt,x0,y0);
          if(M_opt)printf("#M# %.3f %.3f NaN\n",x0,y0);
          break;
        }
        if ( (uv = SQRT(vxi*vxi + vyi*vyi)) <= 0.0 ) {
          log_break_zero(xt,yt,x0,y0);
          if(M_opt)printf("#M# %.3f %.3f NaN\n",x0,y0);
          break;
        }
        dx2 = dir * delta * vxi/uv;
        dy2 = dir * delta * vyi/uv;

        if (verbose > 2) {
          fprintf(stderr,"#\t x2=%15.3f, y2=%15.3f,",xt,yt);
          fprintf(stderr," dx2=%15.5f, dy2=%15.5f\n",dx2,dy2);
        }


        /*
         * RK-STEP 3
         */
        xt = xi + dx2;
        yt = yi + dy2;
        (void)interp2(nx,ny,p_x,p_y,p_vx,p_vy,xt,yt,&vxi,&vyi);

        if(isnan(vxi) || isnan(vyi)) {
          log_break_nan(xt,yt,x0,y0);
          if(M_opt)printf("#M# %.3f %.3f NaN\n",x0,y0);
          break;
        }
        if ( (uv = SQRT(vxi*vxi + vyi*vyi)) <= 0.0 ) {
          log_break_zero(xt,yt,x0,y0);
          if(M_opt)printf("#M# %.3f %.3f NaN\n",x0,y0);
          break;
        }
        dx3 = dir * delta * vxi/uv;
        dy3 = dir * delta * vyi/uv;

        if(verbose>2) {
          fprintf(stderr,"#\t x3=%15.3f, y3=%15.3f,", xt,yt);
          fprintf(stderr," dx3=%15.5f, dy3=%15.5f\n",dx3,dy3);
        }

        /*
         * RK-STEP update
         */
        dx = ( dx0/6.0f + dx1/3.0f + dx2/3.0f + dx3/6.0f );
        dy = ( dy0/6.0f + dy1/3.0f + dy2/3.0f + dy3/6.0f );


        /*
         * check final step 
         */
        if ( (xi+dx > p_x[nx-1]) || (xi+dx < p_x[0]) ) {
          log_break_dx(xi+dx, x0, y0);
          if(M_opt)printf("#M# %.3f %.3f NaN\n",x0,y0);
          break;
        }
        if ( (yi+dy > p_y[ny-1]) || (yi+dy < p_y[0]) ) {
          log_break_dy(yi+dy, x0, y0);
          if(M_opt)printf("#M# %.3f %.3f NaN\n",x0,y0);
          break;
        }


        xi += dx;
        yi += dy;



        /*
         * check if already visited
         *
         * is not working now for backward streamlines! Why?
         */

        /* compute index in blank array */
        ib = (size_t)( ((xi - p_x[0]) / xb_inc));
        jb = (size_t)( ((yi - p_y[0]) / yb_inc));
        blank[jb*nbx + ib] = npoly;

        /* if (verbose) { */
          
        /*   fprintf(stderr," xi=%15.5f yi=%15.5f\n",xi,yi); */
        /*   fprintf(stderr," ib=%d jb=%d\n",ib,jb); */
        /*   fprintf(stderr," xb=%15.5f yb=%15.5f\n",p_xc[ib],p_yc[jb]); */

        /* } */

        
        
        /* if (D_opt) {         */
          
        /*   if (ib <0 || ib > nbx-1) { */
        /*     fprintf(stdout,"# error: %lu %lu %f\n",ib,nbx,(xi - p_x[0])); */
        /*     break; */
        /*   } */
        /*   if (jb <0 || jb > nby-1) { */
        /*     fprintf(stdout,"# error: %lu %lu %f\n",ib,nby,(yi - p_y[0])); */
        /*     break; */
        /*   } */
          
        /*   if ( (blank[jb*nbx + ib] > 0) && (blank[jb*nbx + ib] != npoly)) { */
            
        /*     if (verbose) { */
        /*       fprintf(stdout,"# blank: %lu %lu %d\n", */
        /*               ib,jb, blank[jb*nbx + ib]); */
        /*       fprintf(stdout,"# Stop: x = %.3f, y = %.3f allready visited by %d\n", */
        /*               xi,yi,(int)blank[jb*nbx + ib]); */
              
        /*     } */
        /*     break; */
            
        /*   } else { */
        /*     blank[jb*nbx + ib] = npoly; */
        /*   } */
        /* } */


        if(verbose>1) {
          fprintf(stderr,"# x = %.3f, y = %.3f,", xi,yi);
          fprintf(stderr," dx = %.3f, dy = %.3f\n", dx,dy);
        }

        lim=SQRT(dx*dx+dy*dy);
        if (lim*1000.0 < (MIN(x_inc,y_inc))) {
          log_break_stepsize(lim, x_inc, y_inc, x0, y0);
          break;
        }

      }
      /* end of work for one stream line */
      if (iter == maxiter) {
        log_break_maxiter(x0, y0);
        if(M_opt)printf("#M# %.3f %.3f NaN\n",x0,y0);
        /* if (M_opt) { */
        /*   fprintf(stderr,"#ERROR: maxiter reached for x0 = %.3f, y0 = %.3f without reaching the mask\n", x0,y0); */
        /* } */
      }
      
    }
  }

  /* write blank array for later usage */
#ifdef LATER_USAGE
  grdwrite ("blank.nc", nbx, nby, p_xc, p_yc, (void *)blank, NC_INT);
#endif
  

  if (fp != stdin) fclose(fp);
  if(p_x) free(p_x);
  if(p_y) free(p_y);
  if(p_xc) free(p_xc);
  if(p_yc) free(p_yc);
  if(p_vx) free(p_vx);
  if(p_vy) free(p_vy);

  if(p_mask) free(p_mask);
  if(p_xm) free(p_xm);
  if(p_ym) free(p_ym);
  
  if(blank) free(blank);

#if 0
  log_finalize();
#endif

  return  EXIT_SUCCESS;
} /* main */



/*****************************************************************************/
int interp2(size_t nx, size_t ny, double* p_x, double* p_y,
            double* p_vx, double* p_vy,
            double xi, double yi, double* p_vxi, double* p_vyi){

  size_t ixel=0,iyel=0; /* where we are */
  size_t ixel1=0,iyel1=0; /* the next */
  size_t i00,i01,i10,i11;
  double p1=0.0,p2=0.0,q1=0.0,q2=0.0;

#if TEST_LOCATE
  size_t ixel_dbg = 0;
#endif
  
  debug_printf(DEBUG_INFO,"search xi = %.3f in x[%lu] = %.3f < xi < x[%lu] = %.3f\n",
               xi,0,p_x[0],nx-1,p_x[nx-1]);
  debug_printf(DEBUG_INFO,"search yi = %.3f in y[%lu] = %.3f < yi < y[%lu] = %.3f\n",
               yi,0,p_y[0],ny-1,p_y[ny-1]);

  if (xi >= p_x[nx-1]) {
    /* right margin */
    ixel = nx-1; 
  } else if (xi <= p_x[0]) {
    /* left margin */
    ixel = 0;
  } else {
    locate(p_x, nx, xi, &ixel);
#if TEST_LOCATE
    locate_slow(p_x, nx, xi, &ixel_dbg);
    fprintf(stderr,"x[%lu] = %.3f  <=> x[%lu] = %.3f\n",
            ixel,p_x[ixel],ixel_dbg,p_x[ixel_dbg]);
#endif    
  }
  debug_printf(DEBUG_INFO,"found xi = %.3f -> x[%lu] = %.3f\n",
               xi,ixel,p_x[ixel]);

  if (yi >= p_y[ny-1]) {
    /* top margin */
    iyel = ny-1;
  } else if (yi <= p_y[0]) {
    /* bottom margin */
    iyel = 0;
  } else {
    locate(p_y, ny, yi, &iyel);
  }
  debug_printf(DEBUG_INFO,"found yi = %.3f -> y[%lu] = %.3f\n",
               yi,iyel,p_y[iyel]);

  ixel1 = MIN(ixel+1,nx-1);
  iyel1 = MIN(iyel+1,ny-1);

  i00 = (iyel)*(nx)+ixel;
  i01 = (iyel)*(nx)+ixel1;
  i10 = (iyel1)*(nx)+ixel;
  i11 = (iyel1)*(nx)+ixel1;

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
  
  return EXIT_SUCCESS;
} /* interp2 */


/**
 *
 */
void locate_slow(double* xx, size_t n, double x, size_t *j)
{
  size_t i;
  i=0;
  while (xx[i] <= x) {
    i++;
  }
  *j = i-1;
} /* locate */


/*
 *
 */
void locate(double* xx, size_t n, double x, size_t *i)
{
  /* fast search
   * table lookup */
  size_t iu,im,il;
  int ascnd;

  il=0;
  iu=n-1;
  ascnd=(xx[n-1] >= xx[0]);
  while (iu-il > 1) {
    im=(iu+il) >> 1;
    if (x >= xx[im] == ascnd)
      il=im;
    else
      iu=im;
  }

  if (x == xx[0])
    *i = 0;
  else if(x == xx[n-1])
    *i = n-1;
  else
    *i = il;

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
          "  -r                  report why a streamline stopped to stderr (default: off)\n"
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
