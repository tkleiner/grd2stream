#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>

#include "debug_printf.h"
#include "grdio.h"


/* rank (number of dimensions) for each variable */
#define RANK_x 1
#define RANK_y 1
#define RANK_z 2

static void
check_err(const int stat, const int line, const char *file) {
  if (stat != NC_NOERR) {
    (void) fprintf(stderr, "line %d of %s: %s\n", 
		   line, file, nc_strerror(stat));
    exit(1);
  }
}


int 
grdwrite (const char *filename, const size_t nx, size_t ny, 
				  const double *p_x, const double *p_y, const double *p_z)
{
  int  stat;                   /* return status */
  int  ncid;                   /* netCDF id */
  
  /* dimension ids */
  int x_dim, y_dim;
  /* variable ids */
  int x_id, y_id, z_id;
  
  /* variable shapes */
  int x_dims[RANK_x];
  int y_dims[RANK_y];
  int z_dims[RANK_z];
  
  /* 
   * enter define mode 
   */
  stat = nc_create(filename, NC_CLOBBER, &ncid);
  check_err(stat,__LINE__,__FILE__);

  /* define dimensions */
  stat = nc_def_dim(ncid, "x", nx, &x_dim);
  check_err(stat,__LINE__,__FILE__);
  stat = nc_def_dim(ncid, "y", ny, &y_dim);
  check_err(stat,__LINE__,__FILE__);
  
  /* define variables */
  x_dims[0] = x_dim;
  stat = nc_def_var(ncid, "x", NC_DOUBLE, RANK_x, x_dims, &x_id);
  check_err(stat,__LINE__,__FILE__);
  
  y_dims[0] = y_dim;
  stat = nc_def_var(ncid, "y", NC_DOUBLE, RANK_y, y_dims, &y_id);
  check_err(stat,__LINE__,__FILE__);
   
  z_dims[0] = y_dim;
  z_dims[1] = x_dim;
  stat = nc_def_var(ncid, "z", NC_DOUBLE, RANK_z, z_dims, &z_id);
  check_err(stat,__LINE__,__FILE__);
  
  /* assign attributes */
  stat = nc_put_att_text(ncid, x_id, "long_name", 1, "x");
  check_err(stat,__LINE__,__FILE__);
  stat = nc_put_att_text (ncid, x_id, "cartesian_axis", 1, "x");
  check_err(stat,__LINE__,__FILE__);
  
  stat = nc_put_att_text(ncid, y_id, "long_name", 1, "y");
  check_err(stat,__LINE__,__FILE__);
  stat = nc_put_att_text (ncid, y_id, "cartesian_axis", 1, "y");
  check_err(stat,__LINE__,__FILE__);
  
  stat = nc_put_att_text(ncid, z_id, "long_name", 1, "z");
  check_err(stat,__LINE__,__FILE__);
  stat = nc_put_att_text (ncid, z_id, "cartesian_axis", 1, "z");
  check_err(stat,__LINE__,__FILE__);
  
  stat = nc_put_att_text(ncid, NC_GLOBAL, "Conventions", 6, "COARDS");
  check_err(stat,__LINE__,__FILE__);
  
  /* leave define mode */
  stat = nc_enddef (ncid);
  check_err(stat,__LINE__,__FILE__);
  
  /* 
   * store the data 
   */
  stat = nc_put_var_double(ncid, x_id, p_x);
  check_err(stat,__LINE__,__FILE__);
  
  stat = nc_put_var_double(ncid, y_id, p_y);
  check_err(stat,__LINE__,__FILE__);
  
  stat = nc_put_var_double(ncid, z_id, p_z);
  check_err(stat,__LINE__,__FILE__);
  
  stat = nc_close(ncid);
  check_err(stat,__LINE__,__FILE__);
  return 0;
}

int 
grdread (const char *filename, size_t *p_nx, size_t *p_ny, 
		       double **pp_x, double **pp_y, double **pp_z)
{

   int stat= 0;
   int ncid=0;
   int ndims=0,nvars=0;
   int xid=0,yid=1,zid=2;
   size_t start[] = { 0, 0 };
   size_t count[] = { 0, 0 };

   /* Open the netCDF file.*/
   stat = nc_open(filename, NC_NOWRITE, &ncid);
   check_err(stat,__LINE__,__FILE__);
   
   stat = nc_inq(ncid, &ndims, &nvars, (int *) 0, (int *) 0);
   check_err(stat,__LINE__,__FILE__);
   
   if(ndims != 2 || nvars != 3) {
     return EXIT_FAILURE;
   }

   stat = nc_inq_dimlen(ncid, 0, p_nx); /* xdim */
   check_err(stat,__LINE__,__FILE__);

   stat = nc_inq_dimlen(ncid, 1, p_ny); /* ydim */
   check_err(stat,__LINE__,__FILE__);

   debug_printf(DEBUG_INFO,"%s -> nx = %lu\n",__FILE__,*p_nx);
   debug_printf(DEBUG_INFO,"%s -> ny = %lu\n",__FILE__,*p_ny);

   if(*pp_x == NULL) *pp_x = (double*) malloc((*p_nx) * sizeof(double));
   if(*pp_y == NULL) *pp_y = (double*) malloc((*p_ny) * sizeof(double));
   if(*pp_z == NULL) *pp_z = (double*) malloc((*p_ny)*(*p_nx) * sizeof(double));

#if 0
   count[0] = (size_t)(*p_nx);
   count[1] = (size_t)(*p_ny);
   if (*pp_x == NULL || *pp_x == NULL || *pp_x == NULL ){
     exit (EXIT_FAILURE);
   } else {
     stat = nc_get_vara_double(ncid, xid, &start[0], &count[0], *pp_x);
     check_err(stat,__LINE__,__FILE__);
     stat = nc_get_vara_double(ncid, yid, &start[1], &count[1], *pp_y);
     check_err(stat,__LINE__,__FILE__);
     stat = nc_get_vara_double(ncid, zid, start, count, *pp_z);
     check_err(stat,__LINE__,__FILE__);
   }

#else
   count[1] = (size_t)(*p_nx);
   count[0] = (size_t)(*p_ny);
   if (*pp_x == NULL || *pp_x == NULL || *pp_x == NULL ){
     exit (EXIT_FAILURE);
   } else {
     stat = nc_get_vara_double(ncid, xid, &start[1], &count[1], *pp_x);
     check_err(stat,__LINE__,__FILE__);
     stat = nc_get_vara_double(ncid, yid, &start[0], &count[0], *pp_y);
     check_err(stat,__LINE__,__FILE__);
     stat = nc_get_vara_double(ncid, zid, start, count, *pp_z);
     check_err(stat,__LINE__,__FILE__);
   }
#endif


   
   /* Close the netCDF file.*/
   stat = nc_close(ncid);
   check_err(stat,__LINE__,__FILE__);
   
   return EXIT_SUCCESS;
} /* read_netCDF_file */

