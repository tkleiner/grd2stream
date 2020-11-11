#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* for memcopy */
#include <math.h>
#include <netcdf.h>

#if ENABLE_GMT_API
#include <gmt.h>
#endif


#include "debug_printf.h"
#include "grdio.h"


/* rank (number of dimensions) for each variable */
#define RANK_x 1
#define RANK_y 1
#define RANK_z 2

/**
 * from gmt_macros.h 
 */
#define GMT_memcpy(to,from,n,type) memcpy(to, from, (n)*sizeof(type))


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
				  const double *p_x, const double *p_y, const void *const p_z, const nc_type vtype)
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

	double dfill = (double)IO_NAN;
	int ifill = IO_NAN;
  
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

	switch (vtype) {
	case NC_DOUBLE:
    stat = nc_def_var(ncid, "z", NC_DOUBLE, RANK_z, z_dims, &z_id);
		stat = nc_put_att_double(ncid, z_id, "missing_value", NC_DOUBLE, 1,	&dfill);
		break;
	case NC_INT:
    stat = nc_def_var(ncid, "z", NC_INT, RANK_z, z_dims, &z_id);
		stat = nc_put_att_int(ncid, z_id, "missing_value", NC_INT, 1,	&ifill);
		break;
	default:
		/*type_error(); */
		break;
	}			/* end switch */


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


	switch (vtype) {
	case NC_DOUBLE:
		stat = nc_put_var_double(ncid, z_id, (const double *)p_z);
		break;
	case NC_INT:
		stat = nc_put_var_int(ncid, z_id, (const int *)p_z);
		break;
	default:
		/*type_error(); */
		break;
	}			/* end switch */

  
  /* 
     stat = nc_put_var_double(ncid, z_id, p_z);
  */
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


// Function to reverse elements of an array
static void reverse(double arr[], int n)
{
	for (int low = 0, high = n - 1; low < high; low++, high--)
	{
		double temp = arr[low];
		arr[low] = arr[high];
		arr[high] = temp;
	}
}



int 
grdread_gmt (const char *filename, size_t *p_nx, size_t *p_ny, 
		       double **pp_x, double **pp_y, double **pp_z)
{

  (void) fprintf(stderr, "Reading <%s> done\n", filename);


#if ENABLE_GMT_API
 
  
   int stat= 0;
   int verbose = 1;
   
   unsigned int row, col, k;
    
   //void *API;                        /* The API control structure */
   struct GMTAPI_CTRL *API = NULL;      /* GMT5 API control structure */

   //struct GMTAPI_CTRL *API
     
   struct GMT_DATASET *D = NULL;     /* Structure to hold input dataset */

   
   /* vx, vy */
   struct GMT_GRID *G = NULL;

   /* from grdinfo.c */
   uint64_t ij, n_nan = 0, n = 0;
   double x_min = 0.0, y_min = 0.0, z_min = 0.0, x_max = 0.0, y_max = 0.0, z_max = 0.0, wesn[4];

  
   /* 1. Initializing new GMT session */
   if ((API = GMT_Create_Session ("grd2stream session", 2U, 0U, NULL)) == NULL) {
     return EXIT_FAILURE;
   }

   (void) fprintf(stderr, "GMT_Create_Session for <%s> done\n", filename);

   /* 2. READ HEADER ONLY */
   /* from grdclip (GMT5: GMT_GRID_HEADER_ONLY GMT6: GMT_CONTAINER_ONLY  ) */
#if ENABLE_GMT_API == 5
   if ((G = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, filename, NULL)) == NULL) {
     return EXIT_FAILURE;
   }
#elif ENABLE_GMT_API == 6
   if ((G = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_CONTAINER_ONLY, NULL, filename, NULL)) == NULL) {
     return EXIT_FAILURE;
   }
#else

#endif

   
   /* name in GMT5, title in GMT6 */ 
#if ENABLE_GMT_API == 5
   if (verbose) {
     fprintf(stderr,"GMT_Read_HEADER %s\n", G->header->name);
   }
#elif ENABLE_GMT_API == 6
   if (verbose) {
     fprintf(stderr,"GMT_Read_HEADER %s\n", G->header->title);
   }
#endif


#if 0
    unsigned int mx, my;             /* Actual dimensions of the grid in memory, allowing for the padding */
    size_t       nm;                 /* Number of data items in this grid (n_columns * n_rows) [padding is excluded] */
    size_t       size;               /* Actual number of items (not bytes) required to hold this grid (= mx * my) */
    size_t       n_alloc;            /* Bytes allocated for this grid */
    uint32_t n_columns;                   /* Number of columns */
    uint32_t n_rows;                      /* Number of rows */
#endif
   
   if (verbose) {
     fprintf (stderr, "\t zmin: %g\n", G->header->z_min);
     fprintf (stderr, "\t zmax: %g\n", G->header->z_max);
     fprintf (stderr, "\tunits: %s\n", G->header->z_units);
     fprintf (stderr, "\t   nx: %u\n", G->header->nx);
     fprintf (stderr, "\t   ny: %u\n", G->header->ny);
     fprintf (stderr, "\t   mx: %u\n", G->header->mx);
     fprintf (stderr, "\t   my: %u\n", G->header->my);
     fprintf (stderr, "\t rows: %u\n", G->header->n_rows);
     fprintf (stderr, "\t cols: %u\n", G->header->n_columns);
     fprintf (stderr, "\t size: %u\n", G->header->size);
   }


   
   /* copy region */
   GMT_memcpy (wesn, G->header->wesn, 4, double);

   if (verbose) {
     fprintf (stderr, "\t  reg: -R%f/%f/%f/%f\n", wesn[0], wesn[1], wesn[2], wesn[3]);
   }
   


   /*
    * Read in the data 
    */
#if ENABLE_GMT_API == 5
   if (GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, wesn, filename, G) == NULL) {	/* Get data */
     return (EXIT_FAILURE);
   }
#elif ENABLE_GMT_API == 6
   if (GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_DATA_ONLY, wesn, filename, G) == NULL) {
     return (EXIT_FAILURE);
   }
#endif
   
   /* size of data, size of x,y, */
   *p_nx = G->header->nx;
   *p_ny = G->header->ny;
   
   if(*pp_x == NULL) *pp_x = (double*) malloc((*p_nx) * sizeof(double));
   if(*pp_y == NULL) *pp_y = (double*) malloc((*p_ny) * sizeof(double));
   if(*pp_z == NULL) *pp_z = (double*) malloc((*p_ny)*(*p_nx) * sizeof(double));
     

   /* struct GMT_GRID {                     /\* To hold a GMT float grid and its header in one container *\/ */
   /*  struct GMT_GRID_HEADER *header;      /\* Pointer to full GMT header for the grid *\/ */
   /*  float                  *data;        /\* Pointer to the float grid *\/ */
   /*  double                 *x, *y;       /\* Vector of coordinates *\/ */
   /*  void *hidden;                        /\* ---- Variables "hidden" from the API ---- *\/ */
   /* }; */

   /* x and y are of type double in GMT6 */
   GMT_memcpy ((*pp_x), G->x, G->header->nx, double);
   GMT_memcpy ((*pp_y), G->y, G->header->ny, double); /* wrong orientation, foes max to min, why */
   reverse((*pp_y), G->header->ny);  /* why needed ? */

   
   /* copy 1D GMT array data into 2D  data */

   /* todo: this grid is in the wrong oruientation */
   for (ij=0; ij<(*p_nx)*(*p_ny); ij++){
     (*pp_z)[ij] = (double)G->data[ij];
   }

   
   /* for (uint32_t i = 0; i < G->header->nx; ++i) { */
   /*   for (uint32_t j = 0; j < G->header->ny; ++j) { */
   /*     ij = j*G->header->nx + i; */
   /*     pp_z[i][j] = (double)G->data[ij]; */
   /*   } */
   /* } */


   /* int GMT_Put_Row (void *API, int row_no, struct GMT_GRID *G, float *row); */
   /* for (uint32_t i = 0; i < G->header->n_rows; i++) { */

   /* } */



   
   
   /* Destroy the GMT session */
   GMT_Destroy_Session (API);
   
   return EXIT_SUCCESS;
#else
   return EXIT_FAILURE;
#endif   

} /* read_netCDF_file */


