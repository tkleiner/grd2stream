#include <math.h>
#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* for memcopy */

#if ENABLE_GMT_API
#include <gmt.h>
/* #include <gmt_grdio.h> */
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
#define GMT_memcpy(to, from, n, type) memcpy(to, from, (n) * sizeof(type))

static void check_err(const int stat, const int line, const char *file) {
  if (stat != NC_NOERR) {
    (void)fprintf(stderr, "line %d of %s: %s\n", line, file, nc_strerror(stat));
    exit(1);
  }
}

int grdwrite(const char *filename, const size_t nx, size_t ny,
             const double *p_x, const double *p_y, const void *const p_z,
             const nc_type vtype) {
  int stat; /* return status */
  int ncid; /* netCDF id */

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
  check_err(stat, __LINE__, __FILE__);

  /* define dimensions */
  stat = nc_def_dim(ncid, "x", nx, &x_dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid, "y", ny, &y_dim);
  check_err(stat, __LINE__, __FILE__);

  /* define variables */
  x_dims[0] = x_dim;
  stat = nc_def_var(ncid, "x", NC_DOUBLE, RANK_x, x_dims, &x_id);
  check_err(stat, __LINE__, __FILE__);

  y_dims[0] = y_dim;
  stat = nc_def_var(ncid, "y", NC_DOUBLE, RANK_y, y_dims, &y_id);
  check_err(stat, __LINE__, __FILE__);

  z_dims[0] = y_dim;
  z_dims[1] = x_dim;

  switch (vtype) {
  case NC_DOUBLE:
    stat = nc_def_var(ncid, "z", NC_DOUBLE, RANK_z, z_dims, &z_id);
    stat = nc_put_att_double(ncid, z_id, "missing_value", NC_DOUBLE, 1, &dfill);
    break;
  case NC_INT:
    stat = nc_def_var(ncid, "z", NC_INT, RANK_z, z_dims, &z_id);
    stat = nc_put_att_int(ncid, z_id, "missing_value", NC_INT, 1, &ifill);
    break;
  default:
    /*type_error(); */
    break;
  } /* end switch */

  check_err(stat, __LINE__, __FILE__);

  /* assign attributes */
  stat = nc_put_att_text(ncid, x_id, "long_name", 1, "x");
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, x_id, "cartesian_axis", 1, "x");
  check_err(stat, __LINE__, __FILE__);

  stat = nc_put_att_text(ncid, y_id, "long_name", 1, "y");
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_text(ncid, y_id, "cartesian_axis", 1, "y");
  check_err(stat, __LINE__, __FILE__);

  stat = nc_put_att_text(ncid, z_id, "long_name", 1, "z");
  check_err(stat, __LINE__, __FILE__);

  stat = nc_put_att_text(ncid, NC_GLOBAL, "Conventions", 6, "COARDS");
  check_err(stat, __LINE__, __FILE__);

  /* leave define mode */
  stat = nc_enddef(ncid);
  check_err(stat, __LINE__, __FILE__);

  /*
   * store the data
   */
  stat = nc_put_var_double(ncid, x_id, p_x);
  check_err(stat, __LINE__, __FILE__);

  stat = nc_put_var_double(ncid, y_id, p_y);
  check_err(stat, __LINE__, __FILE__);

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
  } /* end switch */

  /*
     stat = nc_put_var_double(ncid, z_id, p_z);
  */
  check_err(stat, __LINE__, __FILE__);

  stat = nc_close(ncid);
  check_err(stat, __LINE__, __FILE__);
  return 0;
}

int grdread(const char *filename, size_t *p_nx, size_t *p_ny, double **pp_x,
            double **pp_y, double **pp_z) {

  int stat = 0;
  int ncid = 0;
  int ndims = 0, nvars = 0;
  int xid = 0, yid = 1, zid = 2;
  size_t start[] = {0, 0};
  size_t count[] = {0, 0};

  /* Open the netCDF file.*/
  stat = nc_open(filename, NC_NOWRITE, &ncid);
  check_err(stat, __LINE__, __FILE__);

  stat = nc_inq(ncid, &ndims, &nvars, (int *)0, (int *)0);
  check_err(stat, __LINE__, __FILE__);

  if (ndims != 2 || nvars != 3) {
    return EXIT_FAILURE;
  }

  stat = nc_inq_dimlen(ncid, 0, p_nx); /* xdim */
  check_err(stat, __LINE__, __FILE__);

  stat = nc_inq_dimlen(ncid, 1, p_ny); /* ydim */
  check_err(stat, __LINE__, __FILE__);

  debug_printf(DEBUG_INFO, "%s -> nx = %lu\n", __FILE__, *p_nx);
  debug_printf(DEBUG_INFO, "%s -> ny = %lu\n", __FILE__, *p_ny);

  if (*pp_x == NULL)
    *pp_x = (double *)malloc((*p_nx) * sizeof(double));
  if (*pp_y == NULL)
    *pp_y = (double *)malloc((*p_ny) * sizeof(double));
  if (*pp_z == NULL)
    *pp_z = (double *)malloc((*p_ny) * (*p_nx) * sizeof(double));

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
  if (*pp_x == NULL || *pp_x == NULL || *pp_x == NULL) {
    exit(EXIT_FAILURE);
  } else {
    stat = nc_get_vara_double(ncid, xid, &start[1], &count[1], *pp_x);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_get_vara_double(ncid, yid, &start[0], &count[0], *pp_y);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_get_vara_double(ncid, zid, start, count, *pp_z);
    check_err(stat, __LINE__, __FILE__);
  }
#endif

  /* Close the netCDF file.*/
  stat = nc_close(ncid);
  check_err(stat, __LINE__, __FILE__);

  return EXIT_SUCCESS;
} /* read_netCDF_file */

// Function to reverse elements of an array
static void reverse(double arr[], int n) {
  for (int low = 0, high = n - 1; low < high; low++, high--) {
    double temp = arr[low];
    arr[low] = arr[high];
    arr[high] = temp;
  }
}

/* GMT6 */
//struct GMT_GRID_HEADER {
//    uint32_t n_columns;                     /* Number of columns */
//    uint32_t n_rows;                        /* Number of rows */
//    uint32_t registration;                  /* GMT_GRID_NODE_REG (0) for node grids,
//                                               GMT_GRID_PIXEL_REG (1) for pixel grids */
//    double wesn[4];                         /* Min/max x and y coordinates */
//    double z_min;                           /* Minimum z value */
//    double z_max;                           /* Maximum z value */
//    double inc[2];                          /* The x and y increments */
//    double z_scale_factor;                  /* Grid values must be multiplied by this factor */
//    double z_add_offset;                    /* After scaling, add this */
//    char   x_units[GMT_GRID_UNIT_LEN80];    /* Units in x-direction */
//    char   y_units[GMT_GRID_UNIT_LEN80];    /* Units in y-direction */
//    char   z_units[GMT_GRID_UNIT_LEN80];    /* Grid value units */
//    char   title[GMT_GRID_TITLE_LEN80];     /* Name of data set */
//    char   command[GMT_GRID_COMMAND_LEN320];/* Name of generating command */
//    char   remark[GMT_GRID_REMARK_LEN160];  /* Comments regarding this data set */
//};


int grdread_gmt(const char *filename, size_t *p_nx, size_t *p_ny, double **pp_x,
                double **pp_y, double **pp_z) {

#if ENABLE_GMT_API

  int stat = 0;
  unsigned int row, col, k;
  struct GMTAPI_CTRL *API = NULL; /* GMT5/GMT6 API control structure */
  struct GMT_DATASET *D = NULL;   /* Structure to hold input dataset */
  struct GMT_GRID *G = NULL;      /* vx, vy */
  /* from grdinfo.c */
  uint64_t ij, node, n_nan = 0, n = 0;
  double x_min = 0.0, y_min = 0.0, z_min = 0.0;
  double x_max = 0.0, y_max = 0.0, z_max = 0.0, wesn[4];

  /* 1. Initializing new GMT session */
  if ((API = GMT_Create_Session("grd2stream session", 2U, 0U, NULL)) == NULL) {
    return EXIT_FAILURE;
  }
  debug_printf(DEBUG_INFO, "%s: GMT_Create_Session for <%s> done\n", __FILE__,
               filename);

  /* 2. READ HEADER ONLY */
  /* from grdclip (GMT5: GMT_GRID_HEADER_ONLY GMT6: GMT_CONTAINER_ONLY  ) */
#if ENABLE_GMT_API == 5
  if ((G = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE,
                         GMT_GRID_HEADER_ONLY, NULL, filename, NULL)) == NULL) {
    return EXIT_FAILURE;
  }
#elif ENABLE_GMT_API == 6
  if ((G = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE,
                         GMT_CONTAINER_ONLY, NULL, filename, NULL)) == NULL) {
    return EXIT_FAILURE;
  }
#else

#endif

  /* name in GMT5, title in GMT6 */
#if ENABLE_GMT_API == 5
  debug_printf(DEBUG_INFO, "%s: GMT_Read_HEADER %s\n", __FILE__,
               G->header->name);
#elif ENABLE_GMT_API == 6
  debug_printf(DEBUG_INFO, "%s: GMT_Read_HEADER %s\n", __FILE__,
               G->header->title);
#endif

  uint32_t ny = 0; /* Number of columns */
  uint32_t nx = 0; /* Number of rows */

  /* copy region */
  GMT_memcpy(wesn, G->header->wesn, 4, double);

#if ENABLE_GMT_API == 5
  nx = G->header->nx;
  ny = G->header->ny;
#elif ENABLE_GMT_API == 6
  nx = G->header->n_columns;
  ny = G->header->n_rows;
#endif
  debug_printf(DEBUG_INFO, "%s:  zmin: %g\n", __FILE__, G->header->z_min);
  debug_printf(DEBUG_INFO, "%s:  zmax: %g\n", __FILE__, G->header->z_max);
  debug_printf(DEBUG_INFO, "%s: units: %s\n", __FILE__, G->header->z_units);
  debug_printf(DEBUG_INFO, "%s:    nx: %u\n", __FILE__, nx);
  debug_printf(DEBUG_INFO, "%s:    ny: %u\n", __FILE__, ny);
  /* mx > nx, e.g 605 > 601 */
  debug_printf(DEBUG_INFO, "%s:    mx: %u\n", __FILE__, G->header->mx);
  /* my > ny, e.g 605 > 601 */
  debug_printf(DEBUG_INFO, "%s:    my: %u\n", __FILE__, G->header->my);
  /* Number of items (not bytes) required to hold the grid (= mx * my) */
  debug_printf(DEBUG_INFO, "%s:  size: %u\n", __FILE__, G->header->size);
  debug_printf(DEBUG_INFO, "%s:   reg:-R%f/%f/%f/%f\n", __FILE__, wesn[0],
               wesn[1], wesn[2], wesn[3]);

  /*
   * Read in the data
   */
#if ENABLE_GMT_API == 5
  if (GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE,
                    GMT_GRID_DATA_ONLY, wesn, filename,
                    G) == NULL) { /* Get data */
    return (EXIT_FAILURE);
  }
#elif ENABLE_GMT_API == 6
  if (GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE,
                    GMT_DATA_ONLY, wesn, filename, G) == NULL) {
    return (EXIT_FAILURE);
  }
#endif

  /* size of data, size of x,y, */
  *p_nx = nx;
  *p_ny = ny;

  if (*pp_x == NULL)
    *pp_x = (double *)malloc(nx * sizeof(double));
  if (*pp_y == NULL)
    *pp_y = (double *)malloc(ny * sizeof(double));
  if (*pp_z == NULL)
    *pp_z = (double *)malloc(nx * ny * sizeof(double));

  /* x and y are of type double in GMT6 */
  GMT_memcpy((*pp_x), G->x, nx, double);
  GMT_memcpy((*pp_y), G->y, ny, double);
  reverse((*pp_y), ny); /* needed because y-axis is inverted */

  /* copy 1D GMT array data into 2D  data */
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      ij = (ny - j - 1) * nx + i;
      node = GMT_Get_Index(API, G->header, j, i); // row,col
      (*pp_z)[ij] = (double)G->data[node]; /* data is float or double  in gmt */
    }
  }

  /* Destroy the GMT session */
  GMT_Destroy_Session(API);

  return EXIT_SUCCESS;
#else
  return EXIT_FAILURE;
#endif

} /* read_netCDF_file */
