#ifndef NCIO_H
#define NCIO_H

#define IO_NAN (-9999)


/**
 * generic GMT grd file reader
 */
int grdread (const char *filename, size_t *p_nx, size_t *p_ny, 
		      double **pp_x, double **pp_y, double **pp_z);


/**
 * wrapper to call gmt api for netcdf/gdal/... file io
 */
int grdread_gmt (const char *filename, size_t *p_nx, size_t *p_ny, 
		      double **pp_x, double **pp_y, double **pp_z);


/**
 *
 */
#if 0
int grdwrite (const char *filename, const size_t nx, size_t ny, 
				  const double *p_x, const double *p_y, const double *p_z);
#endif

int grdwrite (const char *filename, const size_t nx, size_t ny, 
              const double *p_x, const double *p_y, const void *const p_z, const nc_type vtype);


#endif
