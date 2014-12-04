#ifndef NCIO_H
#define NCIO_H

/**
 *
 */
int grdread (const char *filename, size_t *p_nx, size_t *p_ny, 
		      double **pp_x, double **pp_y, double **pp_z);

/**
 *
 */
int grdwrite (const char *filename, const size_t nx, size_t ny, 
				  const double *p_x, const double *p_y, const double *p_z);

#endif
