#ifndef NCIO_H
#define NCIO_H

/**
 *
 */
int grdread (const char *filename, size_t *p_nx, size_t *p_ny, 
		      float **pp_x, float **pp_y, float **pp_z);

/**
 *
 */
int grdwrite (const char *filename, const size_t nx, size_t ny, 
				  const float *p_x, const float *p_y, const float *p_z);

#endif
