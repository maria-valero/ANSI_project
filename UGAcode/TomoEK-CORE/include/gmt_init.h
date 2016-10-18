/*
 * gmt_init.h
 *
 *  Created on: Feb 16, 2015
 *      Author: nishita
 */

#ifndef GMT_INIT_H_
#define GMT_INIT_H_

//enum GMT_enum_coord {GMT_GEOGRAPHIC = 0,	/* Means coordinates are lon,lat : compute spherical distances */
	//GMT_CARTESIAN,	/* Means coordinates are Cartesian x,y : compute Cartesian distances */
//	GMT_GEO2CART,	/* Means coordinates are lon,lat but must be mapped to (x,y) : compute Cartesian distances */
	//GMT_CART2GEO};	

//enum GMT_enum_sph {GMT_DIST_M = 10,	/* 2-D lon, lat data, convert distance to meter */
	//GMT_DIST_DEG = 20,	/* 2-D lon, lat data, convert distance to spherical degree */
	//GMT_DIST_COS = 30};	/* 2-D lon, lat data, convert distance to cos of spherical degree */

//enum GMT_enum_mdist {GMT_FLATEARTH = 1,	/* Compute Flat Earth distances */
//	GMT_GREATCIRCLE,	/* Compute great circle distances */
	//GMT_GEODESIC,		/* Compute geodesic distances */
	//GMT_LOXODROME};	


int gmt_parse_R_option (struct GMT_CTRL *GMT, char *item);
void GMT_set_grddim (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *h);
void GMT_set_pad (struct GMT_CTRL *GMT, unsigned int pad);
int GMTAPI_init_grid (struct GMT_CTRL *GMT, struct GMT_OPTION *opt, double *range, double *inc, int registration, unsigned int mode, struct GMT_GRID *G);
int GMT_BC_init (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *h);
struct GMT_OPTION * GMT_Create_Options (void *V_API, int n_args_in, void *in);
enum GMT_enum_units GMT_get_unit_number (struct GMT_CTRL *GMT, char unit);
struct GMT_CTRL * GMT_begin_module (struct GMTAPI_CTRL *API, const char *lib_name, const char *mod_name, struct GMT_CTRL **Ccopy);
int GMT_default_error (struct GMT_CTRL *GMT, char option);


#endif /* GMT_INIT_H_ */
