#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <assert.h>
#include "surface.h"
#include "memory.h"
#include "gmt_macros.h"
#include "gmt_resources.h"
#include "string_op.h"
#include "gmt_grd.h"
#include "gmt_unit.h"
#include "common_math.h"
#include "gmt_type.h"
#include "gmt_parce.h"
#include "grd_io.h"
#include "gmt_map.h"
#include "gmt_grdio.h"

#include "gmt_io.h"
#include "gmt_api.h"
#include <time.h>

#define PATH_MAX 1024
#define GMT_NO_PROJ		-1	/* Projection not specified (initial value) */
#define GMT_CONV_LIMIT	1.0e-8		/* Fairly tight convergence limit or "close to zero" limit */
#define GMT_SMALL	1.0e-4		/* Needed when results aren't exactly zero but close */
#define GMTAPI_MAX_ID 999999	/* Largest integer that will fit in the %06d format */

#define METERS_IN_A_FOOT		0.3048			/* 2.54 * 12 / 100 */
#define METERS_IN_A_SURVEY_FOOT		(1200.0/3937.0)		/* ~0.3048006096 m */
#define METERS_IN_A_KM			1000.0
#define METERS_IN_A_MILE		1609.433	/* meters in statute mile */
#define METERS_IN_A_NAUTICAL_MILE	1852.0


static inline struct GMTAPI_CTRL * gmt_get_api_ptr (struct GMTAPI_CTRL *ptr) {return (ptr);} //nishita


/* Conic projections tagged 200-299 */
#define GMT_IS_CONICAL(C) (C->current.proj.projection / 100 == 2)

//extern static inline uint64_t gmt_n_cols_needed_for_gaps (struct GMT_CTRL *GMT, uint64_t n);
//extern static inline void gmt_update_prev_rec (struct GMT_CTRL *GMT, uint64_t n_use);

extern int GMT_nc_read_grd_info (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *header);
extern int GMT_nc_update_grd_info (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *header);
extern int GMT_nc_write_grd_info (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *header);
extern int GMT_nc_write_grd (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *header, float *grid, double wesn[], unsigned int *pad, unsigned int complex_mode);
extern int GMT_nc_read_grd (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *header, float *grid, double wesn[], unsigned int *pad, unsigned int complex_mode);

void GMT_set_pad (struct GMT_CTRL *GMT, unsigned int pad)
{	/* Changes the 4 GMT default pad values to given isotropic pad */
	GMT->current.io.pad[XLO] = GMT->current.io.pad[XHI] = GMT->current.io.pad[YLO] = GMT->current.io.pad[YHI] = pad;
}

void GMT_set_geographic (struct GMT_CTRL *GMT, unsigned int dir)
{
	/* Eliminate lots of repeated statements to do this: */
	GMT->current.io.col_type[dir][GMT_X] = GMT_IS_LON;
	GMT->current.io.col_type[dir][GMT_Y] = GMT_IS_LAT;
}

int gmt_get_grdtype (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *h)
{	/* Determine if grid is Cartesian or geographic, and if so if longitude range is <360, ==360, or >360 */
	if (GMT_x_is_lon (GMT, GMT_IN)) {	/* Data set is geographic with x = longitudes */
		if (fabs (h->wesn[XHI] - h->wesn[XLO] - 360.0) < GMT_SMALL) {
			printf ( "Geographic grid, longitudes span exactly 360\n");
			/* If w/e is 360 and gridline reg then we have a repeat entry for 360.  For pixel there are never repeat pixels */
			return ((h->registration == GMT_GRID_NODE_REG) ? GMT_GRID_GEOGRAPHIC_EXACT360_REPEAT : GMT_GRID_GEOGRAPHIC_EXACT360_NOREPEAT);
		}
		else if (fabs (h->nx * h->inc[GMT_X] - 360.0) < GMT_SMALL) {
			printf ( "Geographic grid, longitude cells span exactly 360\n");
			/* If n*xinc = 360 and previous test failed then we do not have a repeat node */
			return (GMT_GRID_GEOGRAPHIC_EXACT360_NOREPEAT);
		}
		else if ((h->wesn[XHI] - h->wesn[XLO]) > 360.0) {
			printf ( "Geographic grid, longitudes span more than 360\n");
			return (GMT_GRID_GEOGRAPHIC_MORE360);
		}
		else {
			printf ("Geographic grid, longitudes span less than 360\n");
			return (GMT_GRID_GEOGRAPHIC_LESS360);
		}
	}
	else if (h->wesn[YLO] >= -90.0 && h->wesn[YHI] <= 90.0) {	/* Here we simply advice the user if grid looks like geographic but is not set as such */
		if (fabs (h->wesn[XHI] - h->wesn[XLO] - 360.0) < GMT_SMALL) {
			printf ("Cartesian grid, yet x spans exactly 360 and -90 <= y <= 90.\n");
			printf ("     To make sure the grid is recognized as geographical and global, use the -fg option\n");
			return (GMT_GRID_CARTESIAN);
		}
		else if (fabs (h->nx * h->inc[GMT_X] - 360.0) < GMT_SMALL) {
			printf ( "Cartesian grid, yet x cells span exactly 360 and -90 <= y <= 90.\n");
			printf ( "     To make sure the grid is recognized as geographical and global, use the -fg option\n");
			return (GMT_GRID_CARTESIAN);
		}
	}
	printf ( "Grid is Cartesian\n");
	return (GMT_GRID_CARTESIAN);
}


void GMT_grd_init (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *header, struct GMT_OPTION *options, bool update)
{	/* GMT_grd_init initializes a grd header to default values and copies the
	 * options to the header variable command.
	 * update = true if we only want to update command line */

	int i;

	if (update)	/* Only clean the command history */
		GMT_memset (header->command, GMT_GRID_COMMAND_LEN320, char);
	else {		/* Wipe the slate clean */
		GMT_memset (header, 1, struct GMT_GRID_HEADER);

		/* Set the variables that are not initialized to 0/false/NULL */
		header->z_scale_factor = 1.0;
		header->row_order      = k_nc_start_south; /* S->N */
		header->z_id           = -1;
		header->n_bands        = 1; /* Grids have at least one band but images may have 3 (RGB) or 4 (RGBA) */
		header->z_min          = GMT->session.d_NaN;
		header->z_max          = GMT->session.d_NaN;
		header->nan_value      = GMT->session.f_NaN;
		if (GMT_is_geographic (GMT, GMT_OUT)) {
			strcpy (header->x_units, "longitude [degrees_east]");
			strcpy (header->y_units, "latitude [degrees_north]");
		}
		else {
			strcpy (header->x_units, "x");
			strcpy (header->y_units, "y");
		}
		strcpy (header->z_units, "z");
		GMT_grd_setpad (GMT, header, GMT->current.io.pad); /* Assign default pad */
	}

	/* Always update command line history, if given */

	if (options) {
		size_t len;
		struct GMTAPI_CTRL *API = GMT->parent;
		int argc = 0; char **argv = NULL;

		if ((argv = GMT_Create_Args (API, &argc, options)) == NULL) {
			printf ("Error: Could not create argc, argv from linked structure options!\n");
			return;
		}
		strncpy (header->command, GMT->init.module_name, GMT_GRID_COMMAND_LEN320);
		len = strlen (header->command);
		for (i = 0; len < GMT_GRID_COMMAND_LEN320 && i < argc; i++) {
			len += strlen (argv[i]) + 1;
			if (len >= GMT_GRID_COMMAND_LEN320) continue;
			strcat (header->command, " ");
			strcat (header->command, argv[i]);
		}
		header->command[len] = 0;
		sprintf (header->title, "Produced by %s", GMT->init.module_name);
		GMT_Destroy_Args (API, argc, &argv);
	}
}


void GMT_RI_prepare (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *h)
{
	/* This routine adjusts the grid header. It computes the correct nx, ny, x_inc and y_inc,
	   based on user input of the -I option and the current settings of x_min, x_max, y_min, y_max and registration.
	   On output the grid boundaries are always gridline or pixel oriented, depending on registration.
	   The routine is not run when nx and ny are already set.
	*/
	unsigned int one_or_zero;
	double s;

	one_or_zero = !h->registration;
	h->xy_off = 0.5 * h->registration;	/* Use to calculate mean location of block */

	/* XINC AND XMIN/XMAX CHECK FIRST */

	/* Adjust x_inc */

	if (GMT->current.io.inc_code[GMT_X] & GMT_INC_IS_NNODES) {	/* Got nx */
		h->inc[GMT_X] = GMT_get_inc (GMT, h->wesn[XLO], h->wesn[XHI], lrint(h->inc[GMT_X]), h->registration);
	}
	else if (GMT->current.io.inc_code[GMT_X] & GMT_INC_UNITS) {	/* Got funny units */
		switch (GMT->current.io.inc_code[GMT_X] & GMT_INC_UNITS) {
			case GMT_INC_IS_FEET:	/* foot */
				s = METERS_IN_A_FOOT;
				break;
			case GMT_INC_IS_KM:	/* km */
				s = METERS_IN_A_KM;
				break;
			case GMT_INC_IS_MILES:	/* Statute mile */
				s = METERS_IN_A_MILE;
				break;
			case GMT_INC_IS_NMILES:	/* Nautical mile */
				s = METERS_IN_A_NAUTICAL_MILE;
				break;
			case GMT_INC_IS_SURVEY_FEET:	/* US survey foot */
				s = METERS_IN_A_SURVEY_FOOT;
				break;
			case GMT_INC_IS_M:	/* Meter */
			default:
				s = 1.0;
				break;
		}
		h->inc[GMT_X] *= s / (GMT->current.proj.DIST_M_PR_DEG * cosd (0.5 * (h->wesn[YLO] + h->wesn[YHI])));	/* Latitude scaling of E-W distances */
		printf ("Distance to degree conversion implies x_inc = %g\n", h->inc[GMT_X]);
	}
	if (!(GMT->current.io.inc_code[GMT_X] & (GMT_INC_IS_NNODES | GMT_INC_IS_EXACT))) {	/* Adjust x_inc to exactly fit west/east */
		s = h->wesn[XHI] - h->wesn[XLO];
		h->nx = urint (s / h->inc[GMT_X]);
		s /= h->nx;
		h->nx += one_or_zero;
		if (fabs (s - h->inc[GMT_X]) > 0.0) {
			h->inc[GMT_X] = s;
			printf ( "Given domain implies x_inc = %g\n", h->inc[GMT_X]);
		}
	}

	/* Determine nx */

	h->nx = GMT_grd_get_nx (GMT, h);

	if (GMT->current.io.inc_code[GMT_X] & GMT_INC_IS_EXACT) {	/* Want to keep x_inc exactly as given; adjust x_max accordingly */
		s = (h->wesn[XHI] - h->wesn[XLO]) - h->inc[GMT_X] * (h->nx - one_or_zero);
		if (fabs (s) > 0.0) {
			h->wesn[XHI] -= s;
			printf ("x_max adjusted to %g\n", h->wesn[XHI]);
		}
	}

	/* YINC AND YMIN/YMAX CHECK SECOND */

	/* Adjust y_inc */

	if (GMT->current.io.inc_code[GMT_Y] & GMT_INC_IS_NNODES) {	/* Got ny */
		h->inc[GMT_Y] = GMT_get_inc (GMT, h->wesn[YLO], h->wesn[YHI], lrint(h->inc[GMT_Y]), h->registration);
		printf ("Given ny implies y_inc = %g\n", h->inc[GMT_Y]);
	}
	else if (GMT->current.io.inc_code[GMT_Y] & GMT_INC_UNITS) {	/* Got funny units */
		switch (GMT->current.io.inc_code[GMT_Y] & GMT_INC_UNITS) {
			case GMT_INC_IS_FEET:	/* feet */
				s = METERS_IN_A_FOOT;
				break;
			case GMT_INC_IS_KM:	/* km */
				s = METERS_IN_A_KM;
				break;
			case GMT_INC_IS_MILES:	/* miles */
				s = METERS_IN_A_MILE;
				break;
			case GMT_INC_IS_NMILES:	/* nmiles */
				s = METERS_IN_A_NAUTICAL_MILE;
				break;
			case GMT_INC_IS_SURVEY_FEET:	/* US survey feet */
				s = METERS_IN_A_SURVEY_FOOT;
				break;
			case GMT_INC_IS_M:	/* Meter */
			default:
				s = 1.0;
				break;
		}
		h->inc[GMT_Y] = (h->inc[GMT_Y] == 0.0) ? h->inc[GMT_X] : h->inc[GMT_Y] * s / GMT->current.proj.DIST_M_PR_DEG;
		printf ( "Distance to degree conversion implies y_inc = %g\n", h->inc[GMT_Y]);
	}
	if (!(GMT->current.io.inc_code[GMT_Y] & (GMT_INC_IS_NNODES | GMT_INC_IS_EXACT))) {	/* Adjust y_inc to exactly fit south/north */
		s = h->wesn[YHI] - h->wesn[YLO];
		h->ny = urint (s / h->inc[GMT_Y]);
		s /= h->ny;
		h->ny += one_or_zero;
		if (fabs (s - h->inc[GMT_Y]) > 0.0) {
			h->inc[GMT_Y] = s;
			printf ("Given domain implies y_inc = %g\n", h->inc[GMT_Y]);
		}
	}

	/* Determine ny */

	h->ny = GMT_grd_get_ny (GMT, h);

	if (GMT->current.io.inc_code[GMT_Y] & GMT_INC_IS_EXACT) {	/* Want to keep y_inc exactly as given; adjust y_max accordingly */
		s = (h->wesn[YHI] - h->wesn[YLO]) - h->inc[GMT_Y] * (h->ny - one_or_zero);
		if (fabs (s) > 0.0) {
			h->wesn[YHI] -= s;
			printf ("y_max adjusted to %g\n", h->wesn[YHI]);
		}
	}

	h->r_inc[GMT_X] = 1.0 / h->inc[GMT_X];
	h->r_inc[GMT_Y] = 1.0 / h->inc[GMT_Y];
}

void GMT_set_grdinc (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *h)
{
	/* Update grid increments based on w/e/s/n, nx/ny, and registration */
	h->inc[GMT_X] = GMT_get_inc (GMT, h->wesn[XLO], h->wesn[XHI], h->nx, h->registration);
	h->inc[GMT_Y] = GMT_get_inc (GMT, h->wesn[YLO], h->wesn[YHI], h->ny, h->registration);
	h->r_inc[GMT_X] = 1.0 / h->inc[GMT_X];	/* Get inverse increments to avoid divisions later */
	h->r_inc[GMT_Y] = 1.0 / h->inc[GMT_Y];
}

void GMT_set_grddim (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *h)
{	/* Assumes pad is set and then computes nx, ny, mx, my, nm, size, xy_off based on w/e/s/n.  */
	h->nx = GMT_grd_get_nx (GMT, h);		/* Set nx, ny based on w/e/s/n and offset */
	h->ny = GMT_grd_get_ny (GMT, h);
	
	h->mx = gmt_grd_get_nxpad (h, h->pad);	/* Set mx, my based on h->{nx,ny} and the current pad */
	h->my = gmt_grd_get_nypad (h, h->pad);
	h->nm = gmt_grd_get_nm (h);		/* Sets the number of actual data items */
	h->size = gmt_grd_get_size (h);		/* Sets the number of items (not bytes!) needed to hold this array, which includes the padding (size >= nm) */
	h->xy_off = 0.5 * h->registration;
	GMT_set_grdinc (GMT, h);
}

int GMT_BC_init (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *h)
{	/* Initialize grid boundary conditions based on grid header and -n settings */
	int i = 0, type;
	bool same;
	char *kind[5] = {"not set", "natural", "periodic", "geographic", "extended data"};

	if (h->no_BC) return (GMT_NOERROR);	/* Told not to deal with BC stuff */

	if (GMT->common.n.bc_set) {	/* Override BCs via -n+<BC> */
		while (GMT->common.n.BC[i]) {
			switch (GMT->common.n.BC[i]) {
				case 'g':	/* Geographic sets everything */
					h->gn = h->gs = true;
					h->BC[0] = h->BC[1] = h->BC[2] = h->BC[3] = GMT_BC_IS_GEO;
					break;
				case 'n':	/* Natural BCs */
					if (GMT->common.n.BC[i+1] == 'x') { h->BC[0] = h->BC[1] = GMT_BC_IS_NATURAL; i++; }
					else if (GMT->common.n.BC[i+1] == 'y') { h->BC[2] = h->BC[3] = GMT_BC_IS_NATURAL; i++; }
					else h->BC[0] = h->BC[1] = h->BC[2] = h->BC[3] = GMT_BC_IS_NATURAL;
					break;
				case 'p':	/* Periodic BCs */
					if (GMT->common.n.BC[i+1] == 'x') { h->BC[0] = h->BC[1] = GMT_BC_IS_PERIODIC; h->nxp = -1; i++; }
					else if (GMT->common.n.BC[i+1] == 'y') { h->BC[2] = h->BC[3] = GMT_BC_IS_PERIODIC; h->nyp = -1; i++; }
					else { h->BC[0] = h->BC[1] = h->BC[2] = h->BC[3] = GMT_BC_IS_PERIODIC; h->nxp = h->nyp = -1; }
					break;
				default:
					printf ("Error: Cannot parse boundary condition %s\n", GMT->common.n.BC);
					return (-1);
					break;
			}
			i++;
		}
		if (h->gn && !(h->BC[0] == GMT_BC_IS_GEO && h->BC[1] == GMT_BC_IS_GEO && h->BC[2] == GMT_BC_IS_GEO && h->BC[3] == GMT_BC_IS_GEO)) {
			printf ("Warning: GMT boundary condition g overrides n[x|y] or p[x|y]\n");
			h->BC[0] = h->BC[1] = h->BC[2] = h->BC[3] = GMT_BC_IS_GEO;
		}
	}
	else {	/* Determine BC based on whether grid is geographic or not */
		type = (GMT_x_is_lon (GMT, GMT_IN)) ? GMT_BC_IS_GEO : GMT_BC_IS_NATURAL;
		for (i = 0; i < 4; i++) if (h->BC[i] == GMT_BC_IS_NOTSET) h->BC[i] = type;
	}

	/* Check if geographic conditions can be used with this grid */
	if (h->gn && !GMT_grd_is_global (GMT, h)) {
		/* User has requested geographical conditions, but grid is not global */
		printf ( "Warning: longitude range too small; geographic boundary condition changed to natural.\n");
		h->nxp = h->nyp = 0;
		h->gn  = h->gs = false;
		for (i = 0; i < 4; i++) if (h->BC[i] == GMT_BC_IS_NOTSET) h->BC[i] = GMT_BC_IS_NATURAL;
	}
	else if (GMT_grd_is_global (GMT, h)) {	/* Grid is truly global */
		double xtest = fmod (180.0, h->inc[GMT_X]) * h->r_inc[GMT_X];
		/* xtest should be within GMT_SMALL of zero or of one.  */
		if (xtest > GMT_SMALL && xtest < (1.0 - GMT_SMALL) ) {
			/* Error.  We need it to divide into 180 so we can phase-shift at poles.  */
			printf ( "Warning: x_inc does not divide 180; geographic boundary condition changed to natural.\n");
			h->nxp = h->nyp = 0;
			h->gn  = h->gs = false;
			for (i = 0; i < 4; i++) if (h->BC[i] == GMT_BC_IS_NOTSET) h->BC[i] = GMT_BC_IS_NATURAL;
		}
		else {
			h->nxp = urint (360.0 * h->r_inc[GMT_X]);
			h->nyp = 0;
			h->gn = ((fabs(h->wesn[YHI] - 90.0)) < (GMT_SMALL * h->inc[GMT_Y]));
			h->gs = ((fabs(h->wesn[YLO] + 90.0)) < (GMT_SMALL * h->inc[GMT_Y]));
			if (!h->gs) h->BC[2] = GMT_BC_IS_NATURAL;
			if (!h->gn) h->BC[3] = GMT_BC_IS_NATURAL;
		}
	}
	else {	/* Either periodic or natural */
		if (h->nxp != 0) h->nxp = (h->registration == GMT_GRID_PIXEL_REG) ? h->nx : h->nx - 1;
		if (h->nyp != 0) h->nyp = (h->registration == GMT_GRID_PIXEL_REG) ? h->ny : h->ny - 1;
	}

	for (i = 1, same = true; same && i < 4; i++) if (h->BC[i] != h->BC[i-1]) same = false;

	if (same)
		printf ( "Chosen boundary condition for all edges: %s\n", kind[h->BC[XLO]]);
	else {
		if (h->BC[XLO] == h->BC[XHI])
			printf ("Chosen boundary condition for left and right edges: %s\n", kind[h->BC[XLO]]);
		else {
			printf ( "Chosen boundary condition for left   edge: %s\n", kind[h->BC[XLO]]);
			printf ("Chosen boundary condition for right  edge: %s\n", kind[h->BC[XHI]]);
		}
		if (h->BC[YLO] == h->BC[YHI])
			printf ("Chosen boundary condition for bottom and top edges: %s\n", kind[h->BC[YLO]]);
		else {
			printf ("Chosen boundary condition for bottom edge: %s\n", kind[h->BC[YLO]]);
			printf ( "Chosen boundary condition for top    edge: %s\n", kind[h->BC[YHI]]);
		}
	}
	/* Set this grid's interpolation parameters */

	h->bcr_interpolant = GMT->common.n.interpolant;
	h->bcr_threshold = GMT->common.n.threshold;
	h->bcr_n = (h->bcr_interpolant == BCR_NEARNEIGHBOR) ? 1 : ((h->bcr_interpolant == BCR_BILINEAR) ? 2 : 4);

	return (GMT_NOERROR);
}

int gmt_init_grdheader (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *header, struct GMT_OPTION *options, double wesn[], double inc[], unsigned int registration, unsigned int mode)
{	/* Convenient way of setting a header struct wesn, inc, and registartion, then compute dimensions, etc. */
	double wesn_dup[4], inc_dup[2];
	//printf("*******************> %s %s %d\n",__FILE__ ,__func__,__LINE__);
	if ((mode & GMT_VIA_OUTPUT) && wesn == NULL && inc == NULL) return (GMT_NOERROR);	/* OK for creating blank container for output */
	if (wesn == NULL) {	/* Must select -R setting */
		if (!GMT->common.R.active) {
			printf ("No wesn given and no -R in effect.  Cannot initialize new grid\n");
			//GMT_exit (GMT, EXIT_FAILURE);
			return EXIT_FAILURE;
		}
	}
	else	/* In case user is passing header->wesn etc we must save them first as GMT_grd_init will clobber them */
		GMT_memcpy (wesn_dup, wesn, 4, double);
	if (inc == NULL) {	/* Must select -I setting */
		if (!GMT->common.API_I.active) {
			printf ( "No inc given and no -I in effect.  Cannot initialize new grid\n");
			//GMT_exit (GMT, EXIT_FAILURE);
			return EXIT_FAILURE;
		}
	}
	else	/* In case user is passing header->inc etc we must save them first as GMT_grd_init will clobber them */
		GMT_memcpy (inc_dup,  inc,  2, double);
	/* Clobber header and reset */
	GMT_grd_init (GMT, header, options, false);	/* This is for new grids only so update is always false */
	if (wesn == NULL)
		GMT_memcpy (header->wesn, GMT->common.R.wesn, 4, double);
	else
		GMT_memcpy (header->wesn, wesn_dup, 4, double);
	if (inc == NULL)
		GMT_memcpy (header->inc, GMT->common.API_I.inc, 2, double);
	else
		GMT_memcpy (header->inc, inc_dup, 2, double);
	/* registration may contain complex mode information */
	if (registration & GMT_GRID_DEFAULT_REG) registration |= GMT->common.r.registration;	/* Set the default registration */
	header->registration = (registration & 1);
	header->complex_mode = (registration & GMT_GRID_IS_COMPLEX_MASK);
	header->grdtype = gmt_get_grdtype (GMT, header);
	GMT_RI_prepare (GMT, header);	/* Ensure -R -I consistency and set nx, ny in case of meter units etc. */


	
	//GMT_err_pass (GMT, GMT_grd_RI_verify (GMT, header, 1), "");

	
	GMT_grd_setpad (GMT, header, GMT->current.io.pad);	/* Assign default GMT pad */
	GMT_set_grddim (GMT, header);	/* Set all dimensions before returning */
	GMT_BC_init (GMT, header);	/* Initialize grid interpolation and boundary condition parameters */
	return (GMT_NOERROR);
}

int GMTAPI_init_grid (struct GMT_CTRL *GMT, struct GMT_OPTION *opt, double *range, double *inc, int registration, unsigned int mode, struct GMT_GRID *G)
{
	gmt_init_grdheader (GMT, G->header, opt, range, inc, registration, mode);
	return (GMT_OK);
}

int GMT_Encode_ID (void *V_API, char *filename, int object_ID)
{
	/* Creates a filename with the embedded GMTAPI Object ID.  Space must exist */

	if (V_API == NULL) return_error (V_API, GMT_NOT_A_SESSION);	/* GMT_Create_Session has not been called */
	if (!filename) return_error (V_API, GMT_MEMORY_ERROR);		/* Oops, not allocated space */
	if (object_ID == GMT_NOTSET) return_error (V_API, GMT_NOT_A_VALID_ID);	/* ID is nont set yet */
	if (object_ID > GMTAPI_MAX_ID) return_error (V_API, GMT_ID_TOO_LARGE);	/* ID is too large to fit in %06d format below */

	sprintf (filename, "@GMTAPI@-%06d", object_ID);	/* Place the object ID in the special GMT API format */
	return (GMT_OK);	/* No error encountered */
}

bool GMT_UTMzone_to_wesn (struct GMT_CTRL *GMT, unsigned int zone_x, char zone_y, int hemi, double wesn[])
{	/* Given the full UTM zone specification, return w/e/s/n */

	bool error = false;

	wesn[XHI] = -180.0 + 6.0 * zone_x;	wesn[XLO] = wesn[XHI] - 6.0;

	if (zone_y == 0) {	/* Latitude zone is not specified */
		if (hemi == -1) {
			wesn[YLO] = -80.0;	wesn[YHI] = 0.0;
		}
		else if (hemi == +1) {
			wesn[YLO] = 0.0;	wesn[YHI] = 84.0;
		}
		else
			error = true;
		return (error);
	}
	else if (zone_y < 'A' || zone_y > 'Z')
		error = true;
	else if (zone_y <= 'B') {
		wesn[YLO] = -90.0;	wesn[YHI] = -80.0;
		wesn[XHI] = 180.0 * (zone_y - 'A');
		wesn[XLO] = wesn[XHI] - 180.0;
	}
	else if (zone_y <= 'I') {	/* I will behave as J */
		wesn[YLO] = -80.0 + 8.0 * (zone_y - 'C');	wesn[YHI] = wesn[YLO] + 8.0;
	}
	else if (zone_y <= 'O') {	/* O will behave as P */
		wesn[YLO] = -80.0 + 8.0 * (zone_y - 'D');	wesn[YHI] = wesn[YLO] + 8.0;
	}
	else if (zone_y <= 'W') {
		wesn[YLO] = -80.0 + 8.0 * (zone_y - 'E');	wesn[YHI] = wesn[YLO] + 8.0;
		if (zone_y == 'V' && zone_x == 31) wesn[XHI] = 3.0;
		if (zone_y == 'V' && zone_x == 32) wesn[XLO] = 3.0;
	}
	else if (zone_y == 'X') {
		wesn[YLO] = 72.0;	wesn[YHI] = 84.0;
		if (zone_x == 31) wesn[XHI] = 9.0;
		if (zone_x == 33) {wesn[XLO] = 9.0; wesn[XHI] = 21.0;}
		if (zone_x == 35) {wesn[XLO] = 21.0; wesn[XHI] = 33.0;}
		if (zone_x == 37) wesn[XLO] = 33.0;
		if (zone_x == 32 || zone_x == 34 || zone_x == 36) error = true;
	}
	else {	/* Y or Z */
		wesn[YLO] = 84.0;	wesn[YHI] = 90.0;
		wesn[XHI] = 180.0 * (zone_y - 'Y');
		wesn[XLO] = wesn[XHI] - 180.0;
	}

	return (error);
}
int GMT_rectR_to_geoR (struct GMT_CTRL *GMT, char unit, double rect[], double out_wesn[], bool get_R)
{
	/* If user gives -Re|f|k|M|n<xmin>/<xmax>/<ymin>/<ymax>[/<zmin>/<zmax>][r] then we must
	 * call GMT_mapproject to convert this to geographic degrees.
	 * get_R is true when this is done to obtain the -R setting.  */

	int object_ID, proj_class;
	uint64_t dim[4] = {1, 1, 2, 2};	/* Just a single data table with one segment with two 2-column records */
	bool was_R, was_J;
	double wesn[4];
	char buffer[GMT_BUFSIZ] = {""}, in_string[GMT_STR16] = {""}, out_string[GMT_STR16] = {""};
	struct GMT_DATASET *In = NULL, *Out = NULL;

	//GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Call GMT_rectR_to_geoR to convert projected -R to geo -R\n");
	if (GMT_is_dnan (GMT->current.proj.lon0)) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Central meridian is not known; cannot convert -R<unit>... to geographic corners\n");
		return (GMT_MAP_NO_PROJECTION);
	}
	if (GMT_is_dnan (GMT->current.proj.lat0)) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Projection standard latitude is not known; cannot convert -R<unit>... to geographic corners\n");
		return (GMT_MAP_NO_PROJECTION);
	}
	/* Create dataset to hold the rect coordinates */
	if ((In = (struct GMT_DATASET *) GMT_Create_Data (GMT->parent, GMT_IS_DATASET, GMT_IS_POINT, 0, dim, NULL, NULL, 0, 0, NULL)) == NULL) return (GMT_MEMORY_ERROR);

	In->table[0]->segment[0]->coord[GMT_X][0] = rect[XLO];
	In->table[0]->segment[0]->coord[GMT_Y][0] = rect[YLO];
	In->table[0]->segment[0]->coord[GMT_X][1] = rect[XHI];
	In->table[0]->segment[0]->coord[GMT_Y][1] = rect[YHI];

	/* Set up machinery to call mapproject */

	/* Register In as input source via ref (this just returns the ID associated with In sinc already registered by GMT_Create_Data) */
	if ((object_ID = GMT_Register_IO (GMT->parent, GMT_IS_DATASET, GMT_IS_REFERENCE, GMT_IS_POINT, GMT_IN, NULL, In)) == GMT_NOTSET) {
		return (GMT->parent->error);
	}
	if (GMT_Encode_ID (GMT->parent, in_string, object_ID) != GMT_OK) {	/* Make filename with embedded object ID */
		return (GMT->parent->error);
	}
	if ((object_ID = GMT_Register_IO (GMT->parent, GMT_IS_DATASET, GMT_IS_DUPLICATE, GMT_IS_POINT, GMT_OUT, NULL, NULL)) == GMT_NOTSET) {
		return (GMT->parent->error);
	}
	if (GMT_Encode_ID (GMT->parent, out_string, object_ID)) {
		return (GMT->parent->error);	/* Make filename with embedded object ID */
	}
	was_R = GMT->common.R.active;	was_J = GMT->common.J.active;
	GMT->common.R.active = GMT->common.J.active = false;	/* To allow new entries */

	/* Determine suitable -R setting for this projection */

	/* Default w/e/s/n is small patch centered on projection center - this may change below */
	wesn[XLO] = GMT->current.proj.lon0 - 1.0;		wesn[XHI] = GMT->current.proj.lon0 + 1.0;
	wesn[YLO] = MAX (GMT->current.proj.lat0 -1.0, -90.0);	wesn[YHI] = MIN (GMT->current.proj.lat0 + 1.0, 90.0);

	proj_class = GMT->current.proj.projection / 100;	/* 1-4 for valid projections */
	if (GMT->current.proj.projection == GMT_AZ_EQDIST) proj_class = 4;	/* Make -JE use global region */
	switch (proj_class) {
		case 1:	/* Cylindrical: pick small equatorial patch centered on central meridian */
			if (GMT->current.proj.projection == GMT_UTM && GMT_UTMzone_to_wesn (GMT, GMT->current.proj.utm_zonex, GMT->current.proj.utm_zoney, GMT->current.proj.utm_hemisphere, wesn))
			{
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: UTM projection insufficiently specified to auto-determine geographic region\n");
				return (GMT_MAP_NO_PROJECTION);
			}
			else {
				wesn[YLO] = -1.0;	wesn[YHI] = 1.0;
			}
			break;
		case 2: /* Conical: Use default patch */
			break;
		case 3: /* Azimuthal: Use default patch, or hemisphere for polar projections */
			wesn[XLO] = GMT->current.proj.lon0 - 180.0;	wesn[XHI] = GMT->current.proj.lon0 + 180.0;
			if (doubleAlmostEqualZero (GMT->current.proj.lat0, 90.0)) {
				wesn[YLO] = 0.0;	wesn[YHI] = 90.0;
			}
			else if (doubleAlmostEqualZero (GMT->current.proj.lat0, -90.0)) {
				wesn[YLO] = -90.0;	wesn[YHI] = 0.0;
			}
			break;
		case 4: /* Global: Give global region */
			wesn[XLO] = 0.0;	wesn[XHI] = 360.0;	wesn[YLO] = -90.0;	wesn[YHI] = 90.0;
			break;
		default:
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: No map projection specified to auto-determine geographic region\n");
			break;
	}
	sprintf (buffer, "-R%g/%g/%g/%g -J%s -I -F%c -C -bi2d -bo2d -<%s ->%s", wesn[XLO], wesn[XHI], wesn[YLO], wesn[YHI], GMT->common.J.string, unit, in_string, out_string);
	//if (get_R) GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Obtaining geographic corner coordinates via mapproject %s\n", buffer);
	//if (GMT_Call_Module (GMT->parent, "mapproject", GMT_MODULE_CMD, buffer) != GMT_OK) {	/* Get the corners in degrees via mapproject */
	//	return (GMT->parent->error);
	//}
	GMT->common.R.active = was_R;	GMT->common.J.active = was_J;
	if ((Out = (struct GMT_DATASET *)GMT_Retrieve_Data (GMT->parent, object_ID)) == NULL) {
		return (GMT->parent->error);
	}
	out_wesn[XLO] = Out->table[0]->segment[0]->coord[GMT_X][0];
	out_wesn[YLO] = Out->table[0]->segment[0]->coord[GMT_Y][0];
	out_wesn[XHI] = Out->table[0]->segment[0]->coord[GMT_X][1];
	out_wesn[YHI] = Out->table[0]->segment[0]->coord[GMT_Y][1];

	//if (get_R) GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "Region selection -R%s is replaced by the equivalent geographic region -R%.12g/%.12g/%.12g/%.12gr\n", GMT->common.R.string, out_wesn[XLO], out_wesn[YLO], out_wesn[XHI], out_wesn[YHI]);

	if (GMT_Destroy_Data (GMT->parent, &In) != GMT_OK) {
		return (GMT->parent->error);
	}
	if (GMT_Destroy_Data (GMT->parent, &Out) != GMT_OK) {
		return (GMT->parent->error);
	}

	return (GMT_NOERROR);
}

enum GMT_enum_units GMT_get_unit_number (struct GMT_CTRL *GMT, char unit) {
	/* Converts character unit (e.g., 'k') to unit number (e.g., GMT_IS_KM) */
	enum GMT_enum_units mode;

	switch (unit) {
		case '\0':
		case 'e':
			mode = GMT_IS_METER;
			break;
		case 'k':
			mode = GMT_IS_KM;
			break;
		case 'M':
			mode = GMT_IS_MILE;
			break;
		case 'n':
			mode = GMT_IS_NAUTICAL_MILE;
			break;
		case 'i':
			mode = GMT_IS_INCH;
			break;
		case 'c':
			mode = GMT_IS_CM;
			break;
		case 'p':
			mode = GMT_IS_PT;
			break;
		case 'f':
			mode = GMT_IS_FOOT;
			break;
		case 'u':
			mode = GMT_IS_SURVEY_FOOT;
			break;
		default:
			mode = GMT_IS_NOUNIT;
	}

	return (mode);
}

int GMT_init_scales (struct GMT_CTRL *GMT, unsigned int unit, double *fwd_scale, double *inv_scale, double *inch_to_unit, double *unit_to_inch, char *unit_name) {
	/* unit is 0-8 (see gmt_project.h for enums) and stands for m, km, mile, nautical mile, inch, cm, point, foot, or (US) survey foot */
	/* fwd_scale is used to convert user distance units to meter */
	/* inv_scale is used to convert meters to user distance units */
	/* inch_to_unit is used to convert internal inches to users units (c, i, p) */
	/* unit_to_inch is used to convert users units (c, i, p) to internal inches */
	/* unit_name (unless NULL) is set to the name of the user's map measure unit (cm/inch/point) */

	if (unit >= GMT_N_UNITS) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "GMT ERROR: Unit id must be 0-%d\n", GMT_N_UNITS-1);
		GMT_exit (GMT, EXIT_FAILURE); return EXIT_FAILURE;
	}

	/* These scales are used when 1:1 is not set to ensure that the
	 * output (or input with -I) is given (taken) in the units set
	 * by PROJ_LENGTH_UNIT */

	switch (GMT->current.setting.proj_length_unit) {
		case GMT_CM:
			*inch_to_unit = 2.54;
			if (unit_name) strcpy (unit_name, "cm");
			break;
		case GMT_INCH:
			*inch_to_unit = 1.0;
			if (unit_name) strcpy (unit_name, "inch");
			break;
		case GMT_PT:
			*inch_to_unit = 72.0;
			if (unit_name) strcpy (unit_name, "point");
			break;
		case GMT_M:
			if (GMT_compat_check (GMT, 4)) {
				*inch_to_unit = 0.0254;
				if (unit_name) strcpy (unit_name, "m");
			}
			break;
	}
	*unit_to_inch = 1.0 / (*inch_to_unit);
	*fwd_scale = 1.0 / GMT->current.proj.m_per_unit[unit];
	*inv_scale = GMT->current.proj.m_per_unit[unit];
	return GMT_OK;
}

bool GMT_check_region (struct GMT_CTRL *GMT, double wesn[])
{	/* If region is given then we must have w < e and s < n */
	return ((wesn[XLO] >= wesn[XHI] || wesn[YLO] >= wesn[YHI]));
}

int gmt_parse_R_option (struct GMT_CTRL *GMT, char *item) {
	unsigned int i, icol, pos, error = 0, n_slash = 0;
	int got, col_type[2], expect_to_read;
	size_t length;
	bool inv_project = false, scale_coord = false;
	char text[GMT_BUFSIZ] = {""}, string[GMT_BUFSIZ] = {""}, r_unit = 0;
	double p[6];
	//printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
	if (!item || !item[0]) return (GMT_PARSE_ERROR);	/* Got nothing */
	/* Parse the -R option.  Full syntax: -R<grdfile> or -Rg or -Rd or -R[g|d]w/e/s/n[/z0/z1][r] */
	length = strlen (item) - 1;
	for (i = 0; i < length; i++) if (item[i] == '/') n_slash++;
	
	strncpy (GMT->common.R.string, item, GMT_LEN256);	/* Verbatim copy */
	if ((item[0] == 'g' || item[0] == 'd') && item[1] == '\0') {	/* Check -Rd|g separately in case user has files called d or g */
		if (item[0] == 'g') {	/* -Rg is shorthand for -R0/360/-90/90 */
			GMT->common.R.wesn[XLO] = 0.0, GMT->common.R.wesn[XHI] = 360.0;
			GMT->current.io.geo.range = GMT_IS_0_TO_P360_RANGE;
		}
		else {			/* -Rd is shorthand for -R-180/+180/-90/90 */
			GMT->common.R.wesn[XLO] = -180.0, GMT->common.R.wesn[XHI] = 180.0;
			GMT->current.io.geo.range = GMT_IS_M180_TO_P180_RANGE;
		}
		GMT->common.R.wesn[YLO] = -90.0;	GMT->common.R.wesn[YHI] = +90.0;
		GMT_set_geographic (GMT, GMT_IN);
		return (GMT_NOERROR);
	}
	if (!GMT_access (GMT, item, R_OK)) {	/* Gave a readable file, presumably a grid */
		struct GMT_GRID *G = NULL;
		//printf("\n++++++++++++++++++++++++++++++>file : %s line : %d func: %s \n",__FILE__,__LINE__,__func__);
		if ((G = (struct GMT_GRID *) GMT_Read_Data (GMT->parent, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, item, NULL)) == NULL) {	/* Read header */
		//	printf("\n++++++++++++++++++++++++++++++>file : %s line : %d func: %s \n",__FILE__,__LINE__,__func__);
			return (GMT->parent->error);
		}
		//printf("\n++++++++++++++++++++++++++++++>file : %s line : %d func: %s \n",__FILE__,__LINE__,__func__);
		GMT_memcpy (&(GMT->current.io.grd_info.grd), G->header, 1, struct GMT_GRID_HEADER);
		if (GMT_Destroy_Data (GMT->parent, &G) != GMT_OK) {
			return (GMT->parent->error);
		}
		if ((GMT->current.proj.projection == GMT_UTM || GMT->current.proj.projection == GMT_TM || GMT->current.proj.projection == GMT_STEREO)) {	/* Perhaps we got an [U]TM or stereographic grid? */
			if (fabs (GMT->current.io.grd_info.grd.wesn[XLO]) > 360.0 || fabs (GMT->current.io.grd_info.grd.wesn[XHI]) > 360.0 \
			  || fabs (GMT->current.io.grd_info.grd.wesn[YLO]) > 90.0 || fabs (GMT->current.io.grd_info.grd.wesn[YHI]) > 90.0) {	/* Yes we probably did, but cannot be sure */
				inv_project = true;
				r_unit = 'e';	/* Must specify the "missing" leading e for meter */
				sprintf (string, "%.16g/%.16g/%.16g/%.16g", GMT->current.io.grd_info.grd.wesn[XLO], GMT->current.io.grd_info.grd.wesn[XHI], GMT->current.io.grd_info.grd.wesn[YLO], GMT->current.io.grd_info.grd.wesn[YHI]);
			}
		}
		if (!inv_project) {	/* Got grid with degrees */
		//	printf("\n++++++++++++++++++++++++++++++>file : %s line : %d func: %s \n",__FILE__,__LINE__,__func__);
			GMT_memcpy (GMT->common.R.wesn, GMT->current.io.grd_info.grd.wesn, 4, double);
			if (GMT->current.io.grd_info.grd.registration == GMT_GRID_NODE_REG && doubleAlmostEqualZero (GMT->common.R.wesn[XHI] - GMT->common.R.wesn[XLO] + GMT->current.io.grd_info.grd.inc[GMT_X], 360.0)) {
				/* Geographic grid with gridline registration that does not contain the repeating column, but is still 360 range */
				//GMT_Report (GMT->parent, GMT_MSG_DEBUG, "-R<file> with gridline registration and non-repeating column detected; return full 360 degree range for -R\n");
				if (GMT_IS_ZERO (GMT->common.R.wesn[XLO]) || doubleAlmostEqualZero (GMT->common.R.wesn[XLO], -180.0))
					GMT->common.R.wesn[XHI] = GMT->common.R.wesn[XLO] + 360.0;
				else
					GMT->common.R.wesn[XLO] = GMT->common.R.wesn[XHI] - 360.0;
			}
			GMT->common.R.wesn[ZLO] = GMT->current.io.grd_info.grd.z_min;	GMT->common.R.wesn[ZHI] = GMT->current.io.grd_info.grd.z_max;
			GMT->current.io.grd_info.active = true;
			return (GMT_NOERROR);
		}
	}
	else if ((item[0] == 'g' || item[0] == 'd') && n_slash == 3) {	/* Here we have a region appended to -Rd|g */
		GMT_set_geographic (GMT, GMT_IN);
		strncpy (string, &item[1], GMT_BUFSIZ);
		GMT->current.io.geo.range = (item[0] == 'g') ? GMT_IS_0_TO_P360_RANGE : GMT_IS_M180_TO_P180_RANGE;
	}
	else if (strchr (GMT_LEN_UNITS2, item[0])) {	/* Specified min/max in projected distance units */
		strncpy (string, &item[1], GMT_BUFSIZ);
		r_unit = item[0];	/* The leading unit */
		if (GMT_IS_LINEAR (GMT))	/* Just scale up the values */
			scale_coord = true;
		else
			inv_project = true;
	}
	else if (item[length] != 'r' && (GMT->current.proj.projection == GMT_UTM || GMT->current.proj.projection == GMT_TM || GMT->current.proj.projection == GMT_STEREO)) {	/* Just _might_ be getting -R in meters, better check */
		double rect[4];
		strncpy (string, item, GMT_BUFSIZ);
		sscanf (string, "%lg/%lg/%lg/%lg", &rect[XLO], &rect[XHI], &rect[YLO], &rect[YHI]);
		if (fabs (rect[XLO]) > 360.0 || fabs (rect[XHI]) > 360.0 || fabs (rect[YLO]) > 90.0 || fabs (rect[YHI]) > 90.0) {	/* Oh, yeah... */
			inv_project = true;
			r_unit = 'e';	/* Must specify the "missing" leading e for meter */
		}
	}
	else
		{/* Plain old -Rw/e/s/n */
		
		strncpy (string, item, GMT_BUFSIZ);
		//printf("file : %s line : %d func: %s Plain old -Rw  :%s\n",__FILE__,__LINE__,__func__, string);
		}

	/* Now decode the string */

	length = strlen (string) - 1;
	col_type[0] = col_type[1] = 0;
	if (string[length] == 'r') {
		GMT->common.R.oblique = true;
		string[strlen(string)-1] = '\0';	/* Remove the trailing r so GMT_scanf will work */
	}
	else
		GMT->common.R.oblique = false;
	i = pos = 0;
	GMT_memset (p, 6, double);

	//GMT->current.io.col_type[GMT_IN][icol] = GMT_IS_UNKNOWN; //added by nishita

	
	while ((GMT_strtok (string, "/", &pos, text))) {
		if (i > 5) {
			error++;
			return (error);		/* Have to break out here to avoid segv on p[6]  */
		}
		/* Figure out what column corresponds to a token to get col_type[GMT_IN] flag  */
		if (i > 3)
			icol = 2;
		else if (GMT->common.R.oblique)
			icol = i%2;
		else
			icol = i/2;
		if (icol < 2 && GMT->current.setting.io_lonlat_toggle[GMT_IN]) icol = 1 - icol;	/* col_types were swapped */
		/* If column is either RELTIME or ABSTIME, use ARGTIME; if inv_project then just read floats via atof */
		if (inv_project)	/* input is distance units */
		{
			p[i] = atof (text);
		}
		else if (GMT->current.io.col_type[GMT_IN][icol] == GMT_IS_UNKNOWN) {	/* No -J or -f set, proceed with caution */
			got = GMT_scanf_arg (GMT, text, GMT->current.io.col_type[GMT_IN][icol], &p[i]);
		//	printf("file : %s line : %d func: %s got %d\n",__FILE__,__LINE__,__func__,got);
			if (got & GMT_IS_GEO)
				GMT->current.io.col_type[GMT_IN][icol] = got;
			else if (got & GMT_IS_RATIME)
				GMT->current.io.col_type[GMT_IN][icol] = got, GMT->current.proj.xyz_projection[icol] = GMT_TIME;
		}
		else {	/* Things are set, do or die */
			expect_to_read = (GMT->current.io.col_type[GMT_IN][icol] & GMT_IS_RATIME) ? GMT_IS_ARGTIME : GMT->current.io.col_type[GMT_IN][icol];
			//printf("file : %s line : %d func: %s projection: %d\n",__FILE__,__LINE__,__func__,GMT->current.proj.projection);
			//expect_to_read = GMT_IS_GEO;
			//printf("file : %s line : %d func: %s  %d  %d\n",__FILE__,__LINE__,__func__,expect_to_read,expect_to_read, GMT->current.io.col_type[GMT_IN][icol]);
			error += GMT_verify_expectations (GMT, expect_to_read, GMT_scanf (GMT, text, expect_to_read, &p[i]), text);
			//printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
		}
		if (error) return (error);

		i++;
	}
	if (GMT->common.R.oblique) double_swap (p[2], p[1]);	/* So w/e/s/n makes sense */
	if (inv_project) {	/* Convert rectangular distances to geographic corner coordinates */
		double wesn[4];
		//printf("\n++++++++++++++++++++++++++++++>file : %s line : %d func: %s \n",__FILE__,__LINE__,__func__);
		GMT->common.R.oblique = false;
		error += GMT_rectR_to_geoR (GMT, r_unit, p, wesn, true);
		GMT_memcpy (p, wesn, 4, double);
		GMT->common.R.oblique = true;
	}
	else if (scale_coord) {	/* Just scale x/y coordinates to meters according to given unit */
		//printf("\n++++++++++++++++++++++++++++++>file : %s line : %d func: %s \n",__FILE__,__LINE__,__func__);
		double fwd_scale, inv_scale = 0.0, inch_to_unit, unit_to_inch;
		int k_unit;
		k_unit = GMT_get_unit_number (GMT, item[0]);
		GMT_init_scales (GMT, k_unit, &fwd_scale, &inv_scale, &inch_to_unit, &unit_to_inch, NULL);
		for (pos = 0; pos < 4; pos++) p[pos] *= inv_scale;
	}

	if (GMT_is_geographic (GMT, GMT_IN)) {	/* Arrange so geographic region always has w < e */
		double w = p[0], e = p[1];
		if (p[0] <= -360.0 || p[1] > 360.0) {	/* Arrange so geographic region always has |w,e| <= 360 */
			double shift = (p[0] <= -360.0) ? 360.0 : -360.0;
			p[0] += shift;	p[1] += shift;
			//GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "Warning -R: Given west and east values [%g %g] were adjusted so not exceed multiples of 360 [%g %g]\n", w, e, p[0], p[1]);
		}
#if 0	/* This causes too much trouble: Better to annoy the person wishing this to work vs annoy all those who made an honest error.  We cannot be mind-readers here so we insist on e > w */
		else if (p[0] > p[1]) {	/* Arrange so geographic region always has w < e */
			if (GMT->current.io.geo.range == GMT_IS_M180_TO_P180_RANGE) p[0] -= 360.0; else p[1] += 360.0;
			//GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "Warning -R: Given west and east values [%g %g] were adjusted so west < east [%g %g]\n", w, e, p[0], p[1]);
		}
#endif
	}
	if (i < 4 || i > 6 || ((!GMT->common.R.oblique && GMT_check_region (GMT, p)) || (i == 6 && p[4] >= p[5]))) error++;
	GMT_memcpy (GMT->common.R.wesn, p, 6, double);	/* This will probably be reset by GMT_map_setup */
	error += GMT_check_condition (GMT, i == 6 && !GMT->current.proj.JZ_set, "Error: -R with six parameters requires -Jz|Z\n");

	return (error);
}

struct GMT_OPTION * GMT_Append_Option (void *V_API, struct GMT_OPTION *new_opt, struct GMT_OPTION *head)
{
	/* Append this entry to the end of the linked list */
	struct GMT_OPTION *current = NULL;

	if (V_API == NULL) return_null (V_API, GMT_NOT_A_SESSION);	/* GMT_Create_Session has not been called */
	if (!new_opt) return_null (V_API, GMT_OPTION_IS_NULL);		/* No option was passed */
	if (!new_opt->arg) return_null (V_API, GMT_ARG_IS_NULL);	/* Option argument must not be null pointer */

	if (head == NULL) return (new_opt);			/* No list yet, let new_opt become the start of the list */

	/* Here the list already existed with head != NULL */
	
	if (new_opt->option == GMT_OPT_OUTFILE) {	/* Only allow one output file on command line */
		/* Search for existing output option */
		for (current = head; current->next && current->option != GMT_OPT_OUTFILE; current = current->next);
		if (current->option == GMT_OPT_OUTFILE) return_null (V_API, GMT_ONLY_ONE_ALLOWED);	/* Cannot have > 1 output file */
		/* Here current is at end so no need to loop again */
	}
	else {	/* Not an output file name so just go to end of list */
		for (current = head; current->next; current = current->next);			/* Go to end of list */
	}
	/* Append new_opt to the list */
	current->next = new_opt;
	new_opt->previous = current;
	
	//printf("in append option new_opt %s\n",new_opt->arg);
	
	return (head);		/* Return head of list */
}



void GMT_chop (char *string) {
	/* Chops off any CR or LF and terminates string */
	char *p;
	assert (string != NULL); /* NULL pointer */
  /* if (string == NULL) return; / NULL pointer */
	if ((p = strpbrk (string, "\r\n")))
		/* Overwrite 1st CR or LF with terminate string */
		*p = '\0';
}

struct GMT_OPTION *GMT_Make_Option (void *V_API, char option, char *arg)
{
	/* Create a structure option given the option character and the optional argument arg */
	struct GMT_OPTION *new_opt = NULL;
	struct GMTAPI_CTRL *API = NULL;




	if (V_API == NULL) return_null (V_API, GMT_NOT_A_SESSION);	/* GMT_Create_Session has not been called */
	API = gmt_get_api_ptr (V_API);	/* Cast void pointer to a GMTAPI_CTRL pointer */

	/* Here we have a program-specific option or a file name.  In either case we create a new option structure */

	new_opt = GMT_memory (API->GMT, NULL, 1, struct GMT_OPTION);	/* Allocate one option structure */

	new_opt->option = option;		/* Assign which option character was used */
	if (!arg)				/* If arg is a NULL pointer: */
		new_opt->arg = strdup ("");	/* Copy an empty string, such that new_opt->arg[0] = '\0', which avoids */
						/* segfaults later on since few functions check for NULL pointers  */
	else {					/* If arg is set to something (may be an empty string): */
		new_opt->arg = strdup (arg);	/* Allocate space for the argument and duplicate it in the option structure */
		GMT_chop (new_opt->arg);	/* Get rid of any trailing \n \r from cross-binary use in Cygwin/Windows */
	}

	return (new_opt);		/* Pass back the pointer to the allocated option structure */
}

bool GMT_not_numeric (struct GMT_CTRL *GMT, char *text)
{
	/* true if text cannot represent a valid number  However,
	 * false does not therefore mean we have a valid number because
	 * <date>T<clock> representations may use all kinds
	 * of punctuations or letters according to the various format
	 * settings in gmt.conf.  Here we just rule out things
	 * that we are sure of. */

	unsigned int i, k, n_digits = 0, n_period = 0, period = 0, n_plus = 0, n_minus = 0;
	static char *valid = "0123456789-+.:WESNT" GMT_LEN_UNITS GMT_DIM_UNITS;
	if (!text) return (true);		/* NULL pointer */
	if (!strlen (text)) return (true);	/* Blank string */
	if (isalpha ((int)text[0])) return (true);	/* Numbers cannot start with letters */
	if (!(text[0] == '+' || text[0] == '-' || text[0] == '.' || isdigit ((int)text[0]))) return (true);	/* Numbers must be [+|-][.][<digits>] */
	for (i = 0; text[i]; i++) {	/* Check each character */
		/* First check for ASCII values that should never appear in any number */
		if (!strchr (valid, text[i])) return (true);	/* Found a char not among valid letters */
		if (isdigit ((int)text[i])) n_digits++;
		if (text[i] == '.') {
			n_period++;
			period = i;
		}
		if (text[i] == '+') n_plus++;
		if (text[i] == '-') n_minus++;
	}
	if (n_digits == 0 || n_period > 1 || (n_plus + n_minus) > 2) return (true);
	if (n_period) {	/* Check if we have filename.ext with ext having no numbers */
		for (i = period + 1, n_digits = k = 0; text[i]; i++, k++) if (isdigit ((int)text[i])) n_digits++;
		if (k > 0 && n_digits == 0) return (true);	/* Probably a file */
	}
	return (false);	/* This may in fact be numeric */
}


struct GMT_OPTION * GMT_Create_Options (void *V_API, int n_args_in, void *in)
{
	/* This function will loop over the n_args_in supplied command line arguments (in) and
	 * returns a linked list of GMT_OPTION structures for each program option.
	 * These will in turn will be processed by the program-specific option parsers.
	 * What actually happens is controlled by n_args_in.  There are these cases:
	 * n_args_in < 0 : This means that in is already a linked list and we just return it.
	 * n_args_in == 0: in is a single text string with multiple options (e.g., "-R0/2/0/5 -Jx1 -O -m > file")
	 *		   and we must first break this command string into separate words.
	 * n_args_in > 0 : in is an array of text strings (e.g., argv[]).
	 *
	 * Note: 1. If argv[0] is the calling program name, make sure to pass argc-1, args+1 instead.
	 *	 2. in == NULL is allowed only for n_args_in <= 0 (an empty list of options).
	 */

	int error = GMT_OK;
	unsigned int arg, first_char, n_args;
	char option, **args = NULL, **new_args = NULL;
	struct GMT_OPTION *head = NULL, *new_opt = NULL;
	struct GMT_CTRL *G = NULL;
	struct GMTAPI_CTRL *API = NULL;

	//printf("file : %s line : %d func: %s number of arg : %s\n",__FILE__,__LINE__,__func__, argv[0]);

	if (V_API == NULL) return_null (V_API, GMT_NOT_A_SESSION);		/* GMT_Create_Session has not been called */
	if (in == NULL && n_args_in) return_null (V_API, GMT_ARGV_LIST_NULL);	/* Gave no argument pointer but said we had at least 1 */
	if (in == NULL) return (NULL);	/* Gave no argument pointer so a null struct is returned */
	if (n_args_in < 0) return (in);	/* Already converted to linked list */
	API = gmt_get_api_ptr (V_API);	/* Convert API to a GMTAPI_CTRL pointer */
	G = API->GMT;	/* GMT control structure */
	if (n_args_in == 0) {	/* Check if a single command line, if so break into tokens */
		unsigned int pos = 0, new_n_args = 0, k;
		bool quoted;
		size_t n_alloc = GMT_SMALL_CHUNK;
		char p[GMT_BUFSIZ] = {""}, *txt_in = in;	/* Passed a single text string */
		new_args = GMT_memory (G, NULL, n_alloc, char *);
		/* txt_in can contain options that take multi-word text strings, e.g., -B+t"My title".  We avoid the problem of splitting
		 * these items by temporarily replacing spaces inside quoted strings with ASCII 31 US (Unit Separator), do the strtok on
		 * space, and then replace all ASCII 31 with space at the end (we do the same for tab using ASCII 29 GS (group separator) */
		for (k = 0, quoted = false; txt_in[k]; k++) {
			if (txt_in[k] == '\"') quoted = !quoted;	/* Initially false, becomes true at start of quote, then false when exit quote */
			else if (quoted && txt_in[k] == '\t') txt_in[k] = 29;
			else if (quoted && txt_in[k] == ' ')  txt_in[k] = 31;
		}
		while ((GMT_strtok (txt_in, " ", &pos, p))) {	/* Break up string into separate words, and strip off double quotes */
			unsigned int i, o;
			for (k = 0; p[k]; k++) if (p[k] == 29) p[k] = '\t'; else if (p[k] == 31) p[k] = ' ';	/* Replace spaces and tabs masked above */
			for (i = o = 0; p[i]; i++) if (p[i] != '\"') p[o++] = p[i];	/* Ignore double quotes */
			p[o] = '\0';
			new_args[new_n_args++] = strdup (p);
			if (new_n_args == n_alloc) {
				n_alloc += GMT_SMALL_CHUNK;
				new_args = GMT_memory (G, new_args, n_alloc, char *);
			}
		}
		for (k = 0; txt_in[k]; k++) if (txt_in[k] == 29) txt_in[k] = '\t'; else if (txt_in[k] == 31) txt_in[k] = ' ';	/* Replace spaces and tabs masked above */
		args = new_args;
		n_args = new_n_args;
	}
	else {
		args = (char **)in;	/* Gave an argv[] argument */
		n_args = n_args_in;
	}
	if (args == NULL && n_args) return_null (API, GMT_ARGV_LIST_NULL);	/* Conflict between # of args and args being NULL */
	
	for (arg = 0; arg < n_args; arg++) {	/* Loop over all command arguments */
		
		if (!args[arg]) continue;	/* Skip any NULL arguments quietly */
		
		
		if (args[arg][0] == '<' && !args[arg][1] && (arg+1) < n_args && args[arg+1][0] != '-')	/* string command with "< file" for input */
			first_char = 0, option = GMT_OPT_INFILE, arg++;
		else if (args[arg][0] == '>' /*&& !args[arg][1]/* && (arg+1) < n_args && args[arg+1][0] != '-'*/)	/* string command with "> file" for output */
			{
			
			first_char = 1, option = GMT_OPT_OUTFILE/*, arg++*/; //first_char = 1 prebiously
			}
		else if (args[arg][0] == '+' && !args[arg][1] && n_args == 1)	/* extended synopsis + */
			first_char = 1, option = GMT_OPT_USAGE, G->common.synopsis.extended = true;
		else if (args[arg][0] != '-')	/* Probably a file (could also be a gmt/grdmath OPERATOR or number, to be handled later by GMT_Make_Option) */
			first_char = 0, option = GMT_OPT_INFILE;
		else if (!args[arg][1])	/* Found the special synopsis option "-" */
			first_char = 1, option = GMT_OPT_SYNOPSIS;
		else if (!strcmp(args[arg], "--help"))	/* Translate '--help' to '-?' */
			first_char = 6, option = GMT_OPT_USAGE;
		else if ((isdigit ((int)args[arg][1]) || args[arg][1] == '.') && !GMT_not_numeric (API->GMT, args[arg])) /* A negative number, most likely; convert to "file" for now */
				first_char = 0, option = GMT_OPT_INFILE;
		else	/* Most likely found a regular option flag (e.g., -D45.0/3) */
			first_char = 2, option = args[arg][1];

		if ((new_opt = GMT_Make_Option (API, option, &args[arg][first_char])) == NULL) return_null (API, error);	/* Create the new option structure given the args, or return the error */
		head = GMT_Append_Option (API, new_opt, head);		/* Hook new option to the end of the list (or initiate list if head == NULL) */
	}
	if (n_args_in == 0) {	/* Free up temporary arg list */
		for (arg = 0; arg < n_args; arg++) free (new_args[arg]);
		GMT_free (G, new_args);
	}

	return (head);		/* We return the linked list */
}
void gmt_init_unit_conversion (struct GMT_CTRL *GMT) {
	/* Loads the m_per_unit array with the scaling factors that converts various units to meters.
	 * Also sets all the names for the units.
	 * See gmt_project.h for enums that can be used as array indices) */

	GMT->current.proj.m_per_unit[GMT_IS_METER]		= 1.0;				/* m in m */
	GMT->current.proj.m_per_unit[GMT_IS_KM]			= METERS_IN_A_KM;		/* m in km */
	GMT->current.proj.m_per_unit[GMT_IS_MILE]		= METERS_IN_A_MILE;		/* m in miles */
	GMT->current.proj.m_per_unit[GMT_IS_NAUTICAL_MILE]	= METERS_IN_A_NAUTICAL_MILE;	/* m in nautical mile */
	GMT->current.proj.m_per_unit[GMT_IS_INCH]		= 0.0254;			/* m in inch */
	GMT->current.proj.m_per_unit[GMT_IS_CM]			= 0.01;				/* m in cm */
	GMT->current.proj.m_per_unit[GMT_IS_PT]			= 0.0254 / 72.0;		/* m in point */
	GMT->current.proj.m_per_unit[GMT_IS_FOOT]		= METERS_IN_A_FOOT;		/* m in foot */
	GMT->current.proj.m_per_unit[GMT_IS_SURVEY_FOOT]	= METERS_IN_A_SURVEY_FOOT;	/* m in US Survey foot */
	
	strcpy (GMT->current.proj.unit_name[GMT_IS_METER],		"m");
	strcpy (GMT->current.proj.unit_name[GMT_IS_KM],		 	"km");
	strcpy (GMT->current.proj.unit_name[GMT_IS_MILE],		"mile");
	strcpy (GMT->current.proj.unit_name[GMT_IS_NAUTICAL_MILE], 	"nautical mile");
	strcpy (GMT->current.proj.unit_name[GMT_IS_INCH],		"inch");
	strcpy (GMT->current.proj.unit_name[GMT_IS_CM],		 	"cm");
	strcpy (GMT->current.proj.unit_name[GMT_IS_PT],		 	"point");
	strcpy (GMT->current.proj.unit_name[GMT_IS_FOOT],		"foot");
	strcpy (GMT->current.proj.unit_name[GMT_IS_SURVEY_FOOT],	"survey foot");
}

void gmt_set_today (struct GMT_CTRL *GMT)
{	/* Gets the rata die of today */
	time_t right_now = time (NULL);			/* Unix time right now */
	struct tm *timeinfo = localtime ( &right_now );
	GMT->current.time.today_rata_die = GMT_rd_from_gymd (GMT, 1900 + timeinfo->tm_year, timeinfo->tm_mon + 1, timeinfo->tm_mday);
}

int GMT_dummy_grd_info (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *header) {
	//printf("file : %s line : %d func: %s \n",__FILE__,__LINE__,__func__);
	return (GMT_GRDIO_UNKNOWN_FORMAT);
}

int GMT_dummy_grd_read (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *header, float *grid, double wesn[], unsigned int *pad, unsigned int complex_mode) {
	return (GMT_GRDIO_UNKNOWN_FORMAT);
}
void GMT_grdio_init (struct GMT_CTRL *GMT) {
	unsigned int id;
	//printf("file : %s line : %d func: %s \n",__FILE__,__LINE__,__func__);

	
	id                        = k_grd_unknown_fmt;
	GMT->session.grdformat[id]  = "Unknown grid format";
	GMT->session.readinfo[id]   = &GMT_dummy_grd_info;
	GMT->session.updateinfo[id] = &GMT_dummy_grd_info;
	GMT->session.writeinfo[id]  = &GMT_dummy_grd_info;
	GMT->session.readgrd[id]    = &GMT_dummy_grd_read;
	GMT->session.writegrd[id]   = &GMT_dummy_grd_read;

	/* FORMAT: GMT netCDF-based (byte) grdio (COARDS compliant) */

	id                        = GMT_GRID_IS_NB;
	GMT->session.grdformat[id]  = "nb = GMT netCDF format (8-bit integer), " GMT_NC_CONVENTION;
	GMT->session.readinfo[id]   = &GMT_nc_read_grd_info;
	GMT->session.updateinfo[id] = &GMT_nc_update_grd_info;
	GMT->session.writeinfo[id]  = &GMT_nc_write_grd_info;
	GMT->session.readgrd[id]    = &GMT_nc_read_grd;
	GMT->session.writegrd[id]   = &GMT_nc_write_grd;

	/* FORMAT: GMT netCDF-based (short) grdio (COARDS compliant) */

	id                        = GMT_GRID_IS_NS;
	GMT->session.grdformat[id]  = "ns = GMT netCDF format (16-bit integer), " GMT_NC_CONVENTION;
	GMT->session.readinfo[id]   = &GMT_nc_read_grd_info;
	GMT->session.updateinfo[id] = &GMT_nc_update_grd_info;
	GMT->session.writeinfo[id]  = &GMT_nc_write_grd_info;
	GMT->session.readgrd[id]    = &GMT_nc_read_grd;
	GMT->session.writegrd[id]   = &GMT_nc_write_grd;

	/* FORMAT: GMT netCDF-based (int) grdio (COARDS compliant) */

	id                        = GMT_GRID_IS_NI;
	GMT->session.grdformat[id]  = "ni = GMT netCDF format (32-bit integer), " GMT_NC_CONVENTION;
	GMT->session.readinfo[id]   = &GMT_nc_read_grd_info;
	GMT->session.updateinfo[id] = &GMT_nc_update_grd_info;
	GMT->session.writeinfo[id]  = &GMT_nc_write_grd_info;
	GMT->session.readgrd[id]    = &GMT_nc_read_grd;
	GMT->session.writegrd[id]   = &GMT_nc_write_grd;

	/* FORMAT: GMT netCDF-based (float) grdio (COARDS compliant) */

	id                        = GMT_GRID_IS_NF;
	GMT->session.grdformat[id]  = "nf = GMT netCDF format (32-bit float), " GMT_NC_CONVENTION;
	GMT->session.readinfo[id]   = &GMT_nc_read_grd_info;
	GMT->session.updateinfo[id] = &GMT_nc_update_grd_info;
	GMT->session.writeinfo[id]  = &GMT_nc_write_grd_info;
	GMT->session.readgrd[id]    = &GMT_nc_read_grd;
	GMT->session.writegrd[id]   = &GMT_nc_write_grd;

	/* FORMAT: GMT netCDF-based (double) grdio (COARDS compliant) */

	id                        = GMT_GRID_IS_ND;
	GMT->session.grdformat[id]  = "nd = GMT netCDF format (64-bit float), " GMT_NC_CONVENTION;
	GMT->session.readinfo[id]   = &GMT_nc_read_grd_info;
	GMT->session.updateinfo[id] = &GMT_nc_update_grd_info;
	GMT->session.writeinfo[id]  = &GMT_nc_write_grd_info;
	GMT->session.readgrd[id]    = &GMT_nc_read_grd;
	GMT->session.writegrd[id]   = &GMT_nc_write_grd;

	/* FORMAT: GMT netCDF-based (byte) grdio */

	id                        = GMT_GRID_IS_CB;
	GMT->session.grdformat[id]  = "cb = GMT netCDF format (8-bit integer, deprecated)";
	GMT->session.readinfo[id]   = &GMT_dummy_grd_info;
	GMT->session.updateinfo[id] = &GMT_dummy_grd_info;
	GMT->session.writeinfo[id]  = &GMT_dummy_grd_info;
	GMT->session.readgrd[id]    = &GMT_dummy_grd_info;
	GMT->session.writegrd[id]   = &GMT_dummy_grd_info;

	/* FORMAT: GMT netCDF-based (short) grdio */

	id                        = GMT_GRID_IS_CS;
	GMT->session.grdformat[id]  = "cs = GMT netCDF format (16-bit integer, deprecated)";
	GMT->session.readinfo[id]   = &GMT_dummy_grd_info;
	GMT->session.updateinfo[id] = &GMT_dummy_grd_info;
	GMT->session.writeinfo[id]  = &GMT_dummy_grd_info;
	GMT->session.readgrd[id]    = &GMT_dummy_grd_info;
	GMT->session.writegrd[id]   = &GMT_dummy_grd_info;

	/* FORMAT: GMT netCDF-based (int) grdio */

	id                        = GMT_GRID_IS_CI;
	GMT->session.grdformat[id]  = "ci = GMT netCDF format (32-bit integer, deprecated)";
	GMT->session.readinfo[id]   = &GMT_dummy_grd_info;
	GMT->session.updateinfo[id] = &GMT_dummy_grd_info;
	GMT->session.writeinfo[id]  = &GMT_dummy_grd_info;
	GMT->session.readgrd[id]    = &GMT_dummy_grd_info;
	GMT->session.writegrd[id]   = &GMT_dummy_grd_info;

	/* FORMAT: GMT netCDF-based (float) grdio */

	id                        = GMT_GRID_IS_CF;
	GMT->session.grdformat[id]  = "cf = GMT netCDF format (32-bit float, deprecated)";
	GMT->session.readinfo[id]   = &GMT_dummy_grd_info;
	GMT->session.updateinfo[id] = &GMT_dummy_grd_info;
	GMT->session.writeinfo[id]  = &GMT_dummy_grd_info;
	GMT->session.readgrd[id]    = &GMT_dummy_grd_info;
	GMT->session.writegrd[id]   = &GMT_dummy_grd_info;

	/* FORMAT: GMT netCDF-based (double) grdio */

	id                        = GMT_GRID_IS_CD;
	GMT->session.grdformat[id]  = "cd = GMT netCDF format (64-bit float, deprecated)";
	GMT->session.readinfo[id]   = &GMT_dummy_grd_info;
	GMT->session.updateinfo[id] = &GMT_dummy_grd_info;
	GMT->session.writeinfo[id]  = &GMT_dummy_grd_info;
	GMT->session.readgrd[id]    = &GMT_dummy_grd_info;
	GMT->session.writegrd[id]   = &GMT_dummy_grd_info;

	/* FORMAT: GMT native binary (bit) grdio */

	id                        = GMT_GRID_IS_BM;
	GMT->session.grdformat[id]  = "bm = GMT native, C-binary format (bit-mask)";
	GMT->session.readinfo[id]   = &GMT_dummy_grd_info;
	GMT->session.updateinfo[id] = &GMT_dummy_grd_info;
	GMT->session.writeinfo[id]  = &GMT_dummy_grd_info;
	GMT->session.readgrd[id]    = &GMT_dummy_grd_info;
	GMT->session.writegrd[id]   = &GMT_dummy_grd_info;

	/* FORMAT: GMT native binary (byte) grdio */

	id                        = GMT_GRID_IS_BB;
	GMT->session.grdformat[id]  = "bb = GMT native, C-binary format (8-bit integer)";
	GMT->session.readinfo[id]   = &GMT_dummy_grd_info;
	GMT->session.updateinfo[id] = &GMT_dummy_grd_info;
	GMT->session.writeinfo[id]  = &GMT_dummy_grd_info;
	GMT->session.readgrd[id]    = &GMT_dummy_grd_info;
	GMT->session.writegrd[id]   = &GMT_dummy_grd_info;

	/* FORMAT: GMT native binary (short) grdio */

	id                        = GMT_GRID_IS_BS;
	GMT->session.grdformat[id]  = "bs = GMT native, C-binary format (16-bit integer)";
	GMT->session.readinfo[id]   = &GMT_dummy_grd_info;
	GMT->session.updateinfo[id] = &GMT_dummy_grd_info;
	GMT->session.writeinfo[id]  = &GMT_dummy_grd_info;
	GMT->session.readgrd[id]    = &GMT_dummy_grd_info;
	GMT->session.writegrd[id]   = &GMT_dummy_grd_info;

	/* FORMAT: GMT native binary (int) grdio */

	id                        = GMT_GRID_IS_BI;
	GMT->session.grdformat[id]  = "bi = GMT native, C-binary format (32-bit integer)";
	GMT->session.readinfo[id]   = &GMT_dummy_grd_info;
	GMT->session.updateinfo[id] = &GMT_dummy_grd_info;
	GMT->session.writeinfo[id]  = &GMT_dummy_grd_info;
	GMT->session.readgrd[id]    = &GMT_dummy_grd_info;
	GMT->session.writegrd[id]   = &GMT_dummy_grd_info;

	/* FORMAT: GMT native binary (float) grdio */

	id                        = GMT_GRID_IS_BF;
	GMT->session.grdformat[id]  = "bf = GMT native, C-binary format (32-bit float)";
	GMT->session.readinfo[id]   = &GMT_dummy_grd_info;
	GMT->session.updateinfo[id] = &GMT_dummy_grd_info;
	GMT->session.writeinfo[id]  = &GMT_dummy_grd_info;
	GMT->session.readgrd[id]    = &GMT_dummy_grd_info;
	GMT->session.writegrd[id]   = &GMT_dummy_grd_info;

	/* FORMAT: GMT native binary (double) grdio */

	id                        = GMT_GRID_IS_BD;
	GMT->session.grdformat[id]  = "bd = GMT native, C-binary format (64-bit float)";
	GMT->session.readinfo[id]   = &GMT_dummy_grd_info;
	GMT->session.updateinfo[id] = &GMT_dummy_grd_info;
	GMT->session.writeinfo[id]  = &GMT_dummy_grd_info;
	GMT->session.readgrd[id]    = &GMT_dummy_grd_info;
	GMT->session.writegrd[id]   = &GMT_dummy_grd_info;

	/* FORMAT: SUN 8-bit standard rasterfile grdio */

	id                        = GMT_GRID_IS_RB;
	GMT->session.grdformat[id]  = "rb = SUN rasterfile format (8-bit standard)";
	GMT->session.readinfo[id]   = &GMT_dummy_grd_info;
	GMT->session.updateinfo[id] = &GMT_dummy_grd_info;
	GMT->session.writeinfo[id]  = &GMT_dummy_grd_info;
	GMT->session.readgrd[id]    = &GMT_dummy_grd_info;
	GMT->session.writegrd[id]   = &GMT_dummy_grd_info;

	/* FORMAT: NOAA NGDC MGG grid format */

	id                        = GMT_GRID_IS_RF;
	GMT->session.grdformat[id]  = "rf = GEODAS grid format GRD98 (NGDC)";
	GMT->session.readinfo[id]   = &GMT_dummy_grd_info;
	GMT->session.updateinfo[id] = &GMT_dummy_grd_info;
	GMT->session.writeinfo[id]  = &GMT_dummy_grd_info;
	GMT->session.readgrd[id]    = &GMT_dummy_grd_info;
	GMT->session.writegrd[id]   = &GMT_dummy_grd_info;

	/* FORMAT: GMT native binary (float) grdio (Surfer format) */

	id                        = GMT_GRID_IS_SF;
	GMT->session.grdformat[id]  = "sf = Golden Software Surfer format 6 (32-bit float)";
	GMT->session.readinfo[id]   = &GMT_dummy_grd_info;
	GMT->session.updateinfo[id] = &GMT_dummy_grd_info;
	GMT->session.writeinfo[id]  = &GMT_dummy_grd_info;
	GMT->session.readgrd[id]    = &GMT_dummy_grd_info;
	GMT->session.writegrd[id]   = &GMT_dummy_grd_info;

	/* FORMAT: GMT native binary (double) grdio (Surfer format) */

	id                        = GMT_GRID_IS_SD;
	GMT->session.grdformat[id]  = "sd = Golden Software Surfer format 7 (64-bit float, read-only)";
	GMT->session.readinfo[id]   = &GMT_dummy_grd_info;
	GMT->session.updateinfo[id] = &GMT_dummy_grd_info;
	GMT->session.writeinfo[id]  = &GMT_dummy_grd_info;
	GMT->session.readgrd[id]    = &GMT_dummy_grd_info;
	GMT->session.writegrd[id]   = &GMT_dummy_grd_info;

	/* FORMAT: GMT native binary (float) grdio (AGC format) */

	id                        = GMT_GRID_IS_AF;
	GMT->session.grdformat[id]  = "af = Atlantic Geoscience Center format AGC (32-bit float)";
	GMT->session.readinfo[id]   = &GMT_dummy_grd_info;
	GMT->session.updateinfo[id] = &GMT_dummy_grd_info;
	GMT->session.writeinfo[id]  = &GMT_dummy_grd_info;
	GMT->session.readgrd[id]    = &GMT_dummy_grd_info;
	GMT->session.writegrd[id]   = &GMT_dummy_grd_info;

	/* FORMAT: ESRI Arc/Info ASCII Interchange Grid format (integer) */

	id                        = GMT_GRID_IS_EI;
	GMT->session.grdformat[id]  = "ei = ESRI Arc/Info ASCII Grid Interchange format (ASCII integer)";
	GMT->session.readinfo[id]   = &GMT_dummy_grd_info;
	GMT->session.updateinfo[id] = &GMT_dummy_grd_info;
	GMT->session.writeinfo[id]  = &GMT_dummy_grd_info;
	GMT->session.readgrd[id]    = &GMT_dummy_grd_info;
	GMT->session.writegrd[id]   = &GMT_dummy_grd_info;

	/* FORMAT: ESRI Arc/Info ASCII Interchange Grid format (float) */

	id                        = GMT_GRID_IS_EF;
	GMT->session.grdformat[id]  = "ef = ESRI Arc/Info ASCII Grid Interchange format (ASCII float)";
	GMT->session.readinfo[id]   = &GMT_dummy_grd_info;
	GMT->session.updateinfo[id] = &GMT_dummy_grd_info;
	GMT->session.writeinfo[id]  = &GMT_dummy_grd_info;
	GMT->session.readgrd[id]    = &GMT_dummy_grd_info;
	GMT->session.writegrd[id]   = &GMT_dummy_grd_info;

	/* FORMAT: Import via the GDAL interface */

	id                        = GMT_GRID_IS_GD;
#ifdef HAVE_GDAL
	GMT->session.grdformat[id]  = "gd = Import/export through GDAL";
	GMT->session.readinfo[id]   = &GMT_gdal_read_grd_info;
	GMT->session.updateinfo[id] = &GMT_gdal_write_grd_info;
	GMT->session.writeinfo[id]  = &GMT_gdal_write_grd_info;
	GMT->session.readgrd[id]    = &GMT_gdal_read_grd;
	GMT->session.writegrd[id]   = &GMT_gdal_write_grd;
#else
	GMT->session.grdformat[id]  = "gd = Import/export through GDAL (not supported)";
	GMT->session.readinfo[id]   = &GMT_dummy_grd_info;
	GMT->session.updateinfo[id] = &GMT_dummy_grd_info;
	GMT->session.writeinfo[id]  = &GMT_dummy_grd_info;
	GMT->session.readgrd[id]    = &GMT_dummy_grd_read;
	GMT->session.writegrd[id]   = &GMT_dummy_grd_read;
#endif

	/* ----------------------------------------------
	 * ADD CUSTOM FORMATS BELOW AS THEY ARE NEEDED */

}

int GMT_set_env (struct GMT_CTRL *GMT)
{
	char *this_c = NULL, path[PATH_MAX+1];

#ifdef SUPPORT_EXEC_IN_BINARY_DIR
	/* If SUPPORT_EXEC_IN_BINARY_DIR is defined we try to set the share dir to
	 * ${GMT_SOURCE_DIR}/share and the user dir to ${GMT_BINARY_DIR}/share in
	 * order to simplify debugging and running in GMT_BINARY_DIR, e.g., when
	 * debugging with Xcode or Visual Studio. This saves us from setting the
	 * env variables GMT_SHAREDIR and GMT_USERDIR and we do not have to install
	 * src/share in its destination dir. */

	/* Only true, when we are running in a subdir of GMT_BINARY_DIR_SRC_DEBUG: */
	bool running_in_bindir_src = !strncmp (GMT->init.runtime_bindir, GMT_BINARY_DIR_SRC_DEBUG, strlen(GMT_BINARY_DIR_SRC_DEBUG));
#endif

	/* Determine GMT->session.SHAREDIR (directory containing coast, cpt, etc. subdirectories) */

	//if ((this_c = getenv ("GMT5_SHAREDIR")) != NULL
	//		&& GMT_verify_sharedir_version (this_c) )
		/* GMT5_SHAREDIR was set */
		//GMT->session.SHAREDIR = strdup (this_c);
		
		GMT->session.SHAREDIR = strdup ("data"); //nishita
		
//	else if ((this_c = getenv ("GMT_SHAREDIR")) != NULL
		//	&& GMT_verify_sharedir_version (this_c) )
		/* GMT_SHAREDIR was set */
	//	GMT->session.SHAREDIR = strdup (this_c);
#ifdef SUPPORT_EXEC_IN_BINARY_DIR
	//else if ( running_in_bindir_src && GMT_verify_sharedir_version (GMT_SHARE_DIR_DEBUG) )
		/* Use ${GMT_SOURCE_DIR}/share to simplify debugging and running in GMT_BINARY_DIR */
	//	GMT->session.SHAREDIR = strdup (GMT_SHARE_DIR_DEBUG);
#endif
	//else if ( GMT_verify_sharedir_version (GMT_SHARE_DIR) )
		/* Found in hardcoded GMT_SHARE_DIR */
		//GMT->session.SHAREDIR = strdup (GMT_SHARE_DIR);
	//else {
		/* SHAREDIR still not found, make a smart guess based on runpath: */
		//if ( GMT_guess_sharedir (path, GMT->init.runtime_bindir) )
		//	GMT->session.SHAREDIR = strdup (path);
		//else {
			/* Still not found */
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL,
			//	"Error: Could not locate share directory that matches the current GMT version %s.\n",
			//	GMT_PACKAGE_VERSION_WITH_SVN_REVISION);
			//GMT_exit (GMT, EXIT_FAILURE); return EXIT_FAILURE;
		//}
	//}
	//DOS_path_fix (GMT->session.SHAREDIR);

	/* Determine HOMEDIR (user home directory) */

	//if ((this_c = getenv ("HOME")) != NULL)
		/* HOME was set */

		char path1[GMT_BUFSIZ] ;//= "/home/nishita/buildGMTUsingMake/data";
				if (getcwd(path1, sizeof(path1)) == NULL)
				{
					printf("Error \n\n\n");
					return -1;		
				}
				/*********Added by Maria**********/
				FILE *f;
				char value[100],value2[100];
				strcat(path1,"/datapath.dat");
				printf("PATH BEFORE: %s\n",path1);
				f=fopen(path1,"r");
				fscanf(f,"%s",value);
				fclose(f);
				value[strlen(value)-6]='\0';
				//strcat(value,'\0');
				strcpy(path,value);
				/*********************************/

		GMT->session.HOMEDIR = strdup (path); //nishita
#ifdef WIN32
	//else if ((this_c = getenv ("HOMEPATH")) != NULL)
		/* HOMEPATH was set */
		//GMT->session.HOMEDIR = strdup (this_c);
#endif
	//else {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Error: Could not determine home directory!\n");
		//GMT_exit (GMT, EXIT_FAILURE); return EXIT_FAILURE;
//	}
	//DOS_path_fix (GMT->session.HOMEDIR);

	/* Determine GMT_USERDIR (directory containing user replacements contents in GMT_SHAREDIR) */

	if ((this_c = getenv ("GMT_USERDIR")) != NULL)
		/* GMT_USERDIR was set */
		GMT->session.USERDIR = strdup (this_c);
#ifdef SUPPORT_EXEC_IN_BINARY_DIR
	//else if ( running_in_bindir_src && access (GMT_USER_DIR_DEBUG, R_OK|X_OK) == 0 )
		/* Use ${GMT_BINARY_DIR}/share to simplify debugging and running in GMT_BINARY_DIR */
		//GMT->session.USERDIR = strdup (GMT_USER_DIR_DEBUG);
#endif
	else {
		/* Use default path for GMT_USERDIR (~/.gmt) */
		sprintf (path, "%s/%s", GMT->session.HOMEDIR, "data");//nishita
		//printf("HOMEDIR HEREEEEEEEEEE %s\n",path);
		GMT->session.USERDIR = strdup (path);
	}
	//DOS_path_fix (GMT->session.USERDIR);
	if (GMT->session.USERDIR != NULL && access (GMT->session.USERDIR, R_OK)) {
		/* If we cannot access this dir then we won't use it */
		free (GMT->session.USERDIR);
		GMT->session.USERDIR = NULL;
	}

	//if (GMT_compat_check (GMT, 4)) {
		/* Check if obsolete GMT_CPTDIR was specified */

		//if ((this_c = getenv ("GMT_CPTDIR")) != NULL) {
			/* GMT_CPTDIR was set */
		//	GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: Environment variable GMT_CPTDIR was set but is no longer used by GMT.\n");
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: System-wide color tables are in %s/cpt.\n", GMT->session.SHAREDIR);
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: Use GMT_USERDIR (%s) instead and place user-defined color tables there.\n", GMT->session.USERDIR);
		//}
	//}

	/* Determine GMT_DATADIR (data directories) */

	//if ((this_c = getenv ("GMT_DATADIR")) != NULL) {
		/* GMT_DATADIR was set */
		//if (strchr (this_c, PATH_SEPARATOR) || access (this_c, R_OK) == 0) {
			/* A list of directories or a single directory that is accessible */
			//GMT->session.DATADIR = strdup (this_c);
			GMT->session.DATADIR = strdup ("data"); //nishita
			//DOS_path_fix (GMT->session.DATADIR);
		//}
	//}

	/* Determine GMT_TMPDIR (for isolation mode). Needs to exist use it. */

	//if ((this_c = getenv ("GMT_TMPDIR")) != NULL) 
		//{
		/* GMT_TMPDIR was set */
		//if (access (this_c, R_OK|W_OK|X_OK)) {
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: Environment variable GMT_TMPDIR was set to %s, but directory is not accessible.\n", this_c);
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: GMT_TMPDIR needs to have mode rwx. Isolation mode switched off.\n");
			//GMT->session.TMPDIR = NULL;
		//}
		//else
			GMT->session.TMPDIR = strdup ("data");//nishita
		//DOS_path_fix (GMT->session.TMPDIR);
	//}

	return GMT_OK;
}

struct GMT_CTRL *New_GMT_Ctrl (char *session, unsigned int pad) {
	/* Allocate and initialize a new common control structure */
	int i;
	char path[PATH_MAX+1];
	char *unit_name[4] = {"cm", "inch", "m", "point"};
	double u2u[4][4] = {	/* Local initialization of unit conversion factors */
		{   1.00,    1.0/2.54,    0.01,         72.0/2.54 },
		{   2.54,    1.0,         0.0254,       72.0 },
		{ 100.00,    100.0/2.54,  1.0,          72.0/0.0254 },
		{ 2.54/72.0, 1.0/72.0,    0.0254/72.0,  1.0 }
	};
	struct GMT_PROJ4 GMT_proj4[GMT_N_PROJ4] = {
	{ "aea" 	 , GMT_ALBERS },
	{ "aeqd"	 , GMT_AZ_EQDIST },
	{ "cyl_stere", GMT_CYL_STEREO },
	{ "cass"	 , GMT_CASSINI },
	{ "cea" 	 , GMT_CYL_EQ },
	{ "eck4"	 , GMT_ECKERT4 },
	{ "eck6"	 , GMT_ECKERT6 },
	{ "eqc" 	 , GMT_CYL_EQDIST },
	{ "eqdc"	 , GMT_ECONIC },
	{ "gnom"	 , GMT_GNOMONIC },
	{ "hammer"	 , GMT_HAMMER },
	{ "laea"	 , GMT_LAMB_AZ_EQ },
	{ "lcc" 	 , GMT_LAMBERT },
	{ "merc"	 , GMT_MERCATOR },
	{ "mill"	 , GMT_MILLER },
	{ "moll"	 , GMT_MOLLWEIDE },
	{ "nsper"	 , GMT_GENPER },
	{ "omerc"	 , GMT_OBLIQUE_MERC },
	{ "omercp"	 , GMT_OBLIQUE_MERC_POLE },
	{ "ortho"	 , GMT_ORTHO },
	{ "polar"	 , GMT_POLAR },
	{ "poly"	 , GMT_POLYCONIC },
	{ "robin"	 , GMT_ROBINSON },
	{ "sinu"	 , GMT_SINUSOIDAL },
	{ "stere"	 , GMT_STEREO },
	{ "tmerc"	 , GMT_TM },
	{ "utm" 	 , GMT_UTM },
	{ "vandg"	 , GMT_VANGRINTEN },
	{ "wintri"	 , GMT_WINKEL },
	{ "xy"		 , GMT_LINEAR },
	{ "z"		 , GMT_ZAXIS }
};
	struct GMT_CTRL *GMT = NULL;
	GMT = calloc (1U, sizeof (struct GMT_CTRL));
	/* Assign the three std* pointers */

	GMT->session.std[GMT_IN]  = stdin;
	GMT->session.std[GMT_OUT] = stdout;
	GMT->session.std[GMT_ERR] = stderr;

	/* Set default verbosity level */
	GMT->current.setting.verbose = GMT_MSG_COMPAT;

	/* We don't know the module or library names yet */
	GMT->init.module_name = GMT->init.module_lib = NULL;

	
	GMT_set_env (GMT);	/* Get GMT_SHAREDIR and other environment path parameters */

	//GMT_init_fonts (GMT);	/* Load in available font names */
	
	//printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
	
	GMT_init_distaz (GMT, GMT_MAP_DIST_UNIT, GMT_GREATCIRCLE, GMT_MAP_DIST);	/* Default spherical distance calculations in m */

	GMT->current.map.n_lon_nodes = 360;
	GMT->current.map.n_lat_nodes = 180;
	GMT->current.map.frame.check_side = false;
	GMT->current.map.frame.horizontal = 0;
	GMT->current.map.dlon = (GMT->common.R.wesn[XHI] - GMT->common.R.wesn[XLO]) / GMT->current.map.n_lon_nodes;
	GMT->current.map.dlat = (GMT->common.R.wesn[YHI] - GMT->common.R.wesn[YLO]) / GMT->current.map.n_lat_nodes;
	
	/* PLOT settings */

	GMT->current.plot.mode_3D = 3;	/* Draw both fore and/or back 3-D box lines [1 + 2] */

	/* PROJ settings */
	
	GMT->current.proj.projection = GMT_NO_PROJ;
	/* We need some defaults here for the cases where we do not actually set these with GMT_map_setup */
	/* z_level will be updated in GMT_init_three_D, but if it doesn't, it does not matter,
	 * because by default, z_scale = 0.0 */
	GMT->current.proj.z_level = DBL_MAX;
		GMT->current.proj.xyz_pos[GMT_X] = GMT->current.proj.xyz_pos[GMT_Y] = GMT->current.proj.xyz_pos[GMT_Z] = true;
	GMT->current.proj.z_project.view_azimuth = 180.0;
	GMT->current.proj.z_project.view_elevation = 90.0;
	GMT->current.proj.z_project.plane = -1;	/* Initialize no perspective projection */
	GMT->current.proj.z_project.level = 0.0;
	for (i = 0; i < 4; i++) GMT->current.proj.edge[i] = true;
	
	GMT_grdio_init (GMT);
	//printf("file : %s line : %d func: %s \n",__FILE__,__LINE__,__func__);
	GMT_set_pad (GMT, pad); /* Sets default number of rows/cols for boundary padding in this session */
	
	GMT->current.proj.f_horizon = 90.0;
	GMT->current.proj.proj4 = GMT_memory (GMT, NULL, GMT_N_PROJ4, struct GMT_PROJ4);
	//printf("file : %s line : %d func: %s \n",__FILE__,__LINE__,__func__);
	for (i = 0; i < GMT_N_PROJ4; i++) {	/* Load up proj4 structure once and for all */
		GMT->current.proj.proj4[i].name = strdup (GMT_proj4[i].name);
		GMT->current.proj.proj4[i].id = GMT_proj4[i].id;
	}
	/* TIME_SYSTEM settings */
	//printf("file : %s line : %d func: %s \n",__FILE__,__LINE__,__func__);
	strcpy (GMT->current.setting.time_system.epoch, "2000-01-01T12:00:00");
	GMT->current.setting.time_system.unit = 'd';
	//printf("file : %s line : %d func: %s \n",__FILE__,__LINE__,__func__);
	/* INIT settings */
	//printf("file : %s line : %d func: %s \n",__FILE__,__LINE__,__func__);
	GMT_memcpy (GMT->session.u2u, u2u, 1, u2u);
	for (i = 0; i < 4; i++) strncpy (GMT->session.unit_name[i], unit_name[i], 8U);
	GMT_make_fnan (GMT->session.f_NaN);
	GMT_make_dnan (GMT->session.d_NaN);
	for (i = 0; i < 3; i++) GMT->session.no_rgb[i] = -1.0;
	//printf("file : %s line : %d func: %s \n",__FILE__,__LINE__,__func__);
		return (GMT);
}

void GMT_init_io_columns (struct GMT_CTRL *GMT, unsigned int dir)
{
	/* Initialize (reset) information per column which may have changed due to -i -o */
	unsigned int i;
	for (i = 0; i < GMT_MAX_COLUMNS; i++) GMT->current.io.col[dir][i].col = GMT->current.io.col[dir][i].order = i;	/* Default order */
	if (dir == GMT_OUT) return;
	for (i = 0; i < GMT_MAX_COLUMNS; i++) GMT->current.io.col_skip[i] = false;	/* Consider all input columns */
}


void GMT_io_init (struct GMT_CTRL *GMT)
{
	/* No need to memset the structure to NULL as this is done in gmt_init.h directory.
	 * The assignments here are done once per GMT session as GMT_io_init is called
	 * from GMT_begin.  Some variables may change later due to --PAR=value parsing.
	 * GMT_io_init must be called before parsing of defaults. */

	unsigned int i;

	/* Pointer assignment for default ASCII input functions */

	//printf("file : %s line : %d func: %s inc: %d\n",__FILE__,__LINE__,__func__,GMT->current.io.inc_code[GMT_Y]);

	GMT->current.io.input  = GMT->session.input_ascii = &gmt_ascii_input;
	GMT->current.io.output = &GMT_ascii_output;

	GMT->current.io.ogr_parser = &gmt_ogr_header_parser;		/* Parse OGR header records to start with */

	/* Assign non-zero/NULL initial values */

	GMT->current.io.give_report = true;
	GMT->current.io.seg_no = GMT->current.io.rec_no = GMT->current.io.rec_in_tbl_no = 0;	/* These gets incremented so 1 means 1st record */
	GMT->current.io.warn_geo_as_cartesion = true;	/* Not yet read geographic data while in Cartesian mode so we want to warn if we find it */
	GMT->current.setting.io_seg_marker[GMT_IN] = GMT->current.setting.io_seg_marker[GMT_OUT] = '>';
	strcpy (GMT->current.io.r_mode, "r");
	strcpy (GMT->current.io.w_mode, "w");
	strcpy (GMT->current.io.a_mode, "a+");
	for (i = 0; i < 4; i++) {
		GMT->current.io.date_input.item_order[i] = GMT->current.io.date_input.item_pos[i] = -1;
		GMT->current.io.date_output.item_order[i] = GMT->current.io.date_output.item_pos[i] = -1;
	}
	for (i = 0; i < 3; i++) {
		GMT->current.io.clock_input.order[i] = GMT->current.io.clock_output.order[i] = GMT->current.io.geo.order[i] = -1;
	}
	strcpy (GMT->current.io.clock_input.ampm_suffix[0],  "am");
	strcpy (GMT->current.io.clock_output.ampm_suffix[0], "am");
	strcpy (GMT->current.io.clock_input.ampm_suffix[1],  "pm");
	strcpy (GMT->current.io.clock_output.ampm_suffix[1], "pm");

	GMT_init_io_columns (GMT, GMT_IN);	/* Set default input column order */
	GMT_init_io_columns (GMT, GMT_OUT);	/* Set default output column order */
	for (i = 0; i < 2; i++) GMT->current.io.skip_if_NaN[i] = true;								/* x/y must be non-NaN */
	for (i = 0; i < 2; i++) GMT->current.io.col_type[GMT_IN][i] = GMT->current.io.col_type[GMT_OUT][i] = GMT_IS_UNKNOWN;	/* Must be told [or find out] what x/y are */
	for (i = 2; i < GMT_MAX_COLUMNS; i++) GMT->current.io.col_type[GMT_IN][i] = GMT->current.io.col_type[GMT_OUT][i] = GMT_IS_FLOAT;	/* Other columns default to floats */
	GMT_memset (GMT->current.io.curr_rec, GMT_MAX_COLUMNS, double);	/* Initialize current and previous records to zero */
	GMT_memset (GMT->current.io.prev_rec, GMT_MAX_COLUMNS, double);
}

struct GMT_CTRL *GMT_begin (struct GMTAPI_CTRL *API, char *session, unsigned int pad)
{
	/* GMT_begin is called once by GMT_Create_Session and does basic
	 * one-time initialization of GMT before the GMT modules take over.
	 * It will load in the gmt.conf settings from the share dir and
	 * reset them with the user's gmt.conf settings (if any).
	 * It then does final processing of defaults so that all internal
	 * GMT parameters are properly initialized and ready to go. This
	 * means it is possible to write a functioning GMT application that
	 * does not require the use of any GMT modules.  However,
	 * most GMT applications will call various GMT modules and these
	 * may need to process additional --PAR=value arguments. This will
	 * require renewed processing of defaults and takes place in GMT_begin_module
	 * which is called at the start of all GMT modules.  This basically
	 * performs a save/restore operation so that when the GMT module
	 * returns the GMT structure is restored to its original values.
	 *
	 * Note: We do not call GMT_exit here since API is not given and
	 * API->do_not_exit have not been modified by external API yet.
	 */

	char path[GMT_LEN256] = {""};
	struct GMT_CTRL *GMT = NULL;
	
	GMT = New_GMT_Ctrl (session, pad);	/* Allocate and initialize a new common control structure */	
	
	API->GMT = GMT;
	
	GMT->parent = API;	/* So we know who's your daddy */
	
#if 0
	GMT->PSL = New_PSL_Ctrl ("GMT5");		/* Allocate a PSL control structure */
	if (!GMT->PSL) {
		GMT_Message (API, GMT_TIME_NONE, "Error: Could not initialize PSL - Aborting.\n");
		Free_GMT_Ctrl (GMT);	/* Deallocate control structure */
		return NULL;
	}
	GMT->PSL->init.unit = PSL_INCH;					/* We use inches internally in PSL */
	//PSL_beginsession (GMT->PSL, 0, GMT->session.SHAREDIR, GMT->session.USERDIR);	/* Initializes the session and sets a few defaults */
	/* Reset session defaults to the chosen GMT settings; these are fixed for the entire PSL session */
	//PSL_setdefaults (GMT->PSL, GMT->current.setting.ps_magnify, GMT->current.setting.ps_page_rgb, GMT->current.setting.ps_encoding.name);
#endif
	GMT_io_init (GMT);		/* Init the table i/o structure before parsing GMT defaults */

	gmt_init_unit_conversion (GMT);	/* Set conversion factors from various units to meters */
	
	//GMT_hash_init (GMT, keys_hashnode, GMT_keywords, GMT_N_KEYS, GMT_N_KEYS);	/* Initialize hash table for GMT defaults */

	/* Set up hash table for colornames (used to convert <colorname> to <r/g/b>) */

	//GMT_hash_init (GMT, GMT->session.rgb_hashnode, GMT_color_name, GMT_N_COLOR_NAMES, GMT_N_COLOR_NAMES);

	/* Initialize the standard GMT system default settings from the system file */
#if 0
	sprintf (path, "%s/conf/gmt.conf", GMT->session.SHAREDIR);
	if (access (path, R_OK)) {
		/* Not found in SHAREDIR, try USERDIR instead */
		if (GMT_getuserpath (GMT, "conf/gmt.conf", path) == NULL) {
			GMT_Message (API, GMT_TIME_NONE, "Error: Could not find system defaults file %s - Aborting.\n", path);
			Free_GMT_Ctrl (GMT);	/* Deallocate control structure */
			return NULL;
		}
	}
	GMT_loaddefaults (GMT, path);	/* Load GMT system default settings [and PSL settings if selected] */
	GMT_getdefaults (GMT, NULL);	/* Override using local GMT default settings (if any) [and PSL if selected] */

	/* There is no longer a -m option in GMT 5 so multi segments are now always true.
	   However, in GMT_COMPAT mode the -mi and -mo options WILL turn off multi in the other direction. */
	GMT_set_segmentheader (GMT, GMT_IN, true);
	GMT_set_segmentheader (GMT, GMT_OUT, false);	/* Will be turned true when either of two situation arises: */
	/* 1. We read a multisegment header
	   2. The -g option is set which will create gaps and thus multiple segments
	 */


	/* Initialize the output and plot format machinery for ddd:mm:ss[.xxx] strings from the default format strings.
	 * While this is also done in the default parameter loop it is possible that when a decimal plain format has been selected
	 * the format_float_out string has not yet been processed.  We clear that up by processing again here. */

	gmt_geo_C_format (GMT);
	//gmt_plot_C_format (GMT);
#endif

	/* Set default for -n parameters */
	GMT->common.n.antialias = true; GMT->common.n.interpolant = BCR_BICUBIC; GMT->common.n.threshold = 0.5;

	//gmt_get_history (GMT);	/* Process and store command shorthands passed to the application */

	//if (GMT->current.setting.io_gridfile_shorthand) gmt_setshorthand (GMT);	/* Load the short hand mechanism from gmt.io */

	//GMT_fft_initialization (GMT);	/* Determine which FFT algos are available and set pointers */

	gmt_set_today (GMT);	/* Determine today's rata die value */

	return (GMT);
}

struct GMT_CTRL * GMT_begin_module (struct GMTAPI_CTRL *API, const char *lib_name, const char *mod_name, struct GMT_CTRL **Ccopy)
{	/* All GMT modules (i.e. GMT_psxy, GMT_blockmean, ...) must call GMT_begin_module
	 * as their first call and call GMT_end_module as their last call.  This
	 * allows us to capture the GMT control structure so we can reset all
	 * parameters to what they were before exiting the module. Note:
	 * 1. Session items that remain unchanged are not replicated if allocated separately.
	 * 2. Items that may grow through session are not replicated if allocated separately.
	 */

	unsigned int i;
	struct GMT_CTRL *GMT = API->GMT, *Csave = NULL;
	//printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);	

	
	Csave = calloc (1U, sizeof (struct GMT_CTRL));

	GMT_free_tmp_arrays (GMT);			/* Free temp memory for vector io or processing */

	/* First memcpy over everything; this will include pointer addresses we will have to fix below */
	
	GMT_memcpy (Csave, GMT, 1, struct GMT_CTRL);

	/* Increment level uint64_t */
	GMT->hidden.func_level++;		/* This lets us know how deeply we are nested when a GMT module is called */

	/* Now fix things that were allocated separately from the main GMT structure.  These are usually text strings
	 * that were allocated via strdup since the structure only have a pointer allocated. */

	/* GMT_INIT */
	if (GMT->session.n_user_media) {
		Csave->session.n_user_media = GMT->session.n_user_media;
		//Csave->session.user_media = GMT_memory (GMT, NULL, GMT->session.n_user_media, struct GMT_MEDIA);
		Csave->session.user_media_name = GMT_memory (GMT, NULL, GMT->session.n_user_media, char *);
		for (i = 0; i < GMT->session.n_user_media; i++) Csave->session.user_media_name[i] = strdup (GMT->session.user_media_name[i]);
	}

	/* GMT_PLOT */
	if (GMT->current.plot.n_alloc) {
		Csave->current.plot.n_alloc = GMT->current.plot.n_alloc;
		Csave->current.plot.x = GMT_memory (GMT, NULL, GMT->current.plot.n_alloc, double);
		Csave->current.plot.y = GMT_memory (GMT, NULL, GMT->current.plot.n_alloc, double);
		Csave->current.plot.pen = GMT_memory (GMT, NULL, GMT->current.plot.n_alloc, unsigned int);
		GMT_memcpy (Csave->current.plot.x, GMT->current.plot.x, GMT->current.plot.n_alloc, double);
		GMT_memcpy (Csave->current.plot.y, GMT->current.plot.y, GMT->current.plot.n_alloc, double);
		GMT_memcpy (Csave->current.plot.pen, GMT->current.plot.pen, GMT->current.plot.n_alloc, unsigned int);
	}

	/* GMT_IO */
	Csave->current.io.OGR = GMT_duplicate_ogr (GMT, GMT->current.io.OGR);	/* Duplicate OGR struct, if set */
	GMT_free_ogr (GMT, &(GMT->current.io.OGR), 1);		/* Free up the GMT/OGR structure, if used */

	GMT_memset (Csave->current.io.o_format, GMT_MAX_COLUMNS, char *);
	for (i = 0; i < GMT_MAX_COLUMNS; i++)
		if (GMT->current.io.o_format[i]) Csave->current.io.o_format[i] = strdup (GMT->current.io.o_format[i]);

	/* GMT_COMMON */
	if (GMT->common.U.label) Csave->common.U.label = strdup (GMT->common.U.label);
	for (i = 0; i < GMT->common.a.n_aspatial; i++)
		if (GMT->common.a.name[i]) Csave->common.a.name[i] = strdup (GMT->common.a.name[i]);
	if (GMT->common.h.title) Csave->common.h.title = strdup (GMT->common.h.title);
	if (GMT->common.h.remark) Csave->common.h.remark = strdup (GMT->common.h.remark);
	if (GMT->common.h.colnames) Csave->common.h.colnames = strdup (GMT->common.h.colnames);

	/* Reset all the common.?.active settings to false */

	GMT->common.B.active[0] = GMT->common.B.active[1] = GMT->common.K.active = GMT->common.O.active = false;
	GMT->common.P.active = GMT->common.U.active = GMT->common.V.active = false;
	GMT->common.X.active = GMT->common.Y.active = false;
	GMT->common.R.active = GMT->common.J.active = false;
	GMT->common.a.active = GMT->common.b.active[GMT_IN] = GMT->common.b.active[GMT_OUT] = GMT->common.c.active = false;
	GMT->common.f.active[GMT_IN] = GMT->common.f.active[GMT_OUT] = GMT->common.g.active = GMT->common.h.active = false;
	GMT->common.p.active = GMT->common.s.active = GMT->common.t.active = GMT->common.colon.active = false;
	GMT_memset (GMT->common.b.ncol, 2, int);

	*Ccopy = Csave; /* Pass back out for safe-keeping by the module until GMT_end_module is called */

	GMT->init.module_name = mod_name;
	GMT->init.module_lib  = lib_name;

	//printf("file : %s line : %d func: %s, GMT->current.io.inc %d, GMT->current.io.inc_code[GMT_y] :%d\n",__FILE__,__LINE__,__func__,GMT->current.io.inc_code[GMT_X],GMT->current.io.inc_code[GMT_Y]);

	return (GMT);
}

int GMT_default_error (struct GMT_CTRL *GMT, char option)
{
	/* GMT_default_error ignores all the common options that have already been processed and returns
	 * true if the option is not an already processed common option.
	 */

	int error = 0;

	switch (option) {

		case '-': break;	/* Skip indiscriminently */
		case '>': break;	/* Skip indiscriminently since dealt with internally */
		case 'B': error += (GMT->common.B.active[0] == false && GMT->common.B.active[1] == false); break;
		case 'J': error += GMT->common.J.active == false; break;
		case 'K': error += GMT->common.K.active == false; break;
		case 'O': error += GMT->common.O.active == false; break;
		case 'P': error += GMT->common.P.active == false; break;
		case 'R': error += GMT->common.R.active == false; break;
		case 'U': error += GMT->common.U.active == false; break;
		case 'V': error += GMT->common.V.active == false; break;
		case 'X': error += GMT->common.X.active == false; break;
		case 'Y': error += GMT->common.Y.active == false; break;
		case 'a': error += GMT->common.a.active == false; break;
		case 'b': error += (GMT->common.b.active[GMT_IN] == false && GMT->common.b.active[GMT_OUT] == false); break;
		case 'c': error += GMT->common.c.active == false; break;
		case 'f': error += (GMT->common.f.active[GMT_IN] == false &&  GMT->common.f.active[GMT_OUT] == false); break;
		case 'g': error += GMT->common.g.active == false; break;
		case 'H':
			if (GMT_compat_check (GMT, 4)) {
				error += GMT->common.h.active == false;
			}
			else
				error++;
			break;
		case 'h': error += GMT->common.h.active == false; break;
		case 'i': error += GMT->common.i.active == false; break;
		case 'n': error += GMT->common.n.active == false; break;
		case 'o': error += GMT->common.o.active == false; break;
		case 'Z':
			if (!GMT_compat_check (GMT, 4)) error++;
			break;
		case 'E':
			if (GMT_compat_check (GMT, 4))
				error += GMT->common.p.active == false;
			else
				error++;
			break;
		case 'p': error += GMT->common.p.active == false; break;
		case 'm': if (!GMT_compat_check (GMT, 4)) error++; break;
		case 'S': if (!GMT_compat_check (GMT, 4)) error++; break;
		case 'F': if (!GMT_compat_check (GMT, 4)) error++; break;
		case 'r': error += GMT->common.r.active == false; break;
		case 's': error += GMT->common.s.active == false; break;
		case 't': error += GMT->common.t.active == false; break;
		case ':': error += GMT->common.colon.active == false; break;

		default:
			/* Not a processed common options */
			error++;
			break;
	}

	if (error) //GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Error: Unrecognized option -%c\n", option);
		printf( "Error: Unrecognized option -%c\n", option);
	return (error);
}


