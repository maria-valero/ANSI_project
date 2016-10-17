//#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <malloc.h>
#include <inttypes.h>
#include <unistd.h>
//#include <cstdint>
//#include "../include/stdint.h"
#include "../include/gmt_config.h"
#include "../include/common_math.h"
#include "../include/gmt_init.h"
#include "../include/gmt_resources.h"
#include "../include/surface.h"
#include "../include/gmt_type.h"
#include "../include/gmt_macros.h"
//#include "../include/io.h"
#include "../include/memory.h"
#include "../include/gmt_api.h"
#include "../include/gmt_io.h"
//using namespace std;
#define PROGRAM_NAME "gmt"
#define THIS_MODULE_NAME	"surface"
#define THIS_MODULE_LIB		"core"
#define THIS_MODULE_PURPOSE	"Grid table data using adjustable tension continuous curvature splines"
#define DAT_FILE_LIST "list.dat"
#define BUFFSIZE_NEW 4096


static char *GMT_unique_option[GMT_N_UNIQUE] = {	/* The common GMT command-line options [ just the subset that accepts arguments (e.g., -O is not listed) ] */
#include "../include/gmt_unique.h"
};

#ifndef HAVE___FUNC__
#	ifdef HAVE___FUNCTION__
#		define __func__ __FUNCTION__
#	else
#		define __func__ "<unknown>"
#	endif
#endif


#define SURFACE_OUTSIDE LONG_MAX	/* Index number indicating data is outside usable area */

#define MIN(x, y) (((x) < (y)) ? (x) : (y))	/* min and max value macros */

#define MAX(x, y) (((x) > (y)) ? (x) : (y))

/* Convenience macro for GMT_memory_func */

#define GMT_PROG_OPTIONS "-:RVabfhirs" GMT_OPT("FH")

static inline struct GMTAPI_CTRL * gmt_get_api_ptr (struct GMTAPI_CTRL *ptr) {return (ptr);}


/* Check condition and report error if true */
#define GMT_check_condition(C,condition,...) ((condition) ? 1/*+GMT_Report(C->parent,GMT_MSG_NORMAL,__VA_ARGS__) */: 0)

#define GMT_more_than_once(GMT,active) (GMT_check_condition (GMT, active, "Warning: Option -%c given more than once\n", option))

extern int GMT_grd2xyz (void *V_API, int mode, void *args);
/* Macro to simplify call to memcpy when duplicating values and memset when zeroing out */
//#define GMT_memcpy(to,from,n,type) memcpy(to, from, (n)*sizeof(type))
//#define GMT_memset(array,n,type) memset(array, 0, (n)*sizeof(type))
/*
 Macros returns true if the two coordinates are lon/lat; way should be GMT_IN or GMT_OUT
#define GMT_x_is_lon(C,way) (C->current.io.col_type[way][GMT_X] == GMT_IS_LON)
#define GMT_y_is_lat(C,way) (C->current.io.col_type[way][GMT_Y] == GMT_IS_LAT)
#define GMT_is_geographic(C,way) (GMT_x_is_lon(C,way) && GMT_y_is_lat(C,way))
*/


extern int createSlowMap(double period, char *pa, char *datapath);
extern struct ISO_AXIS CalculateISOMap(double min, double max, int N_bin, char* fileName, char *pa, char *datapath) ;
extern void calculateTiimeDiff(char* stationFile);
extern void createArrivalTime(char *arrivalTimefile) ;

struct SURFACE_GLOBAL {		/* Things needed inside compare function must be global for now */
	int block_ny;		/* Number of nodes in y-dir for a given grid factor */
	double grid_xinc, grid_yinc;	/* size of each grid cell for a given grid factor */
	double x_min, y_min;		/* Lower left corner of grid */
} GMT_Surface_Global;

struct ISO_AXIS {
	int x_coor;
	int y_coor;
} axis;


struct GMTAPI_CTRL * GMT_get_API_ptr (struct GMTAPI_CTRL *ptr)
{	/* Clean casting of void to API pointer at start of a module
 	 * If ptr is NULL we are in deep trouble...
	 */
	if (ptr == NULL)
		{
			//printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
			 return_null (NULL, GMT_NOT_A_SESSION);
		}
	return (ptr);
}

uint64_t iterate (struct GMT_CTRL *GMT, struct SURFACE_INFO *C, int mode)
{	/* Main finite difference solver */
	uint64_t ij, briggs_index, ij_v2, iteration_count = 0, ij_sw, ij_se;
	int i, j, k, kase;
	int x_case, y_case, x_w_case, x_e_case, y_s_case, y_n_case;
	char *iu = C->iu;

	double current_limit = C->converge_limit / C->grid;
	double change, max_change = 0.0, busum, sum_ij;
	double b0, b1, b2, b3, b4, b5;
	float *u = C->Grid->data;

	double x_0_const = 4.0 * (1.0 - C->boundary_tension) / (2.0 - C->boundary_tension);
	double x_1_const = (3 * C->boundary_tension - 2.0) / (2.0 - C->boundary_tension);
	double y_denom = 2 * C->l_epsilon * (1.0 - C->boundary_tension) + C->boundary_tension;
	double y_0_const = 4 * C->l_epsilon * (1.0 - C->boundary_tension) / y_denom;
	double y_1_const = (C->boundary_tension - 2 * C->l_epsilon * (1.0 - C->boundary_tension) ) / y_denom;

	//sprintf (C->format,"%%4ld\t%%c\t%%8" PRIu64 "\t%s\t%s\t%%10" PRIu64 "\n", GMT->current.setting.format_float_out, GMT->current.setting.format_float_out);

	do {
		briggs_index = 0;	/* Reset the constraint table stack pointer  */

		max_change = -1.0;

		/* Fill in auxiliary boundary values (in new way) */

		/* First set d2[]/dn2 = 0 along edges */
		/* New experiment : (1-T)d2[]/dn2 + Td[]/dn = 0  */

		for (i = 0; i < C->nx; i += C->grid) {
			/* set d2[]/dy2 = 0 on south side */
			ij = C->ij_sw_corner + i * C->my;
			/* u[ij - 1] = 2 * u[ij] - u[ij + grid];  */
			u[ij - 1] = (float)(y_0_const * u[ij] + y_1_const * u[ij + C->grid]);
			/* set d2[]/dy2 = 0 on north side */
			ij = C->ij_nw_corner + i * C->my;
			/* u[ij + 1] = 2 * u[ij] - u[ij - grid];  */
			u[ij + 1] = (float)(y_0_const * u[ij] + y_1_const * u[ij - C->grid]);

		}
		if (C->periodic) {	/* Set periodic boundary conditions in longitude */
			for (j = 0; j < C->ny; j += C->grid) {
				ij_sw = C->ij_sw_corner + j;
				ij_se = C->ij_se_corner + j;
				u[ij_sw+C->offset[0][5]]  = u[ij_se+C->offset[20][5]];
				u[ij_se+C->offset[20][6]] = u[ij_sw+C->offset[0][6]];
				u[ij_se] = u[ij_sw] = 0.5f * (u[ij_se] + u[ij_sw]);	/* Set to average of east and west */
			}
		}
		else {	/* Regular natural BC */
			for (j = 0; j < C->ny; j += C->grid) {
				/* set d2[]/dx2 = 0 on west side */
				ij = C->ij_sw_corner + j;
				/* u[ij - my] = 2 * u[ij] - u[ij + grid_east];  */
				u[ij - C->my] = (float)(x_1_const * u[ij + C->grid_east] + x_0_const * u[ij]);
				/* set d2[]/dx2 = 0 on east side */
				ij = C->ij_se_corner + j;
				/* u[ij + my] = 2 * u[ij] - u[ij - grid_east];  */
				u[ij + C->my] = (float)(x_1_const * u[ij - C->grid_east] + x_0_const * u[ij]);
			}
		}

		/* Now set d2[]/dxdy = 0 at each corner */

		ij = C->ij_sw_corner;
		u[ij - C->my - 1] = u[ij + C->grid_east - 1] + u[ij - C->my + C->grid] - u[ij + C->grid_east + C->grid];

		ij = C->ij_nw_corner;
		u[ij - C->my + 1] = u[ij + C->grid_east + 1] + u[ij - C->my - C->grid] - u[ij + C->grid_east - C->grid];

		ij = C->ij_se_corner;
		u[ij + C->my - 1] = u[ij - C->grid_east - 1] + u[ij + C->my + C->grid] - u[ij - C->grid_east + C->grid];

		ij = C->ij_ne_corner;
		u[ij + C->my + 1] = u[ij - C->grid_east + 1] + u[ij + C->my - C->grid] - u[ij - C->grid_east - C->grid];

		/* Now set (1-T)dC/dn + Tdu/dn = 0 at each edge */
		/* New experiment: only dC/dn = 0  */

		x_w_case = 0;
		x_e_case = C->block_nx - 1;
		for (i = 0; i < C->nx; i += C->grid, x_w_case++, x_e_case--) {

			if(x_w_case < 2)
				x_case = x_w_case;
			else if(x_e_case < 2)
				x_case = 4 - x_e_case;
			else
				x_case = 2;

			/* South side */
			kase = x_case * 5;
			ij = C->ij_sw_corner + i * C->my;
			u[ij + C->offset[kase][11]] =
				(float)(u[ij + C->offset[kase][0]] + C->eps_m2*(u[ij + C->offset[kase][1]] + u[ij + C->offset[kase][3]]
					- u[ij + C->offset[kase][8]] - u[ij + C->offset[kase][10]])
					+ C->two_plus_em2 * (u[ij + C->offset[kase][9]] - u[ij + C->offset[kase][2]]) );
				/*  + tense * C->eps_m2 * (u[ij + C->offset[kase][2]] - u[ij + C->offset[kase][9]]) / (1.0 - tense);  */
			/* North side */
			kase = x_case * 5 + 4;
			ij = C->ij_nw_corner + i * C->my;
			u[ij + C->offset[kase][0]] =
				-(float)(-u[ij + C->offset[kase][11]] + C->eps_m2 * (u[ij + C->offset[kase][1]] + u[ij + C->offset[kase][3]]
					- u[ij + C->offset[kase][8]] - u[ij + C->offset[kase][10]])
					+ C->two_plus_em2 * (u[ij + C->offset[kase][9]] - u[ij + C->offset[kase][2]]) );
				/*  - tense * C->eps_m2 * (u[ij + C->offset[kase][2]] - u[ij + C->offset[kase][9]]) / (1.0 - tense);  */
		}

		y_s_case = 0;
		y_n_case = C->block_ny - 1;
		for (j = 0; j < C->ny; j += C->grid, y_s_case++, y_n_case--) {

			if(y_s_case < 2)
				y_case = y_s_case;
			else if(y_n_case < 2)
				y_case = 4 - y_n_case;
			else
				y_case = 2;

			if (C->periodic) {	/* Set periodic boundary conditions in longitude */
				/* West side */
				kase = y_case;
				ij_sw = C->ij_sw_corner + j;
				ij_se = C->ij_se_corner + j;
				u[ij_sw+C->offset[kase][4]] = u[ij_se+C->offset[20+kase][4]];
				/* East side */
				kase = 20 + y_case;
				u[ij_se + C->offset[kase][7]] = u[ij_sw+C->offset[y_case][7]];
			}
			else {	/* Natural BCs */
				/* West side */
				kase = y_case;
				ij = C->ij_sw_corner + j;
				u[ij+C->offset[kase][4]] =
					u[ij + C->offset[kase][7]] + (float)(C->eps_p2 * (u[ij + C->offset[kase][3]] + u[ij + C->offset[kase][10]]
					-u[ij + C->offset[kase][1]] - u[ij + C->offset[kase][8]])
					+ C->two_plus_ep2 * (u[ij + C->offset[kase][5]] - u[ij + C->offset[kase][6]]));
					/*  + tense * (u[ij + C->offset[kase][6]] - u[ij + C->offset[kase][5]]) / (1.0 - tense);  */
				/* East side */
				kase = 20 + y_case;
				ij = C->ij_se_corner + j;
				u[ij + C->offset[kase][7]] =
					- (float)(-u[ij + C->offset[kase][4]] + C->eps_p2 * (u[ij + C->offset[kase][3]] + u[ij + C->offset[kase][10]]
					- u[ij + C->offset[kase][1]] - u[ij + C->offset[kase][8]])
					+ C->two_plus_ep2 * (u[ij + C->offset[kase][5]] - u[ij + C->offset[kase][6]]) );
					/*  - tense * (u[ij + C->offset[kase][6]] - u[ij + C->offset[kase][5]]) / (1.0 - tense);  */
			}
		}



		/* That's it for the boundary points.  Now loop over all data  */

		x_w_case = 0;
		x_e_case = C->block_nx - 1;
		for (i = 0; i < C->nx; i += C->grid, x_w_case++, x_e_case--) {

			if(x_w_case < 2)
				x_case = x_w_case;
			else if(x_e_case < 2)
				x_case = 4 - x_e_case;
			else
				x_case = 2;

			y_s_case = 0;
			y_n_case = C->block_ny - 1;

			ij = C->ij_sw_corner + i * C->my;

			for (j = 0; j < C->ny; j += C->grid, ij += C->grid, y_s_case++, y_n_case--) {

				if (iu[ij] == 5) continue;	/* Point is fixed  */

				if(y_s_case < 2)
					y_case = y_s_case;
				else if(y_n_case < 2)
					y_case = 4 - y_n_case;
				else
					y_case = 2;

				kase = x_case * 5 + y_case;
				sum_ij = 0.0;

				if (iu[ij] == 0) {		/* Point is unconstrained  */
					for (k = 0; k < 12; k++) {
						sum_ij += (u[ij + C->offset[kase][k]] * C->coeff[0][k]);
					}
				}
				else {				/* Point is constrained  */

					b0 = C->briggs[briggs_index].b[0];
					b1 = C->briggs[briggs_index].b[1];
					b2 = C->briggs[briggs_index].b[2];
					b3 = C->briggs[briggs_index].b[3];
					b4 = C->briggs[briggs_index].b[4];
					b5 = C->briggs[briggs_index].b[5];
					briggs_index++;
					if (iu[ij] < 3) {
						if (iu[ij] == 1) {	/* Point is in quadrant 1  */
							busum = b0 * u[ij + C->offset[kase][10]]
								+ b1 * u[ij + C->offset[kase][9]]
								+ b2 * u[ij + C->offset[kase][5]]
								+ b3 * u[ij + C->offset[kase][1]];
						}
						else {			/* Point is in quadrant 2  */
							busum = b0 * u[ij + C->offset[kase][8]]
								+ b1 * u[ij + C->offset[kase][9]]
								+ b2 * u[ij + C->offset[kase][6]]
								+ b3 * u[ij + C->offset[kase][3]];
						}
					}
					else {
						if (iu[ij] == 3) {	/* Point is in quadrant 3  */
							busum = b0 * u[ij + C->offset[kase][1]]
								+ b1 * u[ij + C->offset[kase][2]]
								+ b2 * u[ij + C->offset[kase][6]]
								+ b3 * u[ij + C->offset[kase][10]];
						}
						else {		/* Point is in quadrant 4  */
							busum = b0 * u[ij + C->offset[kase][3]]
								+ b1 * u[ij + C->offset[kase][2]]
								+ b2 * u[ij + C->offset[kase][5]]
								+ b3 * u[ij + C->offset[kase][8]];
						}
					}
					for (k = 0; k < 12; k++) {
						sum_ij += (u[ij + C->offset[kase][k]] * C->coeff[1][k]);
					}
					sum_ij = (sum_ij + C->a0_const_2 * (busum + b5))
						/ (C->a0_const_1 + C->a0_const_2 * b4);
				}

				/* New relaxation here  */
				sum_ij = u[ij] * C->relax_old + sum_ij * C->relax_new;

				if (C->constrained) {	/* Must check limits.  Note lower/upper is in standard scanline format and need ij_v2! */
					ij_v2 = GMT_IJP (C->Grid->header, C->ny - j - 1, i);
					if (C->set_low && !GMT_is_fnan (C->Low->data[ij_v2]) && sum_ij < C->Low->data[ij_v2])
						sum_ij = C->Low->data[ij_v2];
					else if (C->set_high && !GMT_is_fnan (C->High->data[ij_v2]) && sum_ij > C->High->data[ij_v2])
						sum_ij = C->High->data[ij_v2];
				}

				change = fabs (sum_ij - u[ij]);
				u[ij] = (float)sum_ij;
				if (change > max_change) max_change = change;
			}
		}
		iteration_count++;
		C->total_iterations++;
		max_change *= C->z_scale;	/* Put max_change into z units  */
	//	GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, C->format,
			//C->grid, C->mode_type[mode], iteration_count, max_change, current_limit, C->total_iterations);

	} while (max_change > current_limit && iteration_count < C->max_iterations);

	//GMT_Report (GMT->parent, GMT_MSG_VERBOSE, C->format,
		//C->grid, C->mode_type[mode], iteration_count, max_change, current_limit, C->total_iterations);

	return (iteration_count);
}

void *New_surface_Ctrl (struct GMT_CTRL *GMT) {	/* Allocate and initialize a new control structure */
	struct SURFACE_CTRL *C;

	C = GMT_memory (GMT,NULL, 1, struct SURFACE_CTRL);

	/* Initialize values whose defaults are not 0/false/NULL */
	C->N.value = 250;
	C->A.value = 1.0;
	C->Z.value = 1.4;

	return (C);
}

void find_nearest_point (struct SURFACE_INFO *C)
{	/* Determines the nearest data point epr bin and sets the
	 * Briggs parameters or, if really close, sets the node value */
	uint64_t ij_v2, k, last_index, iu_index, briggs_index;
	int i, j, block_i, block_j;
	double x0, y0, dx, dy, xys, xy1, btemp, b0, b1, b2, b3, b4, b5;
	float z_at_node, *u = C->Grid->data;
	char *iu = C->iu;
	struct GMT_GRID_HEADER *h = C->Grid->header;

	last_index = UINTMAX_MAX;
	C->small = 0.05 * ((C->grid_xinc < C->grid_yinc) ? C->grid_xinc : C->grid_yinc);

	for (i = 0; i < C->nx; i += C->grid)	/* Reset grid info */
		for (j = 0; j < C->ny; j += C->grid)
			iu[C->ij_sw_corner + i*C->my + j] = 0;

	briggs_index = 0;
	for (k = 0; k < C->npoints; k++) {	/* Find constraining value  */
		if (C->data[k].index != last_index) {
			block_i = (int)C->data[k].index/C->block_ny;
			block_j = (int)C->data[k].index%C->block_ny;
			last_index = C->data[k].index;
	 		iu_index = C->ij_sw_corner + (block_i * C->my + block_j) * C->grid;
	 		x0 = h->wesn[XLO] + block_i*C->grid_xinc;
	 		y0 = h->wesn[YLO] + block_j*C->grid_yinc;
	 		dx = (C->data[k].x - x0)*C->r_grid_xinc;
	 		dy = (C->data[k].y - y0)*C->r_grid_yinc;
	 		if (fabs(dx) < C->small && fabs(dy) < C->small) {	/* Close enough to assign value to node */
	 			iu[iu_index] = 5;
	 			/* v3.3.4: NEW CODE
	 			 * Since point is basically moved from (dx, dy) to (0,0) we must adjust for
	 			 * the C->small change in the planar trend between the two locations, and then
	 			 * possibly clip the range if constraining surfaces were given.  Note that
	 			 * dx, dy is in -1/1 range normalized by (grid * x|y_inc) so to recover the
	 			 * dx,dy in final grid fractions we must scale by grid */

	 			z_at_node = C->data[k].z + (float) (C->r_z_scale * C->grid * (C->plane_c1 * dx + C->plane_c2 * dy));
	 			if (C->constrained) {	/* Must use ij_v2 since constrained grids are in standard scanline format */
					ij_v2 = GMT_IJP (C->Grid->header, C->ny - block_j * C->grid - 1, block_i * C->grid);
					if (C->set_low  && !GMT_is_fnan (C->Low->data[ij_v2]) && z_at_node < C->Low->data[ij_v2])
						z_at_node = C->Low->data[ij_v2];
					else if (C->set_high && !GMT_is_fnan (C->High->data[ij_v2]) && z_at_node > C->High->data[ij_v2])
						z_at_node = C->High->data[ij_v2];
	 			}
	 			u[iu_index] = z_at_node;
	 		}
	 		else {
	 			if (dx >= 0.0) {
	 				if (dy >= 0.0)
	 					iu[iu_index] = 1;
	 				else
	 					iu[iu_index] = 4;
	 			}
	 			else {
	 				if (dy >= 0.0)
	 					iu[iu_index] = 2;
	 				else
	 					iu[iu_index] = 3;
	 			}
	 			dx = fabs(dx);
	 			dy = fabs(dy);
	 			btemp = 2 * C->one_plus_e2 / ( (dx + dy) * (1.0 + dx + dy) );
	 			b0 = 1.0 - 0.5 * (dx + (dx * dx)) * btemp;
	 			b3 = 0.5 * (C->e_2 - (dy + (dy * dy)) * btemp);
	 			xys = 1.0 + dx + dy;
	 			xy1 = 1.0 / xys;
	 			b1 = (C->e_2 * xys - 4 * dy) * xy1;
	 			b2 = 2 * (dy - dx + 1.0) * xy1;
	 			b4 = b0 + b1 + b2 + b3 + btemp;
	 			b5 = btemp * C->data[k].z;
	 			C->briggs[briggs_index].b[0] = b0;
	 			C->briggs[briggs_index].b[1] = b1;
	 			C->briggs[briggs_index].b[2] = b2;
	 			C->briggs[briggs_index].b[3] = b3;
	 			C->briggs[briggs_index].b[4] = b4;
	 			C->briggs[briggs_index].b[5] = b5;
	 			briggs_index++;
	 		}
	 	}
	 }
}

void set_grid_parameters (struct SURFACE_INFO *C)
{	/* Updates the grid space parameters given the new C->grid setting */
	GMT_Surface_Global.block_ny = C->block_ny = (C->ny - 1) / C->grid + 1;
	C->block_nx = (C->nx - 1) / C->grid + 1;
	GMT_Surface_Global.grid_xinc = C->grid_xinc = C->grid * C->Grid->header->inc[GMT_X];
	GMT_Surface_Global.grid_yinc = C->grid_yinc = C->grid * C->Grid->header->inc[GMT_Y];
	C->grid_east = C->grid * C->my;
	C->r_grid_xinc = 1.0 / C->grid_xinc;
	C->r_grid_yinc = 1.0 / C->grid_yinc;
}

void set_coefficients (struct SURFACE_INFO *C)
{	/* These are the coefficients in the finite-difference expressionss */
	double e_4, loose, a0;

	loose = 1.0 - C->interior_tension;
	C->e_2 = C->l_epsilon * C->l_epsilon;
	e_4 = C->e_2 * C->e_2;
	C->eps_p2 = C->e_2;
	C->eps_m2 = 1.0/C->e_2;
	C->one_plus_e2 = 1.0 + C->e_2;
	C->two_plus_ep2 = 2.0 + 2.0*C->eps_p2;
	C->two_plus_em2 = 2.0 + 2.0*C->eps_m2;

	C->x_edge_const = 4 * C->one_plus_e2 - 2 * (C->interior_tension / loose);
	C->e_m2 = 1.0 / C->e_2;
	C->y_edge_const = 4 * (1.0 + C->e_m2) - 2 * (C->interior_tension * C->e_m2 / loose);


	a0 = 1.0 / ( (6 * e_4 * loose + 10 * C->e_2 * loose + 8 * loose - 2 * C->one_plus_e2) + 4*C->interior_tension*C->one_plus_e2);
	C->a0_const_1 = 2 * loose * (1.0 + e_4);
	C->a0_const_2 = 2.0 - C->interior_tension + 2 * loose * C->e_2;

	C->coeff[1][4] = C->coeff[1][7] = -loose;
	C->coeff[1][0] = C->coeff[1][11] = -loose * e_4;
	C->coeff[0][4] = C->coeff[0][7] = -loose * a0;
	C->coeff[0][0] = C->coeff[0][11] = -loose * e_4 * a0;
	C->coeff[1][5] = C->coeff[1][6] = 2 * loose * C->one_plus_e2;
	C->coeff[0][5] = C->coeff[0][6] = (2 * C->coeff[1][5] + C->interior_tension) * a0;
	C->coeff[1][2] = C->coeff[1][9] = C->coeff[1][5] * C->e_2;
	C->coeff[0][2] = C->coeff[0][9] = C->coeff[0][5] * C->e_2;
	C->coeff[1][1] = C->coeff[1][3] = C->coeff[1][8] = C->coeff[1][10] = -2 * loose * C->e_2;
	C->coeff[0][1] = C->coeff[0][3] = C->coeff[0][8] = C->coeff[0][10] = C->coeff[1][1] * a0;

	C->e_2 *= 2;		/* We will need these in boundary conditions  */
	C->e_m2 *= 2;

	C->ij_sw_corner = 2 * C->my + 2;			/*  Corners of array of actual data  */
	C->ij_se_corner = C->ij_sw_corner + (C->nx - 1) * C->my;
	C->ij_nw_corner = C->ij_sw_corner + C->ny - 1;
	C->ij_ne_corner = C->ij_se_corner + C->ny - 1;

}

void initialize_grid (/*struct GMT_CTRL *GMT,*/ struct SURFACE_INFO *C)
{	/*
	 * For the initial gridsize, compute weighted averages of data inside the search radius
	 * and assign the values to u[i,j] where i,j are multiples of gridsize.
	 */
	uint64_t index_1, index_2, k, k_index;
	int irad, jrad, i, j, imin, imax, jmin, jmax, ki, kj;
	double r, rfact, sum_w, sum_zw, weight, x0, y0;
	float *u = C->Grid->data;
	struct GMT_GRID_HEADER *h = C->Grid->header;

	 irad = irint (ceil(C->radius/C->grid_xinc));
	 jrad = irint (ceil(C->radius/C->grid_yinc));
	 rfact = -4.5/(C->radius*C->radius);
	 for (i = 0; i < C->block_nx; i ++ ) {
	 	x0 = h->wesn[XLO] + i*C->grid_xinc;
	 	for (j = 0; j < C->block_ny; j ++ ) {
	 		y0 = h->wesn[YLO] + j*C->grid_yinc;
	 		imin = i - irad;
	 		if (imin < 0) imin = 0;
	 		imax = i + irad;
	 		if (imax >= C->block_nx) imax = C->block_nx - 1;
	 		jmin = j - jrad;
	 		if (jmin < 0) jmin = 0;
	 		jmax = j + jrad;
	 		if (jmax >= C->block_ny) jmax = C->block_ny - 1;
	 		index_1 = imin*C->block_ny + jmin;
	 		index_2 = imax*C->block_ny + jmax + 1;
	 		sum_w = sum_zw = 0.0;
	 		k = 0;
	 		while (k < C->npoints && C->data[k].index < index_1) k++;
	 		for (ki = imin; k < C->npoints && ki <= imax && C->data[k].index < index_2; ki++) {
	 			for (kj = jmin; k < C->npoints && kj <= jmax && C->data[k].index < index_2; kj++) {
	 				k_index = ki*C->block_ny + kj;
	 				while (k < C->npoints && C->data[k].index < k_index) k++;
	 				while (k < C->npoints && C->data[k].index == k_index) {
	 					r = (C->data[k].x-x0)*(C->data[k].x-x0) + (C->data[k].y-y0)*(C->data[k].y-y0);
	 					weight = exp (rfact*r);
	 					sum_w += weight;
	 					sum_zw += weight*C->data[k].z;
	 					k++;
	 				}
	 			}
	 		}
	 		if (sum_w == 0.0) {

	 			//sprintf (C->format, "Warning: no data inside search radius at: %s %s\n", GMT->current.setting.format_float_out, GMT->current.setting.format_float_out);
	 			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, C->format, x0, y0);
	 			u[C->ij_sw_corner + (i * C->my + j) * C->grid] = (float)C->z_mean;
	 		}
	 		else {
	 			u[C->ij_sw_corner + (i*C->my+j)*C->grid] = (float)(sum_zw/sum_w);
	 		}
		}
	}
}


void smart_divide (struct SURFACE_INFO *C)
{	/* Divide grid by its largest prime factor */
	C->grid /= C->factors[C->n_fact - 1];
	C->n_fact--;
}

void set_offset (struct SURFACE_INFO *C)
{	/* Because of the multigrid approach the distance from the central node to
	 * its neighbors needed in the finite difference expressions varies.  E.g.,
	 * when C->grid is 8 then the next neighbor to the right is 8 columns over.
	 * But at the boundaries the spacing is always 1.  Thus, as the current node
	 * moves over the interior of the grid the distances change as we get close
	 * to any of the 4 boundaries.  The offset array is used to determine what
	 * the offset in rows and columns are relative to the current point, given
	 * what "kase" we are examining.  kase is a combination of x and y information
	 * related to how close we are to the left/right or top/bottom boundary.
	 */
	int add_w[5], add_e[5], add_s[5], add_n[5], add_w2[5], add_e2[5], add_s2[5], add_n2[5];
	unsigned int i, j, kase;

	add_w[0] = -C->my; add_w[1] = add_w[2] = add_w[3] = add_w[4] = -C->grid_east;
	add_w2[0] = -2 * C->my;  add_w2[1] = -C->my - C->grid_east;  add_w2[2] = add_w2[3] = add_w2[4] = -2 * C->grid_east;
	add_e[4] = C->my; add_e[0] = add_e[1] = add_e[2] = add_e[3] = C->grid_east;
	add_e2[4] = 2 * C->my;  add_e2[3] = C->my + C->grid_east;  add_e2[2] = add_e2[1] = add_e2[0] = 2 * C->grid_east;

	add_n[4] = 1; add_n[3] = add_n[2] = add_n[1] = add_n[0] = C->grid;
	add_n2[4] = 2;  add_n2[3] = C->grid + 1;  add_n2[2] = add_n2[1] = add_n2[0] = 2 * C->grid;
	add_s[0] = -1; add_s[1] = add_s[2] = add_s[3] = add_s[4] = -C->grid;
	add_s2[0] = -2;  add_s2[1] = -C->grid - 1;  add_s2[2] = add_s2[3] = add_s2[4] = -2 * C->grid;

	for (i = 0, kase = 0; i < 5; i++) {
		for (j = 0; j < 5; j++, kase++) {
			C->offset[kase][0] = add_n2[j];
			C->offset[kase][1] = add_n[j] + add_w[i];
			C->offset[kase][2] = add_n[j];
			C->offset[kase][3] = add_n[j] + add_e[i];
			C->offset[kase][4] = add_w2[i];
			C->offset[kase][5] = add_w[i];
			C->offset[kase][6] = add_e[i];
			C->offset[kase][7] = add_e2[i];
			C->offset[kase][8] = add_s[j] + add_w[i];
			C->offset[kase][9] = add_s[j];
			C->offset[kase][10] = add_s[j] + add_e[i];
			C->offset[kase][11] = add_s2[j];
		}
	}
}
int compare_points (const void *point_1v, const void *point_2v)
{
		/*  Routine for qsort to sort data structure for fast access to data by node location.
		    Sorts on index first, then on radius to node corresponding to index, so that index
		    goes from low to high, and so does radius.
		*/
	uint64_t block_i, block_j, index_1, index_2;
	double x0, y0, dist_1, dist_2;
	const struct SURFACE_DATA *point_1 = point_1v, *point_2 = point_2v;

	index_1 = point_1->index;
	index_2 = point_2->index;
	if (index_1 < index_2) return (-1);
	if (index_1 > index_2) return (1);
	if (index_1 == SURFACE_OUTSIDE) return (0);
	/* Points are in same grid cell, find the one who is nearest to grid point */
	block_i = point_1->index/GMT_Surface_Global.block_ny;
	block_j = point_1->index%GMT_Surface_Global.block_ny;
	x0 = GMT_Surface_Global.x_min + block_i * GMT_Surface_Global.grid_xinc;
	y0 = GMT_Surface_Global.y_min + block_j * GMT_Surface_Global.grid_yinc;
	dist_1 = (point_1->x - x0) * (point_1->x - x0) + (point_1->y - y0) * (point_1->y - y0);
	dist_2 = (point_2->x - x0) * (point_2->x - x0) + (point_2->y - y0) * (point_2->y - y0);
	if (dist_1 < dist_2) return (-1);
	if (dist_1 > dist_2) return (1);
	return (0);
}

void set_index (struct SURFACE_INFO *C)
{	/* recomputes data[k].index for new value of grid,
	   sorts data on index and radii, and throws away
	   data which are now outside the usable limits. */
	int i, j;
	uint64_t k, k_skipped = 0;
	struct GMT_GRID_HEADER *h = C->Grid->header;

	for (k = 0; k < C->npoints; k++) {
		i = irint (floor(((C->data[k].x-h->wesn[XLO])*C->r_grid_xinc) + 0.5));
		j = irint (floor(((C->data[k].y-h->wesn[YLO])*C->r_grid_yinc) + 0.5));
		if (i < 0 || i >= C->block_nx || j < 0 || j >= C->block_ny) {
			C->data[k].index = SURFACE_OUTSIDE;
			k_skipped++;
		}
		else
			C->data[k].index = i * C->block_ny + j;
	}

	qsort (C->data, C->npoints, sizeof (struct SURFACE_DATA), compare_points);

	C->npoints -= k_skipped;

}


void fill_in_forecast (struct SURFACE_INFO *C) {

	/* Fills in bilinear estimates into new node locations
	   after grid is divided.
	 */

	uint64_t index_0, index_1, index_2, index_3, index_new;
	int ii, jj, i, j;
	char *iu = C->iu;
	double delta_x, delta_y, a0, a1, a2, a3, old_size;
	float *u = C->Grid->data;

	old_size = 1.0 / (double)C->old_grid;

	/* first do from southwest corner */

	for (i = 0; i < (C->nx-1); i += C->old_grid) {

		for (j = 0; j < (C->ny-1); j += C->old_grid) {

			/* get indices of bilinear square */
			index_0 = C->ij_sw_corner + i * C->my + j;
			index_1 = index_0 + C->old_grid * C->my;
			index_2 = index_1 + C->old_grid;
			index_3 = index_0 + C->old_grid;

			/* get coefficients */
			a0 = u[index_0];
			a1 = u[index_1] - a0;
			a2 = u[index_3] - a0;
			a3 = u[index_2] - a0 - a1 - a2;

			/* find all possible new fill ins */

			for (ii = i;  ii < (i + C->old_grid); ii += C->grid) {
				delta_x = (ii - i) * old_size;
				for (jj = j;  jj < (j + C->old_grid); jj += C->grid) {
					index_new = C->ij_sw_corner + ii * C->my + jj;
					if (index_new == index_0) continue;
					delta_y = (jj - j) * old_size;
					u[index_new] = (float)(a0 + a1 * delta_x + delta_y * ( a2 + a3 * delta_x));
					iu[index_new] = 0;
				}
			}
			iu[index_0] = 5;
		}
	}

	/* now do linear guess along east edge */

	for (j = 0; j < (C->ny-1); j += C->old_grid) {
		index_0 = C->ij_se_corner + j;
		index_3 = index_0 + C->old_grid;
		for (jj = j;  jj < j + C->old_grid; jj += C->grid) {
			index_new = C->ij_se_corner + jj;
			delta_y = (jj - j) * old_size;
			u[index_new] = u[index_0] + (float)(delta_y * (u[index_3] - u[index_0]));
			iu[index_new] = 0;
		}
		iu[index_0] = 5;
	}
	/* now do linear guess along north edge */
	for (i = 0; i < (C->nx-1); i += C->old_grid) {
		index_0 = C->ij_nw_corner + i * C->my;
		index_1 = index_0 + C->old_grid * C->my;
		for (ii = i;  ii < i + C->old_grid; ii += C->grid) {
			index_new = C->ij_nw_corner + ii * C->my;
			delta_x = (ii - i) * old_size;
			u[index_new] = u[index_0] + (float)(delta_x * (u[index_1] - u[index_0]));
			iu[index_new] = 0;
		}
		iu[index_0] = 5;
	}
	/* now set northeast corner to fixed and we're done */
	iu[C->ij_ne_corner] = 5;
}

void replace_planar_trend (struct SURFACE_INFO *C)
{	/* Restore the LS plan we removed */
	int i, j;
	uint64_t ij;
	float *u = C->Grid->data;

	 for (i = 0; i < C->nx; i++) {
	 	for (j = 0; j < C->ny; j++) {
	 		ij = C->ij_sw_corner + i * C->my + j;
	 		u[ij] = (float)((u[ij] * C->z_scale) + (C->plane_c0 + C->plane_c1 * i + C->plane_c2 * j));
		}
	}
}

int write_output_surface (struct GMT_CTRL *GMT, struct SURFACE_INFO *C, char *grdfile)
{	/* Uses v.2.0 netCDF grd format - hence need to transpose original grid to be GMT compatible.  This will be rewritten, eventually */
	uint64_t index, k;
	int i, j, err;
	float *u = C->Grid->data, *v2 = NULL;
	if ((err = load_constraints (GMT, C, false))) return (err);	/* Reload constraints but this time do not transform data */

	strcpy (C->Grid->header->title, "Data gridded with continuous surface splines in tension");

	v2 = GMT_memory_aligned (GMT, NULL, C->Grid->header->size, float);
	index = C->ij_sw_corner;
	if (GMT->common.r.active) {	/* Pixel registration request. Reset limits to the original extents */
		GMT_memcpy (C->Grid->header->wesn, C->wesn_orig, 4, double);
		C->Grid->header->registration = GMT->common.r.registration;
		/* Must reduce nx,ny by 1 to exclude the extra padding for pixel grids */
		C->Grid->header->nx--;	C->nx--;
		C->Grid->header->ny--;	C->ny--;
	}
	for (i = 0; i < C->nx; i++, index += C->my) {
		for (j = 0; j < C->ny; j++) {
			k = GMT_IJP (C->Grid->header, j, i);
			v2[k] = u[index + C->ny - j - 1];
			//printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
			if (C->set_low  && !GMT_is_fnan (C->Low->data[k]) && v2[k] < C->Low->data[k]) v2[k] = C->Low->data[k];
			if (C->set_high && !GMT_is_fnan (C->High->data[k]) && v2[k] > C->High->data[k]) v2[k] = C->High->data[k];
		}
	}
	if (C->periodic) {	/* Ensure periodicity of E-W boundaries */
		for (j = 0; j < C->ny; j++) {
			k = GMT_IJP (C->Grid->header, j, 0);
			v2[k] = v2[k+C->nx-1] = (float)(0.5 * (v2[k] + v2[k+C->nx-1]));	/* Set these to the same as their average */
			//printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
		}
	}
	GMT_free_aligned (GMT, C->Grid->data);	/* Free original column-oriented grid */
	C->Grid->data = v2;			/* Hook in new scanline-oriented grid */
	if (GMT_Write_Data (GMT->parent, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, grdfile, C->Grid) != GMT_OK) {
		return (GMT->parent->error);
	}
	//if ((C->set_low  > 0 && C->set_low  < 3) && GMT_Destroy_Data (GMT->parent, &C->Low) != GMT_OK) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Failed to free C->Low\n");
	//}
	//if ((C->set_high > 0 && C->set_high < 3) && GMT_Destroy_Data (GMT->parent, &C->High) != GMT_OK) {
	//	GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Failed to free C->High\n");
	//}
	return (0);
}

/* Various functions from surface that are now used elsewhere as well */

unsigned int GMT_gcd_euclid (unsigned int a, unsigned int b)
{
	/* Returns the greatest common divisor of u and v by Euclid's method.
	 * I have experimented also with Stein's method, which involves only
	 * subtraction and left/right shifting; Euclid is faster, both for
	 * integers of size 0 - 1024 and also for random integers of a size
	 * which fits in a long integer.  Stein's algorithm might be better
	 * when the integers are HUGE, but for our purposes, Euclid is fine.
	 *
	 * Walter H. F. Smith, 25 Feb 1992, after D. E. Knuth, vol. II  */

	unsigned int u, v, r;

	u = MAX (a, b);
	v = MIN (a, b);

	while (v > 0) {
		r = u % v;	/* Knuth notes that u < 2v 40% of the time;  */
		u = v;		/* thus we could have tried a subtraction  */
		v = r;		/* followed by an if test to do r = u%v  */
	}
	return (u);
}


unsigned int GMT_get_prime_factors (/*struct GMT_CTRL *GMT,*/ uint64_t n, unsigned int *f)
{
	/* Fills the integer array f with the prime factors of n.
	 * Returns the number of locations filled in f, which is
	 * one if n is prime.
	 *
	 * f[] should have been malloc'ed to enough space before
	 * calling prime_factors().  We can be certain that f[32]
	 * is enough space, for if n fits in a long, then n < 2**32,
	 * and so it must have fewer than 32 prime factors.  I think
	 * that in general, ceil(log2((double)n)) is enough storage
	 * space for f[].
	 *
	 * Tries 2,3,5 explicitly; then alternately adds 2 or 4
	 * to the previously tried factor to obtain the next trial
	 * factor.  This is done with the variable two_four_toggle.
	 * With this method we try 7,11,13,17,19,23,25,29,31,35,...
	 * up to a maximum of sqrt(n).  This shortened list results
	 * in 1/3 fewer divisions than if we simply tried all integers
	 * between 5 and sqrt(n).  We can reduce the size of the list
	 * of trials by an additional 20% by removing the multiples
	 * of 5, which are equal to 30m +/- 5, where m >= 1.  Starting
	 * from 25, these are found by alternately adding 10 or 20.
	 * To do this, we use the variable ten_twenty_toggle.
	 *
	 * W. H. F. Smith, 26 Feb 1992, after D.E. Knuth, vol. II  */

	unsigned int current_factor = 0;	/* The factor currently being tried  */
	unsigned int max_factor;		/* Don't try any factors bigger than this  */
	unsigned int n_factors = 0;		/* Returned; one if n is prime  */
	unsigned int two_four_toggle = 0;	/* Used to add 2 or 4 to get next trial factor  */
	unsigned int ten_twenty_toggle = 0;	/* Used to add 10 or 20 to skip_five  */
	unsigned int skip_five = 25;	/* Used to skip multiples of 5 in the list  */
	unsigned int base_factor[3] = {2, 3, 5};	/* Standard factors to try */
	uint64_t m = n;			/* Used to keep a working copy of n  */
	unsigned int k;			/* counter */

	/* Initialize max_factor  */

	if (m < 2) return (0);
	max_factor = urint (floor(sqrt((double)m)));

	/* First find the 2s, 3s, and 5s */
	for (k = 0; k < 3; k++) {
		current_factor = base_factor[k];
		while (!(m % current_factor)) {
			m /= current_factor;
			f[n_factors++] = current_factor;
		}
		if (m == 1) return (n_factors);
	}

	/* Unless we have already returned we now try all the rest  */

	while (m > 1 && current_factor <= max_factor) {

		/* Current factor is either 2 or 4 more than previous value  */

		if (two_four_toggle) {
			current_factor += 4;
			two_four_toggle = 0;
		}
		else {
			current_factor += 2;
			two_four_toggle = 1;
		}

		/* If current factor is a multiple of 5, skip it.  But first,
			set next value of skip_five according to 10/20 toggle:  */

		if (current_factor == skip_five) {
			if (ten_twenty_toggle) {
				skip_five += 20;
				ten_twenty_toggle = 0;
			}
			else {
				skip_five += 10;
				ten_twenty_toggle = 1;
			}
			continue;
		}

		/* Get here when current_factor is not a multiple of 2,3 or 5:  */

		while (!(m % current_factor)) {
			m /= current_factor;
			f[n_factors++] = current_factor;
		}
	}

	/* Get here when all factors up to floor(sqrt(n)) have been tried.  */

	if (m > 1) f[n_factors++] = (unsigned int)m;	/* m is an additional prime factor of n  */

	return (n_factors);
}

void load_parameters_surface (struct SURFACE_INFO *C, struct SURFACE_CTRL *Ctrl)
{	/* Place program options into the surface struct.  This was done this way
	 * since surface.c relied heavily on global variables which are a no-no
	 * in GMT5.  The simplest solution was to collect all those variables into
	 * a single structure and pass a pointer to that structure to functions.
	 */
	if (Ctrl->S.active) {
		if (Ctrl->S.unit == 'm') Ctrl->S.radius /= 60.0;
		if (Ctrl->S.unit == 's') Ctrl->S.radius /= 3600.0;
	}
	C->radius = Ctrl->S.radius;
	C->relax_new = Ctrl->Z.value;
	C->max_iterations = Ctrl->N.value;
	C->radius = Ctrl->S.radius;
	C->low_file = Ctrl->L.low;
	C->high_file = Ctrl->L.high;
	C->set_low = Ctrl->L.lmode;
	C->low_limit = Ctrl->L.min;
	C->set_high = Ctrl->L.hmode;
	C->high_limit = Ctrl->L.max;
	C->boundary_tension = Ctrl->T.b_tension;
	C->interior_tension = Ctrl->T.i_tension;
	C->l_epsilon = Ctrl->A.value;
	C->converge_limit = Ctrl->C.value;
}


int GMT_Complete_Options (struct GMT_CTRL *GMT, struct GMT_OPTION *options)
{
	/* Go through the given arguments and look for shorthands,
	 * i.e., -B, -J, -R, -X, -x, -Y, -y, -c, -p. given without arguments.
	 * If found, see if we have a matching command line history and then
	 * update that entry in the option list.
	 * Finally, keep the option arguments in the history list.
	 * However, when func_level > 1, do update the entry, but do not
	 * remember it in history. Note, there are two special cases here:
	 * -J is special since we also need to deal with the sub-species
	 *    like -JM, -JX, etc.  So these all have separate entries.
	 * -B is special because the option is repeatable for different
	 *    aspects of the basemap.  We concatenate all of them to store
	 *    in the history file and use RS = ASCII 30 as separator.
	 *    Also, a single -B in the options may expand to several
	 *    separate -B<args> so the linked options list may grow.
	 */

	int id = 0, k, n_B = 0, B_id;
	unsigned int pos = 0, B_replace = 1;
	bool update, remember;
	struct GMT_OPTION *opt = NULL, *opt2 = NULL, *B_next = NULL;
	char str[3] = {""}, B_string[BUFFSIZE_NEW] = {""}, p[BUFFSIZE_NEW] = {""}, B_delim[2] = {30, 0};	/* Use ASCII 30 RS Record Separator between -B strings */

	remember = (GMT->hidden.func_level == 1);	/* Only update the history for top level function */

	for (opt = options; opt; opt = opt->next) if (opt->option == 'B') {	/* Do some initial counting of how many -B options and determine if there is just one with no args */
		if (n_B > 0 || opt->arg[0]) B_replace = 0;
		n_B++;
	}
	for (k = 0, B_id = -1; k < GMT_N_UNIQUE && B_id == -1; k++) if (!strcmp (GMT_unique_option[k], "B")) B_id = k;	/* B_id === 0 but just in case this changes we do this search anyway */

	for (opt = options; opt; opt = opt->next) {
		if (!strchr (GMT_SHORTHAND_OPTIONS, opt->option)) continue;	/* Not one of the shorthand options */
		update = false;
		//GMT_Report (GMT->parent, GMT_MSG_DEBUG, "History: Process -%c%s.\n", opt->option, opt->arg);

		str[0] = opt->option; str[1] = str[2] = '\0';
		if (opt->option == 'J') {	/* -J is special since it can be -J or -J<code> */
			/* Always look up "J" first. It comes before "J?" and tells what the last -J was */
			for (k = 0, id = -1; k < GMT_N_UNIQUE && id == -1; k++) if (!strcmp (GMT_unique_option[k], str)) id = k;
			if (id < 0) return;
			if (opt->arg && opt->arg[0]) {	/* Gave -J<code>[<args>] so we either use or update history and continue */
				str[1] = opt->arg[0];
				/* Remember this last -J<code> for later use as -J, but do not remember it when -Jz|Z */
				if (str[1] != 'Z' && str[1] != 'z' && remember) {
					if (GMT->init.history[id]) free (GMT->init.history[id]);
					GMT->init.history[id] = strdup (&str[1]);
				}
				if (opt->arg[1]) update = true;	/* Gave -J<code><args> so we want to update history and continue */
			}
			else {
				if (!GMT->init.history[id]) return;
				str[1] = GMT->init.history[id][0];
			}
			/* Continue looking for -J<code> */
			for (k = id + 1, id = -1; k < GMT_N_UNIQUE && id == -1; k++) if (!strcmp (GMT_unique_option[k], str)) id = k;
			if (id < 0) return;
		}
		else if (opt->option == 'B') {	/* -B is also special since there may be many of these, or just -B */
			if (B_replace) {	/* Only -B is given and we want to use the history */
				if (B_replace == 2) continue;	/* Already done this */
				if (!GMT->init.history[B_id]) return;
				opt2 = opt;			/* Since we dont want to change the opt loop avove */
				B_next = opt->next;		/* Pointer to option following the -B option */
				if (opt2->arg) free (opt2->arg);	/* Free previous pointer to arg */
				GMT_strtok (GMT->init.history[B_id], B_delim, &pos, p);	/* Get the first argument */
				opt2->arg = strdup (p);		/* Update arg */
				while (GMT_strtok (GMT->init.history[B_id], B_delim, &pos, p)) {	/* Parse any additional |<component> statements */
					opt2->next = GMT_Make_Option (GMT->parent, 'B', p);	/* Create new struct */
					opt2->next->previous = opt2;
					opt2 = opt2->next;
				}
				opt2->next = B_next;	/* Hook back onto main option list */
				B_replace = 2;	/* Flag to let us know we are done with -B */
			}
			else {	/* One of possibly several -B<arg> options; concatenate and separate by RS */
				if (B_string[0]) strcat (B_string, B_delim);	/* Add RS separator between args */
				strcat (B_string, opt->arg);
			}
		}
		else {	/* Gave -R[<args>], -V[<args>] etc., so we either use or update the history and continue */
			for (k = 0, id = -1; k < GMT_N_UNIQUE && id == -1; k++) if (!strcmp (GMT_unique_option[k], str)) id = k;	/* Find entry in history array */
			if (id < 0) return;	/* Error: user gave shorthand option but there is no record in the history */
			if (opt->arg && opt->arg[0]) update = true;	/* Gave -R<args>, -V<args> etc. so we we want to update history and continue */
		}
		if (opt->option != 'B') {	/* Do -B separately again after the loop so skip it here */
			if (update) {	/* Gave -J<code><args>, -R<args>, -V<args> etc. so we update history and continue */
				if (remember) {
					if (GMT->init.history[id]) free (GMT->init.history[id]);
					GMT->init.history[id] = strdup (opt->arg);
				}
			}
			else {	/* Gave -J<code>, -R, -J etc. so we complete the option and continue */
				if (!GMT->init.history[id]) return;
				if (opt->arg) free (opt->arg);	/* Free previous pointer to arg */
				opt->arg = strdup (GMT->init.history[id]);
			}
		}
	}

	if (B_string[0]) {	/* Got a concatenated string with one or more individual -B args, now separated by the RS character (ascii 30) */
		if (GMT->init.history[B_id]) free (GMT->init.history[B_id]);
		GMT->init.history[B_id] = strdup (B_string);
	}

	return (GMT_NOERROR);
}

int gmt_compare_cols (const void *point_1, const void *point_2)
{
	/* Sorts cols into ascending order  */
	if (((struct GMT_COL_INFO *)point_1)->col < ((struct GMT_COL_INFO *)point_2)->col) return (-1);
	if (((struct GMT_COL_INFO *)point_1)->col > ((struct GMT_COL_INFO *)point_2)->col) return (+1);
	return (0);
}

int64_t gmt_parse_range (struct GMT_CTRL *GMT, char *p, int64_t *start, int64_t *stop)
{	/* Parses p looking for range or columns or individual columns.
	 * If neither then we just increment both start and stop. */
	int64_t inc = 1;
	int got;
	char *c = NULL;
	if ((c = strchr (p, '-'))) {	/* Range of columns given. e.g., 7-9 */
		got = sscanf (p, "%" PRIu64 "-%" PRIu64, start, stop);
		if (got != 2) inc = 0L;	/* Error flag */
	}
	else if ((c = strchr (p, ':'))) {	/* Range of columns given. e.g., 7:9 or 1:2:5 */
		got = sscanf (p, "%" PRIu64 ":%" PRIu64 ":%" PRIu64, start, &inc, stop);
		if (got == 2) { *stop = inc; inc = 1L;}	/* Just got start:stop with implicit inc = 1 */
		else if (got != 3 || inc < 1) inc = 0L;	/* Error flag */
	}
	else if (isdigit ((int)p[0]))	/* Just a single column, e.g., 3 */
		*start = *stop = atol (p);
	else				/* Just assume it goes column by column */
		(*start)++, (*stop)++;
	if ((*stop) < (*start)) inc = 0L;	/* Not good */
	if (inc == 0)
		printf("inc == 0 \n");
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Bad range [%s]: col, start-stop, start:stop, or start:step:stop must yield monotonically increasing positive selections\n", p);
	return (inc);	/* Either > 0 or 0 for error */
}


int gmt_parse_i_option (struct GMT_CTRL *GMT, char *arg)
{
	/* Routine will decode the -i<col>|<colrange>[l][s<scale>][o<offset>],... arguments */

	char copy[BUFFSIZE_NEW] = {""}, p[BUFFSIZE_NEW] = {""}, *c = NULL;
	char txt_a[GMT_LEN256] = {""}, txt_b[GMT_LEN256] = {""};
	unsigned int k = 0, pos = 0;
	int64_t i, start = -1, stop = -1, inc;
	int convert;
	double scale, offset;

	if (!arg || !arg[0]) return (GMT_PARSE_ERROR);	/* -i requires an argument */

	strncpy (copy, arg, BUFFSIZE_NEW);
	for (i = 0; i < GMT_MAX_COLUMNS; i++) GMT->current.io.col_skip[i] = true;	/* Initially, no input column is requested */

	while ((GMT_strtok (copy, ",", &pos, p))) {	/* While it is not empty, process it */
		convert = 0, scale = 1.0, offset = 0.0;

		if ((c = strchr (p, 'o'))) {	/* Look for offset */
			c[0] = '\0';	/* Wipe out the 'o' so that next scan terminates there */
			convert |= 1;
			offset = atof (&c[1]);
		}
		if ((c = strchr (p, 's'))) {	/* Look for scale factor */
			c[0] = '\0';	/* Wipe out the 's' so that next scan terminates there */
			if (GMT_compat_check (GMT, 4)) {	/* GMT4 */
				i = (int)strlen (p) - 1;
				convert = (p[i] == 'l') ? 2 : 1;
				i = sscanf (&c[1], "%[^/]/%[^l]", txt_a, txt_b);
				//if (i == 0) GMT_Report (GMT->parent, GMT_MSG_NORMAL, "-i...s contains bad scale info\n");
				scale = atof (txt_a);
				if (i == 2) offset = atof (txt_b);
			} else {
				convert |= 1;
				scale = atof (&c[1]);
			}
		}
		if ((c = strchr (p, 'l'))) {	/* Look for log indicator */
			c[0] = '\0';	/* Wipe out the 's' so that next scan terminates there */
			convert = 2;
		}
		if ((inc = gmt_parse_range (GMT, p, &start, &stop)) == 0) return (GMT_PARSE_ERROR);

		/* Now set the code for these columns */

		for (i = start; i <= stop; i += inc, k++) {
			GMT->current.io.col_skip[i] = false;	/* Do not skip these */
			GMT->current.io.col[GMT_IN][k].col = (unsigned int)i;	/* Requested order of columns */
			GMT->current.io.col[GMT_IN][k].order = k;		/* Requested order of columns */
			GMT->current.io.col[GMT_IN][k].convert = convert;
			GMT->current.io.col[GMT_IN][k].scale = scale;
			GMT->current.io.col[GMT_IN][k].offset = offset;
		}
	}
	qsort (GMT->current.io.col[GMT_IN], k, sizeof (struct GMT_COL_INFO), gmt_compare_cols);
	GMT->common.i.n_cols = k;

	return (GMT_NOERROR);
}

int gmt_parse_colon_option (struct GMT_CTRL *GMT, char *item) {
	int error = 0, way, off = 0;
	bool ok[2] = {false, false};
	static char *mode[4] = {"i", "o", "", ""}, *dir[2] = {"input", "output"};
	char kase = (item) ? item[0] : '\0';
	/* Parse the -: option.  Full syntax: -:[i|o].
	 * We know that if -f was given it has already been parsed due to the parsing order imposed.
	 * Must check that -: does not conflict with -f */

	switch (kase) {
		case 'i':	/* Toggle on input data only */
			ok[GMT_IN] = true;
			break;
		case 'o':	/* Toggle on output data only */
			ok[GMT_OUT] = true;
			break;
		case '\0':	/* Toggle both input and output data */
			ok[GMT_IN] = ok[GMT_OUT] = true;
			off = 2;
			break;
		default:
			error++;	/* Error */
			break;
	}
	for (way = 0; !error && way < 2; way++) if (ok[way]) {
		if (GMT->current.io.col_type[way][GMT_X] == GMT_IS_UNKNOWN && GMT->current.io.col_type[way][GMT_Y] == GMT_IS_UNKNOWN)	/* Dont know what x/y is yet */
			GMT->current.setting.io_lonlat_toggle[way] = true;
		else if (GMT->current.io.col_type[way][GMT_X] == GMT_IS_FLOAT && GMT->current.io.col_type[way][GMT_Y] == GMT_IS_FLOAT)	/* Cartesian x/y vs y/x cannot be identified */
			GMT->current.setting.io_lonlat_toggle[way] = true;
		else if (GMT_is_geographic (GMT, way))	/* Lon/lat becomes lat/lon */
			GMT->current.setting.io_lonlat_toggle[way] = true;
		//else if (GMT->current.io.col_type[way][GMT_X] == GMT_IS_LAT && GMT->current.io.col_type[way][GMT_Y] == GMT_IS_LON)	/* Already lat/lon! */
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: -:%s given but %s order already set by -f; -:%s ignored.\n", mode[way+off], dir[way], mode[way+off]);
		else {
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Error: -:%s given but %s first two columns do not hold x/y or lon/lat\n", mode[way+off], dir[way]);
			error++;
		}
	}
	if (error) GMT->current.setting.io_lonlat_toggle[GMT_IN] = GMT->current.setting.io_lonlat_toggle[GMT_OUT] = false;	/* Leave in case we had errors */
	return (error);
}

int gmt_parse_f_option (struct GMT_CTRL *GMT, char *arg)
{
	/* Routine will decode the -f[i|o]<col>|<colrange>[t|T|g],... arguments */

	char copy[BUFFSIZE_NEW] = {""}, p[BUFFSIZE_NEW] = {""};
	unsigned int k = 1, ic, pos = 0, code, *col = NULL;
	size_t len;
	int64_t i, start = -1, stop = -1, inc;
	enum GMT_enum_units unit = GMT_IS_METER;
	bool both_i_and_o = false;

	if (!arg || !arg[0]) return (GMT_PARSE_ERROR);	/* -f requires an argument */

	if (arg[0] == 'i')	/* Apply to input columns only */
		col = GMT->current.io.col_type[GMT_IN];
	else if (arg[0] == 'o')	/* Apply to output columns only */
		col = GMT->current.io.col_type[GMT_OUT];
	else {			/* Apply to both input and output columns */
		both_i_and_o = true;
		k = 0;
	}

	strncpy (copy, &arg[k], BUFFSIZE_NEW);	/* arg should NOT have a leading i|o part */



	if (copy[0] == 'g' || copy[0] == 'p') {	/* Got -f[i|o]g which is shorthand for -f[i|o]0x,1y, or -fp[<unit>] (see below) */
		if (both_i_and_o) {
			//printf(">>>>>>>>>>file : %s line : %d func: %s \n",__FILE__,__LINE__,__func__);
			GMT_set_geographic (GMT, GMT_IN);
			GMT_set_geographic (GMT, GMT_OUT);
			//printf(">>>>>>>>>>file : %s line : %d func: %s \n",__FILE__,__LINE__,__func__);
		}
		else {
			col[GMT_X] = GMT_IS_LON;
			col[GMT_Y] = GMT_IS_LAT;
		}
		pos = 1;
		start = stop = 1;
	}
	if (copy[0] == 'p') {	/* Got -f[i|o]p[<unit>] for projected floating point map coordinates (e.g., UTM meters) */
		if (copy[1] && strchr (GMT_LEN_UNITS2, copy[1])) {	/* Given a unit via -fp<unit>*/
			if ((unit = GMT_get_unit_number (GMT, copy[1])) == GMT_IS_NOUNIT) {
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Error: Malformed -f argument [%s] - bad projected unit\n", arg);
				return 1;
			}
			pos++;
		}
		GMT->current.proj.inv_coordinates = true;
		GMT->current.proj.inv_coord_unit = unit;
	}



	while ((GMT_strtok (copy, ",", &pos, p))) {	/* While it is not empty, process it */

		if ((inc = gmt_parse_range (GMT, p, &start, &stop)) == 0) return (GMT_PARSE_ERROR);
		//printf(">>>>>>>>>>file : %s line : %d func: %s \n",__FILE__,__LINE__,__func__);
		len = strlen (p);	/* Length of the string p */
		ic = (int) p[len-1];	/* Last char in p is the potential code T, t, x, y, or f. */
		switch (ic) {
			case 'T':	/* Absolute calendar time */
				code = GMT_IS_ABSTIME;
				break;
			case 't':	/* Relative time (units since epoch) */
				code = GMT_IS_RELTIME;
				break;
			case 'x':	/* Longitude coordinates */
				code = GMT_IS_LON;
				break;
			case 'y':	/* Latitude coordinates */
				code = GMT_IS_LAT;
				break;
			case 'f':	/* Plain floating point coordinates */
				code = GMT_IS_FLOAT;
				break;
			default:	/* No suffix, consider it an error */
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Error: Malformed -f argument [%s]\n", arg);
				return 1;
				break;
		}

		/* Now set the code for these columns */

		if (both_i_and_o)
			for (i = start; i <= stop; i += inc) GMT->current.io.col_type[GMT_IN][i] = GMT->current.io.col_type[GMT_OUT][i] = code;
		else
			for (i = start; i <= stop; i += inc) col[i] = code;
	}
	//printf(">>>>>>>>>>file : %s line : %d func: %s \n",__FILE__,__LINE__,__func__);
	return (GMT_NOERROR);
}

int GMT_parse_common_options (struct GMT_CTRL *GMT, char *list, char option, char *item)
{
	/* GMT_parse_common_options interprets the command line for the common, unique options
	 * -B, -J, -K, -O, -P, -R, -U, -V, -X, -Y, -b, -c, -f, -g, -h, -i, -n, -o, -p, -r, -s, -t, -:, -- and -^.
	 * The list passes all of these that we should consider.
	 * The API will also consider -I for grid increments.
	 */

	int error = 0, i = 0;	/* The i and i+= GMT_more_than_once are there to avoid compiler warnings... */
	//printf("=====> option %c file : %s line : %d func: %s\n",option,__FILE__,__LINE__,__func__);
	if (!list || !strchr (list, option)) return (0);	/* Not a common option we accept */

	//if (GMT_compat_check (GMT, 4)) {
		/* Translate some GMT4 options */
		//switch (option) {
			//case 'E': GMT_COMPAT_OPT ('p'); break;
			//case 'F': GMT_COMPAT_OPT ('r'); break;
			//case 'H': GMT_COMPAT_OPT ('h'); break;
		//}
//	}
	switch (option) {	/* Handle parsing of this option, if allowed here */

		//case 'B':
			//switch (item[0]) {	/* Check for -B[p] and -Bs */
				//case 's': GMT->common.B.active[1] = true; break;
				//default:  GMT->common.B.active[0] = true; break;
			//}
			//if (!error) error = gmt_parse_B_option (GMT, item);
			//break;

		case 'I':
			if (GMT->hidden.func_level > 0) return (0);	/* Just skip if we are inside a GMT module. -I is an API common option only */

			if (GMT_getinc (GMT, item, GMT->common.API_I.inc)) {
				//GMT_inc_syntax (GMT, 'I', 1);
				error++;
			}
			GMT->common.API_I.active = true;
			break;

		case 'J':
			//if (item && (item[0] == 'Z' || item[0] == 'z')) {	/* -JZ or -Jz */
				//error += (GMT_check_condition (GMT, GMT->common.J.zactive, "Warning: Option -JZ|z given more than once\n") || gmt_parse_J_option (GMT, item));
				//GMT->common.J.zactive = true;
			//}
			//else {	/* Horizontal map projection */
				//error += (GMT_check_condition (GMT, GMT->common.J.active, "Warning: Option -J given more than once\n") || gmt_parse_J_option (GMT, item));
				//GMT->common.J.active = true;
			//}
			break;

		case 'K':
			//i += GMT_more_than_once (GMT, GMT->common.K.active);
			//GMT->common.K.active = true;
			break;

		case 'O':
			//i += GMT_more_than_once (GMT, GMT->common.O.active);
			//GMT->common.O.active = true;
			break;

		case 'P':
			//i += GMT_more_than_once (GMT, GMT->common.P.active);
			//GMT->common.P.active = true;
			break;

		case 'Q':
		case 'S':
			/*if (GMT_compat_check (GMT, 4)) {
				GMT_Report (GMT->parent, GMT_MSG_COMPAT, "Warning: Option -%c is deprecated. Use -n instead.\n" GMT_COMPAT_INFO, option);
				error += backwards_SQ_parsing (GMT, option, item);
			}
			else {
				GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Option -%c is not a recognized common option\n", option);
				return (1);
			}
			break;*/

		case 'R':
			//printf("\nDDDDDDDDDDDDDDDDDDDDDDDDDDDfile : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
			error += (GMT_more_than_once (GMT, GMT->common.R.active) || gmt_parse_R_option (GMT, item));
			GMT->common.R.active = true;
			break;

		case 'U':
			/*error += (GMT_more_than_once (GMT, GMT->common.U.active) || gmt_parse_U_option (GMT, item));
			GMT->common.U.active = true;
			break;
			*/
		case 'V':
			//i += GMT_more_than_once (GMT, GMT->common.V.active);
			//GMT->common.V.active = true;
			//if (item && item[0]) {	/* Specified a verbosity level */
				//if (gmt_parse_V_option (GMT, item[0])) {
					//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Error: Unknown argument to -V option, -V%c\n", item[0]);
					//error++;
				//}
			//}
			//else
				//GMT->current.setting.verbose = GMT_MSG_VERBOSE;
			//break;

		case 'X':
			/*error += (GMT_more_than_once (GMT, GMT->common.X.active) || gmt_parse_XY_option (GMT, GMT_X, item));
			GMT->common.X.active = true;
			break;*/

		case 'Y':
			/*error += (GMT_more_than_once (GMT, GMT->common.Y.active) || gmt_parse_XY_option (GMT, GMT_Y, item));
			GMT->common.Y.active = true;
			break;*/

		case 'Z':	/* GMT4 Backwards compatibility */
			//if (GMT_compat_check (GMT, 4)) {
				//GMT_Report (GMT->parent, GMT_MSG_COMPAT, "Warning: Option -Z[<zlevel>] is deprecated. Use -p<azim>/<elev>[/<zlevel>] instead.\n" GMT_COMPAT_INFO);
				//if (item && item[0]) GMT->current.proj.z_level = atof (item);
			//}
			//else {
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Option -%c is not a recognized common option\n", option);
				//return (1);
			//}
			//break;

		case 'a':
			/*error += (GMT_more_than_once (GMT, GMT->common.a.active) || gmt_parse_a_option (GMT, item));
			GMT->common.a.active = true;
			break;*/

		case 'b':
			/*switch (item[0]) {
				case 'i':
					error += GMT_check_condition (GMT, GMT->common.b.active[GMT_IN], "Warning Option -bi given more than once\n");
					GMT->common.b.active[GMT_IN] = true;
					break;
				case 'o':
					error += GMT_check_condition (GMT, GMT->common.b.active[GMT_OUT], "Warning Option -bo given more than once\n");
					GMT->common.b.active[GMT_OUT] = true;
					break;
				default:
					error += GMT_check_condition (GMT, GMT->common.b.active[GMT_IN] + GMT->common.b.active[GMT_OUT], "Warning Option -b given more than once\n");
					GMT->common.b.active[GMT_IN] = GMT->common.b.active[GMT_OUT] = true;
					break;
			}
			error += gmt_parse_b_option (GMT, item);
			break;*/

		case 'c':
			/*error += (GMT_more_than_once (GMT, GMT->common.c.active) || gmt_parse_c_option (GMT, item));
			GMT->common.c.active = true;
			break;*/

		case 'f':
			switch (item[0]) {
				case 'i':
					//error += GMT_check_condition (GMT, GMT->common.f.active[GMT_IN], "Warning Option -fi given more than once\n");
					GMT->common.f.active[GMT_IN] = true;
					break;
				case 'o':
					//error += GMT_check_condition (GMT, GMT->common.f.active[GMT_OUT], "Warning Option -fo given more than once\n");
					GMT->common.f.active[GMT_OUT] = true;
					break;
				default:
					//error += GMT_check_condition (GMT, GMT->common.f.active[GMT_IN] | GMT->common.f.active[GMT_OUT], "Warning Option -f given more than once\n");
					GMT->common.f.active[GMT_IN] = GMT->common.f.active[GMT_OUT] = true;
					break;
			}
			error += gmt_parse_f_option (GMT, item);
			break;
		case 'g':
			/*error += gmt_parse_g_option (GMT, item);
			GMT->common.g.active = true;
			break;*/

		case 'h':
			/*error += (GMT_more_than_once (GMT, GMT->common.h.active) || gmt_parse_h_option (GMT, item));
			GMT->common.h.active = true;
			break;*/

		case 'i':

			error += (GMT_more_than_once (GMT, GMT->common.i.active) || gmt_parse_i_option (GMT, item));
			GMT->common.i.active = true;
			break;

		case 'M':	/* Backwards compatibility */
		case 'm':
			/*if (GMT_compat_check (GMT, 4)) {
				GMT_Report (GMT->parent, GMT_MSG_COMPAT, "Warning: Option -%c is deprecated. Segment headers are automatically identified.\n", option);
			}
			else {
				GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Option -%c is not a recognized common option\n", option);
				return (1);
			}
			break;*/

		case 'n':
			/*error += (GMT_more_than_once (GMT, GMT->common.n.active) || gmt_parse_n_option (GMT, item));
			GMT->common.n.active = true;
			break;*/

		case 'o':
			/*error += (GMT_more_than_once (GMT, GMT->common.o.active) || gmt_parse_o_option (GMT, item));
			GMT->common.o.active = true;
			break;*/

		case 'p':
			/*error += (GMT_more_than_once (GMT, GMT->common.p.active) || gmt_parse_p_option (GMT, item));
			GMT->common.p.active = true;
			break;*/

		case 'r':
			//if (GMT->current.io.grd_info.active) GMT->common.r.active = false;	/* OK to override registration given via -Rfile */
			//error += GMT_more_than_once (GMT, GMT->common.r.active);
			//GMT->common.r.active = true;
			//GMT->common.r.registration = GMT_GRID_PIXEL_REG;
			//break;

		case 's':
			//error += (GMT_more_than_once (GMT, GMT->common.s.active) || gmt_parse_s_option (GMT, item));
			//GMT->common.s.active = true;
			//break;

		case 't':
			//error += GMT_more_than_once (GMT, GMT->common.t.active);
			//if (item[0]) {
				//GMT->common.t.active = true;
				//GMT->common.t.value = atof (item);
			//}
			//else {
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Option -t was not given any value (please add transparency in 0-100%% range)!\n");
				//error++;
			//}
			//break;

		case ':':
		//	printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
			error += (GMT_more_than_once (GMT, GMT->common.colon.active) || gmt_parse_colon_option (GMT, item));
			GMT->common.colon.active = true;
			break;

		case '^':
			//if (GMT->common.synopsis.active)
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: Option - given more than once\n");
			GMT->common.synopsis.active = true;
			break;

		case '-':
			//error += GMT_parse_dash_option (GMT, item);
			break;

		case '>':	/* Registered output file; nothing to do here */
			//printf("printing optionssssssssss>     ........\n");
			break;

		default:	/* Here we end up if an unrecognized option is passed (should not happen, though) */
		//	printf ("Option -%c is not a recognized common option\n", option);
			return (1);
			break;
	}

	/* On error, give syntax message */

	//if (error) GMT_syntax (GMT, option);
	return (error);
}

int GMT_Parse_Common (void *V_API, char *given_options, struct GMT_OPTION *options)
{
	/* GMT_Parse_Common parses the option list for a program and detects the GMT common options.
	 * These are processed in the order required by GMT regardless of order given.
	 * The settings will override values set previously by other commands.
	 * It ignores filenames and only return errors if GMT common options are misused.
	 * Note: GMT_CRITICAL_OPT_ORDER is defined in gmt_common.h
	 */

	struct GMT_OPTION *opt = NULL;
	char list[2] = {0, 0}, critical_opt_order[] = GMT_CRITICAL_OPT_ORDER;
	unsigned int i, n_errors = 0;
	struct GMTAPI_CTRL *API = NULL;
	//printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
	if (V_API == NULL)
		{
		printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
		return_error (V_API, GMT_NOT_A_SESSION);	/* GMT_Create_Session has not been called */
		}
	/* Check if there are short-hand commands present (e.g., -J with no arguments); if so complete these to full options
	 * by consulting the current GMT history machinery.  If not possible then we have an error to report */

	API = gmt_get_api_ptr (V_API);	/* Cast void pointer to a GMTAPI_CTRL pointer *//* Cast void pointer to a GMTAPI_CTRL pointer */
	//printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
	if(API->GMT==NULL)
	{
		printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
		return ;
	}

	//fflush(stdout);
	if (GMT_Complete_Options (API->GMT, options))
		{
			printf("checking error");
			return_error (API, GMT_OPTION_HISTORY_ERROR);	/* Replace shorthand failed */
		}
	//if (API->GMT->common.B.mode == 0) API->GMT->common.B.mode = gmt_check_b_options (API->GMT, options);	/* Determine the syntax of the -B option(s) */
	//printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
	/* First parse the common options in the order they appear in GMT_CRITICAL_OPT_ORDER */
	for (i = 0; i < strlen (critical_opt_order); i++) {	/* These are the GMT options that must be parsed in this particular order, if present */
		if (strchr (given_options, critical_opt_order[i]) == NULL) continue;	/* Not a selected option */
		list[0] = critical_opt_order[i];	/* Just look for this particular option in the linked opt list */
		for (opt = options; opt; opt = opt->next) {
			if (opt->option != list[0]) continue;

			n_errors += GMT_parse_common_options (API->GMT, list, opt->option, opt->arg);
		}
	}
	//printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
	/* Now any critical option given has been parsed.  Next we parse anything NOT a critical GMT option */
	for (i = 0; given_options[i]; i++) {	/* Loop over all options given */
		if (strchr (critical_opt_order, given_options[i])) continue;	/* Skip critical option */
		list[0] = given_options[i];	/* Just look for this particular option */
		for (opt = options; opt; opt = opt->next) {
			//printf(">>>>>>>>>>file : %s line : %d func: %s opt->arg: %s\n",__FILE__,__LINE__,__func__,opt->arg);
			n_errors += GMT_parse_common_options (API->GMT, list, opt->option, opt->arg);
			if (opt->option != list[0]) continue;
		}
	}

	/* Update [read-only] pointer to the current option list */
	API->GMT->current.options = options;
	if (n_errors) return_error (API, GMT_PARSE_ERROR);	/* One or more options failed to parse */
	return (GMT_OK);
}

bool GMT_check_filearg (struct GMT_CTRL *GMT, char option, char *file, unsigned int direction)
{	/* Return true if a file arg was given and, if direction is GMT_IN, check that the file
	 * exists and is readable. Otherwise wre return false. */
	unsigned int k = 0;
	char message[GMT_LEN16] = {""};
	if (option == GMT_OPT_INFILE)
		sprintf (message, "for input file");
	else if (option == GMT_OPT_OUTFILE)
		sprintf (message, "for output file");
	else
		sprintf (message, "option -%c", option);

	if (!file || file[0] == '\0') {
		printf ( "Error %s: No filename provided\n");
		return false;	/* No file given */
	}
	if (direction == GMT_OUT) return true;		/* Cannot check any further */
	if (file[0] == '=') k = 1;	/* Gave a list of files with =<filelist> mechanism in x2sys */





	if (GMT_access (GMT, &file[k], F_OK)) {	/* Cannot find the file anywhere GMT looks */
		printf ( "Error %s: No such file (%s)\n", message, &file[k]);
		return false;	/* Could not find this file */
	}
	if (GMT_access (GMT, &file[k], R_OK)) {	/* Cannot read this file (permissions?) */
		printf ( "Error %s: Cannot read file (%s) - check permissions\n", message, &file[k]);
		return false;	/* Could not find this file */
	}
	return true;	/* Seems OK */
}

int GMT_check_binary_io (struct GMT_CTRL *GMT, uint64_t n_req) {
	int n_errors = 0;

	/* Check the binary options that are used with most GMT programs.
	 * GMT is the pointer to the GMT structure.
	 * n_req is the number of required columns. If 0 then it relies on
	 *    GMT->common.b.ncol[GMT_IN] to be non-zero.
	 * Return value is the number of errors that where found.
	 */

	if (!GMT->common.b.active[GMT_IN]) return (GMT_NOERROR);	/* Let machinery figure out input cols for ascii */

	/* These are specific tests for binary input */

	if (GMT->common.b.ncol[GMT_IN] == 0) GMT->common.b.ncol[GMT_IN] = n_req;
	if (GMT->common.b.ncol[GMT_IN] == 0) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Syntax error: Must specify number of columns in binary input data (-bi)\n");
		n_errors++;
	}
	else if (n_req > GMT->common.b.ncol[GMT_IN]) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Syntax error: Binary input data (-bi) provides %d but must have at least %d columns\n", GMT->common.b.ncol[GMT_IN], n_req);
		n_errors++;
	}
	if (GMT->common.b.ncol[GMT_IN] < GMT->common.i.n_cols) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Syntax error: Binary input data (-bi) provides %d but column selection (-i) asks for %d columns\n", GMT->common.b.ncol[GMT_IN], GMT->common.i.n_cols);
		n_errors++;
	}
	if (GMT->common.b.active[GMT_OUT] && GMT->common.b.ncol[GMT_OUT] && (GMT->common.b.ncol[GMT_OUT] < GMT->common.o.n_cols)) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Syntax error: Binary output data (-bo) provides %d but column selection (-o) asks for %d columns\n", GMT->common.b.ncol[GMT_OUT], GMT->common.o.n_cols);
		n_errors++;
	}

	//GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "Provides %d, expects %d-column binary data\n", GMT->common.b.ncol[GMT_IN], n_req);

	return (n_errors);
}

int GMT_surface_parse (struct GMT_CTRL *GMT, struct SURFACE_CTRL *Ctrl, struct GMT_OPTION *options)
{
	/* This parses the options provided to surface and sets parameters in CTRL.
	 * Any GMT common options will override values set previously by other commands.
	 * It also replaces any file names specified as input or output with the data ID
	 * returned when registering these sources/destinations with the API.
	 */
	//printf("file : %s line : %d func: %s opt->option :%s\n",__FILE__,__LINE__,__func__);
	unsigned int n_errors = 0, k;
	char modifier;
	struct GMT_OPTION *opt = NULL;
	struct GMTAPI_CTRL *API = GMT->parent;

	//printf("file : %s line : %d func: %s, GMT->current.io.inc %d, GMT->current.io.inc_code[GMT_y] :%d\n",__FILE__,__LINE__,__func__,GMT->current.io.inc_code[GMT_X],GMT->current.io.inc_code[GMT_Y]);
//return;
	for (opt = options; opt; opt = opt->next) {

		switch (opt->option) {

			case '<':	/* Skip input files */

				if (!GMT_check_filearg (GMT, '<', opt->arg, GMT_IN)) n_errors++;

		//continue;;

				break;

			/* Processes program-specific parameters */

			case 'A':
				Ctrl->A.active = true;
				Ctrl->A.value = atof (opt->arg);
				break;
			case 'C':
				Ctrl->C.active = true;
				Ctrl->C.value = atof (opt->arg);
				break;
			case 'D':
				if ((Ctrl->D.active = GMT_check_filearg (GMT, 'D', opt->arg, GMT_IN)))
					Ctrl->D.file = strdup (opt->arg);
				else
					n_errors++;
				break;
			case 'G':

				if ((Ctrl->G.active = GMT_check_filearg (GMT, 'G', opt->arg, GMT_OUT)))
					Ctrl->G.file = strdup (opt->arg);
				else
					n_errors++;
				break;
			case 'I':

				Ctrl->I.active = true;
				if (GMT_getinc (GMT, opt->arg, Ctrl->I.inc)) {
					//GMT_inc_syntax (GMT, 'I', 1);
					n_errors++;
				}
				break;
			case 'L':	/* Set limits */
				Ctrl->L.active = true;
				switch (opt->arg[0]) {
					case 'l':	/* Lower limit  */
						n_errors += GMT_check_condition (GMT, opt->arg[1] == 0, "Syntax error -Ll option: No argument given\n");
						Ctrl->L.low = strdup (&opt->arg[1]);
						if (!GMT_access (GMT, Ctrl->L.low, F_OK))	/* file exists */
							Ctrl->L.lmode = 3;
						else if (Ctrl->L.low[0] == 'd')		/* Use data minimum */
							Ctrl->L.lmode = 1;
						else {
							Ctrl->L.lmode = 2;		/* Use given value */
							Ctrl->L.min = atof (&opt->arg[1]);
						}
						break;
					case 'u':	/* Upper limit  */
						n_errors += GMT_check_condition (GMT, opt->arg[1] == 0, "Syntax error -Lu option: No argument given\n");
						Ctrl->L.high = strdup (&opt->arg[1]);
						if (!GMT_access (GMT, Ctrl->L.high, F_OK))	/* file exists */
							Ctrl->L.hmode = 3;
						else if (Ctrl->L.high[0] == 'd')	/* Use data maximum */
							Ctrl->L.hmode = 1;
						else {
							Ctrl->L.hmode = 2;		/* Use given value */
							Ctrl->L.max = atof (&opt->arg[1]);
						}
						break;
					default:
						n_errors++;
						break;
				}
				break;
			case 'N':
				Ctrl->N.active = true;
				Ctrl->N.value = atoi (opt->arg);
				break;
			case 'Q':
				Ctrl->Q.active = true;
				break;
			case 'S':
				//Ctrl->S.active = true;
				//Ctrl->S.radius = atof (opt->arg);
				//Ctrl->S.unit = opt->arg[strlen(opt->arg)-1];
				//if (Ctrl->S.unit == 'c' && GMT_compat_check (GMT, 4)) {
				//	GMT_Report (API, GMT_MSG_COMPAT, "Warning: Unit c is deprecated; use s instead.\n");
				//	Ctrl->S.unit = 's';
				//}
				//if (!strchr ("sm ", Ctrl->S.unit)) {
				//	GMT_Report (API, GMT_MSG_NORMAL, "Syntax error -S option: Unrecognized unit %c\n", Ctrl->S.unit);
				//	n_errors++;
				//}
				//break;
			case 'T':
				Ctrl->T.active = true;
				k = 0;
				//if (GMT_compat_check (GMT, 4)) {	/* GMT4 syntax allowed for upper case */
					//modifier = opt->arg[strlen(opt->arg)-1];
					//if (modifier == 'B') modifier = 'b';
					//else if (modifier == 'I') modifier = 'i';
					//if (!(modifier == 'b' || modifier == 'i'))
					//	modifier = opt->arg[0], k = 1;
				//}
				//else {
					modifier = opt->arg[0];
					k = 1;
				//}
				if (modifier == 'b') {
					Ctrl->T.b_tension = atof (&opt->arg[k]);
				}
				else if (modifier == 'i') {
					Ctrl->T.i_tension = atof (&opt->arg[k]);
				}
				else if (modifier == '.' || (modifier >= '0' && modifier <= '9')) {
					Ctrl->T.i_tension = Ctrl->T.b_tension = atof (opt->arg);
				}
				else {
					printf ("Syntax error -T option: Unrecognized modifier %c\n", modifier);
					n_errors++;
				}
				break;
			case 'Z':
				Ctrl->Z.active = true;
				Ctrl->Z.value = atof (opt->arg);
				break;

			default:	/* Report bad options */
				//n_errors += GMT_default_error (GMT, opt->option);
				break;
		}

	}

	//return;
	GMT_check_lattice (GMT, Ctrl->I.inc, NULL, &Ctrl->I.active);

	n_errors += GMT_check_condition (GMT, !GMT->common.R.active, "Syntax error: Must specify -R option\n");
	n_errors += GMT_check_condition (GMT, Ctrl->I.inc[GMT_X] <= 0.0 || Ctrl->I.inc[GMT_Y] <= 0.0, "Syntax error -I option: Must specify positive increment(s)\n");
	n_errors += GMT_check_condition (GMT, Ctrl->N.value < 1, "Syntax error -N option: Max iterations must be nonzero\n");
	n_errors += GMT_check_condition (GMT, Ctrl->Z.value < 1.0 || Ctrl->Z.value > 2.0, "Syntax error -Z option: Relaxation value must be 1 <= z <= 2\n");
	n_errors += GMT_check_condition (GMT, !Ctrl->G.file, "Syntax error option -G: Must specify output file\n");
	n_errors += GMT_check_binary_io (GMT, 3);

	return (n_errors ? GMT_PARSE_ERROR : GMT_OK);
}

struct GMT_TEXTSET * GMT_create_textset (struct GMT_CTRL *GMT, uint64_t n_tables, uint64_t n_segments, uint64_t n_rows, bool alloc_only)
{	/* Create an empty text set structure with the required number of empty tables, all set to hold n_segments with n_rows */
	/* Allocate the new textset structure given the specified dimensions.
	 * IF alloc_only is true then we do NOT set the corresponding counters (i.e., n_segments).  */
	uint64_t tbl, seg;
	struct GMT_TEXTTABLE *T = NULL;
	struct GMT_TEXTSET *D = NULL;

	D = GMT_memory (GMT, NULL, 1, struct GMT_TEXTSET);
	D->table = GMT_memory (GMT, NULL, n_tables, struct GMT_TEXTTABLE *);
	D->n_alloc = D->n_tables = n_tables;
	if (!alloc_only) D->n_segments = n_tables * n_segments;
	for (tbl = 0; tbl < n_tables; tbl++) {
		D->table[tbl] = GMT_memory (GMT, NULL, 1, struct GMT_TEXTTABLE);
		T = D->table[tbl];
		T->n_alloc = n_segments;
		T->segment = GMT_memory (GMT, NULL, T->n_alloc, struct GMT_TEXTSEGMENT *);
		if (!alloc_only) T->n_segments = n_segments;
		for (seg = 0; seg < T->n_segments; seg++) {
			T->segment[seg] = GMT_memory (GMT, NULL, 1, struct GMT_TEXTSEGMENT);
			T->segment[seg]->record = GMT_memory (GMT, NULL, n_rows, char *);
			T->segment[seg]->n_alloc = n_rows;
			//T->segment[seg]->n_rows = n_rows;

		}
	}
	D->alloc_mode = GMT_ALLOCATED_BY_GMT;		/* Memory can be freed by GMT. */
	D->alloc_level = GMT->hidden.func_level;	/* Must be freed at this level. */
	D->id = GMT->parent->unique_var_ID++;		/* Give unique identifier */

	return (D);
}

bool GMT_x_is_outside (struct GMT_CTRL *GMT, double *x, double left, double right)
{
	/* Determines if this x is inside the effective x-domain.  This is normally
	 * west to east, but when gridding is concerned it can be extended by +-0.5 * dx
	 * for gridline-registered grids.  Also, if x is longitude we must check for
	 * wrap-arounds by 360 degrees, and x may be modified accordingly.
	 */
	if (GMT_is_dnan (*x)) return (true);
	if (GMT_x_is_lon (GMT, GMT_IN)) {	/* Periodic longitude test */
		while ((*x) > left) (*x) -= 360.0;	/* Make sure we start west or west */
		while ((*x) < left) (*x) += 360.0;	/* See if we are outside east */
		return (((*x) > right) ? true : false);
	}
	else	/* Cartesian test */
		return (((*x) < left || (*x) > right) ? true : false);
}


int read_data_surface (struct GMT_CTRL *GMT, struct SURFACE_INFO *C, struct GMT_OPTION *options)
{	/* Procdss input data into data structure */
	int i, j, error;
	uint64_t k, kmax = 0, kmin = 0, n_dup = 0;
	double *in, half_dx, zmin = DBL_MAX, zmax = -DBL_MAX, wesn_lim[4];
	struct GMT_GRID_HEADER *h = C->Grid->header;


	C->data = GMT_memory (GMT, NULL, C->n_alloc, struct SURFACE_DATA);

	/* Read in xyz data and computes index no and store it in a structure */

	if ((error = GMT_set_cols (GMT, GMT_IN, 3)) != GMT_OK) {
		return (error);
	}

	if (GMT_Init_IO (GMT->parent, GMT_IS_DATASET, GMT_IS_POINT, GMT_IN, GMT_ADD_DEFAULT, 0, options) != GMT_OK) {	/* Establishes data input */
		return (GMT->parent->error);
	}

	k = 0;
	C->z_mean = 0.0;
	/* Initially allow points to be within 1 grid spacing of the grid */
	wesn_lim[XLO] = h->wesn[XLO] - C->grid_xinc;	wesn_lim[XHI] = h->wesn[XHI] + C->grid_xinc;
	wesn_lim[YLO] = h->wesn[YLO] - C->grid_yinc;	wesn_lim[YHI] = h->wesn[YHI] + C->grid_yinc;
	half_dx = 0.5 * C->grid_xinc;

	if (GMT_Begin_IO (GMT->parent, GMT_IS_DATASET, GMT_IN, GMT_HEADER_ON) != GMT_OK) {	/* Enables data input and sets access mode */

		return (GMT->parent->error);
	}

	do {	/* Keep returning records until we reach EOF */

		if ((in = GMT_Get_Record (GMT->parent, GMT_READ_DOUBLE, NULL)) == NULL) {	/* Read next record, get NULL if special case */

			if (GMT_REC_IS_ERROR (GMT)) 		/* Bail if there are any read errors */
				return (GMT_RUNTIME_ERROR);
			if (GMT_REC_IS_ANY_HEADER (GMT)) 	/* Skip all headers */
				continue;
			if (GMT_REC_IS_EOF (GMT)) 		/* Reached end of file */
				break;
		}

		/* Data record to process */

		if (GMT_is_dnan (in[GMT_Z])) continue;
		if (GMT_y_is_outside (GMT, in[GMT_Y], wesn_lim[YLO], wesn_lim[YHI])) continue;	/* Outside y-range */
		if (GMT_x_is_outside (GMT, &in[GMT_X], wesn_lim[XLO], wesn_lim[XHI])) continue;	/* Outside x-range (or longitude) */
		if (C->periodic && (h->wesn[XHI]-in[GMT_X] < half_dx)) {	/* Push all values to the western nodes */
			in[GMT_X] -= 360.0;	/* Make this point be constraining the western node value */
			i = 0;
		}
		else
			i = irint (floor(((in[GMT_X]-h->wesn[XLO])*C->r_grid_xinc) + 0.5));

		if (i < 0 || i >= C->block_nx) continue;
		j = irint (floor(((in[GMT_Y]-h->wesn[YLO])*C->r_grid_yinc) + 0.5));
		if (j < 0 || j >= C->block_ny) continue;

		C->data[k].index = i * C->block_ny + j;
		C->data[k].x = (float)in[GMT_X];
		C->data[k].y = (float)in[GMT_Y];
		C->data[k].z = (float)in[GMT_Z];
		if (zmin > in[GMT_Z]) zmin = in[GMT_Z], kmin = k;
		if (zmax < in[GMT_Z]) zmax = in[GMT_Z], kmax = k;
		k++;
		C->z_mean += in[GMT_Z];
		if (k == C->n_alloc) {
			C->n_alloc <<= 1;
			C->data = GMT_memory (GMT, C->data, C->n_alloc, struct SURFACE_DATA);
		}
		if (C->periodic && i == 0) {	/* Replicate information to eastern boundary */
			i = C->block_nx - 1;
			C->data[k].index = i * C->block_ny + j;
			C->data[k].x = (float)(in[GMT_X] + 360.0);
			C->data[k].y = (float)in[GMT_Y];
			C->data[k].z = (float)in[GMT_Z];
			if (zmin > in[GMT_Z]) zmin = in[GMT_Z], kmin = k;
			if (zmax < in[GMT_Z]) zmax = in[GMT_Z], kmax = k;
			k++;
			C->z_mean += in[GMT_Z];
			if (k == C->n_alloc) {
				C->n_alloc <<= 1;
				C->data = GMT_memory (GMT, C->data, C->n_alloc, struct SURFACE_DATA);
			}
			n_dup++;
		}
	} while (true);


	if (GMT_End_IO (GMT->parent, GMT_IN, 0) != GMT_OK) {	/* Disables further data input */
		return (GMT->parent->error);
	}

	C->npoints = k;

	if (C->npoints == 0) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, " No datapoints inside region, aborts\n");
		return (EXIT_FAILURE);
	}

	C->z_mean /= k;
/*	if (GMT_is_verbose (GMT, GMT_MSG_VERBOSE)) {
		sprintf (C->format, "%s %s %s\n", GMT->current.setting.format_float_out, GMT->current.setting.format_float_out, GMT->current.setting.format_float_out);
		GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "Minimum value of your dataset x,y,z at: ");
		GMT_Report (GMT->parent, GMT_MSG_VERBOSE, C->format, (double)C->data[kmin].x, (double)C->data[kmin].y, (double)C->data[kmin].z);
		GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "Maximum value of your dataset x,y,z at: ");
		GMT_Report (GMT->parent, GMT_MSG_VERBOSE, C->format, (double)C->data[kmax].x, (double)C->data[kmax].y, (double)C->data[kmax].z);
		if (C->periodic && n_dup) GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "Number of input values shared between repeating west and east column nodes: %" PRIu64 "\n", n_dup);
	}*/
	C->data = GMT_memory (GMT, C->data, C->npoints, struct SURFACE_DATA);

	if (C->set_low == 1)
		C->low_limit = C->data[kmin].z;
	else if (C->set_low == 2 && C->low_limit > C->data[kmin].z) {
	/*	C->low_limit = data[kmin].z;	*/
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: Your lower value is > than min data value.\n");
	}
	if (C->set_high == 1)
		C->high_limit = C->data[kmax].z;
	else if (C->set_high == 2 && C->high_limit < C->data[kmax].z) {
	/*	C->high_limit = data[kmax].z;	*/
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: Your upper value is < than max data value.\n");
	}
	return (0);
}

void interp_breakline (struct GMT_CTRL *GMT, struct SURFACE_INFO *C, struct GMT_DATATABLE *xyzline)
{	/* Add constraints from breaklines */

	uint64_t n_tot = 0, this_ini = 0, this_end = 0, n_int = 0;
	uint64_t k = 0, n, kmax = 0, kmin = 0, row, seg;
	int srow, scol;
	size_t n_alloc;
	double *x = NULL, *y = NULL, *z = NULL, dx, dy, dz, r_dx, r_dy, zmin = DBL_MAX, zmax = -DBL_MAX;

	n_alloc = GMT_CHUNK;
	x = GMT_memory (GMT, NULL, n_alloc, double);
	y = GMT_memory (GMT, NULL, n_alloc, double);
	z = GMT_memory (GMT, NULL, n_alloc, double);

	r_dx = 1.0 / C->grid_xinc;
	r_dy = 1.0 / C->grid_yinc;
	for (seg = 0; seg < xyzline->n_segments; seg++) {
		for (row = 0; row < xyzline->segment[seg]->n_rows - 1; row++) {
			dx = xyzline->segment[seg]->coord[GMT_X][row+1] - xyzline->segment[seg]->coord[GMT_X][row];
			dy = xyzline->segment[seg]->coord[GMT_Y][row+1] - xyzline->segment[seg]->coord[GMT_Y][row];
			dz = xyzline->segment[seg]->coord[GMT_Z][row+1] - xyzline->segment[seg]->coord[GMT_Z][row];
			n_int = lrint (MAX (fabs(dx) * r_dx, fabs(dy) * r_dy ) ) + 1;
			this_end += n_int;

			if (n_alloc >= this_end) {
				n_alloc += MAX (GMT_CHUNK, n_int);
				x = GMT_memory (GMT, x, n_alloc, double);
				y = GMT_memory (GMT, x, n_alloc, double);
				z = GMT_memory (GMT, x, n_alloc, double);
			}

			dx /= (floor((double)n_int) - 1);
			dy /= (floor((double)n_int) - 1);
			dz /= (floor((double)n_int) - 1);
			for (k = this_ini, n = 0; k < this_end - 1; k++, n++) {
				x[k] = xyzline->segment[seg]->coord[GMT_X][row] + n * dx;
				y[k] = xyzline->segment[seg]->coord[GMT_Y][row] + n * dy;
				z[k] = xyzline->segment[seg]->coord[GMT_Z][row] + n * dz;
			}
			x[this_end - 1] = xyzline->segment[seg]->coord[GMT_X][row+1];
			y[this_end - 1] = xyzline->segment[seg]->coord[GMT_Y][row+1];
			z[this_end - 1] = xyzline->segment[seg]->coord[GMT_Z][row+1];

			this_ini += n_int;
		}

		n_tot += this_end;
	}

	/* Now add the interpolated breakline to the C structure */

	k = C->npoints;
	C->data = GMT_memory (GMT, C->data, k+n_tot, struct SURFACE_DATA);
	C->z_mean *= k;		/* It was already computed, reset it to sum */
	if (C->set_low == 1)
		zmin = C->low_limit;
	if (C->set_high == 1)
		zmax = C->high_limit;

	for (n = 0; n < n_tot; n++) {

		if (GMT_is_dnan (z[n])) continue;

		scol = irint (floor (((x[n] - C->Grid->header->wesn[XLO]) * C->r_grid_xinc) + 0.5));
		if (scol < 0 || scol >= C->block_nx) continue;
		srow = irint (floor (((y[n] - C->Grid->header->wesn[YLO]) * C->r_grid_yinc) + 0.5));
		if (srow < 0 || srow >= C->block_ny) continue;

		C->data[k].index = scol * C->block_ny + srow;
		C->data[k].x = (float)x[n];
		C->data[k].y = (float)y[n];
		C->data[k].z = (float)z[n];
		if (zmin > z[n]) zmin = z[n], kmin = k;
		if (zmax < z[n]) zmax = z[n], kmax = k;
		k++;
		C->z_mean += z[n];
	}

	if (k != (C->npoints + n_tot))		/* We had some NaNs */
		C->data = GMT_memory (GMT, C->data, k, struct SURFACE_DATA);

	C->npoints = k;
	C->z_mean /= k;

	if (C->set_low == 1)
		C->low_limit = C->data[kmin].z;
	if (C->set_high == 1)
		C->high_limit = C->data[kmax].z;

	GMT_free (GMT, x);
	GMT_free (GMT, y);
	GMT_free (GMT, z);
}

void throw_away_unusables (struct GMT_CTRL *GMT, struct SURFACE_INFO *C)
{
	/* This is a new routine to eliminate data which will become
		unusable on the final iteration, when grid = 1.
		It assumes grid = 1 and set_grid_parameters has been
		called.  We sort, mark redundant data as SURFACE_OUTSIDE, and
		sort again, chopping off the excess.

		Experimental modification 5 Dec 1988 by Smith, as part
		of a new implementation using core memory for b[6]
		coefficients, eliminating calls to temp file.
	*/

	uint64_t last_index, n_outside, k;

	/* Sort the data  */

	qsort (C->data, C->npoints, sizeof (struct SURFACE_DATA), compare_points);

	/* If more than one datum is indexed to same node, only the first should be kept.
		Mark the additional ones as SURFACE_OUTSIDE
	*/
	last_index = UINTMAX_MAX;
	n_outside = 0;
	for (k = 0; k < C->npoints; k++) {
		if (C->data[k].index == last_index) {
			C->data[k].index = SURFACE_OUTSIDE;
			n_outside++;
		}
		else
			last_index = C->data[k].index;
	}

	if (n_outside) {	/* Sort again; this time the SURFACE_OUTSIDE points will be thrown away  */
		qsort (C->data, C->npoints, sizeof (struct SURFACE_DATA), compare_points);
		C->npoints -= n_outside;
		C->data = GMT_memory (GMT, C->data, C->npoints, struct SURFACE_DATA);
		//GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "%" PRIu64 " unusable points were supplied; these will be ignored.\n", n_outside);
		//GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "You should have pre-processed the data with block-mean, -median, or -mode.\n");
	}
}

void remove_planar_trend (struct SURFACE_INFO *C)
{	/* Fit LS plane and remove from data; we restore before output */
	uint64_t i;
	double a, b, c, d, xx, yy, zz;
	double sx, sy, sz, sxx, sxy, sxz, syy, syz;
	struct GMT_GRID_HEADER *h = C->Grid->header;

	sx = sy = sz = sxx = sxy = sxz = syy = syz = 0.0;

	for (i = 0; i < C->npoints; i++) {

		xx = (C->data[i].x - h->wesn[XLO]) * h->r_inc[GMT_X];
		yy = (C->data[i].y - h->wesn[YLO]) * h->r_inc[GMT_Y];
		zz = C->data[i].z;

		sx += xx;
		sy += yy;
		sz += zz;
		sxx += (xx * xx);
		sxy += (xx * yy);
		sxz += (xx * zz);
		syy += (yy * yy);
		syz += (yy * zz);
	}

	d = C->npoints*sxx*syy + 2*sx*sy*sxy - C->npoints*sxy*sxy - sx*sx*syy - sy*sy*sxx;

	if (d == 0.0) {
		C->plane_c0 = C->plane_c1 = C->plane_c2 = 0.0;
		return;
	}

	a = sz*sxx*syy + sx*sxy*syz + sy*sxy*sxz - sz*sxy*sxy - sx*sxz*syy - sy*syz*sxx;
	b = C->npoints*sxz*syy + sz*sy*sxy + sy*sx*syz - C->npoints*sxy*syz - sz*sx*syy - sy*sy*sxz;
	c = C->npoints*sxx*syz + sx*sy*sxz + sz*sx*sxy - C->npoints*sxy*sxz - sx*sx*syz - sz*sy*sxx;

	C->plane_c0 = a / d;
	C->plane_c1 = b / d;
	C->plane_c2 = c / d;
	if (C->periodic) C->plane_c1 = 0.0;	/* Cannot have x-trend for periodic geographic data */

	for (i = 0; i < C->npoints; i++) {
		xx = (C->data[i].x - h->wesn[XLO]) * h->r_inc[GMT_X];
		yy = (C->data[i].y - h->wesn[YLO]) * h->r_inc[GMT_Y];
		C->data[i].z -= (float)(C->plane_c0 + C->plane_c1 * xx + C->plane_c2 * yy);
	}
}

int rescale_z_values (struct GMT_CTRL *GMT, struct SURFACE_INFO *C)
{	/* Find and normalize data by its rms */
	uint64_t i;
	double ssz = 0.0;

	for (i = 0; i < C->npoints; i++) ssz += (C->data[i].z * C->data[i].z);

	/* Set z_scale = rms(z) */

	C->z_scale = sqrt (ssz / C->npoints);

	if (C->z_scale < GMT_CONV_LIMIT) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: Input data lie exactly on a plane.\n");
		C->r_z_scale = C->z_scale = 1.0;
		return (1);	/* Flag to tell the main to just write out the plane */
	}
	else
		C->r_z_scale = 1.0 / C->z_scale;

	for (i = 0; i < C->npoints; i++) C->data[i].z *= (float)C->r_z_scale;

	if (C->converge_limit == 0.0) C->converge_limit = 0.001 * C->z_scale; /* i.e., 1 ppt of L2 scale */
	return (0);
}

int load_constraints (struct GMT_CTRL *GMT, struct SURFACE_INFO *C, int transform)
{	/* Deal with the constants or grids supplied via -L */
	unsigned int i, j;
	uint64_t ij;
	double yy;
	struct GMTAPI_CTRL *API = GMT->parent;

	/* Load lower/upper limits, verify range, deplane, and rescale */

	if (C->set_low > 0) {
		if (C->set_low < 3) {
			if ((C->Low = GMT_Duplicate_Data (API, GMT_IS_GRID, GMT_DUPLICATE_ALLOC, C->Grid)) == NULL) return (API->error);
			for (ij = 0; ij < C->mxmy; ij++) C->Low->data[ij] = (float)C->low_limit;
		}
		else {

			if ((C->Low = GMT_Read_Data (GMT->parent, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, C->low_file, NULL)) == NULL) return (API->error);	/* Get header only */
			if (C->Low->header->nx != C->Grid->header->nx || C->Low->header->ny != C->Grid->header->ny) {
				//GMT_Report (API, GMT_MSG_NORMAL, "Lower limit file not of proper dimension!\n");
				return (EXIT_FAILURE);
			}
			if (GMT_Read_Data (GMT->parent, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, C->low_file, C->Low) == NULL) return (API->error);
		}
		if (transform) {
			for (j = 0; j < C->Grid->header->ny; j++) {
				yy = (double)(C->Grid->header->ny - j - 1);
				for (i = 0; i < C->Grid->header->nx; i++) {
					ij = GMT_IJP (C->Grid->header, j, i);
					if (GMT_is_fnan (C->Low->data[ij])) continue;
					C->Low->data[ij] -= (float)(C->plane_c0 + C->plane_c1 * i + C->plane_c2 * yy);
					C->Low->data[ij] *= (float)C->r_z_scale;
				}
			}
		}
		C->constrained = true;
	}
	if (C->set_high > 0) {
		if (C->set_high < 3) {
			if ((C->High = GMT_Duplicate_Data (API, GMT_IS_GRID, GMT_DUPLICATE_ALLOC, C->Grid)) == NULL) return (API->error);
			for (ij = 0; ij < C->mxmy; ij++) C->High->data[ij] = (float)C->high_limit;
		}
		else {
			if ((C->High = GMT_Read_Data (GMT->parent, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, C->high_file, NULL)) == NULL) return (API->error);	/* Get header only */
			if (C->High->header->nx != C->Grid->header->nx || C->High->header->ny != C->Grid->header->ny) {
				//GMT_Report (API, GMT_MSG_NORMAL, "Upper limit file not of proper dimension!\n");
				return (EXIT_FAILURE);
			}
			if (GMT_Read_Data (GMT->parent, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, C->high_file, C->High) == NULL) return (API->error);
		}
		if (transform) {
			for (j = 0; j < C->Grid->header->ny; j++) {
				yy = (double)(C->ny - j - 1);
				for (i = 0; i < C->Grid->header->nx; i++) {
					ij = GMT_IJP (C->Grid->header, j, i);
					if (GMT_is_fnan (C->High->data[ij])) continue;
					C->High->data[ij] -= (float)(C->plane_c0 + C->plane_c1 * i + C->plane_c2 * yy);
					C->High->data[ij] *= (float)C->r_z_scale;
				}
			}
		}
		C->constrained = true;
	}

	return (0);
}

int GMT_surface (void *V_API, int mode, void *args)
{
	int error = 0, key, one = 1;
	double wesn[4];

	struct GMT_DATATABLE *xyzline = NULL;
	struct GMT_DATASET *Lin = NULL;
	struct SURFACE_INFO C;
	struct SURFACE_CTRL *Ctrl = NULL;
	struct GMT_CTRL *GMT = NULL, *GMT_cpy = NULL;
	struct GMT_OPTION *options = NULL;
	struct GMTAPI_CTRL *API = GMT_get_API_ptr (V_API);	/* Cast from void to GMTAPI_CTRL pointer */



	//confusing part start
	options = GMT_Create_Options (API, mode, args);	if (API->error) return (API->error);	/* Set or get option list */
	/*----------------------- Standard module initialization and parsing ----------------------*/
	/* Parse the command-line arguments */
	GMT = GMT_begin_module (API, THIS_MODULE_LIB, THIS_MODULE_NAME, &GMT_cpy); /* Save current state */

	if (GMT_Parse_Common (API, GMT_PROG_OPTIONS, options)) return (API->error);

	//parts end
	//printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);


	Ctrl = New_surface_Ctrl (GMT);	/* Allocate and initialize a new control structure */

	if (Ctrl == NULL)
		printf("something wrong\n");
	if ((error = GMT_surface_parse (GMT, Ctrl, options)))
		return (error);
	//return;

	/*---------------------------- This is the surface main code ----------------------------*/

	GMT_memset (&C, 1, struct SURFACE_INFO);
	GMT_memset (&GMT_Surface_Global, 1, struct SURFACE_GLOBAL);
	C.n_alloc = GMT_CHUNK;
	C.z_scale = C.r_z_scale = 1.0;
	C.mode_type[0] = 'I';
	C.mode_type[1] = 'D';	/* D means include data points when iterating */

	GMT_memcpy (C.wesn_orig, GMT->common.R.wesn, 4, double);	/* Save original region in case of -r */
	GMT_memcpy (wesn, GMT->common.R.wesn, 4, double);		/* Specified region */
	C.periodic = (GMT_is_geographic (GMT, GMT_IN) && GMT_360_RANGE (wesn[XLO], wesn[XHI]));
	if (GMT->common.r.active) {		/* Pixel registration request. Use the trick of offsetting area by x_inc(y_inc) / 2 */
		wesn[XLO] += Ctrl->I.inc[GMT_X] / 2.0;	wesn[XHI] += Ctrl->I.inc[GMT_X] / 2.0;
		wesn[YLO] += Ctrl->I.inc[GMT_Y] / 2.0;	wesn[YHI] += Ctrl->I.inc[GMT_Y] / 2.0;
		one++;	/* Just so we can report correct nx,ny for the grid; internally it is the same until output */
		/* nx,ny remains the same for now but nodes are in "pixel" position.  Must reset to original wesn and reduce nx,ny by 1 when we write result */
	}

	if ((C.Grid = GMT_Create_Data (API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, wesn, Ctrl->I.inc, \
		GMT_GRID_NODE_REG, 2, Ctrl->G.file)) == NULL) return (API->error);	/* Need a pad of 2 */



	if (C.Grid->header->nx < 4 || C.Grid->header->ny < 4) {
		printf ( "Error: Grid must have at least 4 nodes in each direction (you have %d by %d) - abort.\n", C.Grid->header->nx, C.Grid->header->ny);
		return (EXIT_FAILURE);
		//return 0; //changed by nishita
		}
	load_parameters_surface (&C, Ctrl);	/* Pass parameters from parsing control to surface INFO structure */

	C.relax_old = 1.0 - C.relax_new;

	C.nx = C.Grid->header->nx;
	C.ny = C.Grid->header->ny;
	C.nxny = C.Grid->header->nm;

	C.mx = C.Grid->header->mx;
	C.my = C.Grid->header->my;
	C.mxmy = C.Grid->header->size;
	GMT_Surface_Global.x_min = C.Grid->header->wesn[XLO];
	GMT_Surface_Global.y_min = C.Grid->header->wesn[YLO];

	/* New stuff here for v4.3: Check out the grid dimensions */
	C.grid = GMT_gcd_euclid (C.nx-1, C.ny-1);




	//if (GMT_is_verbose (GMT, GMT_MSG_VERBOSE) || Ctrl->Q.active) {
	//	sprintf (C.format, "Grid domain: W: %s E: %s S: %s N: %s nx: %%d ny: %%d [", GMT->current.setting.format_float_out, GMT->current.setting.format_float_out, GMT->current.setting.format_float_out, GMT->current.setting.format_float_out);
	//	(GMT->common.r.active) ? strcat (C.format, "pixel registration]\n") : strcat (C.format, "gridline registration]\n");
	///	GMT_Report (API, GMT_MSG_VERBOSE, C.format, C.wesn_orig[XLO], C.wesn_orig[XHI], C.wesn_orig[YLO], C.wesn_orig[YHI], C.nx-one, C.ny-one);
	//}
	if (C.grid == 1)
		//cout << "Warning: Your grid dimensions are mutually prime" <<endl;
	//if ((C.grid == 1 && GMT_is_verbose (GMT, GMT_MSG_VERBOSE)) || Ctrl->Q.active) suggest_sizes_for_surface (GMT, C.factors, C.nx-1, C.ny-1);
	if (Ctrl->Q.active) return (EXIT_SUCCESS);

	//printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);


	/* New idea: set grid = 1, read data, setting index.  Then throw
		away data that can't be used in end game, constraining
		size of briggs->b[6] structure.  */

	C.grid = 1;
	set_grid_parameters (&C);
	if (read_data_surface (GMT, &C, options)) return (EXIT_FAILURE);
	if (Ctrl->D.active) {
		if ((Lin = GMT_Read_Data (API, GMT_IS_DATASET, GMT_IS_FILE, GMT_IS_LINE, GMT_READ_NORMAL, NULL, Ctrl->D.file, NULL)) == NULL)
			return (API->error);
		xyzline = Lin->table[0];			/* Can only be one table since we read a single file */

		interp_breakline (GMT, &C, xyzline);
	}
	throw_away_unusables (GMT, &C);
	remove_planar_trend (&C);
	key = rescale_z_values (GMT, &C);
	if (key == 1) {	/* Data lies exactly on a plane; just return the plane grid */
		GMT_free (GMT, C.data);
		if (GMT_Create_Data (API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, NULL, NULL, 0, 0, C.Grid) == NULL)
			return (API->error);
		C.ij_sw_corner = 2 * C.my + 2;			/*  Corners of array of actual data  */
		replace_planar_trend (&C);
		if ((error = write_output_surface (GMT, &C, Ctrl->G.file))) return (error);
		return (EXIT_SUCCESS);
	}


	load_constraints (GMT, &C, true);





	/* Set up factors and reset grid to first value  */

	C.grid = GMT_gcd_euclid (C.nx-1, C.ny-1);
	C.n_fact = GMT_get_prime_factors (/*GMT, */C.grid, C.factors);
	set_grid_parameters (&C);
	while (C.block_nx < 4 || C.block_ny < 4) {
		smart_divide (&C);
		set_grid_parameters (&C);
	}
	set_offset (&C);
	set_index (&C);






	/* Now the data are ready to go for the first iteration.  */

	/* Allocate more space  */

	C.briggs = GMT_memory (GMT, NULL, C.npoints, struct SURFACE_BRIGGS);
	C.iu = GMT_memory (GMT, NULL, C.mxmy, char);
	if (GMT_Create_Data (API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, NULL, NULL, 0, 0, C.Grid) == NULL)
		return (API->error);





	if (C.radius > 0) initialize_grid (/*GMT,*/ &C); /* Fill in nodes with a weighted avg in a search radius  */

	//cout << "Grid\tMode\tIteration\tMax Change\tConv Limit\tTotal Iterations" <<endl;
	set_coefficients (&C);
	C.old_grid = C.grid;
	find_nearest_point (&C);
	iterate (GMT, &C, 1);
	while (C.grid > 1) {
		smart_divide (&C);
		set_grid_parameters (&C);
		set_offset (&C);
		set_index (&C);
		fill_in_forecast (&C);
		iterate (GMT, &C, 0);
		C.old_grid = C.grid;
		find_nearest_point (&C);
		iterate (GMT, &C, 1);
		//iterate ( &C, 1); // remove GMT
	}



	//if (GMT_is_verbose (GMT, GMT_MSG_VERBOSE)) check_errors (GMT, &C);




	replace_planar_trend (&C);





	GMT_free (GMT, C.data);
	GMT_free (GMT, C.briggs);
	GMT_free (GMT, C.iu);


	/*
	if ((C.set_low  > 0 && C.set_low < 3) && GMT_Destroy_Data (API, &C.Low) != GMT_OK) {
		GMT_Report (API, GMT_MSG_NORMAL, "Failed to free C.Low\n");
	}
	if ((C.set_high > 0 && C.set_high < 3) && GMT_Destroy_Data (API, &C.High) != GMT_OK) {
		GMT_Report (API, GMT_MSG_NORMAL, "Failed to free C.High\n");
	}
	*/



	if ((error = write_output_surface (GMT, &C, Ctrl->G.file))) return (error);


	return 0;
}

static int GMTAPI_session_counter = 0;	/* Keeps track of the ID of new sessions for multi-session programs */

void * GMT_Create_Session (char *session, unsigned int pad, unsigned int mode, int (*print_func) (FILE *, const char *))
{
	/* Initializes the GMT API for a new session. This is typically called once in a program,
	 * but programs that manage many threads might call it several times to create as many
	 * sessions as needed. [Note: There is of yet no thread support built into the GMT API
	 * but you could still manage many sessions at once].
	 * The session argument is a textstring used when reporting errors or messages from activity
	 *   originating within this session.
	 * Pad sets the default number or rows/cols used for grid padding.  GMT uses 2; users of
	 *   the API may wish to use 0 if they have no need for BCs, etc.
	 * The mode argument is a bitflag that controls a few things:
	 *   bit 1 means call return and not exit
	 *   bit 2 means we are called by an external API (e.g., Matlab, Python).
	 * We return the pointer to the allocated API structure.
	 * If any error occurs we report the error, set the error code via API->error, and return NULL.
	 * We terminate each session with a call to GMT_Destroy_Session.
	 */

	struct GMTAPI_CTRL *API = NULL;
	static char *unknown = "unknown";

	if ((API = calloc (1, sizeof (struct GMTAPI_CTRL))) == NULL) return_null (NULL, GMT_MEMORY_ERROR);	/* Failed to allocate the structure */
	API->pad = pad;		/* Preserve the default pad value for this session */
	//API->print_func = (print_func == NULL) ? gmt_print_func : print_func;	/* Pointer to the print function to use in GMT_Message|Report */
	API->do_not_exit = mode & 1;	/* if set, then API_exit & GMT_exit are simply a return; otherwise they call exit */
	API->mode = mode & 2;		/* if false|0 then we dont list read and write as modules */
	if (API->internal) API->leave_grid_scaled = 1;	/* Do NOT undo grid scaling after write since modules do not reuse grids we same some CPU */
	/* GMT_begin initializes, among onther things, the settings in the user's (or the system's) gmt.conf file */
	//GMT_begin (API, session, pad);

	if (GMT_begin (API, session, pad) == NULL) {		/* Initializing GMT and PSL machinery failed */
		free (API);	/* Free API */

		return_null (API, GMT_MEMORY_ERROR);
	}

	//GMT_Report (API, GMT_MSG_DEBUG, "GMT_Create_Session initialized GMT structure\n");

	/* Allocate memory to keep track of registered data resources */

	API->n_objects_alloc = GMT_SMALL_CHUNK;	/* Start small; this may grow as more resources are registered */
	API->object = GMT_memory (API->GMT, NULL, API->n_objects_alloc, struct GMTAPI_DATA_OBJECT *);

	/* Set the unique Session parameters */

	API->session_ID = GMTAPI_session_counter++;		/* Guarantees each session ID will be unique and sequential from 0 up */
	if (session) {
		API->session_tag = strdup ("surface");	/* Only used in reporting and error messages */	//GMT_basename (session)
		API->GMT->init.module_name = API->session_tag;	/* So non-modules can report name of program, */
	}
	else
		API->GMT->init.module_name = unknown; /* or unknown */

	//GMTAPI_init_sharedlibs (API);				/* Count how many shared libraries we should know about, and get their names and paths */

	return (API);	/* Pass the structure back out */
}

int GMT_Destroy_Session (void *V_API)
{
	/* GMT_Destroy_Session terminates the information for the specified session and frees all memory.
	 * Returns false if all is well and true if there were errors. */

	unsigned int i;
	struct GMTAPI_CTRL *API = gmt_get_api_ptr (V_API);

	if (API == NULL) return_error (API, GMT_NOT_A_SESSION);	/* GMT_Create_Session has not been called */

	//GMT_Report (API, GMT_MSG_DEBUG, "Entering GMT_Destroy_Session\n");
	//GMT_Garbage_Collection (API, GMT_NOTSET);	/* Free any remaining memory from data registration during the session */
	//GMTAPI_free_sharedlibs (API);			/* Close shared libraries and free list */

	/* Deallocate all remaining objects associated with NULL pointers (e.g., rec-by-rec i/o) */
	for (i = 0; i < API->n_objects; i++) GMTAPI_Unregister_IO (API, API->object[i]->ID, GMT_NOTSET);
	GMT_free (API->GMT, API->object);
	//GMT_end (API->GMT);	/* Terminate GMT machinery */
	free(API->GMT);
	if (API->session_tag) free (API->session_tag);
	GMT_memset (API, 1U, struct GMTAPI_CTRL);	/* Wipe it clean first */
 	free (API);	/* Not GMT_free since this item was allocated before GMT was initialized */

	return (GMT_OK);
}
void call_surface(char *pa)
{
	//path
	printf("Path in surface function: %s\n",pa);
	int status = GMT_NOT_A_VALID_MODULE;	/* Default status code */
	bool gmt_main = false;		/* Set to true if no module specified */
	unsigned int modulename_arg_n = 0;	/* Argument number in script[] that contains module name */

	struct GMTAPI_CTRL *api_ctrl = NULL;	// GMT API control structure 
	char gmt_module[GMT_LEN32] = "gmt";
	char *progname = "surface";			// Last component from the pathname 
	char *module = "surface";			// Module name 
	char* command1 = "surface";
	int argc = 6;
	char *script[6];
	char staEvent[500]={""};
	char grdfile [500]={""};
	char outfile [500]={""};

	//char staEvent[GMT_LEN16];
	//char grdfile [GMT_LEN32];
	//char outfile [GMT_LEN64];
	FILE *fp;
	char path[BUFFSIZE_NEW]={""};
	int i,count ;
	int round = 2; //use 0.0 and 0.2 interpretation

	//if (getcwd(path, sizeof(path)) == NULL)
	//	{
		//	printf("Error \n\n\n");
		//	return -1;
		//}

		strcpy(path,pa);
		//strcpy(path,"/home/labm/.core/TomoEK/data/");
		//strcat(path,DAT_FILE_LIST);


	fp = fopen(path,"r");
	if (fp == NULL)
		printf("File %s Not found!!\n",DAT_FILE_LIST);
	fscanf (fp, "%s", &staEvent);
	while (!feof (fp))
	{
		for (i = 0; i < round; i++)
			{
			if (i  == 0)
				{
				count = 0;
				argc = 5;
				//char *script[5];
				script[0] = "surface";
				script[1] = staEvent;
				sprintf(grdfile,"-G%s_0.%d.grd",staEvent,count);
				script[2] = grdfile;
				script[3] = "-I0.2";
				script[4] = "-R235/255/25/50";
				script[5] ="\0";
printf("I am an if\n");
				}
			else
				{
				argc = 6;
				//char *script[6];
				count = round;
				script[0] = "surface";
				script[1] = staEvent;
				sprintf(grdfile,"-G%s_0.%d.grd",staEvent,count);
				script[2] = "-T0.2";
				script[3] = grdfile;
				script[4] = "-I0.2";
				script[5] = "-R235/255/25/50";
				}
		//sprintf(grdfile,"%s_0.%d.grd",staEvent,count);


		//script[2] = "-R235/255/25/50";
		//sprintf(outfile,">%s.HD.0.%d",staEvent,count);
		//script[3] = outfile;


	// Initialize new GMT session 
	printf("GMT Script[0] = %s\n",script[0]);
	printf("GMT Script[1] = %s\n",script[1]);
	printf("GMT Script[2] = %s\n",script[2]);
	printf("GMT Script[3] = %s\n",script[3]);
	printf("GMT Script[4] = %s\n",script[4]);
	printf("GMT Script[5] = %s\n",script[5]);
	if ((api_ctrl = GMT_Create_Session (script[0], 2U, 0U, NULL)) == NULL)
		return EXIT_FAILURE;

	api_ctrl->internal = true;	// This is a proper GMT internal session (external programs will default to false) 


	// Test if script[0] contains a module name: 
	module = progname;	// Try this module name unless it equals PROGRAM_NAME in which case we just enter the test if argc > 1 
	gmt_main = !strcmp (module, PROGRAM_NAME);	// true if running the main program, false otherwise 

	if(api_ctrl==NULL)
		printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
	status = GMT_surface (api_ctrl, argc-1-modulename_arg_n, script+1+modulename_arg_n);
			}
		fscanf (fp, "%s", &staEvent);
	}
	return;
exit:
	if (progname) free (progname);
	// Destroy GMT session 
	if (GMT_Destroy_Session (api_ctrl))
		return EXIT_FAILURE;

	return status; // Return the status from the module 

}
void call_grd2xyz(char *pa, char *datapath)
{
	int status = GMT_NOT_A_VALID_MODULE;	/* Default status code */
	bool gmt_main = false;		/* Set to true if no module specified */
	unsigned int modulename_arg_n = 0;	/* Argument number in script[] that contains module name */
	struct GMTAPI_CTRL *api_ctrl = NULL;	/* GMT API control structure */
	char gmt_module[GMT_LEN32] = "gmt";
	char *progname = "grd2xyz";			/* Last component from the pathname */
	char *module = "grd2xyz";			/* Module name */
	char* command1 = "grd2xyz";
	int argc = 4; //3//
	char *script[4];

	char staEvent[500]={""};
	char grdfile [500]={""};
	char outfile [500]={""};

	//char staEvent[GMT_LEN16];
	//char grdfile [GMT_LEN32];
	//char outfile [GMT_LEN64];
	FILE *fp;
	char path[BUFFSIZE_NEW]={""};
	int i,count ;
	int round = 2; //use 0.0 and 0.2 interpretation

	//if (getcwd(path, sizeof(path)) == NULL)
	//	{
		//	printf("Error \n\n\n");
		//	return -1;
		//}

		strcpy(path,pa);		
		//strcpy(path,"/home/labm/.core/TomoEK/data/");
		//strcat(path,DAT_FILE_LIST);

	fp = fopen(path,"r");
	if (fp == NULL)
		printf("File %s Not found!!\n",DAT_FILE_LIST);
	fscanf (fp, "%s", &staEvent);
	while (!feof (fp))
	{
		for (i = 0; i < round; i++)
			{
			if (i  == 0)
				count = 0;
			else
				count = round;
		//sprintf(grdfile,"%s_0.%d.grd",staEvent,count);
		sprintf(grdfile,"%s%s_0.%d.grd",datapath,staEvent,count);
printf("GRDFILE NAME: -------------->>>>>> %s\n",grdfile);
		script[0] = "grd2xyz";
		script[1] = grdfile; //"travel_time_114A.txt_0.2.grd_test";
	//script[1] ="travel_time_114A.txt_0.0.grd_test";
		script[2] = "-R235/255/25/50";
		sprintf(outfile,">%s.HD_0.%d",staEvent,count);
		script[3] = outfile;
		//script[3] = ">travel_time_114A.txt_v1.HD.0.2_test";
	//script[3] = ">travel_time_114A.txt_v1.HD_test";
	//script[3] = "-I0.2";
	//script[4] = "travel_time_114A.txt_v1.HD";


	/* Initialize new GMT session */
	if ((api_ctrl = GMT_Create_Session (script[0], 2U, 0U, NULL)) == NULL)
		return EXIT_FAILURE;

	api_ctrl->internal = true;	/* This is a proper GMT internal session (external programs will default to false) */


	/* Test if script[0] contains a module name: */
	module = progname;	/* Try this module name unless it equals PROGRAM_NAME in which case we just enter the test if argc > 1 */
	gmt_main = !strcmp (module, PROGRAM_NAME);	/* true if running the main program, false otherwise */

	if(api_ctrl==NULL)
		printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
	status = GMT_grd2xyz (api_ctrl, argc-1-modulename_arg_n, script+1+modulename_arg_n);
	}
	fscanf (fp, "%s", &staEvent);
	}
	return;
exit:
	if (progname) free (progname);
	/* Destroy GMT session */
	if (GMT_Destroy_Session (api_ctrl))
		return EXIT_FAILURE;

	return status; /* Return the status from the module */

}

int writedatapath(char *completepath, char *node){
	FILE *f;
	char newpath[100],writefile[100];
	strcpy(writefile,completepath);
	strcat(writefile,"/data/");
	f=fopen("datapath.dat","w");
	fprintf(f,"%s",writefile);
	fclose(f);
	return 0;
}

double averageaggregate(double v1, double v2, double v3, double v4){
	
     double count1 =0.0;
     double acum,acum2,acum3;
     acum = v1 + v2;
     if(v1!=0 && v2!=0){
	acum = acum/2;
     }
     acum2 = acum + v3;
     if(v3!=0 && (v1!=0 || v2!=0)){
	acum2 = acum2/2;
     }
     acum3 = acum2 + v4;
     if(v4!=0 && (v1!=0 || v2!=0)){
	acum3 = acum3/2;
     }

  return acum3;

}

double simpleaverage(double v1, double v2, double v3, double v4){
	
     double count1 =0.0;
     double acum = v1 + v2 + v3 + v4;

     if(v1!=0) count1++;
     if(v2!=0) count1++;
     if(v3!=0) count1++;
     if(v4!=0) count1++;
     if(count1==0) return 0;
     else return acum/count1;

}

void aggregation(char *datapath){

	char path1[100],path2[100],path3[100],path4[100],pathc1[100],central[100];
	int i;
	strcpy(path1,datapath); strcpy(path2,datapath); strcpy(path3,datapath); strcpy(path4,datapath);
	strcpy(central,datapath);
	strcat(path1,"map1.tmo"); strcat(path2,"map2.tmo"); strcat(path3,"map3.tmo"); strcat(path4,"map4.tmo");
	strcat(central,"map600complete.tmo");
	//strcpy(pathc,datapath); 
	//strcat(pathc,"map600complete.tmo");
	printf("\n\n*****PATHS*******\n");
	printf("Path1: %s\n",path1);
	printf("Path2: %s\n",path2);
	printf("Path3: %s\n",path3);
	printf("Path4: %s\n",path4);
	char path_core[100];
	if (getcwd(path_core, sizeof(path_core)) == NULL)
	{
		printf("Error \n\n\n");
		return -1;		
	}
	strcat(path_core,"/");
	strcpy(pathc1,path_core);
	strcat(pathc1,"map600complete.tmo");
	printf("\n\n*****PATHS CORE*******\n");
	printf("Pathc1: %s\n",pathc1);



    double header,average;
    int dx,dx1,dy,dz,cont=0, real=0;
    double value1,value2,value3,value4,valuecom;
    FILE *f1_249, *f250_450, *f451_650, *f651, *endfile, *cen;
    f1_249 = fopen(path1,"r");
    f250_450 = fopen(path2,"r");
    f451_650 = fopen(path3,"r");
    f651 = fopen(path4,"r");
    cen = fopen(central,"w");
    endfile = fopen(pathc1,"w");

    for(i=0;i<12;i++){
        fscanf(f1_249,"%lf",&header);
        fscanf(f250_450,"%lf",&header);
        fscanf(f451_650,"%lf",&header);
        fscanf(f651,"%lf",&header);
        fprintf(endfile,"%lf\n",header);
	fprintf(cen,"%lf\n",header);
    }
    fscanf(f1_249,"%d %d %d",&dx1,&dy,&dz);
    fscanf(f250_450,"%d %d %d",&dx,&dy,&dz);
    fscanf(f451_650,"%d %d %d",&dx,&dy,&dz);
    fscanf(f651,"%d %d %d",&dx,&dy,&dz);
    fprintf(endfile,"%d\n%d\n%d\n",dx1,dy,dz);
    fprintf(cen,"%d\n%d\n%d\n",dx1,dy,dz);

    for(i=0;i<dx1*dy*dz;i++){
        fscanf(f1_249,"%lf",&value1);
        fscanf(f250_450,"%lf",&value2);
        fscanf(f451_650,"%lf",&value3);
        fscanf(f651,"%lf",&value4);
      
            average=averageaggregate(value1,value2,value3,value4);
            fprintf(endfile,"%lf\n",average);
	    fprintf(cen,"%lf\n",average);
            cont++;
        
    }
    
    fclose(f1_249);
    fclose(f250_450);
    fclose(f451_650);
    fclose(f651);
    fclose(endfile);
    fclose(cen);

}

int main(int argc,char *argv[]) {
	 if(argc != 4)
	    {
	      printf("usage:complete_path list_path node_number");
	      return 0;
	    }
	
	char V[100],pa[100],datapath[100];
	strcpy(pa,argv[1]);
	strcpy(datapath,argv[1]);
	strcat(pa,"/");
	strcat(datapath,"/");
	strcat(pa,argv[2]);
	strcat(datapath,"data/");
	printf("Complete path for list: %s\n",pa);
	printf("Data folder path: %s\n",datapath);
	printf("NODE: %s\n",argv[3]);
	writedatapath(argv[1],argv[3]);
	//axis
	//createArrivalTime(argv[1]);
	//calculateTiimeDiff(argv[2]);
	printf("It is entering in the main\n");
	//call_surface(pa);
        printf("Finishing surface function\n");
	//printf("\n\n\n\n\n\n\n");
	//call_grd2xyz(pa,datapath);
	createSlowMap(24,pa,datapath);
	
	char mapid[100];
	strcpy(mapid,"map");
	strcat(mapid,argv[3]);

	CalculateISOMap(0, 360, 18, mapid,pa,datapath) ;
	aggregation(datapath);
        printf("Other thing\n");
	return 0;

}
