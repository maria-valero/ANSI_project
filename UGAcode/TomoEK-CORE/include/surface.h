#ifndef _SURFACE_H
#define _SURFACE_H

#include "gmt_nan.h"
#include "gmt_resources.h"
#include "gmt_type.h"

//#define GMT_N_UNIQUE		59	/* Lines in gmt_unique.h */

struct SURFACE_DATA {	/* Data point and index to node it currently constrains  */
	float x;
	float y;
	float z;
	uint64_t index;
};

struct SURFACE_BRIGGS {		/* Coefficients in Taylor series for Laplacian(z) a la I. C. Briggs (1974)  */
	double b[6];
};

/*struct GMT_GRID {	 To hold a GMT float grid and its header in one container
	struct GMT_GRID_HEADER *header;	 Pointer to full GMT header for the grid
	float *data;			 Pointer to the float grid
 ---- Variables "hidden" from the API ----
	unsigned int id;		 The internal number of the grid
	unsigned int alloc_level;	 The level it was allocated at
	enum GMT_enum_alloc alloc_mode;	 Allocation mode [GMT_ALLOCATED_BY_GMT]
	void *extra;			 Row-by-row machinery information [NULL]
};*/

struct SURFACE_INFO {	/* Control structure for surface setup and execution */
	char *iu;			/* Pointer to grid info array */
	char mode_type[2];		/* D means include data points when iterating
					 * I means just interpolate from larger grid */
	char format[GMT_BUFSIZ];
	char *low_file, *high_file;	/* Pointers to grids with low and high limits, if selected */
	int grid, old_grid;	/* Node spacings  */
	unsigned int n_fact;		/* Number of factors in common (ny-1, nx-1) */
	unsigned int factors[32];		/* Array of common factors */
	unsigned int set_low;		/* 0 unconstrained,1 = by min data value, 2 = by user value */
	unsigned int set_high;		/* 0 unconstrained,1 = by max data value, 2 = by user value */
	size_t n_alloc;
	uint64_t npoints;			/* Number of data points */
	uint64_t ij_sw_corner, ij_se_corner,ij_nw_corner, ij_ne_corner;
	uint64_t n_empty;		/* No of unconstrained nodes at initialization  */
	int nx;				/* Number of nodes in x-dir. */
	int ny;				/* Number of nodes in y-dir. (Final grid) */
	uint64_t nxny;		/* Total number of grid nodes without boundaries  */
	int mx;
	int my;
	uint64_t mxmy;		/* Total number of grid nodes with boundaries  */
	int block_nx;		/* Number of nodes in x-dir for a given grid factor */
	int block_ny;		/* Number of nodes in y-dir for a given grid factor */
	unsigned int max_iterations;	/* Max iter per call to iterate */
	uint64_t total_iterations;
	bool periodic;		/* true if geographic grid and west-east == 360 */
	int grid_east;
	int offset[25][12];	/* Indices of 12 nearby points in 25 cases of edge conditions  */
	bool constrained;		/* true if set_low or set_high is true */
	double low_limit, high_limit;	/* Constrains on range of solution */
	double grid_xinc, grid_yinc;	/* size of each grid cell for a given grid factor */
	double r_grid_xinc, r_grid_yinc;	/* Reciprocals  */
	double converge_limit;		/* Convergence limit */
	double radius;			/* Search radius for initializing grid  */
	double tension;			/* Tension parameter on the surface  */
	double boundary_tension;
	double interior_tension;
	double a0_const_1, a0_const_2;	/* Constants for off grid point equation  */
	double e_2, e_m2, one_plus_e2;
	double eps_p2, eps_m2, two_plus_ep2, two_plus_em2;
	double x_edge_const, y_edge_const;
	double l_epsilon;
	double z_mean;
	double z_scale;			/* Root mean square range of z after removing planar trend  */
	double r_z_scale;		/* reciprocal of z_scale  */
	double plane_c0, plane_c1, plane_c2;	/* Coefficients of best fitting plane to data  */
	double small;			/* Let data point coincide with node if distance < C->small */
	double coeff[2][12];		/* Coefficients for 12 nearby points, constrained and unconstrained  */
	double relax_old, relax_new;	/* Coefficients for relaxation factor to speed up convergence */
	double wesn_orig[4];		/* Original -R domain as we might have shifted it due to -r */
	struct SURFACE_DATA  *data;
	struct SURFACE_BRIGGS *briggs;
	struct GMT_GRID *Grid;			/* The final grid */
	struct GMT_GRID *Low, *High;		/* arrays for minmax values, if set */
};

struct SURFACE_CTRL {
	struct A {	/* -A<aspect_ratio> */
		bool active;
		double value;
	} A;
	struct C {	/* -C<converge_limit> */
		bool active;
		double value;
	} C;
	struct D {	/* -D<line.xyz> */
		bool active;
		char *file;	/* Name of file with breaklines */
	} D;
	struct G {	/* -G<file> */
		bool active;
		char *file;
	} G;
	struct I {	/* -Idx[/dy] */
		bool active;
		double inc[2];
	} I;
	struct L {	/* -Ll|u<limit> */
		bool active;
		char *low, *high;
		double min, max;
		unsigned int lmode, hmode;
	} L;
	struct N {	/* -N<max_iterations> */
		bool active;
		unsigned int value;
	} N;
	struct Q {	/* -Q */
		bool active;
	} Q;
	struct S {	/* -S<radius>[m|c] */
		bool active;
		double radius;
		char unit;
	} S;
	struct T {	/* -T<tension>[i][b] */
		bool active;
		double b_tension, i_tension;
	} T;
	struct Z {	/* -Z<over_relaxation_parameter> */
		bool active;
		double value;
	} Z;
};




/*struct GMTAPI_CTRL {
	 Master controller which holds all GMT API related information at run-time for a single session.
	 * Users can run several GMT sessions concurrently; each session requires its own structure.
	 * Use GMTAPI_Create_Session to initialize a new session and GMTAPI_Destroy_Session to end it.

	uint64_t current_rec[2];		 Current record number >= 0 in the combined virtual dataset (in and out)
	unsigned int n_objects;			 Number of currently active input and output data objects
	unsigned int unique_ID;			 Used to create unique IDs for duration of session
	unsigned int session_ID;		 ID of this session
	unsigned int unique_var_ID;		 Used to create unique object IDs (grid,dataset, etc) for duration of session
	unsigned int current_item[2];		 Array number of current dataset being processed (in and out)
	unsigned int pad;			 Session default for number of rows/cols padding for grids [2]
	unsigned int mode;			 1 if called via external API (Matlab, Python) [0]
	unsigned int leave_grid_scaled;		 1 if we dont want to unpack a grid after we packed it for writing [0]
	bool registered[2];			 true if at least one source/destination has been registered (in and out)
	bool io_enabled[2];			 true if access has been allowed (in and out)
	size_t n_objects_alloc;			 Allocation counter for data objects
	int error;				 Error code from latest API call [GMT_OK]
	int last_error;				 Error code from previous API call [GMT_OK]
	int shelf;				 Place to pass hidden values within API
	unsigned int io_mode[2];		 1 if access as set, 0 if record-by-record
	struct GMT_CTRL *GMT;			 Key structure with low-level GMT internal parameters
	struct GMTAPI_DATA_OBJECT **object;	 List of registered data objects
	char *session_tag;			 Name tag for this session (or NULL)
	bool internal;				 true if session was initiated by gmt.c
	bool deep_debug;			 temporary for debugging
	int (*print_func) (FILE *, const char *);	 Pointer to fprintf function (may be reset by external APIs like MEX)
	unsigned int do_not_exit;		 0 by default, mieaning it is OK to call exit  (may be reset by external APIs like MEX to call return instead)
	//struct Gmt_libinfo *lib;		/* List of shared libs to consider */
	//unsigned int n_shared_libs;		 How many in lib
//};*/



	/* These are the 5 named option codes */
	enum GMT_enum_opt {
		GMT_OPT_USAGE = 	'?',	/* Command-line option for full usage */
		GMT_OPT_SYNOPSIS =	'^',	/* Command-line option for synopsis */
		GMT_OPT_PARAMETER = '-',	/* Command-line option for GMT defaults parameter */
		GMT_OPT_INFILE =	'<',	/* Command-line option for input file */
		GMT_OPT_OUTFILE =	'>'};	/* Command-line option for output file */


#endif
