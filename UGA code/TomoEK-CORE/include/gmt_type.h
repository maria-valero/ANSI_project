#ifndef _GMT_TYPE_H
#define _GMT_TYPE_H
#include <stdint.h>
#include <stdio.h>
#include "gmt_common.h"
#include "gmt_grd.h"
#include "gmt_time.h"
#include "gmt_map.h"
#include "io.h"
#include "gmt_defaults.h"
#include "gmt_private.h"
#include "gmt_calclock.h"

/* p_to_io_func is used as a pointer to functions such as GMT_read_d in assignments
 * and is used to declare GMT_get_io_ptr in gmt_io.c and gmt_prototypes.h */
typedef int (*p_to_io_func) (struct GMT_CTRL *, FILE *, uint64_t, double *);
	
#define GMT_DEG2MIN_F	60.0
#define GMT_DEG2SEC_F	3600.0

#define GMT_MIN2DEG	(1.0 / GMT_DEG2MIN_F)
#define GMT_SEC2DEG	(1.0 / GMT_DEG2SEC_F)
#define GMT_DAY2SEC_F	86400.0
#define GMT_DAY2SEC_I	86400
#define GMT_HR2SEC_I	3600
#define GMT_MIN2SEC_I	60

#define GMT_compat_check(C,version) (C->current.setting.compatibility <= version)	/* true if this section should be processed with backwards compatibility to given version */
#define GMT_fread(ptr,size,nmemb,stream) fread(ptr,size,nmemb,stream)
#define GMT_fwrite(ptr,size,nmemb,stream) fwrite(ptr,size,nmemb,stream)

#define GMT_N_UNIQUE		69	/* Lines in gmt_unique.h */


struct GMT_IO {				/* Used to process input data records */
	void * (*input) (struct GMT_CTRL *, FILE *, uint64_t *, int *);	/* Pointer to function reading ascii or binary tables */
	int (*output) (struct GMT_CTRL *, FILE *, uint64_t, double *);	/* Pointer to function writing ascii or binary tables */
	int (*read_item) (struct GMT_CTRL *, FILE *, uint64_t, double *);		/* Pointer to function reading 1-col z tables in grd2xyz */
	int (*write_item) (struct GMT_CTRL *, FILE *, uint64_t, double *);		/* Pointer to function writing 1-col z tables in xyz2grd */
	bool (*ogr_parser) (struct GMT_CTRL *, char *);				/* Set to handle either header or data OGR records */

	unsigned int pad[4];		/* pad[0] = west, pad[1] = east, pad[2] = south, pad[3] = north */
	unsigned int inc_code[2];
	double curr_rec[GMT_MAX_COLUMNS];	/* The most recently processed data record */
	double prev_rec[GMT_MAX_COLUMNS];	/* The previous data record */
	struct GMT_GRID_INFO grd_info;

	bool multi_segments[2];	/* true if current Ascii input/output file has multiple segments */
	bool skip_bad_records;	/* true if records where x and/or y are NaN or Inf */
	bool give_report;		/* true if functions should report how many bad records were skipped */
	bool skip_duplicates;	/* true if we should ignore duplicate x,y records */
	bool read_mixed;		/* true if we are reading ascii x y [z] [variable numbers of text] */
	bool need_previous;		/* true if when parsing a record we need access to previous record values (e.g., for gap or duplicate checking) */
	bool warn_geo_as_cartesion;	/* true if we should warn if we read a record with geographic data while the expected format has not been set (i.e., no -J or -fg) */

	uint64_t seg_no;		/* Number of current multi-segment in entire data set */
	uint64_t seg_in_tbl_no;		/* Number of current multi-segment in current table */
	uint64_t n_clean_rec;		/* Number of clean records read (not including skipped records or comments or blanks) */
	uint64_t n_bad_records;		/* Number of bad records encountered during i/o */
	unsigned int tbl_no;		/* Number of current table in entire data set */
	unsigned int io_nan_ncols;	/* Number of columns to consider for -s option */
	enum GMT_ogr_status ogr;	/* Tells us if current input source has OGR/GMT metadata (GMT_OGR_TRUE) or not (GMT_OGR_FALSE) or not set (GMT_OGR_UNKNOWN) */
	unsigned int status;		/* 0	All is ok
					   1	Current record is segment header
					   2	Mismatch between actual and expected fields
					   4	EOF
					   8	NaNs encountered in first 2/3 cols */
	uint64_t rec_no;		/* Number of current records (counts headers etc) in entire data set */
	uint64_t rec_in_tbl_no;		/* Number of current record (counts headers etc) in current table */
	uint64_t pt_no;			/* Number of current valid points in a row  */
	uint64_t curr_pos[2][3];	/* Keep track of current input/output table, segment, and row (for rec-by-rec action) */
	char r_mode[4];			/* Current file opening mode for reading (r or rb) */
	char w_mode[4];			/* Current file opening mode for writing (w or wb) */
	char a_mode[4];			/* Current file append mode for writing (a+ or ab+) */
	char current_record[GMT_BUFSIZ];	/* Current ascii record */
	char segment_header[GMT_BUFSIZ];	/* Current ascii segment header */
	char current_filename[2][GMT_BUFSIZ];	/* Current filenames (or <stdin>/<stdout>) */
	char *o_format[GMT_MAX_COLUMNS];	/* Custom output ascii format to overrule format_float_out */
	int ncid;			/* NetCDF file ID (when opening netCDF file) */
	int nvars;			/* Number of requested variablesin netCDF file */
	uint64_t ncols;			/* Number of total columns in netCDF file */
	size_t t_index[GMT_MAX_COLUMNS][5];		/* Indices for cross-sections (netCDF only) */
	size_t count[GMT_MAX_COLUMNS][5];		/* Count used for cross-sections (netCDF only) */
	size_t ndim;			/* Length of the column dimension */
	size_t nrec;			/* Record count */
	struct GMT_DATE_IO date_input;	/* Has all info on how to decode input dates */
	struct GMT_DATE_IO date_output;	/* Has all info on how to write output dates */
	struct GMT_CLOCK_IO clock_input;	/* Has all info on how to decode input clocks */
	struct GMT_CLOCK_IO clock_output;	/* Has all info on how to write output clocks */
	struct GMT_GEO_IO geo;		/* Has all the info on how to write geographic coordinates */
	bool skip_if_NaN[GMT_MAX_COLUMNS];	/* true if column j cannot be NaN and we must skip the record */
	bool col_skip[GMT_MAX_COLUMNS];	/* true of input column is to be ignored [Default reads all columns, but see -i] */
	unsigned int col_type[2][GMT_MAX_COLUMNS];	/* Type of column on input and output: Time, geographic, etc, see GMT_IS_<TYPE> */
	unsigned int io_nan_col[GMT_MAX_COLUMNS];	/* Array of columns to consider for -s option ir true */
	struct GMT_COL_INFO col[2][GMT_MAX_COLUMNS];	/* Order of columns on input and output unless 0,1,2,3,... */
	struct GMT_COL_TYPE fmt[2][GMT_MAX_COLUMNS];	/* Formatting information for binary data */
	struct GMT_OGR *OGR;		/* Pointer to GMT/OGR info used during reading */
	/* The remainder are just pointers to memory allocated elsewhere */
	int *varid;			/* Array of variable IDs (netCDF only) */
	double *scale_factor;		/* Array of scale factors (netCDF only) */
	double *add_offset;		/* Array of offsets (netCDF only) */
	double *missing_value;		/* Array of missing values (netCDF only) */
};

enum GMT_enum_cyl {GMT_MERCATOR = 100,
	GMT_CYL_EQ,
	GMT_CYL_EQDIST,
	GMT_CYL_STEREO,
	GMT_MILLER,
	GMT_TM,
	GMT_UTM,
	GMT_CASSINI,
	GMT_OBLIQUE_MERC = 150,
	GMT_OBLIQUE_MERC_POLE};

enum GMT_enum_azim {GMT_STEREO = 300,
	GMT_LAMB_AZ_EQ,
	GMT_ORTHO,
	GMT_AZ_EQDIST,
	GMT_GNOMONIC,
	GMT_GENPER,
	GMT_POLAR = 350};

struct GMT_THREE_D {
	double view_azimuth, view_elevation;
	double cos_az, sin_az, cos_el, sin_el;
	double corner_x[4], corner_y[4];
	double xmin, xmax, ymin, ymax;
	double world_x, world_y, world_z;	/* Users coordinates of fixed point */
	double view_x, view_y;			/* Desired projected 2-D coordinates of fixed point */
	double x_off, y_off;			/* Offsets to the final projected coordinates */
	double sign[4];		/* Used to determine direction of tickmarks etc */
	double level;		/* Indicates the last level of the perspective plane (if any) */
	unsigned int view_plane;	/* Determines on which plane needs to be projected */
	int plane;		/* Indicates which last plane was plotted in perspective (-1 = none) */
	unsigned int quadrant;	/* quadrant we're looking from */
	unsigned int z_axis;	/* Which z-axis to draw. */
	unsigned int face[3];	/* Tells if this facet has normal in pos direction */
	bool draw[4];	/* axes to draw */
	bool fixed;		/* true if we want a given point to be fixed in the projection [for animations] */
	bool world_given;	/* true if a fixed world point was given in -E ..+glon/lat/z */
	bool view_given;	/* true if a fixed projected point was given in -E ..+cx0/y0 */
};

struct GMT_PROJ4 {	/* Used to assign proj4 projections from GMT projections */
	char *name;
	unsigned int id;
};

#define GMT_N_PROJ4 31

#define GMT_ZAXIS		50

enum GMT_enum_conic {GMT_ALBERS = 200,
	GMT_ECONIC,
	GMT_POLYCONIC,
	GMT_LAMBERT = 250};
enum GMT_enum_misc {GMT_MOLLWEIDE = 400,
	GMT_HAMMER,
	GMT_SINUSOIDAL,
	GMT_VANGRINTEN,
	GMT_ROBINSON,
	GMT_ECKERT4,
	GMT_ECKERT6,
	GMT_WINKEL};

struct GMT_LATSWAP_CONSTS {
	double  c[GMT_LATSWAP_N][4];	/* Coefficients in 4-term series  */
	double	ra;			/* Authalic   radius (sphere for equal-area)  */
	double	rm;			/* Meridional radius (sphere for N-S distance)  */
	bool spherical;		/* True if no conversions need to be done.  */
};

struct GMT_PROJ {

	struct GMT_THREE_D z_project;
	//struct GMT_DATUM_CONV datum;	/* For datum conversions */
	struct GMT_PROJ4 *proj4;	/* A read-only resource we allocate once and pass pointer around */
	//void (*fwd) (struct GMT_CTRL *, double, double, double *, double *);/* Pointers to the selected forward mapping function */
	void (*inv) (struct GMT_CTRL *, double *, double *, double, double);/* Pointers to the selected inverse mapping function */
	//void (*fwd_x) (struct GMT_CTRL *, double, double *);	/* Pointers to the selected linear x forward function */
	//void (*fwd_y) (struct GMT_CTRL *, double, double *);	/* Pointers to the selected linear y forward function */
	//void (*fwd_z) (struct GMT_CTRL *, double, double *);	/* Pointers to the selected linear z forward function */
	//void (*inv_x) (struct GMT_CTRL *, double *, double);	/* Pointers to the selected linear x inverse function */
	//void (*inv_y) (struct GMT_CTRL *, double *, double);	/* Pointers to the selected linear y inverse function */
	//void (*inv_z) (struct GMT_CTRL *, double *, double);	/* Pointers to the selected linear z inverse function */

	double pars[10];		/* Raw unprocessed map-projection parameters as passed on command line */
	double z_pars[2];		/* Raw unprocessed z-projection parameters as passed on command line */

	/* Common projection parameters */

	int projection;		/* Gives the id number for the projection used (-1 if not set) */

	bool units_pr_degree;	/* true if scale is given as inch (or cm)/degree.  false for 1:xxxxx */
	bool north_pole;		/* true if projection is on northern hemisphere, false on southern */
	bool edge[4];		/* true if the edge is a map boundary */
	bool three_D;		/* Parameters for 3-D projections */
	bool JZ_set;		/* true if -Jz|Z was set */
	bool GMT_convert_latitudes;	/* true if using spherical code with authalic/conformal latitudes */
	bool inv_coordinates;	/* true if -fp[unit] was given and we must first recover lon,lat during reading */
	unsigned int n_antipoles;	/* Number of antipole coordinates so far [used for -JE only] */
	struct GMT_LATSWAP_CONSTS GMT_lat_swap_vals;

	enum GMT_enum_units inv_coord_unit;		/* Index to scale that converts input map coordinates to meter before inverting for lon,lat */
	char unit_name[GMT_N_UNITS][GMT_LEN16];	/* Names of the various distance units */
	double m_per_unit[GMT_N_UNITS];	/* Meters in various units.  Use to scale units to meters */
	double origin[3];		/* Projected values of the logical origin for the projection (x, y, z) */
	double rect[4], zmin, zmax;	/* Extreme projected values */
	double rect_m[4];		/* Extreme projected original meter values */
	double scale[3];		/* Scaling for meters to map-distance (typically inch) conversion (x, y, z) */
	double i_scale[3];		/* Inverse Scaling for meters to map-distance (typically inch) conversion (x, y, z) */
	double z_level;			/* Level at which to draw basemap [0] */
	double unit;			/* Gives meters pr plot unit (0.01 or 0.0254) */
	double central_meridian;	/* Central meridian for projection [NaN] */
	double lon0, lat0;		/* Projection center [NaN/NaN if not specified in -J] */
	double pole;			/* +90 pr -90, depending on which pole */
	double mean_radius;		/* Mean radius given the PROJ_* settings */
	double EQ_RAD, i_EQ_RAD;	/* Current ellipsoid parameters */
	double ECC, ECC2, ECC4, ECC6;	/* Powers of eccentricity */
	double M_PR_DEG, KM_PR_DEG;	/* Current spherical approximations to convert degrees to dist */
	double DIST_M_PR_DEG;		/* Current spherical approximations to convert degrees to m even if -J was not set */
	double DIST_KM_PR_DEG;		/* Current spherical approximations to convert degrees to km even if -J was not set */
	double half_ECC, i_half_ECC;	/* 0.5 * ECC and 0.5 / ECC */
	double one_m_ECC2, i_one_m_ECC2; /* 1.0 - ECC2 and inverse */
	unsigned int gave_map_width;	/* nonzero if map width (1), height (2), max dim (3) or min dim (4) is given instead of scale.  0 for 1:xxxxx */

	uint64_t n_geodesic_calls;	/* Number of calls for geodesics in this session */
	uint64_t n_geodesic_approx;	/* Number of calls for geodesics in this session that exceeded iteration limit */

	double f_horizon, rho_max;	/* Azimuthal horizon (deg) and in plot coordinates */

	/* Linear plot parameters */

	unsigned int xyz_projection[3];	/* For linear projection, 0 = linear, 1 = log10, 2 = pow */
	bool xyz_pos[3];		/* true if x,y,z-axis increases in normal positive direction */
	bool compute_scale[3];	/* true if axes lengths were set rather than scales */
	double xyz_pow[3];		/* For GMT_POW projection */
	double xyz_ipow[3];

	/* Center of radii for all conic projections */

	double c_x0, c_y0;

	/* Lambert conformal conic parameters. */

	double l_N, l_i_N, l_Nr, l_i_Nr;
	double l_F, l_rF, l_i_rF;
	double l_rho0;

	/* Oblique Mercator Projection (Spherical version )*/

	double o_sin_pole_lat, o_cos_pole_lat;	/* Pole of rotation */
	double o_pole_lon, o_pole_lat;	/* In degrees */
	double o_beta;			/* lon' = beta for central_meridian (degrees) */
	double o_FP[3], o_FC[3], o_IP[3], o_IC[3];

	/* TM and UTM Projections */

	double t_lat0;
	double t_e2, t_M0;
	double t_c1, t_c2, t_c3, t_c4;
	double t_i1, t_i2, t_i3, t_i4, t_i5;
	double t_r, t_ir;		/* Short for GMT->current.proj.EQ_RAD * GMT->current.setting.proj_scale_factor and its inverse */
	int utm_hemisphere;	/* -1 for S, +1 for N, 0 if to be set by -R */
	unsigned int utm_zonex;	/* The longitude component 1-60 */
	char utm_zoney;			/* The latitude component A-Z */

	/* Lambert Azimuthal Equal-Area Projection */

	double sinp;
	double cosp;
	double Dx, Dy, iDx, iDy;	/* Fudge factors for projections w/authalic lats */

	/* Stereographic Projection */

	double s_c, s_ic;
	double r;		/* Radius of projected sphere in plot units (inch or cm) */
	bool polar;		/* True if projection pole coincides with S or N pole */

	/* Mollweide, Hammer-Aitoff and Winkel Projection */

	double w_x, w_y, w_iy, w_r;

	/* Winkel Tripel Projection */

	double r_cosphi1;	/* = cos (50.467) */

	/* Robinson Projection */

	double n_cx, n_cy;	/* = = 0.8487R, 1.3523R */
	double n_i_cy;
	//double n_phi[GMT_N_ROBINSON], n_X[GMT_N_ROBINSON], n_Y[GMT_N_ROBINSON];
	//double n_x_coeff[3*GMT_N_ROBINSON], n_y_coeff[3*GMT_N_ROBINSON], n_iy_coeff[3*GMT_N_ROBINSON];

	/* Eckert IV Projection */

	double k4_x, k4_y, k4_ix, k4_iy;

	/* Eckert VI Projection */

	double k6_r, k6_ir;

	/* Cassini Projection */

	double c_M0, c_c1, c_c2, c_c3, c_c4;
	double c_i1, c_i2, c_i3, c_i4, c_i5, c_p;

	/* All Cylindrical Projections */

	double j_x, j_y, j_ix, j_iy;

	/* Albers Equal-area conic parameters. */

	double a_n, a_i_n;
	double a_C, a_n2ir2, a_test, a_Cin;
	double a_rho0;

	/* Equidistant conic parameters. */

	double d_n, d_i_n;
	double d_G, d_rho0;

	/* Van der Grinten parameters. */

	double v_r, v_ir;

        /* General Perspective parameters */
        double g_H, g_R;
        double g_P, g_P_inverse;
        double g_lon0;
        double g_sphi1, g_cphi1;
        double g_phig, g_sphig, g_cphig;
        double g_sdphi, g_cdphi;
        double g_B, g_D, g_L, g_G, g_J;
        double g_BLH, g_DG, g_BJ, g_DHJ, g_LH2, g_HJ;
        double g_sin_tilt, g_cos_tilt;
        double g_azimuth, g_sin_azimuth, g_cos_azimuth;
        double g_sin_twist, g_cos_twist;
        double g_width;
        double g_yoffset;
        double g_rmax;
        double g_max_yt;
        double g_xmin, g_xmax;
        double g_ymin, g_ymax;

        unsigned int g_debug;
        int g_box, g_outside, g_longlat_set, g_sphere, g_radius, g_auto_twist;

	/* Polar (cylindrical) projection */

	double p_base_angle;
	bool got_azimuths, got_elevations, z_down;

};

struct GMT_TIME_CONV {		/* Holds all time-related parameters */
	struct GMT_TRUNCATE_TIME truncate;
	struct GMT_Y2K_FIX Y2K_fix;		/* Used to convert 2-digit years to 4-digit years */
	struct GMT_TIME_LANGUAGE language;	/* For time axis */
	time_t tic;				/* Last system time marker */
	int64_t today_rata_die;			/* The rata die of current day at start of program */
};

struct GMT_PLOT_CALCLOCK {
	struct GMT_DATE_IO date;
	struct GMT_CLOCK_IO clock;
	struct GMT_GEO_IO geo;
};

struct GMT_PLOT {		/* Holds all plotting-related parameters */
	uint64_t n;			/* Number of such points */
	size_t n_alloc;			/* Size of allocated plot arrays */
	bool r_theta_annot;		/* true for special r-theta map annotation (see GMT_get_annot_label) */
	unsigned int mode_3D;		/* Determines if we draw fore and/or back 3-D box lines [Default is both] */
	unsigned int *pen;		/* Pen (PSL_MOVE = up, PSL_DRAW = down) for these points */
	struct GMT_PLOT_CALCLOCK calclock;
	/* The rest of the struct contains pointers that may point to memory not included by this struct */
	double *x;			/* Holds the x/y (inches) of a line to be plotted */
	double *y;
	char format[3][2][GMT_LEN256];	/* Keeps the 6 formats for dd:mm:ss plot output */
};

/* This struct is used to pass program options in/out of GMT modules */

struct GMT_OPTION { 			 /* Structure for a single GMT command option */
	char option;				 /* 1-char command line -<option> (e.g. D in -D) identifying the option (* if file) */
	char *arg;					 /* If not NULL, contains the argument for this option */
	struct GMT_OPTION *next;	 /* Pointer to next option in a linked list */
	struct GMT_OPTION *previous; /* Pointer to previous option in a linked list */
};

struct GMT_PLOT_AXIS_ITEM {		/* Information for one type of tick/annotation */
	double interval;		/* Distance between ticks in user units */
	unsigned int parent;		/* Id of axis this item belongs to (0,1,2) */
	bool active;			/* true if we want to use this item */
	bool special;		/* true if custom interval annotations */
	unsigned int flavor;		/* Index into month/day name abbreviation array (0-2) */
	bool upper_case;		/* true if we want upper case text (used with flavor) */
	char type;			/* One of a, A, i, I, f, F, g, G */
	char unit;			/* User's interval unit (y, M, u, d, h, m, c) */
};

struct GMT_PLOT_AXIS {		/* Information for one time axis */
	unsigned int id;		/* 0 (x), 1(y), or 2(z) */
	unsigned int type;		/* GMT_LINEAR, GMT_LOG10, GMT_POW, GMT_TIME */
	unsigned int special;		/* 0, GMT_CUSTOM, GMT_CPT */
	struct GMT_PLOT_AXIS_ITEM item[6];	/* see above defines for which is which */
	double phase;			/* Phase offset for strides: (knot-phase)%interval = 0  */
	char label[GMT_LEN256];	/* Label of the axis */
	char unit[GMT_LEN64];	/* Axis unit appended to annotations */
	char prefix[GMT_LEN64];	/* Axis prefix starting all annotations */
	char *file_custom;		/* File with custom annotations */
};

struct GMT_PLOT_FRAME {		/* Various parameters for plotting of time axis boundaries */
	struct GMT_PLOT_AXIS axis[3];	/* One each for x, y, and z */
	char header[GMT_LEN256];	/* Plot title */
	struct GMT_FILL fill;		/* Fill for the basemap inside, if paint == true */
	bool plotted_header;		/* true if header has been plotted */
	bool init;			/* true if -B was used */
	bool draw;			/* true if -B<int> was used, even -B0, as sign to draw axes */
	bool paint;			/* true if -B +g<fill> was used */
	bool draw_box;			/* true is a 3-D Z-box is desired */
	bool check_side;		/* true if lon and lat annotations should be on x and y axis only */
	bool primary;			/* true if current axis is primary, false if secondary */
	bool slash;			/* true if slashes were used in the -B argument */
	bool obl_grid;			/* true if +o was given to draw oblique gridlines */
	unsigned int set_frame[2];	/* 1 if a -B<WESNframe> setting was given */
	unsigned int horizontal;	/* 1 is S/N annotations should be parallel to axes, 2 if forced */
	unsigned int side[5];		/* Which sides (0-3 in plane; 4 = z) to plot. 2 is annot/draw, 1 is draw, 0 is not */
	unsigned int z_axis[4];		/* Which axes to use for the 3-D z-axis [auto] */
};

struct GMT_MAP {		/* Holds all map-related parameters */
	struct GMT_PLOT_FRAME frame;		/* Everything about the frame parameters */
	int this_x_status;			/* Tells us what quadrant old and new points are in (-4/4) */
	int this_y_status;
	int prev_x_status;
	int prev_y_status;
	int corner;			/* Tells us which corner 1-4 or -1 if not a corner */
	bool on_border_is_outside;		/* true if a point exactly on the map border shoud be considered outside the map */
	bool is_world;			/* true if map has 360 degrees of longitude range */
	bool is_world_tm;			/* true if GMT_TM map is global? */
	bool lon_wrap;			/* true when longitude wrapping over 360 degrees is allowed */
	bool z_periodic;			/* true if grid values are 0-360 degrees (phases etc) */
	bool loxodrome;				/* true if we are computing loxodrome distances */
	unsigned int meridian_straight;		/* 1 if meridians plot as straight lines, 2 for special case */
	unsigned int parallel_straight;		/* 1 if parallels plot as straight lines, 2 for special case */
	unsigned int n_lon_nodes;		/* Somewhat arbitrary # of nodes for lines in longitude (may be reset in gmt_map.c) */
	unsigned int n_lat_nodes;		/* Somewhat arbitrary # of nodes for lines in latitude (may be reset in gmt_map.c) */
	unsigned int path_mode;		/* 0 if we should call GMT_fix_up_path to resample across gaps > path_step, 1 to leave alone */
	double width;				/* Full width in inches of this world map */
	double height;				/* Full height in inches of this world map */
	double half_width;			/* Half width in inches of this world map */
	double half_height;			/* Half height of this world map */
	double dlon;				/* Steps taken in longitude along gridlines (gets reset in gmt_init.c) */
	double dlat;				/* Steps taken in latitude along gridlines (gets reset in gmt_init.c) */
	double path_step;			/* Sampling interval if resampling of paths should be done */
	bool (*outside) (struct GMT_CTRL *, double, double);	/* Pointer to function checking if a lon/lat point is outside map */
	bool (*overlap) (struct GMT_CTRL *, double, double, double, double);	/* Pointer to function checking for overlap between 2 regions */
	bool (*will_it_wrap) (struct GMT_CTRL *, double *, double *, uint64_t, uint64_t *);	/* true if consecutive points indicate wrap */
	int (*jump) (struct GMT_CTRL *, double, double, double, double);	/* true if we jump in x or y */
	unsigned int (*crossing) (struct GMT_CTRL *, double, double, double, double, double *, double *, double *, double *, unsigned int *);	/* Pointer to functions returning crossover point at boundary */
	uint64_t (*clip) (struct GMT_CTRL *, double *, double *, uint64_t, double **, double **, uint64_t *);	/* Pointer to functions that clip a polygon to fit inside map */
	double (*left_edge) (struct GMT_CTRL *, double);	/* Pointers to functions that return left edge of map */
	double (*right_edge) (struct GMT_CTRL *, double);	/* Pointers to functions that return right edge of map */
	struct GMT_DIST dist[3];		/* struct with pointers to functions/scales returning distance between two points points */
	bool (*near_lines_func) (struct GMT_CTRL *, double, double, struct GMT_DATATABLE *, unsigned int, double *, double *, double *);	/* Pointer to function returning distance to nearest line among a set of lines */
	bool (*near_a_line_func) (struct GMT_CTRL *, double, double, uint64_t, struct GMT_DATASEGMENT *, unsigned int, double *, double *, double *);	/* Pointer to function returning distance to line */
	bool (*near_point_func) (struct GMT_CTRL *, double, double, struct GMT_DATATABLE *, double);	/* Pointer to function returning distance to nearest point */	
	unsigned int (*wrap_around_check) (struct GMT_CTRL *, double *, double, double, double, double, double *, double *, unsigned int *);	/* Does x or y wrap checks */
	double (*azimuth_func) (struct GMT_CTRL *, double, double, double, double, bool);	/* Pointer to function returning azimuth between two points points */
	void (*get_crossings) (struct GMT_CTRL *, double *, double *, double, double, double, double);	/* Returns map crossings in x or y */
};

struct GMT_CURRENT {
	/* These are internal parameters that need to be passed around between
	 * many GMT functions.  These values may change by user interaction. */
	struct GMT_DEFAULTS setting;	/* Holds all GMT defaults parameters */
	struct GMT_IO io;		/* Holds all i/o-related parameters */
	struct GMT_PROJ proj;		/* Holds all projection-related parameters */
	struct GMT_MAP map;		/* Holds all projection-related parameters */
	struct GMT_PLOT plot;		/* Holds all plotting-related parameters */
	struct GMT_TIME_CONV time;	/* Holds all time-related parameters */
	//struct GMT_PS ps;		/* Hold parameters related to PS setup */
	struct GMT_OPTION *options;	/* Pointer to current program's options */
	//struct GMT_FFT_HIDDEN fft;	/* Structure with info that must survive between FFT calls */
};
struct GMT_INTERNAL {
	/* These are internal parameters that need to be passed around between
	 * many GMT functions.  These may change during execution but are not
	 * modified directly by user interaction. */
	unsigned int func_level;	/* Keeps track of what level in a nested GMT_func calling GMT_func etc we are.  0 is top function */
	size_t mem_cols;		/* Current number of allocated columns for temp memory */
	size_t mem_rows;		/* Current number of allocated rows for temp memory */
	double **mem_coord;		/* Columns of temp memory */
//#ifdef MEMDEBUG
	//struct MEMORY_TRACKER *mem_keeper;
//#endif
};

#define GMT_N_GRD_FORMATS 25 /* Number of formats above plus 1 */

struct GMT_SHORTHAND {	/* Holds information for each grid extension shorthand read from the user's .gmtio file */
	char *suffix; /* suffix of file */
	char *format; /* format: ff/scale/offset/invalid */
};

struct GMT_SESSION {
	/* These are parameters that is set once at the start of a GMT session and
	 * are essentially read-only constants for the duration of the session */
	FILE *std[3];			/* Pointers for standard input, output, and error */
	void * (*input_ascii) (struct GMT_CTRL *, FILE *, uint64_t *, int *);	/* Pointer to function reading ascii tables only */
	int (*output_ascii) (struct GMT_CTRL *, FILE *, uint64_t, double *);	/* Pointer to function writing ascii tables only */
	unsigned int n_fonts;		/* Total number of fonts returned by GMT_init_fonts */
	unsigned int n_user_media;	/* Total number of user media returned by gmt_load_user_media */
	size_t min_meminc;		/* with -DMEMDEBUG, sets min/max memory increments */
	size_t max_meminc;
	float f_NaN;			/* Holds the IEEE NaN for floats */
	double d_NaN;			/* Holds the IEEE NaN for doubles */
	double no_rgb[4];		/* To hold {-1, -1, -1, 0} when needed */
	double u2u[4][4];		/* u2u is the 4x4 conversion matrix for cm, inch, m, pt */
	char unit_name[4][8];		/* Full name of the 4 units cm, inch, m, pt */
	//struct GMT_HASH rgb_hashnode[GMT_N_COLOR_NAMES];/* Used to translate colornames to r/g/b */
	bool rgb_hashnode_init;		/* true once the rgb_hashnode array has been loaded; false otherwise */
	unsigned int n_shorthands;			/* Length of arrray with shorthand information */
	char *grdformat[GMT_N_GRD_FORMATS];	/* Type and description of grid format */
	int (*readinfo[GMT_N_GRD_FORMATS]) (struct GMT_CTRL *, struct GMT_GRID_HEADER *);	/* Pointers to grid read header functions */
	int (*updateinfo[GMT_N_GRD_FORMATS]) (struct GMT_CTRL *, struct GMT_GRID_HEADER *);	/* Pointers to grid update header functions */
	int (*writeinfo[GMT_N_GRD_FORMATS]) (struct GMT_CTRL *, struct GMT_GRID_HEADER *);	/* Pointers to grid write header functions */
	int (*readgrd[GMT_N_GRD_FORMATS]) (struct GMT_CTRL *, struct GMT_GRID_HEADER *, float *, double *, unsigned int *, unsigned int);	/* Pointers to grid read functions */
	int (*writegrd[GMT_N_GRD_FORMATS]) (struct GMT_CTRL *, struct GMT_GRID_HEADER *, float *, double *, unsigned int *, unsigned int);	/* Pointers to grid read functions */
	//int (*fft1d[k_n_fft_algorithms]) (struct GMT_CTRL *, float *, unsigned int, int, unsigned int);	/* Pointers to available 1-D FFT functions (or NULL if not configured) */
	//int (*fft2d[k_n_fft_algorithms]) (struct GMT_CTRL *, float *, unsigned int, unsigned int, int, unsigned int);	/* Pointers to available 2-D FFT functions (or NULL if not configured) */
	/* This part contains pointers that may point to additional memory outside this struct */
	char *DCWDIR;			/* Path to the DCW directory */
	char *GSHHGDIR;			/* Path to the GSHHG directory */
	char *SHAREDIR;			/* Path to the GMT share directory */
	char *HOMEDIR;			/* Path to the user's home directory */
	char *USERDIR;			/* Path to the user's GMT settings directory */
	char *DATADIR;			/* Path to one or more directories with data sets */
	char *TMPDIR;			/* Path to the directory directory for isolation mode */
	char *CUSTOM_LIBS;		/* Names of one or more comma-separated GMT-compatible shared libraries */
	char **user_media_name;		/* Length of array with custom media dimensions */
	//struct GMT_FONTSPEC *font;		/* Array with font names and height specification */
	//struct GMT_MEDIA *user_media;		/* Array with custom media dimensions */
	struct GMT_SHORTHAND *shorthand;	/* Array with info about shorthand file extension magic */
};
struct GMT_INIT { /* Holds misc run-time parameters */
	unsigned int n_custom_symbols;
	const char *module_name;      /* Name of current module or NULL if not set */
	const char *module_lib;       /* Name of current shared library or NULL if not set */
	/* The rest of the struct contains pointers that may point to memory not included by this struct */
	char *runtime_bindir;         /* Directory that contains the main exe at run-time */
	char *runtime_libdir;         /* Directory that contains the main shared lib at run-time */
	char *history[GMT_N_UNIQUE];  /* The internal gmt.history information */
	//struct GMT_CUSTOM_SYMBOL **custom_symbol; /* For custom symbol plotting in psxy[z]. */
};

struct GMT_CTRL {
	/* Master structure for a GMT invokation.  All internal settings for GMT is accessed here */
	struct GMT_SESSION session;	/* Structure with all values that do not change throughout a session */
	struct GMT_INIT init;		/* Structure with all values that do not change in a GMT_func call */
	struct GMT_COMMON common;	/* Structure with all the common GMT command settings (-R -J ..) */
	struct GMT_CURRENT current;	/* Structure with all the GMT items that can change during execution, such as defaults settings (pens, colors, fonts.. ) */
	struct GMT_INTERNAL hidden;	/* Internal global variables that are not to be changed directly by users */
	//struct PSL_CTRL *PSL;		/* Pointer to the PSL structure [or NULL] */
	struct GMTAPI_CTRL *parent;	/* Owner of this structure [or NULL]; gives access to the API from functions being passed *GMT only */
};

/*#define bswap16 gnuc_bswap16

static inline uint16_t gnuc_bswap16(uint16_t x) {
		if (__builtin_constant_p(x))
			x = inline_bswap16(x);
		else {
#		ifdef __x86_64__
			__asm__("xchgb %h0, %b0" : "+Q" (x));
#		elif defined __i386__
			__asm__("xchgb %h0, %b0" : "+q" (x));
#		endif
		}
		return x;
	}*/

/* If GMT is not set or no_not_exit is false then we call system exit, else we move along */
static inline void GMT_exit (struct GMT_CTRL *GMT, int code) {
	if (GMT == NULL || GMT->parent == NULL || GMT->parent->do_not_exit == false)
		exit (code);
}

#define GMT_fdopen(handle, mode) fdopen(handle, mode)
#endif

