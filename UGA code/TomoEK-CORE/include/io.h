#ifndef _IO_H
#define _IO_H
#include <stdio.h>
/* Low-level structures used internally */
struct GMT_COL_INFO {	/* Used by -i and input parsing */
	unsigned int col;	/* The column number in the order requested via -i */
	unsigned int order;	/* The initial order (0,1,...) but this will be sorted on col */
	unsigned int convert;	/* 2 if we must convert the data by log10, 1 if scale, offset */
	double scale;		/* Multiplier for raw in value */
	double offset;		/* Offset applied after multiplier */
};

/* THere are three GMT/OGR status values */
enum GMT_ogr_status {
	GMT_OGR_UNKNOWN = -1,	/* We have not parsed enough records to know yet */
	GMT_OGR_FALSE,		/* This is NOT a GMT/OGR file */
	GMT_OGR_TRUE};		/* This is a GMT/OGR file */

struct GMT_COL_TYPE {	/* Used by -b for binary formatting */
	unsigned int type;	/* Data type e.g., GMT_FLOAT */
	off_t skip;		/* Rather than read/write an item, jump |skip| bytes before (-ve) or after (+ve) read/write */
	int (*io) (struct GMT_CTRL *, FILE *, uint64_t, double *);	/* Pointer to the correct read or write function given type/swab */
};

struct GMT_QUAD {	/* Counting parameters needed to determine proper longitude min/max range */
	uint64_t quad[4];		/* Keeps track if a longitude fell in these quadrants */
	unsigned int range[2];	/* The format for reporting longitude */
	double min[2], max[2];		/* Min/max values in either -180/180 or 0/360 counting */
};

struct GMT_GRID_ROWBYROW {	/* Holds book-keeping information needed for row-by-row actions */
	size_t size;		/* Bytes per item [4 for float, 1 for byte, etc] */
	size_t n_byte;		/* Number of bytes for row */
	unsigned int row;	/* Current row */
	bool open;		/* true if we have already opened the file */
	bool check;		/* true if we must replace NaNs with another representation on i/o */
	bool auto_advance;	/* true if we want to read file sequentially */

	int fid;		/* NetCDF file number [netcdf files only] */
	size_t edge[2];		/* Dimension arrays [netcdf files only] */
	size_t start[2];	/* Position arrays [netcdf files only] */

	FILE *fp;		/* File pointer [for native files] */

	void *v_row;		/* Void Row pointer for any data format */
};

/* Types of possible column entries in a file: */

enum GMT_col_enum {
	GMT_IS_NAN   =   0,	/* Returned by GMT_scanf routines when read fails */
	GMT_IS_FLOAT		=   1,	/* Generic (double) data type, no special format */
	GMT_IS_LAT		=   2,
	GMT_IS_LON		=   4,
	GMT_IS_GEO		=   6,	/* data type is either Lat or Lon */
	GMT_IS_RELTIME		=   8,	/* For I/O of data in user units */
	GMT_IS_ABSTIME		=  16,	/* For I/O of data in calendar+clock units */
	GMT_IS_RATIME		=  24,	/* To see if time is either Relative or Absolute */
	GMT_IS_ARGTIME		=  32,	/* To invoke GMT_scanf_argtime()  */
	GMT_IS_DIMENSION	=  64,	/* A float with [optional] unit suffix, e.g., 7.5c, 0.4i; convert to inch  */
	GMT_IS_GEOANGLE		= 128,	/* An angle to be converted via map projection to angle on map  */
	GMT_IS_STRING		= 256,	/* An text argument [internally used, not via -f]  */
	GMT_IS_UNKNOWN		= 512};	/* Input type is not knowable without -f */

struct GMT_DATE_IO {
	bool skip;			/* Only true if a format string was pass as NULL */
	unsigned int T_pos;		/* String position of the expected 'T' marker (INPUT only) */
	int item_order[4];		/* The sequence year, month, day, day-of-year in input calendar string (-ve if unused) */
	int item_pos[4];		/* Which position year, month, day, day-of-year has in calendar string (-ve if unused) */
	bool Y2K_year;		/* true if we have 2-digit years */
	bool truncated_cal_is_ok;	/* true if we have YMD or YJ order so smallest unit is to the right */
	bool iso_calendar;		/* true if we do ISO week calendar */
	bool day_of_year;		/* true if we do day-of-year rather than month/day */
	bool mw_text;		/* true if we must plot the month name or Week rather than a numeral */
	bool compact;		/* true if we do not want leading zeros in items (e.g., 03) */
	char format[GMT_LEN64];	/* Actual C format used to input/output date */
	char delimiter[2][2];		/* Delimiter strings in date, e.g. "-" */
};

struct GMT_CLOCK_IO {
	bool skip;			/* Only true if a format string was pass as NULL */
	double f_sec_to_int;		/* Scale to convert 0.xxx seconds to integer xxx (used for formatting) */
	int order[3];		/* The relative order of hour, mn, sec in input clock string (-ve if unused) */
	unsigned int n_sec_decimals;	/* Number of digits in decimal seconds (0 for whole seconds) */
	bool compact;		/* true if we do not want leading zeros in items (e.g., 03) */
	bool twelve_hr_clock;	/* true if we are doing am/pm on output */
	char ampm_suffix[2][8];		/* Holds the strings to append am or pm */
	char format[GMT_LEN64];	/* Actual C format used to output clock */
	char delimiter[2][2];		/* Delimiter strings in clock, e.g. ":" */
};

struct GMT_GEO_IO {			/* For geographic output and plotting */
	double f_sec_to_int;		/* Scale to convert 0.xxx seconds to integer xxx (used for formatting) */
	unsigned int n_sec_decimals;	/* Number of digits in decimal seconds (0 for whole seconds) */
	unsigned int range;		/* 0 for 0/360, 1 for -360/0, 2 for -180/+180 */
	unsigned int wesn;		/* 1 if we want sign encoded with suffix W, E, S, N, 2 if also want space before letter */
	int order[3];			/* The relative order of degree, minute, seconds in form (-ve if unused) */
	bool decimal;			/* true if we want to use the D_FORMAT for decimal degrees only */
	bool no_sign;			/* true if we want absolute values (plot only) */
	char x_format[GMT_LEN64];	/* Actual C format used to plot/output longitude */
	char y_format[GMT_LEN64];	/* Actual C format used to plot/output latitude */
	char delimiter[2][2];		/* Delimiter strings in date, e.g. "-" */
};

/* Various ways to report longitudes */

enum GMT_lon_enum {
	GMT_IS_GIVEN_RANGE 			= 0,	/* Report lon as is */
	GMT_IS_0_TO_P360_RANGE			= 1,	/* Report 0 <= lon <= 360 */
	GMT_IS_0_TO_P360			= 2,	/* Report 0 <= lon < 360 */
	GMT_IS_M360_TO_0_RANGE			= 3,	/* Report -360 <= lon <= 0 */
	GMT_IS_M360_TO_0			= 4,	/* Report -360 < lon <= 0 */
	GMT_IS_M180_TO_P180_RANGE		= 5,	/* Report -180 <= lon <= +180 */
	GMT_IS_M180_TO_P180			= 6,	/* Report -180 <= lon < +180 */
	GMT_IS_M180_TO_P270_RANGE		= 7};	/* Report -180 <= lon < +270 [GSHHG only] */
#endif

