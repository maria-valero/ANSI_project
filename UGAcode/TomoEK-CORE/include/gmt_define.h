/*--------------------------------------------------------------------
 *	$Id: gmt_define.h 12822 2014-01-31 23:39:56Z remko $
 *
 *	Copyright (c) 1991-2014 by P. Wessel, W. H. F. Smith, R. Scharroo, J. Luis and F. Wobbe
 *	See LICENSE.TXT file for copying and redistribution conditions.
 *
 *	This program is free software; you can redistribute it and/or modify
 *	it under the terms of the GNU Lesser General Public License as published by
 *	the Free Software Foundation; version 3 or any later version.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU Lesser General Public License for more details.
 *
 *	Contact info: gmt.soest.hawaii.edu
 *--------------------------------------------------------------------*/

#ifndef _GMT_DEFINE_H
#define _GMT_DEFINE_H
#include "gmt_constants.h"

#ifndef isfinite
#	ifdef HAVE__FINITE /* only WIN32 */
#		define isfinite _finite
#	else
#		define isfinite(x) (!isinf(x) && !isnan(x))
#	endif
#endif

#	define PATH_SEPARATOR ';' /* Win uses ; while Unix uses : */

enum GMT_enum_api {
	GMT_USAGE	= 0,	/* Want to report full program usage message */
	GMT_SYNOPSIS	= 1,	/* Just want the synopsis of usage */
	GMT_STR16	= 16	/* Bytes needed to hold the @GMTAPI@-###### resource names */
};

/* These data primitive identifiers are as follows: */
enum GMT_enum_type {
	GMT_CHAR = 0,  /* int8_t, 1-byte signed integer type */
	GMT_UCHAR,     /* uint8_t, 1-byte unsigned integer type */
	GMT_SHORT,     /* int16_t, 2-byte signed integer type */
	GMT_USHORT,    /* uint16_t, 2-byte unsigned integer type */
	GMT_INT,       /* int32_t, 4-byte signed integer type */
	GMT_UINT,      /* uint32_t, 4-byte unsigned integer type */
	GMT_LONG,      /* int64_t, 8-byte signed integer type */
	GMT_ULONG,     /* uint64_t, 8-byte unsigned integer type */
	GMT_FLOAT,     /* 4-byte data float type */
	GMT_DOUBLE,    /* 8-byte data float type */
	GMT_TEXT,      /* Arbitrarily long text string [OGR/GMT use only] */
	GMT_DATETIME,  /* string with date/time info [OGR/GMT use only] */
	GMT_N_TYPES};  /* The number of supported data types above */

enum GMT_enum_annot {GMT_LINEAR = 0,
	GMT_LOG10,	/* These numbers are only used for GMT->current.proj.xyz_projection[3], */
	GMT_POW,	/* while GMT->current.proj.projection = 0 */
	GMT_TIME,
	GMT_ANNOT_CPT,
	GMT_CUSTOM};
/* Linear projections tagged 0-99 */
#define GMT_IS_LINEAR(C) (C->current.proj.projection / 100 == 0)
#define GMT_IS_ZERO(x) (fabs (x) < GMT_CONV_LIMIT)
/* Macros returns true if the two coordinates are lon/lat; way should be GMT_IN or GMT_OUT */
#define GMT_x_is_lon(C,way) (C->current.io.col_type[way][GMT_X] == GMT_IS_LON)
#define GMT_y_is_lat(C,way) (C->current.io.col_type[way][GMT_Y] == GMT_IS_LAT)
#define GMT_is_geographic(C,way) (GMT_x_is_lon(C,way) && GMT_y_is_lat(C,way))
//#define GMT_Report (V_API,level,format) {  printf (format);}
//#define GMT_Report (V_API,level,format, ...) \
   // printf( "%s ",format)

#endif /* _GMT_DEFINE_H */
