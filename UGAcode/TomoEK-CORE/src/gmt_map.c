/*--------------------------------------------------------------------
 *	$Id: gmt_map.c 12874 2014-02-08 02:29:35Z pwessel $
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
/*
 *			G M T _ M A P . C
 *
 *- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * GMT_map.c contains code related to generic coordinate transformation.
 * For the actual projection functions, see gmt_proj.c
 *
 *	Map Transformation Setup Routines
 *	These routines initializes the selected map transformation
 *	The names and main function are listed below
 *	NB! Note that the transformation function does not check that they are
 *	passed valid lon,lat numbers. I.e asking for log10 scaling using values
 *	<= 0 results in problems.
 *
 * The ellipsoid used is selectable by editing the gmt.conf in your
 * home directory.  If no such file, create one by running gmtdefaults.
 *
 * Usage: Initialize system by calling GMT_map_setup (separate module), and
 * then just use GMT_geo_to_xy() and GMT_xy_to_geo() functions.
 *
 * Author:	Paul Wessel
 * Date:	1-JAN-2010
 * Version:	5
 *
 *
 * PUBLIC GMT Functions include:
 *
 *	GMT_azim_to_angle :	Converts azimuth to angle on the map
 *	GMT_clip_to_map :	Force polygon points to be inside map
 *	GMT_compact_line :	Remove redundant pen movements
 *	GMT_geo_to_xy :		Generic lon/lat to x/y
 *	GMT_geo_to_xy_line :	Same for polygons
 *	GMT_geoz_to_xy :	Generic 3-D lon/lat/z to x/y
 *	GMT_grd_project :	Generalized grid projection with interpolation
 *	GMT_great_circle_dist :	Returns great circle distance in degrees
 *	GMT_img_project :	Generalized image projection with interpolation
 *	GMT_map_outside :	Generic function determines if we're outside map boundary
 *	GMT_map_path :		Return GMT_latpath or GMT_lonpath
 *	GMT_map_setup :		Initialize map projection
 *	GMT_project_init :	Initialize parameters for grid/image transformations
 *	GMT_xy_to_geo :		Generic inverse x/y to lon/lat projection
 *	GMT_xyz_to_xy :		Generic xyz to xy projection
 *	GMT_xyz_to_xy_n :	Same for an array
 *
 * Internal GMT Functions include:
 *
 *	gmt_get_origin :		Find origin of projection based on pole and 2nd point
 *	gmt_get_rotate_pole :		Find rotation pole based on two points on great circle
 *	gmt_ilinearxy :			Inverse linear projection
 *	gmt_init_three_D :		Initializes parameters needed for 3-D plots
 *	gmt_map_crossing :		Generic function finds crossings between line and map boundary
 *	GMT_latpath :			Return path between 2 points of equal latitide
 *	GMT_lonpath :			Return path between 2 points of equal longitude
 *	gmt_radial_crossing :		Determine map crossing in the Lambert azimuthal equal area projection
 *	GMT_left_boundary :		Return left boundary in x-inches
 *	gmt_linearxy :			Linear xy projection
 *	gmt_lon_inside :		Accounts for wrap-around in longitudes and checks for inside
 *	gmt_ellipse_crossing :		Find map crossings in the Mollweide projection
 *	gmt_move_to_rect :		Move an outside point straight in to nearest edge
 *	gmt_polar_outside :		Determines if a point is outside polar projection region
 *	gmt_pole_rotate_forward :	Compute positions from oblique coordinates
 *	gmt_radial_clip :		Clip path outside radial region
 *	gmt_radial_outside :		Determine if point is outside radial region
 *	gmt_radial_overlap :		Determine overlap, always true for his projection
 *	gmt_rect_clip :			Clip to rectangular region
 *	gmt_rect_crossing :		Find crossing between line and rect region
 *	gmt_rect_outside :		Determine if point is outside rect region
 *	gmt_rect_outside2 :		Determine if point is outside rect region (azimuthal proj only)
 *	gmt_rect_overlap :		Determine overlap between rect regions
 *	GMT_right_boundary :		Return x value of right map boundary
 *	gmt_xy_search :			Find xy map boundary
 *	GMT_wesn_clip:			Clip polygon to wesn boundaries
 *	gmt_wesn_crossing :		Find crossing between line and lon/lat rectangle
 *	gmt_wesn_outside :		Determine if a point is outside a lon/lat rectangle
 *	gmt_wesn_overlap :		Determine overlap between lon/lat rectangles
 *	gmt_wesn_search :		Search for extreme coordinates
 *	GMT_wrap_around_check_{x,tm} :	Check if line wraps around due to Greenwich
 *	GMT_x_to_xx :			Generic linear x projection
 *	GMT_xx_to_x :			Generic inverse linear x projection
 *	GMT_y_to_yy :			Generic linear y projection
 *	GMT_yy_to_y :			Generic inverse linear y projection
 *	GMT_z_to_zz :			Generic linear z projection
 *	GMT_zz_to_z :			Generic inverse linear z projection
 */


#include <stdbool.h>
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

#include "gmt_map.h"
#include "gmt_resources.h"
#include "gmt_constants.h"
#include "gmt_define.h"
#include "gmt_type.h"
#include "gmt_macros.h"

#include "io.h"



/* Basic error reporting when things go badly wrong. This Return macro can be
 * used in stead of regular return(code) to print out what the code is before
 * it returns.  We assume the GMT pointer is available in the function!
 */
#define VAR_TO_STR(arg)      #arg
#define Return(err) { GMT_Report(GMT->parent,GMT_MSG_NORMAL,"Internal Error = %s\n",VAR_TO_STR(err)); return (err);}


#define GMT_latg_to_latc(C,lat) GMT_lat_swap_quick (C, lat, C->current.proj.GMT_lat_swap_vals.c[GMT_LATSWAP_G2C])
#define GMT_latg_to_lata(C,lat) GMT_lat_swap_quick (C, lat, C->current.proj.GMT_lat_swap_vals.c[GMT_LATSWAP_G2A])
#define GMT_latc_to_latg(C,lat) GMT_lat_swap_quick (C, lat, C->current.proj.GMT_lat_swap_vals.c[GMT_LATSWAP_C2G])
#define GMT_lata_to_latg(C,lat) GMT_lat_swap_quick (C, lat, C->current.proj.GMT_lat_swap_vals.c[GMT_LATSWAP_A2G])

//#define GMT_x_is_lon(C,way) (C->current.io.col_type[way][GMT_X] == GMT_IS_LON)
//#define GMT_y_is_lat(C,way) (C->current.io.col_type[way][GMT_Y] == GMT_IS_LAT)
//#define GMT_is_geographic(C,way) (GMT_x_is_lon(C,way) && GMT_y_is_lat(C,way))


/* Note by P. Wessel, 18-Oct-2012:
 * In the olden days, GMT only did great circle distances.  In GMT 4 we implemented geodesic
 * distances by Rudoe's formula as given in Bomford [1971].  However, that geodesic is not
 * exactly what we wanted as it is a normal section and do not strictly follow the geodesic.
 * Other candidates are Vincenty [1975], which is widely used and Karney [2012], which is super-
 * accurate.  At this point their differences are in the micro-meter level.  For GMT 5 we have
 * now switched to the Vincenty algorithm as provided by Gerald Evenden, USGS [author of proj4],
 * which is a modified translation of the NGS algorithm and not exactly what is in proj4's geod
 * program (which Evenden thinks is inferior.)  I ran a comparison between many algorithms that
 * either were available via codes or had online calculators.  I sought the geodesic distance
 * from (0,0) to (10,10) on WGS-84; the results were (in meters):
 *
 *	GMT4 (Rudoe):		1565109.099232116
 *	proj4:			1565109.095557918
 *	vdist(0,0,10,10) [0]	1565109.09921775
 *	Karney [1]: 		1565109.09921789
 *	Vincenty [2]:		1565109.099218036
 *	NGS [3]			1565109.0992
 *
 * [0] via Joaquim Luis, supposedly Vincenty [2012]
 * [1] via online calculator at max precision http://geographiclib.sourceforge.net/cgi-bin/Geod
 * [2] downloading, compiling and running http://article.gmane.org/gmane.comp.gis.proj-4.devel/3478.
 *     This is not identical to Vincenty in proj4 but written by Evenden (proj.4 author)
 * [3] via online calculator http://www.ngs.noaa.gov/cgi-bin/Inv_Fwd/inverse2.prl. Their notes says
 *     this is Vincenty; unfortunately I cannot control the output precision.
 *
 * Based on these comparisons we decided to implement the Vincenty [2] code as given.  The older Rudoe
 * code remains in this file for reference.  The define of USE_VINCENTY below selects the new Vincenty code.
 * The choice was based on the readily available C code versus having to reimplement Karney in C.
 */

#define USE_VINCENTY 1	/* New GMT 5 behavior */

double gmt_get_angle (struct GMT_CTRL *GMT, double lon1, double lat1, double lon2, double lat2);

#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif
#ifndef R2D
#define R2D (180.0/M_PI)
#endif

double GMT_lat_swap_quick (struct GMT_CTRL *GMT, double lat, double c[])
{
	/* Return latitude, in degrees, given latitude, in degrees, based on coefficients c */

	double delta, cos2phi, sin2phi;

	/* First deal with trivial cases */

	if (lat >=  90.0) return ( 90.0);
	if (lat <= -90.0) return (-90.0);
	if (GMT_IS_ZERO (lat)) return (0.0);

	sincosd (2.0 * lat, &sin2phi, &cos2phi);

	delta = sin2phi * (c[0] + cos2phi * (c[1] + cos2phi * (c[2] + cos2phi * c[3])));

	return (lat + R2D * delta);
}

double gmt_flatearth_dist_degree (struct GMT_CTRL *GMT, double x0, double y0, double x1, double y1)
{
	/* Calculates the approximate flat earth distance in degrees.
	   If difference in longitudes exceeds 180 we pick the other
	   offset (360 - offset)
	 */
	double dlon;
	
	GMT_set_delta_lon (x0, x1, dlon);
	return (hypot ( dlon * cosd (0.5 * (y1 + y0)), (y1 - y0)));
}

double gmt_flatearth_dist_meter (struct GMT_CTRL *GMT, double x0, double y0, double x1, double y1)
{
	/* Calculates the approximate flat earth distance in km.
	   If difference in longitudes exceeds 180 we pick the other
	   offset (360 - offset)
	 */
	return (gmt_flatearth_dist_degree (GMT, x0, y0, x1, y1) * GMT->current.proj.DIST_M_PR_DEG);
}

double gmt_az_backaz_flatearth (struct GMT_CTRL *GMT, double lonE, double latE, double lonS, double latS, bool baz)
{
	/* Calculate azimuths or backazimuths.  Flat earth code.
	 * First point is considered "Event" and second "Station".
	 * Azimuth is direction from Station to Event.
	 * BackAzimuth is direction from Event to Station */

	double az, dx, dy, dlon;

	if (baz) {	/* exchange point one and two */
		double_swap (lonS, lonE);
		double_swap (latS, latE);
	}
	GMT_set_delta_lon (lonS, lonE, dlon);
	dx = dlon * cosd (0.5 * (latE + latS));
	dy = latE - latS;
	az = (dx == 0.0 && dy == 0.0) ? GMT->session.d_NaN : 90.0 - atan2d (dy, dx);
	if (az < 0.0) az += 360.0;
	return (az);
}

double gmt_geodesic_dist_meter (struct GMT_CTRL *GMT, double lonS, double latS, double lonE, double latE)
{
	/* Compute length of geodesic between locations in meters
	 * We use Rudoe's equation from Bomford.
	 */

	double e1, el, sinthi, costhi, sinthk, costhk, tanthi, tanthk, sina12, cosa12, d_lon;
	double al, a12top, a12bot, a12, e2, e3, c0, c2, c4, v1, v2, z1, z2, x2, y2, dist;
	double e1p1, sqrte1p1, sin_dl, cos_dl, u1bot, u1, u2top, u2bot, u2, b0, du, pdist;


	/* Equations are unstable for latitudes of exactly 0 degrees. */

	if (latE == 0.0) latE = 1.0e-08;
	if (latS == 0.0) latS = 1.0e-08;

	/* Now compute the distance between the two points using Rudoe's
	 * formula given in Bomford's GEODESY, section 2.15(b).
	 * (Unclear if it is 1971 or 1980 edition)
	 * (There is some numerical problem with the following formulae.
	 * If the station is in the southern hemisphere and the event in
	 * in the northern, these equations give the longer, not the
	 * shorter distance between the two locations.  Since the equations
	 * are fairly messy, the simplist solution is to reverse the
	 * meanings of the two locations for this case.)
	 */

	if (latS < 0.0) {	/* Station in southern hemisphere, swap */
		double_swap (lonS, lonE);
		double_swap (latS, latE);
	}
	el = GMT->current.proj.ECC2 / GMT->current.proj.one_m_ECC2;
	e1 = 1.0 + el;
	sincosd (latE, &sinthi, &costhi);
	sincosd (latS, &sinthk, &costhk);
	GMT_set_delta_lon (lonS, lonE, d_lon);
	sincosd (d_lon, &sin_dl, &cos_dl);
	tanthi = sinthi / costhi;
	tanthk = sinthk / costhk;
	al = tanthi / (e1 * tanthk) + GMT->current.proj.ECC2 * sqrt ((e1 + tanthi * tanthi) / (e1 + tanthk * tanthk));
	a12top = sin_dl;
	a12bot = (al - cos_dl) * sinthk;
	a12 = atan2 (a12top,a12bot);
	sincos (a12, &sina12, &cosa12);
	e1 = el * (pow (costhk * cosa12, 2.0) + sinthk * sinthk);
	e2 = e1 * e1;
	e3 = e1 * e2;
	c0 = 1.0 + 0.25 * e1 - (3.0 / 64.0) * e2 + (5.0 / 256.0) * e3;
	c2 = -0.125 * e1 + (1.0 / 32) * e2 - (15.0 / 1024.0) * e3;
	c4 = -(1.0 / 256.0) * e2 + (3.0 / 1024.0) * e3;
	v1 = GMT->current.proj.EQ_RAD / sqrt (1.0 - GMT->current.proj.ECC2 * sinthk * sinthk);
	v2 = GMT->current.proj.EQ_RAD / sqrt (1.0 - GMT->current.proj.ECC2 * sinthi * sinthi);
	z1 = v1 * (1.0 - GMT->current.proj.ECC2) * sinthk;
	z2 = v2 * (1.0 - GMT->current.proj.ECC2) * sinthi;
	x2 = v2 * costhi * cos_dl;
	y2 = v2 * costhi * sin_dl;
	e1p1 = e1 + 1.0;
	sqrte1p1 = sqrt (e1p1);
	u1bot = sqrte1p1 * cosa12;
	u1 = atan2 (tanthk, u1bot);
	u2top = v1 * sinthk + e1p1 * (z2 - z1);
	u2bot = sqrte1p1 * (x2 * cosa12 - y2 * sinthk * sina12);
	u2 = atan2 (u2top, u2bot);
	b0 = v1 * sqrt (1.0 + el * pow (costhk * cosa12, 2.0)) / e1p1;
	du = u2  - u1;
	if (fabs (du) > M_PI) du = copysign (TWO_PI - fabs (du), du);
	pdist = b0 * (c2 * (sin (2.0 * u2) - sin(2.0 * u1)) + c4 * (sin (4.0 * u2) - sin (4.0 * u1)));
	dist = fabs (b0 * c0 * du + pdist);

	return (dist);
}

double GMT_lat_swap (struct GMT_CTRL *GMT, double lat, unsigned int itype)
{
	/* Return latitude, in degrees, given latitude, in degrees, based on itype */

	double delta, cos2phi, sin2phi;

	/* First deal with trivial cases */

	if (lat >=  90.0) return ( 90.0);
	if (lat <= -90.0) return (-90.0);
	if (GMT_IS_ZERO (lat)) return (0.0);

	if (GMT->current.proj.GMT_lat_swap_vals.spherical) return (lat);

	if (itype >= GMT_LATSWAP_N) {
		/* This should never happen -?- or do we want to allow the
			possibility of using itype = -1 to do nothing  */
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "GMT_lat_swap(): Invalid choice, programming bug.\n");
		return(lat);
	}

	sincosd (2.0 * lat, &sin2phi, &cos2phi);

	delta = sin2phi * (GMT->current.proj.GMT_lat_swap_vals.c[itype][0]
		+ cos2phi * (GMT->current.proj.GMT_lat_swap_vals.c[itype][1]
		+ cos2phi * (GMT->current.proj.GMT_lat_swap_vals.c[itype][2]
		+ cos2phi * GMT->current.proj.GMT_lat_swap_vals.c[itype][3])));

	return (lat + R2D * delta);
}

double gmt_haversine (struct GMT_CTRL *GMT, double lon1, double lat1, double lon2, double lat2)
{
	/* Haversine formula for great circle distance.  Intermediate function that returns sin^2 (half_angle).
	 * This avoids problems with short distances where cos(c) is close to 1 and acos is inaccurate.
	 */

	double sx, sy, sin_half_squared;

	if (lat1 == lat2 && lon1 == lon2) return (0.0);

	if (GMT->current.setting.proj_aux_latitude != GMT_LATSWAP_NONE) {	/* Use selected auxiliary latitude */
		lat1 = GMT_lat_swap (GMT, lat1, GMT->current.setting.proj_aux_latitude);
		lat2 = GMT_lat_swap (GMT, lat2, GMT->current.setting.proj_aux_latitude);
	}

	sy = sind (0.5 * (lat2 - lat1));
	sx = sind (0.5 * (lon2 - lon1));	/* If there is a 360 wrap here then the sign of sx is wrong but we only use sx^2 */
	sin_half_squared = sy * sy + cosd (lat2) * cosd (lat1) * sx * sx;

	return (sin_half_squared);
}

double GMT_great_circle_dist_degree (struct GMT_CTRL *GMT, double lon1, double lat1, double lon2, double lat2)
{	/* Great circle distance on a sphere in degrees */

	double sin_half_squared = gmt_haversine (GMT, lon1, lat1, lon2, lat2);
	return (2.0 * d_asind (d_sqrt (sin_half_squared)));
}

double GMT_great_circle_dist_meter (struct GMT_CTRL *GMT, double lon1, double lat1, double lon2, double lat2)
{	/* Calculates the great circle distance in meter */
	return (GMT_great_circle_dist_degree (GMT, lon1, lat1, lon2, lat2) * GMT->current.proj.DIST_M_PR_DEG);
}

double gmt_az_backaz_sphere (struct GMT_CTRL *GMT, double lonE, double latE, double lonS, double latS, bool baz)
{
	/* Calculate azimuths or backazimuths.  Spherical code.
	 * First point is considered "Event" and second "Station".
	 * Azimuth is direction from Station to Event.
	 * BackAzimuth is direction from Event to Station */

	double az, sin_yS, cos_yS, sin_yE, cos_yE, sin_dlon, cos_dlon;

	if (baz) {	/* exchange point one and two */
		double_swap (lonS, lonE);
		double_swap (latS, latE);
	}
	sincosd (latS, &sin_yS, &cos_yS);
	sincosd (latE, &sin_yE, &cos_yE);
	sincosd (lonS - lonE, &sin_dlon, &cos_dlon);
	az = atan2d (cos_yS * sin_dlon, cos_yE * sin_yS - sin_yE * cos_yS * cos_dlon);
	if (az < 0.0) az += 360.0;
	return (az);
}

#define GEOD_TEXT "Rudoe"
double gmt_az_backaz_geodesic (struct GMT_CTRL *GMT, double lonE, double latE, double lonS, double latS, bool baz)
{
	/* Calculate azimuths or backazimuths for geodesics using geocentric latitudes.
	 * First point is considered "Event" and second "Station".
	 * Azimuth is direction from Station to Event.
	 * BackAzimuth is direction from Event to Station */

	double az, a, b, c, d, e, f, g, h, a1, b1, c1, d1, e1, f1, g1, h1, thg, ss, sc;

	/* Equations are unstable for latitudes of exactly 0 degrees. */
	if (latE == 0.0) latE = 1.0e-08;
	if (latS == 0.0) latS = 1.0e-08;

	/* Must convert from geographic to geocentric coordinates in order
	 * to use the spherical trig equations.  This requires a latitude
	 * correction given by: 1-ECC2=1-2*f + f*f = GMT->current.proj.one_m_ECC2
	 */

	thg = atan (GMT->current.proj.one_m_ECC2 * tand (latE));
	sincos (thg, &c, &f);		f = -f;
	sincosd (lonE, &d, &e);		e = -e;
	a = f * e;
	b = -f * d;
	g = -c * e;
	h = c * d;

	/* Calculating some trig constants. */

	thg = atan (GMT->current.proj.one_m_ECC2 * tand (latS));
	sincos (thg, &c1, &f1);		f1 = -f1;
	sincosd (lonS, &d1, &e1);	e1 = -e1;
	a1 = f1 * e1;
	b1 = -f1 * d1;
	g1 = -c1 * e1;
	h1 = c1 * d1;

	/* Spherical trig relationships used to compute angles. */

	if (baz) {	/* Get Backazimuth */
		ss = pow(a-d1,2.0) + pow(b-e1,2.0) + c * c - 2.0;
		sc = pow(a-g1,2.0) + pow(b-h1,2.0) + pow(c-f1,2.0) - 2.0;
	}
	else {		/* Get Azimuth */
		ss = pow(a1-d, 2.0) + pow(b1-e, 2.0) + c1 * c1 - 2.0;
		sc = pow(a1-g, 2.0) + pow(b1-h, 2.0) + pow(c1-f, 2.0) - 2.0;
	}
	az = atan2d (ss,sc);
	if (az < 0.0) az += 360.0;
	return (az);
}

double GMT_great_circle_dist_cos (struct GMT_CTRL *GMT, double lon1, double lat1, double lon2, double lat2)
{
	/* Return cosine of great circle distance */

	double sin_half_squared = gmt_haversine (GMT, lon1, lat1, lon2, lat2);
	return (1.0 - 2.0 * sin_half_squared);	/* Convert sin^2 (half-angle) to cos (angle) */
}

double gmt_geodesic_dist_degree (struct GMT_CTRL *GMT, double lonS, double latS, double lonE, double latE)
{
	/* Compute the great circle arc length in degrees on an ellipsoidal
	 * Earth.  We do this by converting to geocentric coordinates.
	 */

	double a, b, c, d, e, f, a1, b1, c1, d1, e1, f1, thg, sc, sd, dist;

	/* Equations are unstable for latitudes of exactly 0 degrees. */

	if (latE == 0.0) latE = 1.0e-08;
	if (latS == 0.0) latS = 1.0e-08;

	/* Must convert from geographic to geocentric coordinates in order
	 * to use the spherical trig equations.  This requires a latitude
	 * correction given by: 1-ECC2=1-2*F + F*F = GMT->current.proj.one_m_ECC2
	 */

	thg = atan (GMT->current.proj.one_m_ECC2 * tand (latE));
	sincos (thg, &c, &f);		f = -f;
	sincosd (lonE, &d, &e);		e = -e;
	a = f * e;
	b = -f * d;

	/* Calculate some trig constants. */

	thg = atan (GMT->current.proj.one_m_ECC2 * tand (latS));
	sincos (thg, &c1, &f1);		f1 = -f1;
	sincosd (lonS, &d1, &e1);	e1 = -e1;
	a1 = f1 * e1;
	b1 = -f1 * d1;

	/* Spherical trig relationships used to compute angles. */

	sc = a * a1 + b * b1 + c * c1;
	sd = 0.5 * sqrt ((pow (a-a1,2.0) + pow (b-b1,2.0) + pow (c-c1,2.0)) * (pow (a+a1,2.0) + pow (b+b1, 2.0) + pow (c+c1, 2.0)));
	dist = atan2d (sd, sc);
	if (dist < 0.0) dist += 360.0;

	return (dist);
}

double GMT_geodesic_dist_cos (struct GMT_CTRL *GMT, double lonS, double latS, double lonE, double latE)
{	/* Convenience function to get cosine instead */
	return (cosd (gmt_geodesic_dist_degree (GMT, lonS, latS, lonE, latE)));
}

double gmt_loxodrome_dist_degree (struct GMT_CTRL *GMT, double lon1, double lat1, double lon2, double lat2)
{	/* Calculates the distance along the loxodrome, in meter */
	double dist, d_lon;
	GMT_set_delta_lon (lon1, lon2, d_lon);
	if (doubleAlmostEqualZero (lat1, lat2)) {	/* Along parallel */
		if (GMT->current.proj.GMT_convert_latitudes) lat1 = GMT_latg_to_latc (GMT, lat1);
		dist = fabs (d_lon) * cosd (lat1);
	}
	else { /* General case */
		double dx, dy, Az;
		if (GMT->current.proj.GMT_convert_latitudes) {
			lat1 = GMT_latg_to_latc (GMT, lat1);
			lat2 = GMT_latg_to_latc (GMT, lat2);
		}
		dx = D2R * d_lon;
		dy = d_log (GMT, tand (45.0 + 0.5 * lat2)) - d_log (GMT, tand (45.0 + 0.5 * lat1));
		Az = atan2 (dx, dy);
		dist = fabs (dx / cos (Az));
	}
	return (dist);
}

double gmt_loxodrome_dist_meter (struct GMT_CTRL *GMT, double lon1, double lat1, double lon2, double lat2)
{	/* Calculates the loxodrome distance in meter */
	return (gmt_loxodrome_dist_degree (GMT, lon1, lat1, lon2, lat2) * GMT->current.proj.DIST_M_PR_DEG);
}

double gmt_az_backaz_loxodrome (struct GMT_CTRL *GMT, double lonE, double latE, double lonS, double latS, bool baz)
{
	/* Calculate azimuths or backazimuths.  Loxodrome mode.
	 * First point is considered "Event" and second "Station".
	 * Azimuth is direction from Station to Event.
	 * BackAzimuth is direction from Event to Station */

	double az, d_lon;

	if (baz) {	/* exchange point one and two */
		double_swap (lonS, lonE);
		double_swap (latS, latE);
	}
	GMT_set_delta_lon (lonE, lonS, d_lon);
	if (doubleAlmostEqualZero (latS, latE))	/* Along parallel */
		az = (d_lon > 0.0) ? 90 : -90.0;
	else { /* General case */
		double dx, dy;
		if (GMT->current.proj.GMT_convert_latitudes) {
			latS = GMT_latg_to_latc (GMT, latS);
			latE = GMT_latg_to_latc (GMT, latE);
		}
		dx = D2R * d_lon;
		dy = d_log (GMT, tand (45.0 + 0.5 * latS)) - d_log (GMT, tand (45.0 + 0.5 * latE));
		az = atan2d (dx, dy);
		if (az < 0.0) az += 360.0;
	}
	return (az);
}

void GMT_geo_to_cart (struct GMT_CTRL *GMT, double lat, double lon, double *a, bool degrees)
{
	/* Convert geographic latitude and longitude (lat, lon)
	   to a 3-vector of unit length (a). If degrees = true,
	   input coordinates are in degrees, otherwise in radian */

	double clat, clon, slon;

	if (degrees) {
		lat *= D2R;
		lon *= D2R;
	}
	sincos (lat, &a[GMT_Z], &clat);
	sincos (lon, &slon, &clon);
	a[GMT_X] = clat * clon;
	a[GMT_Y] = clat * slon;
}

double GMT_distance_type (struct GMT_CTRL *GMT, double lonS, double latS, double lonE, double latE, unsigned int id)
{	/* Generic function available to programs for contour/label distance calculations */
	return (GMT->current.map.dist[id].scale * GMT->current.map.dist[id].func (GMT, lonS, latS, lonE, latE));
}

double GMT_distance (struct GMT_CTRL *GMT, double lonS, double latS, double lonE, double latE)
{	/* Generic function available to programs */
	return (GMT_distance_type (GMT, lonS, latS, lonE, latE, 0));
}

void GMT_cross3v (struct GMT_CTRL *GMT, double *a, double *b, double *c)
{
	c[GMT_X] = a[GMT_Y] * b[GMT_Z] - a[GMT_Z] * b[GMT_Y];
	c[GMT_Y] = a[GMT_Z] * b[GMT_X] - a[GMT_X] * b[GMT_Z];
	c[GMT_Z] = a[GMT_X] * b[GMT_Y] - a[GMT_Y] * b[GMT_X];
}

double GMT_mag3v (struct GMT_CTRL *GMT, double *a)
{
	return (d_sqrt(a[GMT_X]*a[GMT_X] + a[GMT_Y]*a[GMT_Y] + a[GMT_Z]*a[GMT_Z]));
}

void GMT_normalize3v (struct GMT_CTRL *GMT, double *a)
{
	double r_length;
	r_length = GMT_mag3v (GMT,a);
	if (r_length != 0.0) {
		r_length = 1.0 / r_length;
		a[GMT_X] *= r_length;
		a[GMT_Y] *= r_length;
		a[GMT_Z] *= r_length;
	}
}

double GMT_dot3v (struct GMT_CTRL *GMT, double *a, double *b)
{
	return (a[GMT_X]*b[GMT_X] + a[GMT_Y]*b[GMT_Y] + a[GMT_Z]*b[GMT_Z]);
}

int GMT_great_circle_intersection (struct GMT_CTRL *GMT, double A[], double B[], double C[], double X[], double *CX_dist)
{
	/* A, B, C are 3-D Cartesian unit vectors, i.e., points on the sphere.
	 * Let points A and B define a great circle, and consider a
	 * third point C.  A second great cirle goes through C and
	 * is orthogonal to the first great circle.  Their intersection
	 * X is the point on (A,B) closest to C.  We must test if X is
	 * between A,B or outside.
	 */
	unsigned int i;
	double P[3], E[3], M[3], Xneg[3], cos_AB, cos_MX1, cos_MX2, cos_test;

	GMT_cross3v (GMT, A, B, P);			/* Get pole position of plane through A and B (and origin O) */
	GMT_normalize3v (GMT, P);			/* Make sure P has unit length */
	GMT_cross3v (GMT, C, P, E);			/* Get pole E to plane through C (and origin) but normal to A,B (hence going through P) */
	GMT_normalize3v (GMT, E);			/* Make sure E has unit length */
	GMT_cross3v (GMT, P, E, X);			/* Intersection between the two planes is oriented line*/
	GMT_normalize3v (GMT, X);			/* Make sure X has unit length */
	/* The X we want could be +x or -X; must determine which might be closest to A-B midpoint M */
	for (i = 0; i < 3; i++) {
		M[i] = A[i] + B[i];
		Xneg[i] = -X[i];
	}
	GMT_normalize3v (GMT, M);			/* Make sure M has unit length */
	/* Must first check if X is along the (A,B) segment and not on its extension */

	cos_MX1 = GMT_dot3v (GMT, M, X);		/* Cos of spherical distance between M and +X */
	cos_MX2 = GMT_dot3v (GMT, M, Xneg);	/* Cos of spherical distance between M and -X */
	if (cos_MX2 > cos_MX1) GMT_memcpy (X, Xneg, 3, double);		/* -X is closest to A-B midpoint */
	cos_AB = fabs (GMT_dot3v (GMT, A, B));	/* Cos of spherical distance between A,B */
	cos_test = fabs (GMT_dot3v (GMT, A, X));	/* Cos of spherical distance between A and X */
	if (cos_test < cos_AB) return 1;	/* X must be on the A-B extension if its distance to A exceeds the A-B length */
	cos_test = fabs (GMT_dot3v (GMT, B, X));	/* Cos of spherical distance between B and X */
	if (cos_test < cos_AB) return 1;	/* X must be on the A-B extension if its distance to B exceeds the A-B length */

	/* X is between A and B.  Now calculate distance between C and X */

	*CX_dist = GMT_dot3v (GMT, C, X);		/* Cos of spherical distance between C and X */
	return (0);				/* Return zero if intersection is between A and B */
}

void GMT_cart_to_geo (struct GMT_CTRL *GMT, double *lat, double *lon, double *a, bool degrees)
{
	/* Convert a 3-vector (a) of unit length into geographic
	   coordinates (lat, lon). If degrees = true, the output coordinates
	   are in degrees, otherwise in radian. */

	if (degrees) {
		*lat = d_asind (a[GMT_Z]);
		*lon = d_atan2d (a[GMT_Y], a[GMT_X]);
	}
	else {
		*lat = d_asin (a[GMT_Z]);
		*lon = d_atan2 (a[GMT_Y], a[GMT_X]);
	}
}

bool gmt_near_a_line_spherical (struct GMT_CTRL *P, double lon, double lat, uint64_t seg, struct GMT_DATASEGMENT *S, unsigned int return_mindist, double *dist_min, double *x_near, double *y_near)
{
	bool perpendicular_only = false, interior, within;
	uint64_t row, prev_row;
	double d, A[3], B[3], GMT[3], X[3], xlon, xlat, cx_dist, cos_dist, dist_AB, fraction;

	/* gmt_near_a_line_spherical works in one of two modes, depending on return_mindist.
	   Since return_mindist is composed of two settings we must first set
	   perpendicular_only = (return_mindist >= 10);
	   return_mindist -= 10 * perpendicular_only;
	   That is, if 10 was added it means perpendicular_only is set and then the 10 is
	   removed.  We now consider what is left of return_mindist:
	   (1) return_mindist == 0:
	      We expect each segment to have its dist variable set to a minimum distance,
	      and if the point is within this distance from the line then we return true;
	      otherwise we return false.  If the segments have not set their distances then
	      it will have been initialized at 0 and only a point on the line will return true.
	      If perpendicular_only we ignore a point that is within the distance of the
	      linesegment endpoints but project onto the extension of the line (i.e., it is
	      "outside" the extent of the line).  We return false in that case.
	   (2) return_mindist != 0:
	      Return the minimum distance via dist_min. In addition, if > 1:
	      If == 2 we also return the coordinate of nearest point via x_near, y_near.
	      If == 3 we instead return segment number and point number (fractional) of that point via x_near, y_near.
	      The function will always return true, except if perpendicular_only is set: then we
	      return false if the point projects onto the extension of the line (i.e., it is "outside"
	      the extent of the line).  */

	if (return_mindist >= 10) {	/* Exclude (or flag) circular region surrounding line endpoints */
		perpendicular_only = true;
		return_mindist -= 10;
	}
	GMT_geo_to_cart (P, lat, lon, GMT, true);	/* Our point to test is now GMT */

	if (S->n_rows <= 0) return (false);	/* Empty ; skip */

	/* Find nearest point on this line */

	if (return_mindist) S->dist = 0.0;	/* Explicitly set dist to zero so the shortest distance can be found */

	for (row = 0; row < S->n_rows; row++) {	/* loop over nodes on current line */
		d = GMT_distance (P, lon, lat, S->coord[GMT_X][row], S->coord[GMT_Y][row]);	/* Distance between our point and row'th node on seg'th line */
		if (return_mindist && d < (*dist_min)) {	/* Update minimum distance */
			*dist_min = d;
			if (return_mindist == 2) *x_near = S->coord[GMT_X][row], *y_near = S->coord[GMT_Y][row];	/* Also update (x,y) of nearest point on the line */
			if (return_mindist == 3) *x_near = (double)seg, *y_near = (double)row;	/* Also update (seg, pt) of nearest point on the line */
		}
		interior = (row > 0 && row < (S->n_rows - 1));	/* Only false if we are processing one of the end points */
		if (d <= S->dist && (interior || !perpendicular_only)) return (true);			/* Node inside the critical distance; we are done */
	}

	if (S->n_rows < 2) return (false);	/* 1-point "line" is a point; skip segment check */

	/* If we get here we must check for intermediate points along the great circle lines between segment nodes.*/

	if (return_mindist)		/* Cosine of the great circle distance we are checking for. 2 ensures failure to be closer */
		cos_dist = 2.0;
	else if (P->current.map.dist[GMT_MAP_DIST].arc)	/* Used angular distance measure */
		cos_dist = cosd (S->dist / P->current.map.dist[GMT_MAP_DIST].scale);
	else	/* Used distance units (e.g., meter, km). Conv to meters, then to degrees */
		cos_dist = cosd ((S->dist / P->current.map.dist[GMT_MAP_DIST].scale) / P->current.proj.DIST_M_PR_DEG);
	GMT_geo_to_cart (P, S->coord[GMT_Y][0], S->coord[GMT_X][0], B, true);		/* 3-D vector of end of last segment */

	for (row = 1, within = false; row < S->n_rows; row++) {				/* loop over great circle segments on current line */
		GMT_memcpy (A, B, 3, double);	/* End of last segment is start of new segment */
		GMT_geo_to_cart (P, S->coord[GMT_Y][row], S->coord[GMT_X][row], B, true);	/* 3-D vector of end of this segment */
		if (GMT_great_circle_intersection (P, A, B, GMT, X, &cx_dist)) continue;	/* X not between A and B */
		if (return_mindist) {		/* Get lon, lat of X, calculate distance, and update min_dist if needed */
			GMT_cart_to_geo (P, &xlat, &xlon, X, true);
			d = GMT_distance (P, xlon, xlat, lon, lat);	/* Distance between our point and closest perpendicular point on seg'th line */
			if (d < (*dist_min)) {	/* Update minimum distance */
				*dist_min = d;
				if (return_mindist == 2) { *x_near = xlon; *y_near = xlat;}	/* Also update (x,y) of nearest point on the line */
				else if (return_mindist == 3) {	/* Also update (seg, pt) of nearest point on the line */
					*x_near = (double)seg;
					prev_row = row - 1;
					dist_AB = GMT_distance (P, S->coord[GMT_X][prev_row], S->coord[GMT_Y][prev_row], S->coord[GMT_X][row], S->coord[GMT_Y][row]);
					fraction = (dist_AB > 0.0) ? GMT_distance (P, S->coord[GMT_X][prev_row], S->coord[GMT_Y][prev_row], xlon, xlat) / dist_AB : 0.0;
					*y_near = (double)prev_row + fraction;
				}
				within = true;	/* Found at least one segment with a valid inside distance */
			}
		}
		if (cx_dist >= cos_dist) return (true);	/* X is on the A-B extension AND within specified distance */
	}

	return (within);	/* All tests failed, we are not close to the line(s), or we return a mindist (see comments above) */
}

bool gmt_near_lines_spherical (struct GMT_CTRL *P, double lon, double lat, struct GMT_DATATABLE *T, unsigned int return_mindist, double *dist_min, double *x_near, double *y_near)
{
	uint64_t seg;
	int mode = return_mindist, status;
	bool OK = false;
	if (mode >= 10) mode -= 10;	/* Exclude (or flag) circular region surrounding line endpoints */
	if (mode) *dist_min = DBL_MAX;	/* Want to find the minimum distance so init to huge */

	for (seg = 0; seg < T->n_segments; seg++) {	/* Loop over each line segment */
		status = gmt_near_a_line_spherical (P, lon, lat, seg, T->segment[seg], return_mindist, dist_min, x_near, y_near);
		if (status) {	/* Got a min distance or satisfied the min dist requirement */
			if (!return_mindist) return (true);	/* Done, we are within distance of one of the lines */
			OK = true;
		}
	}
	return (OK);
}

bool gmt_near_a_point_spherical (struct GMT_CTRL *GMT, double x, double y, struct GMT_DATATABLE *T, double dist)
{
	uint64_t row, seg;
	bool each_point_has_distance;
	double d;

	each_point_has_distance = (dist <= 0.0 && T->segment[0]->n_columns > 2);
	for (seg = 0; seg < T->n_segments; seg++) {
		for (row = 0; row < T->segment[seg]->n_rows; row++) {
			d = GMT_distance (GMT, x, y, T->segment[seg]->coord[GMT_X][row], T->segment[seg]->coord[GMT_Y][row]);
			if (each_point_has_distance) dist = T->segment[seg]->coord[GMT_Z][row];
			if (d <= dist) return (true);
		}
	}
	return (false);
}

void gmt_set_distaz (struct GMT_CTRL *GMT, unsigned int mode, unsigned int type)
{	/* Assigns pointers to the chosen distance and azimuth functions */
	char *type_name[3] = {"Map", "Contour", "Contour annotation"};
	char *aux[6] = {"no", "authalic", "conformal", "meridional", "geocentric", "parametric"};
	char *rad[5] = {"mean (R_1)", "authalic (R_2)", "volumetric (R_3)", "meridional", "quadratic"};
	int choice = (GMT->current.setting.proj_aux_latitude == GMT_LATSWAP_NONE) ? 0 : 1 + GMT->current.setting.proj_aux_latitude/2;
	GMT->current.map.dist[type].scale = 1.0;	/* Default scale */


	switch (mode) {	/* Set pointers to distance functions */
		//case GMT_CARTESIAN_DIST:	/* Cartesian 2-D x,y data */
			//GMT->current.map.dist[type].func = &GMT_cartesian_dist;
			//GMT->current.map.azimuth_func = &gmt_az_backaz_cartesian;
			//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "%s distance calculation will be Cartesian\n", type_name[type]);
			//break;
		//case GMT_CARTESIAN_DIST2:	/* Cartesian 2-D x,y data, use r^2 instead of hypot */
			//GMT->current.map.dist[type].func = &GMT_cartesian_dist2;
			//GMT->current.map.azimuth_func = &gmt_az_backaz_cartesian;
			//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "%s distance calculation will be Cartesian\n", type_name[type]);
			//break;
		//case GMT_CARTESIAN_DIST_PROJ:	/* Cartesian distance after projecting 2-D lon,lat data */
			//GMT->current.map.dist[type].func = &GMT_cartesian_dist_proj;
			//GMT->current.map.azimuth_func = &gmt_az_backaz_cartesian_proj;
			//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "%s distance calculation will be Cartesian after first projecting via -J\n", type_name[type]);
			//break;
		//case GMT_CARTESIAN_DIST_PROJ2:	/* Cartesian distance after projecting 2-D lon,lat data, use r^2 instead of hypot  */
			//GMT->current.map.dist[type].func = &GMT_cartesian_dist_proj2;
			//GMT->current.map.azimuth_func = &gmt_az_backaz_cartesian_proj;
			//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "%s distance calculation will be Cartesian after first projecting via -J\n", type_name[type]);
			//break;
		case GMT_DIST_M+GMT_FLATEARTH:	/* 2-D lon, lat data, but scale to Cartesian flat earth in meter */
			GMT->current.map.dist[type].func = &gmt_flatearth_dist_meter;
			GMT->current.map.azimuth_func  = &gmt_az_backaz_flatearth;
			//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "%s distance calculation will be Flat Earth in meters\n", type_name[type]);
			break;
		case GMT_DIST_M+GMT_GREATCIRCLE:	/* 2-D lon, lat data, use spherical distances in meter */

			GMT->current.map.dist[type].func = &GMT_great_circle_dist_meter;
			GMT->current.map.azimuth_func = &gmt_az_backaz_sphere;
			//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "%s distance calculation will be using great circle approximation with %s auxiliary latitudes and %s radius = %.4f m.\n",
				//type_name[type], aux[choice], rad[GMT->current.setting.proj_mean_radius], GMT->current.proj.mean_radius);
			break;
		case GMT_DIST_M+GMT_GEODESIC:	/* 2-D lon, lat data, use geodesic distances in meter */
			GMT->current.map.dist[type].func = &gmt_geodesic_dist_meter;
			GMT->current.map.azimuth_func = &gmt_az_backaz_geodesic;
			//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "%s distance calculation will be using %s geodesics in meters\n", type_name[type], GEOD_TEXT);
			break;
		case GMT_DIST_DEG+GMT_FLATEARTH:	/* 2-D lon, lat data, use Flat Earth distances in degrees */
			//GMT->current.map.dist[type].func = gmt_flatearth_dist_degree;
			//GMT->current.map.azimuth_func = &gmt_az_backaz_flatearth;
			//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "%s distance calculation will be Flat Earth in degrees\n", type_name[type]);
			//break;
		case GMT_DIST_DEG+GMT_GREATCIRCLE:	/* 2-D lon, lat data, use spherical distances in degrees */
			//GMT->current.map.dist[type].func = &GMT_great_circle_dist_degree;
			//GMT->current.map.azimuth_func = &gmt_az_backaz_sphere;
			//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "%s distance calculation will be using great circle approximation with %s auxiliary latitudes and return lengths in degrees.\n",
				//type_name[type], aux[choice]);
			//break;
		case GMT_DIST_DEG+GMT_GEODESIC:	/* 2-D lon, lat data, use geodesic distances in degrees */
			//GMT->current.map.dist[type].func = &gmt_geodesic_dist_degree;
			//GMT->current.map.azimuth_func = &gmt_az_backaz_geodesic;
			//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "%s distance calculation will be using geodesics in degrees\n", type_name[type]);
			//break;
		case GMT_DIST_COS+GMT_GREATCIRCLE:	/* 2-D lon, lat data, and Green's function needs cosine of spherical distance */
			GMT->current.map.dist[type].func = &GMT_great_circle_dist_cos;
			GMT->current.map.azimuth_func = &gmt_az_backaz_sphere;
			//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "%s distance calculation will be using great circle approximation with %s auxiliary latitudes and return cosine of spherical angles.\n",
				//type_name[type], aux[choice]);
			break;
		case GMT_DIST_COS+GMT_GEODESIC:	/* 2-D lon, lat data, and Green's function needs cosine of geodesic distance */
			GMT->current.map.dist[type].func = &GMT_geodesic_dist_cos;
			GMT->current.map.azimuth_func = &gmt_az_backaz_geodesic;
			//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "%s distance calculation will be using cosine of geodesic angle\n", type_name[type]);
			break;
		case GMT_DIST_M+GMT_LOXODROME:	/* 2-D lon, lat data, but measure distance along rhumblines in meter */
			GMT->current.map.dist[type].func = &gmt_loxodrome_dist_meter;
			GMT->current.map.azimuth_func  = &gmt_az_backaz_loxodrome;
			//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "%s distance calculation will be along loxodromes in meters\n", type_name[type]);
			break;
		case GMT_DIST_DEG+GMT_LOXODROME:	/* 2-D lon, lat data, but measure distance along rhumblines in degrees */
			GMT->current.map.dist[type].func = &gmt_loxodrome_dist_degree;
			GMT->current.map.azimuth_func = &gmt_az_backaz_loxodrome;
			//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "%s distance calculation will be along loxodromes with %s auxiliary latitudes and return lengths in degrees.\n",
				//type_name[type], aux[choice]);
			break;
		default:	/* Cannot happen unless we make a bug */
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Mode (=%d) for distance function is unknown. Must be bug.\n", mode);
			exit (EXIT_FAILURE);
			break;
	}
	if (type > 0) return;	/* Contour-related assignemnts end here */

	/* Mapping only */
	if (mode == GMT_CARTESIAN_DIST || mode == GMT_CARTESIAN_DIST2)	{	/* Cartesian data */
		//GMT->current.map.near_lines_func   = &gmt_near_lines_cartesian;
		//GMT->current.map.near_a_line_func  = &gmt_near_a_line_cartesian;
		//GMT->current.map.near_point_func   = &gmt_near_a_point_cartesian;
	}
	else {	/* Geographic data */
		GMT->current.map.near_lines_func   = &gmt_near_lines_spherical;
		GMT->current.map.near_a_line_func  = &gmt_near_a_line_spherical;
		GMT->current.map.near_point_func   = &gmt_near_a_point_spherical;
	}
}

unsigned int GMT_init_distaz (struct GMT_CTRL *GMT, char unit, unsigned int mode, unsigned int type)
{
	/* Initializes distance calcuation given the selected values for:
	 * Distance unit: must be on of the following:
	 *  1) d|e|f|k|m|M|n|s
	 *  2) GMT (Cartesian distance after projecting with -J) | X (Cartesian)
	 *  3) S (cosine distance) | P (cosine after first inverse projecting with -J)
	 * distance-calculation modifier mode: 0 (Cartesian), 1 (flat Earth), 2 (great-circle, 3 (geodesic), 4 (loxodrome)
	 * type: 0 = map distances, 1 = contour distances, 2 = contour annotation distances
	 * We set distance and azimuth functions and scales for this type.
	 * At the moment there is only one azimuth function pointer for all.
	 *
	 * The input args for GMT_init_distaz normally comes from calling GMT_get_distance.
	 */

	unsigned int proj_type = GMT_GEOGRAPHIC;	/* Default is to just use the geographic coordinates as they are */

	if (strchr (GMT_LEN_UNITS, unit) && !GMT_is_geographic (GMT, GMT_IN)) {	/* Want geographic distance units but -fg (or -J) not set */
		GMT_parse_common_options (GMT, "f", 'f', "g");
		//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "Your distance unit (%c) implies geographic data; -fg has been set.\n", unit);
	}
	//printf("file : %s line : %d func: %s unit %c \n",__FILE__,__LINE__,__func__,unit);
	switch (unit) {
			/* First the three arc angular distance units */
			
		case 'd':	/* Arc degrees on spherical body using desired metric mode */
			//gmt_set_distaz (GMT, GMT_DIST_DEG + mode, type);
			//GMT->current.map.dist[type].arc = true;	/* Angular measure */
			//break;
		case 'm':	/* Arc minutes on spherical body using desired metric mode */
			//gmt_set_distaz (GMT, GMT_DIST_DEG + mode, type);
			//GMT->current.map.dist[type].scale = GMT_DEG2MIN_F;
			//GMT->current.map.dist[type].arc = true;	/* Angular measure */
			//break;
		case 's':	/* Arc seconds on spherical body using desired metric mode */
			//gmt_set_distaz (GMT, GMT_DIST_DEG + mode, type);
			//GMT->current.map.dist[type].scale = GMT_DEG2SEC_F;
			//GMT->current.map.dist[type].arc = true;	/* Angular measure */
			//break;
			
			/* Various distance units on the planetary body */
			
		case 'e':	/* Meters on spherical body using desired metric mode */
			gmt_set_distaz (GMT, GMT_DIST_M + mode, type);
			break;
		case 'f':	/* Feet on spherical body using desired metric mode */
			gmt_set_distaz (GMT, GMT_DIST_M + mode, type);
			GMT->current.map.dist[type].scale = 1.0 / METERS_IN_A_FOOT;
			break;
		case 'k':	/* Km on spherical body using desired metric mode */
			gmt_set_distaz (GMT, GMT_DIST_M + mode, type);
			GMT->current.map.dist[type].scale = 1.0 / METERS_IN_A_KM;
			break;
		case 'M':	/* Statute Miles on spherical body using desired metric mode  */
			gmt_set_distaz (GMT, GMT_DIST_M + mode, type);
			GMT->current.map.dist[type].scale = 1.0 / METERS_IN_A_MILE;
			break;
		case 'n':	/* Nautical miles on spherical body using desired metric mode */
			gmt_set_distaz (GMT, GMT_DIST_M + mode, type);
			GMT->current.map.dist[type].scale = 1.0 / METERS_IN_A_NAUTICAL_MILE;
			break;
		case 'u':	/* Survey feet on spherical body using desired metric mode */
			gmt_set_distaz (GMT, GMT_DIST_M + mode, type);
			GMT->current.map.dist[type].scale = 1.0 / METERS_IN_A_SURVEY_FOOT;
			break;
			
			/* Cartesian distances.  Note: The X|C|R|Z|S|P 'units' are only passed internally and are not available as user selections directly */
			
		case 'X':	/* Cartesian distances in user units */
			//proj_type = GMT_CARTESIAN;
			//gmt_set_distaz (GMT, GMT_CARTESIAN_DIST, type);
			break;
		case 'C':	/* Cartesian distances (in PROJ_LENGTH_UNIT) after first projecting input coordinates with -J */
			//gmt_set_distaz (GMT, GMT_CARTESIAN_DIST_PROJ, type);
			//proj_type = GMT_GEO2CART;
			break;
			
		case 'R':	/* Cartesian distances squared in user units */
			//proj_type = GMT_CARTESIAN;
			//gmt_set_distaz (GMT, GMT_CARTESIAN_DIST2, type);
			break;
		case 'Z':	/* Cartesian distances squared (in PROJ_LENGTH_UNIT^2) after first projecting input coordinates with -J */
			//gmt_set_distaz (GMT, GMT_CARTESIAN_DIST_PROJ2, type);
			//proj_type = GMT_GEO2CART;
			break;

			/* Specialized cosine distances used internally only (e.g., greenspline) */
			
		case 'S':	/* Spherical cosine distances (for various gridding functions) */
			//gmt_set_distaz (GMT, GMT_DIST_COS + mode, type);
			break;
		case 'P':	/* Spherical distances after first inversily projecting Cartesian coordinates with -J */
			//gmt_set_distaz (GMT, GMT_CARTESIAN_DIST_PROJ_INV, type);
			//proj_type = GMT_CART2GEO;
			break;
			
		default:
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Syntax error: Distance units must be one of %s\n", GMT_LEN_UNITS_DISPLAY);
			GMT_exit (GMT, EXIT_FAILURE); return EXIT_FAILURE;
			break;
	}
	
	GMT->current.map.dist[type].init = true;	/* OK, we have now initialized the info for this type */
	return (proj_type);
}

