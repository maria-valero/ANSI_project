#include <stddef.h>
#include <stdbool.h>
#include <stdio.h>
#include <inttypes.h>
#include <netcdf.h>
#include "common_math.h"
#include "gmt_private.h"
#include "gmt_type.h"
#include "memory.h"
#include "gmt_macros.h"
#include "gmt_common.h"
#include "gmt_nan.h"
#include "gmt_io.h"
#include "grd_io.h"
#include "netcdf.h"

/* Tiny functions to tell if a value is <, <=, >=, > than the limit */
bool gmt_inside_lower_boundary (double val, double min) {return (val >= min);}
bool gmt_inside_upper_boundary (double val, double max) {return (val <= max);}
bool gmt_outside_lower_boundary (double val, double min) {return (val < min);}
bool gmt_outside_upper_boundary (double val, double max) {return (val > max);}


/* Functions and macros used for new rectangular clipping using the Sutherland/Hodgman algorithm
 * in which we clip the polygon against each of the 4 sides.  To avoid lots of if/switch I have
 * two clip functions (for x and y line) and two in/out functions that tells us if a point is
 * inside the polygon relative to the line.  Then, pointers to these functions are passed to
 * make sure the right functions are used for each side in the loop over sides.  Unless I can
 * figure out a more clever recursive way I need to have 2 temporary arrays to shuffle the
 * intermediate results around.
 *
 * P.Wessel, March 2008
 */

/* This macro calculates the x-coordinates where the line segment crosses the border x = border.
 * By swapping x and y in the call we can use it for finding the y intersection. This macro is
 * never called when (y_prev - y_curr) = 0 so we don't divide by zero.
 */
#define INTERSECTION_COORD(x_curr,y_curr,x_prev,y_prev,border) x_curr + (x_prev - x_curr) * (border - y_curr) / (y_prev - y_curr)
void gmt_alloc_ogr_seg (struct GMT_CTRL *GMT, struct GMT_DATASEGMENT *S, int n_aspatial)
{	/* Allocates the OGR structure for a given segment and copies current values from table OGR segment */
	if (S->ogr) return;	/* Already allocated */
	S->ogr = GMT_memory (GMT, NULL, 1, struct GMT_OGR_SEG);
	S->ogr->n_aspatial = n_aspatial;
	if (n_aspatial) {
		S->ogr->tvalue = GMT_memory (GMT, NULL, n_aspatial, char *);
		S->ogr->dvalue = GMT_memory (GMT, NULL, n_aspatial, double);
	}
}

bool GMT_polygon_is_open (struct GMT_CTRL *GMT, double x[], double y[], uint64_t n)
{	/* Returns true if the first and last point is not identical */
	if (n < 2) return false;	/*	A single point is by definition closed */
	if (!doubleAlmostEqualZero (y[0], y[n-1]))
		return true;	/* y difference exceeds threshold: polygon is OPEN */
	if (!doubleAlmostEqualZero (x[0], x[n-1])) {	/* The x values exceeds threshold, check further if by 360 under geo */
		if (GMT->current.io.col_type[GMT_IN][GMT_X] & GMT_IS_GEO) {	/* Geographical coordinates: Worry about a 360 jump */
			double dlon = fabs (x[0] - x[n-1]);	/* If exactly 360 then we are OK */
			if (!doubleAlmostEqualZero (dlon, 360.0))
				return true;	/* x difference exceeds threshold for an exact 360 offset: polygon is OPEN */
		}
		else	/* Cartesian case */
			return true;	/* x difference exceeds threshold: polygon is OPEN */
	}
	/* Here, first and last are ~identical - to be safe we enforce exact closure */
	x[n-1] = x[0];	y[n-1] = y[0];	/* Note: For geo data, this step may change a 0 to 360 or vice versa */
	return false;	/* Passed the tests so polygon is CLOSED */
}

static inline bool gmt_same_longitude (double a, double b) {
	/* return true if a and b are the same longitude */
	while (a < 0.0)   a += 360.0;
	while (a > 360.0) a -= 360.0;
	while (b < 0.0)   b += 360.0;
	while (b > 360.0) b -= 360.0;
	return doubleAlmostEqualZero (a, b);
}

uint64_t gmt_getprevpoint (double plon, double lon[], uint64_t n, uint64_t this_p)
{	/* Return the previous point that does NOT equal plon */
	uint64_t ip = (this_p == 0) ? n - 2 : this_p - 1;	/* Previous point (-2 because last is a duplicate of first) */
	while (doubleAlmostEqualZero (plon, lon[ip]) || doubleAlmostEqual (fabs(plon - lon[ip]), 360.0)) {	/* Same as plon */
		if (ip == 0)
			ip = n - 2;
		else
			ip--;
	}
	return (ip);
}

int gmt_inonout_sphpol_count (double plon, double plat, const struct GMT_DATASEGMENT *P, unsigned int count[])
{	/* Case of a polar cap */
	uint64_t i, in, ip, prev;
	int cut;
	double W, E, S, N, lon, lon1, lon2, dlon, x_lat;

	/* Draw meridian through P and count all the crossings with the line segments making up the polar cap S */

	GMT_memset (count, 2, unsigned int);	/* Initialize counts to zero */
	for (i = 0; i < P->n_rows - 1; i++) {	/* -1, since we know last point repeats the first */
		/* Here lon1 and lon2 are the end points (in longitude) of the current line segment in S.  There are
		 * four cases to worry about:
		 * 1) lon equals lon1 (i.e., the meridian through lon goes right through lon1)
		 * 2) lon equals lon2 (i.e., the meridian through lon goes right through lon2)
		 * 3) lon lies between lon1 and lon2 and crosses the segment
		 * 4) none of the above
		 * Since we want to obtain either ONE or ZERO intersections per segment we will skip to next
		 * point if case (2) occurs: this avoids counting a crossing twice for consequtive segments.
		 */
		if (gmt_same_longitude (plon, P->coord[GMT_X][i]) && GMT_SAME_LATITUDE (plat, P->coord[GMT_Y][i])) return (1);	/* Point is on the perimeter */
		in = i + 1;			/* Next point index */
		/* First skip segments that have no actual length: consecutive points with both latitudes == -90 or +90 */
		if (fabs (P->coord[GMT_Y][i]) == 90.0 && doubleAlmostEqualZero (P->coord[GMT_Y][i], P->coord[GMT_Y][in]))
			continue;
		/* Next deal with case when the longitude of P goes ~right through the second of the line nodes */
		if (gmt_same_longitude (plon, P->coord[GMT_X][in])) continue;	/* Line goes through the 2nd node - ignore */
		lon1 = P->coord[GMT_X][i];	/* Copy the first of two longitudes since we may need to mess with them */
		lon2 = P->coord[GMT_X][in];	/* Copy the second of two longitudes since we may need to mess with them */
		if (gmt_same_longitude (plon, lon1)) {	/* Line goes through the 1st node */
			/* Must check that the two neighboring points are on either side; otherwise it is just a tangent line */
			ip = gmt_getprevpoint (plon, P->coord[GMT_X], P->n_rows, i);	/* Index of previous point != plon */
			if ((lon2 >= lon1 && P->coord[GMT_X][ip] > lon1) || (lon2 <= lon1 && P->coord[GMT_X][ip] < lon1)) continue;	/* Both on same side */
			cut = (P->coord[GMT_Y][i] > plat) ? 0 : 1;	/* node is north (0) or south (1) of P */
			count[cut]++;
			prev = ip + 1;	/* Always exists because ip is <= n-2 */
			/* If prev < i then we have a vertical segment of 2 or more points; prev points to the other end of the segment.
			 * We must then check if our points plat is within that range, meaning the point lies on the segment */
			if (prev < i && ((plat <= P->coord[GMT_Y][prev] && plat >= P->coord[GMT_Y][i]) || (plat <= P->coord[GMT_Y][i] && plat >= P->coord[GMT_Y][prev]))) return (1);	/* P is on segment boundary; we are done*/
			continue;
		}
		/* OK, not exactly on a node, deal with crossing a line */
		dlon = lon2 - lon1;
		if (dlon > 180.0)		/* Jumped across Greenwich going westward */
			lon2 -= 360.0;
		else if (dlon < -180.0)		/* Jumped across Greenwich going eastward */
			lon1 -= 360.0;
		if (lon1 <= lon2) {	/* Segment goes W to E (or N-S) */
			W = lon1;
			E = lon2;
		}
		else {			/* Segment goes E to W */
			W = lon2;
			E = lon1;
		}
		lon = plon;			/* Local copy of plon, below adjusted given the segment lon range */
		while (lon > W) lon -= 360.0;	/* Make sure we rewind way west for starters */
		while (lon < W) lon += 360.0;	/* Then make sure we wind to inside the lon range or way east */
		if (lon > E) continue;	/* Not crossing this segment */
		if (dlon == 0.0) {	/* Special case of N-S segment: does P lie on it? */
			if (P->coord[GMT_Y][in] < P->coord[GMT_Y][i]) {	/* Get N and S limits for segment */
				S = P->coord[GMT_Y][in];
				N = P->coord[GMT_Y][i];
			}
			else {
				N = P->coord[GMT_Y][in];
				S = P->coord[GMT_Y][i];
			}
			if (plat < S || plat > N) continue;	/* P is not on this segment */
			return (1);	/* P is on segment boundary; we are done*/
		}
		/* Calculate latitude at intersection */
		if (GMT_SAME_LATITUDE (P->coord[GMT_Y][i], P->coord[GMT_Y][in]) && GMT_SAME_LATITUDE (plat, P->coord[GMT_Y][in])) return (1);	/* P is on S boundary */
		x_lat = P->coord[GMT_Y][i] + ((P->coord[GMT_Y][in] - P->coord[GMT_Y][i]) / (lon2 - lon1)) * (lon - lon1);
		if (doubleAlmostEqualZero (x_lat, plat))
			return (1);	/* P is on S boundary */

		cut = (x_lat > plat) ? 0 : 1;	/* Cut is north (0) or south (1) of P */
		count[cut]++;
	}

	return (0);	/* This means no special cases were detected that warranted an immediate return */
}

unsigned int GMT_inonout_sphpol (struct GMT_CTRL *GMT, double plon, double plat, const struct GMT_DATASEGMENT *P)
/* This function is used to see if some point P is located inside, outside, or on the boundary of the
 * spherical polygon S read by GMT_import_table.  Note GMT->current.io.skip_duplicates must be true when the polygon
 * was read so there are NO duplicate (repeated) points.
 * Returns the following values:
 *	0:	P is outside of S
 *	1:	P is inside of S
 *	2:	P is on boundary of S
 */
{
	/* Algorithm:
	 * Case 1: The polygon S contains a geographical pole
	 *	   a) if P is beyond the far latitude then P is outside
	 *	   b) Draw meridian through P and count intersections:
	 *		odd: P is outside; even: P is inside
	 * Case 2: S does not contain a pole
	 *	   a) If P is outside range of latitudes then P is outside
	 *	   c) Draw meridian through P and count intersections:
	 *		odd: P is inside; even: P is outside
	 * In all cases, we check if P is on the outline of S
	 */

	unsigned int count[2];

	if (P->pole) {	/* Case 1 of an enclosed polar cap */
		if (P->pole == +1) {	/* N polar cap */
			if (plat < P->min[GMT_Y]) return (GMT_OUTSIDE);	/* South of a N polar cap */
			if (plat > P->lat_limit) return (GMT_INSIDE);	/* Clearly inside of a N polar cap */
		}
		if (P->pole == -1) {	/* S polar cap */
			if (plat > P->max[GMT_Y]) return (GMT_OUTSIDE);	/* North of a S polar cap */
			if (plat < P->lat_limit) return (GMT_INSIDE);	/* Clearly inside of a S polar cap */
		}

		/* Tally up number of intersections between polygon and meridian through P */

		if (gmt_inonout_sphpol_count (plon, plat, P, count)) return (GMT_ONEDGE);	/* Found P is on S */

		if (P->pole == +1 && count[0] % 2 == 0) return (GMT_INSIDE);
		if (P->pole == -1 && count[1] % 2 == 0) return (GMT_INSIDE);

		return (GMT_OUTSIDE);
	}

	/* Here is Case 2.  First check latitude range */

	if (plat < P->min[GMT_Y] || plat > P->max[GMT_Y]) return (GMT_OUTSIDE);

	/* Longitudes are tricker and are tested with the tallying of intersections */

	if (gmt_inonout_sphpol_count (plon, plat, P, count)) return (GMT_ONEDGE);	/* Found P is on S */

	if (count[0] % 2) return (GMT_INSIDE);

	return (GMT_OUTSIDE);	/* Nothing triggered the tests; we are outside */
}

unsigned int GMT_non_zero_winding (struct GMT_CTRL *GMT, double xp, double yp, double *x, double *y, uint64_t n_path)
{
	/* Routine returns (2) if (xp,yp) is inside the
	   polygon x[n_path], y[n_path], (0) if outside,
	   and (1) if exactly on the path edge.
	   Uses non-zero winding rule in Adobe PostScript
	   Language reference manual, section 4.6 on Painting.
	   Path should have been closed first, so that
	   x[n_path-1] = x[0], and y[n_path-1] = y[0].

	   This is version 2, trying to kill a bug
	   in above routine:  If point is on X edge,
	   fails to discover that it is on edge.

	   We are imagining a ray extending "up" from the
	   point (xp,yp); the ray has equation x = xp, y >= yp.
	   Starting with crossing_count = 0, we add 1 each time
	   the path crosses the ray in the +x direction, and
	   subtract 1 each time the ray crosses in the -x direction.
	   After traversing the entire polygon, if (crossing_count)
	   then point is inside.  We also watch for edge conditions.

	   If two or more points on the path have x[i] == xp, then
	   we have an ambiguous case, and we have to find the points
	   in the path before and after this section, and check them.
	   */

	uint64_t i, j, k, jend, crossing_count;
	bool above;
	double y_sect;

	if (n_path < 2) return (GMT_OUTSIDE);	/* Cannot be inside a null set or a point so default to outside */

	if (GMT_polygon_is_open (GMT, x, y, n_path)) {
		////GMT_Report (GMT->parent, GMT_MSG_NORMAL, "given non-closed polygon\n");
		GMT_exit (GMT, EXIT_FAILURE); return EXIT_FAILURE;
	}

	above = false;
	crossing_count = 0;

	/* First make sure first point in path is not a special case:  */
	j = jend = n_path - 1;
	if (x[j] == xp) {
		/* Trouble already.  We might get lucky:  */
		if (y[j] == yp) return (GMT_ONEDGE);

		/* Go backward down the polygon until x[i] != xp:  */
		if (y[j] > yp) above = true;
		i = j - 1;
		while (x[i] == xp && i > 0) {
			if (y[i] == yp) return (GMT_ONEDGE);
			if (!(above) && y[i] > yp) above = true;
			i--;
		}

		/* Now if i == 0 polygon is degenerate line x=xp;
		   since we know xp,yp is inside bounding box,
		   it must be on edge:  */
		if (i == 0) return (GMT_ONEDGE);

		/* Now we want to mark this as the end, for later:  */
		jend = i;

		/* Now if (j-i)>1 there are some segments the point could be exactly on:  */
		for (k = i+1; k < j; k++) {
			if ((y[k] <= yp && y[k+1] >= yp) || (y[k] >= yp && y[k+1] <= yp)) return (GMT_ONEDGE);
		}


		/* Now we have arrived where i is > 0 and < n_path-1, and x[i] != xp.
			We have been using j = n_path-1.  Now we need to move j forward
			from the origin:  */
		j = 1;
		while (x[j] == xp) {
			if (y[j] == yp) return (GMT_ONEDGE);
			if (!(above) && y[j] > yp) above = true;
			j++;
		}

		/* Now at the worst, j == jstop, and we have a polygon with only 1 vertex
			not at x = xp.  But now it doesn't matter, that would end us at
			the main while below.  Again, if j>=2 there are some segments to check:  */
		for (k = 0; k < j-1; k++) {
			if ((y[k] <= yp && y[k+1] >= yp) || (y[k] >= yp && y[k+1] <= yp)) return (GMT_ONEDGE);
		}


		/* Finally, we have found an i and j with points != xp.  If (above) we may have crossed the ray:  */
		if (above && x[i] < xp && x[j] > xp)
			crossing_count++;
		else if (above && x[i] > xp && x[j] < xp)
			crossing_count--;

		/* End nightmare scenario for x[0] == xp.  */
	}

	else {
		/* Get here when x[0] != xp:  */
		i = 0;
		j = 1;
		while (x[j] == xp && j < jend) {
			if (y[j] == yp) return (GMT_ONEDGE);
			if (!(above) && y[j] > yp) above = true;
			j++;
		}
		/* Again, if j==jend, (i.e., 0) then we have a polygon with only 1 vertex
			not on xp and we will branch out below.  */

		/* if ((j-i)>2) the point could be on intermediate segments:  */
		for (k = i+1; k < j-1; k++) {
			if ((y[k] <= yp && y[k+1] >= yp) || (y[k] >= yp && y[k+1] <= yp)) return (GMT_ONEDGE);
		}

		/* Now we have x[i] != xp and x[j] != xp.
			If (above) and x[i] and x[j] on opposite sides, we are certain to have crossed the ray.
			If not (above) and (j-i)>1, then we have not crossed it.
			If not (above) and j-i == 1, then we have to check the intersection point.  */

		if (x[i] < xp && x[j] > xp) {
			if (above)
				crossing_count++;
			else if ((j-i) == 1) {
				y_sect = y[i] + (y[j] - y[i]) * ((xp - x[i]) / (x[j] - x[i]));
				if (y_sect == yp) return (GMT_ONEDGE);
				if (y_sect > yp) crossing_count++;
			}
		}
		if (x[i] > xp && x[j] < xp) {
			if (above)
				crossing_count--;
			else if ((j-i) == 1) {
				y_sect = y[i] + (y[j] - y[i]) * ((xp - x[i]) / (x[j] - x[i]));
				if (y_sect == yp) return (GMT_ONEDGE);
				if (y_sect > yp) crossing_count--;
			}
		}

		/* End easier case for x[0] != xp  */
	}

	/* Now MAIN WHILE LOOP begins:
		Set i = j, and search for a new j, and do as before.  */

	i = j;
	while (i < jend) {
		above = false;
		j = i+1;
		while (x[j] == xp) {
			if (y[j] == yp) return (GMT_ONEDGE);
			if (!(above) && y[j] > yp) above = true;
			j++;
		}
		/* if ((j-i)>2) the point could be on intermediate segments:  */
		for (k = i+1; k < j-1; k++) {
			if ((y[k] <= yp && y[k+1] >= yp) || (y[k] >= yp && y[k+1] <= yp)) return (GMT_ONEDGE);
		}

		/* Now we have x[i] != xp and x[j] != xp.
			If (above) and x[i] and x[j] on opposite sides, we are certain to have crossed the ray.
			If not (above) and (j-i)>1, then we have not crossed it.
			If not (above) and j-i == 1, then we have to check the intersection point.  */

		if (x[i] < xp && x[j] > xp) {
			if (above)
				crossing_count++;
			else if ((j-i) == 1) {
				y_sect = y[i] + (y[j] - y[i]) * ((xp - x[i]) / (x[j] - x[i]));
				if (y_sect == yp) return (GMT_ONEDGE);
				if (y_sect > yp) crossing_count++;
			}
		}
		if (x[i] > xp && x[j] < xp) {
			if (above)
				crossing_count--;
			else if ((j-i) == 1) {
				y_sect = y[i] + (y[j] - y[i]) * ((xp - x[i]) / (x[j] - x[i]));
				if (y_sect == yp) return (GMT_ONEDGE);
				if (y_sect > yp) crossing_count--;
			}
		}

		/* That's it for this piece.  Advance i:  */

		i = j;
	}

	/* End of MAIN WHILE.  Get here when we have gone all around without landing on edge.  */

	return ((crossing_count) ? GMT_INSIDE: GMT_OUTSIDE);
}

unsigned int gmt_inonout_sub (struct GMT_CTRL *GMT, double x, double y, const struct GMT_DATASEGMENT *S)
{	/* Front end for both spherical and Cartesian in-on-out functions */
	unsigned int side;

	if (GMT_is_geographic (GMT, GMT_IN)) {	/* Assumes these are input polygons */
		if (S->pole)	/* 360-degree polar cap, must check fully */
			side = GMT_inonout_sphpol (GMT, x, y, S);
		else {	/* See if we are outside range of longitudes for polygon */
			while (x > S->min[GMT_X]) x -= 360.0;	/* Wind clear of west */
			while (x < S->min[GMT_X]) x += 360.0;	/* Wind east until inside or beyond east */
			if (x > S->max[GMT_X]) return (GMT_OUTSIDE);	/* Point outside, no need to assign value */
			side = GMT_inonout_sphpol (GMT, x, y, S);
		}
	}
	else {	/* Cartesian case */
		if (x < S->min[GMT_X] || x > S->max[GMT_X]) return (GMT_OUTSIDE);	/* Point outside, no need to assign value */
		side = GMT_non_zero_winding (GMT, x, y, S->coord[GMT_X], S->coord[GMT_Y], S->n_rows);
	}
	return (side);
}

unsigned int GMT_inonout (struct GMT_CTRL *GMT, double x, double y, const struct GMT_DATASEGMENT *S)
{	/* Front end for both spherical and Cartesian in-on-out functions.
 	 * Knows to check for polygons with holes as well. */
	unsigned int side, side_h;
	struct GMT_DATASEGMENT *H = NULL;

	if ((side = gmt_inonout_sub (GMT, x, y, S)) <= GMT_ONEDGE) return (side);	/* Outside polygon or on perimeter, we are done */

	/* Here, point is inside the polygon perimeter. See if there are holes */

	if (GMT->current.io.OGR && (H = S->next)) {	/* Must check for and skip if inside a hole */
		side_h = GMT_OUTSIDE;	/* We are outside a hole until we are found to be inside it */
		while (side_h == GMT_OUTSIDE && H && H->ogr && H->ogr->pol_mode == GMT_IS_HOLE) {	/* Found a hole */
			/* Must check if point is inside this hole polygon */
			side_h = gmt_inonout_sub (GMT, x, y, H);
			H = H->next;	/* Move to next polygon hole */
		}
		if (side_h == GMT_INSIDE) side = GMT_OUTSIDE;	/* Inside one of the holes, hence outside polygon; go to next perimeter polygon */
		if (side_h == GMT_ONEDGE) side = GMT_ONEDGE;	/* On path of one of the holes, hence on polygon path; update side */
	}

	/* Here, point is inside or on edge, we return the value */
	return (side);
}

bool GMT_crossing_dateline (struct GMT_CTRL *GMT, struct GMT_DATASEGMENT *S)
{	/* Return true if this line or polygon feature contains points on either side of the Dateline */
	uint64_t k;
	bool east = false, west = false, cross = false;
	for (k = 0; !cross && k < S->n_rows; k++) {
		if ((S->coord[GMT_X][k] > 180.0 && S->coord[GMT_X][k] < 270.0) || (S->coord[GMT_X][k] > -180.0 && S->coord[GMT_X][k] <  -90.0)) west = true;
		if ((S->coord[GMT_X][k] >  90.0 && S->coord[GMT_X][k] < 180.0) || (S->coord[GMT_X][k] > -270.0 && S->coord[GMT_X][k] < -180.0)) east = true;
		if (east && west) cross = true;
	}
	return (cross);
}

unsigned int gmt_clip_we (double x_prev, double y_prev, double x_curr, double y_curr, double x[], double y[], double border, bool (*inside) (double, double), bool (*outside) (double, double), int *cross)
{	/* Clip against the west or east boundary (i.e., a vertical line with x = border) */
	*cross = 0;
	if (doubleAlmostEqualZero (x_prev, x_curr) && doubleAlmostEqualZero (y_prev, y_curr))
		return (0);	/* Do nothing for duplicates */
	if (outside (x_prev, border)) {	/* Previous point is outside... */
		if (outside (x_curr, border)) return 0;	/* ...as is the current point. Do nothing. */
		/* Here, the line segment intersects the border - return both intersection and inside point */
		x[0] = border;	y[0] = INTERSECTION_COORD (y_curr, x_curr, y_prev, x_prev, border);
		*cross = +1;	/* Crossing to the inside */
		x[1] = x_curr;	y[1] = y_curr;	return (2);
	}
	/* Here x_prev is inside */
	if (inside (x_curr, border)) {	/* Return current point only */
		x[0] = x_curr;	y[0] = y_curr;	return (1);
	}
	/* Segment intersects border - return intersection only */
	*cross = -1;	/* Crossing to the outside */
	x[0] = border;	y[0] = INTERSECTION_COORD (y_curr, x_curr, y_prev, x_prev, border);	return (1);
}

void GMT_duplicate_ogr_seg (struct GMT_CTRL *GMT, struct GMT_DATASEGMENT *S_to, struct GMT_DATASEGMENT *S_from)
{	/* Allocates the OGR structure for a given segment and copies current values from table OGR segment */
	unsigned int col;

	if (!S_from->ogr) return;	/* No data */
	gmt_alloc_ogr_seg (GMT, S_to, S_from->ogr->n_aspatial);
	for (col = 0; col < S_from->ogr->n_aspatial; col++) {
		if (S_from->ogr->tvalue[col]) S_to->ogr->tvalue[col] = strdup (S_from->ogr->tvalue[col]);
		S_to->ogr->dvalue[col] = S_from->ogr->dvalue[col];
	}
	S_to->ogr->pol_mode = S_from->ogr->pol_mode;
}

/* GMT_dateline_clip simply clips a polygon agains the dateline and results in two polygons in L */

unsigned int GMT_split_poly_at_dateline (struct GMT_CTRL *GMT, struct GMT_DATASEGMENT *S, struct GMT_DATASEGMENT ***Lout)
{
	int side, j, np, cross = 0;
	uint64_t row, m;
	size_t n_alloc = 0;
	char label[GMT_BUFSIZ] = {""}, *part = "EW";
	double xx[2], yy[2];
	struct GMT_DATASEGMENT **L = NULL;
	bool (*inside[2]) (double, double);
	bool (*outside[2]) (double, double);


	inside[0] = gmt_inside_upper_boundary;	outside[0] = gmt_outside_upper_boundary;
	inside[1] = gmt_inside_lower_boundary;	outside[1] = gmt_outside_lower_boundary;
	L = GMT_memory (GMT, NULL, 2, struct GMT_DATASEGMENT *);	/* The two polygons */

	for (row = 0; row < S->n_rows; row++) GMT_lon_range_adjust (GMT_IS_0_TO_P360_RANGE, &S->coord[GMT_X][row]);	/* First enforce 0 <= lon < 360 so we dont have to check again */

	for (side = 0; side < 2; side++) {	/* Do it twice to get two truncated polygons */
		L[side] = GMT_memory (GMT, NULL, 1, struct GMT_DATASEGMENT);
		n_alloc = lrint (1.05*S->n_rows+5);	/* Anticipate just a few crossings (5%)+5, allocate more later if needed */
		GMT_alloc_segment (GMT, L[side], n_alloc, S->n_columns, true);	/* Temp segment with twice the number of points as we will add crossings*/
		m = 0;		/* Start with nuthin' */

		/* Must ensure we copy the very first point if it is left of the Dateline */
		if (S->coord[GMT_X][0] < 180.0) { L[side]->coord[GMT_X][0] = S->coord[GMT_X][0]; L[side]->coord[GMT_Y][0] = S->coord[GMT_Y][0]; }	/* First point is inside; add it */
		for (row = 1; row < S->n_rows; row++) {	/* For each line segment */
			np = gmt_clip_we (S->coord[GMT_X][row-1], S->coord[GMT_Y][row-1], S->coord[GMT_X][row], S->coord[GMT_Y][row], xx, yy, 180.0, inside[side], outside[side], &cross);	/* Returns 0, 1, or 2 points */
			for (j = 0; j < np; j++) {	/* Add the np returned points to the new clipped polygon path */
				if (m == n_alloc) GMT_alloc_segment (GMT, L[side], n_alloc << 2, S->n_columns, false);
				L[side]->coord[GMT_X][m] = xx[j]; L[side]->coord[GMT_Y][m] = yy[j]; m++;
			}
		}
		if (GMT_polygon_is_open (GMT, L[side]->coord[GMT_X], L[side]->coord[GMT_Y], m)) {	/* Do we need to explicitly close this clipped polygon? */
			if (m == n_alloc) GMT_alloc_segment (GMT, L[side], n_alloc << 2, S->n_columns, false);
			L[side]->coord[GMT_X][m] = L[side]->coord[GMT_X][0];	L[side]->coord[GMT_Y][m] = L[side]->coord[GMT_Y][0];	m++;	/* Yes. */
		}
		if (m != n_alloc) GMT_alloc_segment (GMT, L[side], m, S->n_columns, false);
		L[side]->n_rows = m;
		if (S->label) {
			sprintf (label, "%s part %c", S->label, part[side]);
			L[side]->label = strdup (label);
		}
		if (S->header) L[side]->header = strdup (S->header);
		if (S->ogr) GMT_duplicate_ogr_seg (GMT, L[side], S);
	}
	L[0]->range = 2;	L[1]->range = 3;
	*Lout = L;
	return (2);
}

bool gmt_straddle_dateline (double x0, double x1) {
	if (fabs (x0 - x1) > 90.0) return (false);	/* Probably Greenwhich crossing with 0/360 discontinuity */
	if ((x0 < 180.0 && x1 > 180.0) || (x0 > 180.0 && x1 < 180.0)) return (true);	/* Crossed Dateline */
	return (false);
}


unsigned int GMT_split_line_at_dateline (struct GMT_CTRL *GMT, struct GMT_DATASEGMENT *S, struct GMT_DATASEGMENT ***Lout)
{	/* Create two or more feature segments by splitting them across the Dateline.
	 * GMT_split_line_at_dateline should ONLY be called when we KNOW we must split. */
	unsigned int n_split;
	uint64_t k, col, seg, row, start, length, *pos = GMT_memory (GMT, NULL, S->n_rows, uint64_t);
	char label[GMT_BUFSIZ] = {""}, *txt = NULL, *feature = "Line";
	double r;
	struct GMT_DATASEGMENT **L = NULL, *Sx = GMT_memory (GMT, NULL, 1, struct GMT_DATASEGMENT);

	for (k = 0; k < S->n_rows; k++) GMT_lon_range_adjust (GMT_IS_0_TO_P360_RANGE, &S->coord[GMT_X][k]);	/* First enforce 0 <= lon < 360 so we dont have to check again */
	GMT_alloc_segment (GMT, Sx, 2*S->n_rows, S->n_columns, true);	/* Temp segment with twice the number of points as we will add crossings*/

	for (k = row = n_split = 0; k < S->n_rows; k++) {	/* Hunt for crossings */
		if (k && gmt_straddle_dateline (S->coord[GMT_X][k-1], S->coord[GMT_X][k])) {	/* Crossed Dateline */
			r = (180.0 - S->coord[GMT_X][k-1]) / (S->coord[GMT_X][k] - S->coord[GMT_X][k-1]);	/* Fractional distance from k-1'th point to 180 crossing */
			Sx->coord[GMT_X][row] = 180.0;	/* Exact longitude is known */
			for (col = 1; col < S->n_columns; col++) Sx->coord[col][row] = S->coord[col][k-1] + r * (S->coord[col][k] - S->coord[col][k-1]);	/* Linear interpolation for other fields */
			pos[n_split++] = row++;		/* Keep track of first point (the crossing) in new section */
		}
		for (col = 0; col < S->n_columns; col++) Sx->coord[col][row] = S->coord[col][k];	/* Append the current point */
		row++;
	}
	Sx->n_rows = row;	/* Number of points in extended feature with explicit crossings */
	if (n_split == 0) {	/* No crossings, should not have been called in the first place */
		////GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "No straddling detected (bug?)\n");
		GMT_free_segment (GMT, &Sx, GMT_ALLOCATED_BY_GMT);
		GMT_free (GMT, pos);
		return 0;
	}
	pos[n_split] = Sx->n_rows - 1;
	n_split++;	/* Now means number of segments */
	L = GMT_memory (GMT, NULL, n_split, struct GMT_DATASEGMENT *);	/* Number of output segments needed are allocated here */
	txt = (S->label) ? S->label : feature;	/* What to label the features */
	start = 0;
	for (seg = 0; seg < n_split; seg++) {	/* Populate the output segment coordinate arrays */
		L[seg] = GMT_memory (GMT, NULL, 1, struct GMT_DATASEGMENT);		/* Allocate space for one segment */
		length = pos[seg] - start + 1;	/* Length of new segment */
		GMT_alloc_segment (GMT, L[seg], length, S->n_columns, true);		/* Allocate array space for coordinates */
		for (col = 0; col < S->n_columns; col++) GMT_memcpy (L[seg]->coord[col], &(Sx->coord[col][start]), length, double);	/* Copy coordinates */
		L[seg]->range = (L[seg]->coord[GMT_X][length/2] > 180.0) ? GMT_IS_M180_TO_P180 : GMT_IS_M180_TO_P180_RANGE;	/* Formatting ID to enable special -180 and +180 formatting on outout */
		/* Modify label to part number */
		sprintf (label, "%s part %" PRIu64, txt, seg);
		L[seg]->label = strdup (label);
		if (S->header) L[seg]->header = strdup (S->header);
		if (S->ogr) GMT_duplicate_ogr_seg (GMT, L[seg], S);
		start = pos[seg];
	}
	GMT_free_segment (GMT, &Sx, GMT_ALLOCATED_BY_GMT);
	GMT_free (GMT, pos);

	*Lout = L;		/* Pass pointer to the array of segments */

	return (n_split);	/* Return how many segments was made */
}

int gmt_prep_ogr_output (struct GMT_CTRL *GMT, struct GMT_DATASET *D) {

	int object_ID, col, stop, n_reg, item, error = 0;
	uint64_t row, seg, seg1, seg2, k;
	char buffer[GMT_BUFSIZ] = {""}, in_string[GMT_STR16] = {""}, out_string[GMT_STR16] = {""};
	struct GMT_DATATABLE *T = NULL;
	struct GMT_DATASET *M = NULL;
	struct GMT_DATASEGMENT *S = NULL;
	struct GMTAPI_DATA_OBJECT O;

	/* When this functions is called we have already registered the output destination.  This will normally
	 * prevent us from register the data set separately in order to call GMT_gmtinfo.  We must temporarily
	 * unregister the output, do our thing, then reregister again. */

	// nishita, check arg geometry in the following func call
	n_reg = GMTAPI_count_objects (GMT->parent, (enum GMT_enum_family) GMT_IS_DATASET, D->geometry, GMT_OUT, &object_ID);	/* Are there outputs registered already? */
	if (n_reg == 1) {	/* Yes, must save and unregister, then reregister later */
		if ((item = GMTAPI_Validate_ID (GMT->parent, GMT_IS_DATASET, object_ID, GMT_OUT)) == GMT_NOTSET) return (GMTAPI_report_error (GMT->parent, error));
		GMT_memcpy (&O, GMT->parent->object[item], 1, struct GMTAPI_DATA_OBJECT);
		GMTAPI_Unregister_IO (GMT->parent, object_ID, GMT_OUT);
	}

	/* Determine w/e/s/n via GMT_gmtinfo */

	/* Create option list, register D as input source via ref */
	if ((object_ID = GMT_Register_IO (GMT->parent, GMT_IS_DATASET, GMT_IS_REFERENCE, GMT_IS_POINT, GMT_IN, NULL, D)) == GMT_NOTSET) {
		return (GMT->parent->error);
	}
	if (GMT_Encode_ID (GMT->parent, in_string, object_ID) != GMT_OK) {
		return (GMT->parent->error);	/* Make filename with embedded object ID */
	}
	if ((object_ID = GMT_Register_IO (GMT->parent, GMT_IS_DATASET, GMT_IS_DUPLICATE, GMT_IS_POINT, GMT_OUT, NULL, NULL)) == GMT_NOTSET) {
		return (GMT->parent->error);
	}
	if (GMT_Encode_ID (GMT->parent, out_string, object_ID)) {
		return (GMT->parent->error);	/* Make filename with embedded object ID */
	}
	sprintf (buffer, "-C -fg -<%s ->%s", in_string, out_string);
#if 0
	if (GMT_Call_Module (GMT->parent, "gmtinfo", GMT_MODULE_CMD, buffer) != GMT_OK) {	/* Get the extent via gmtinfo */
		return (GMT->parent->error);
	}
#endif
	if ((M = GMT_Retrieve_Data (GMT->parent, object_ID)) == NULL) {
		return (GMT->parent->error);
	}

	/* Time to reregister the original destination */

	if ((object_ID = GMT_Register_IO (GMT->parent, GMT_IS_DATASET, GMT_IS_FILE, GMT_IS_POINT, GMT_OUT, NULL, D)) == GMT_NOTSET) {
		return (GMT->parent->error);
	}
	if ((item = GMTAPI_Validate_ID (GMT->parent, GMT_IS_DATASET, object_ID, GMT_OUT)) == GMT_NOTSET) {
		return (GMTAPI_report_error (GMT->parent, error));
	}
	GMT_memcpy (GMT->parent->object[item], &O, 1, struct GMTAPI_DATA_OBJECT);	/* Restore what we had before */

	T = D->table[0];
	T->ogr = GMT_memory (GMT, NULL, 1, struct GMT_OGR);
	sprintf (buffer, "%.8g/%.8g/%.8g/%.8g", M->table[0]->segment[0]->coord[0][0], M->table[0]->segment[0]->coord[1][0], M->table[0]->segment[0]->coord[2][0], M->table[0]->segment[0]->coord[3][0]);
	if (GMT_Destroy_Data (GMT->parent, &M) != GMT_OK) {
		return (GMT->parent->error);
	}
	T->ogr->region = strdup (buffer);
	T->ogr->proj[1] = strdup ("\"-Jx1d --PROJ_ELLIPSOID=WGS84\"");
	T->ogr->proj[2] = strdup ("\"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs\"");
	T->ogr->proj[3] = strdup ("\"GEOGCS[\"GCS_WGS_1984\",DATUM[\"WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]]\"");
	T->ogr->geometry = GMT->common.a.geometry;
	T->ogr->n_aspatial = GMT->common.a.n_aspatial;
	if (T->ogr->n_aspatial) {	/* Copy over the command-line settings */
		T->ogr->name = GMT_memory (GMT, NULL, T->ogr->n_aspatial, char *);
		T->ogr->type = GMT_memory (GMT, NULL, T->ogr->n_aspatial, unsigned int);
		T->ogr->dvalue = GMT_memory (GMT, NULL, T->ogr->n_aspatial, double);
		for (k = 0; k < T->ogr->n_aspatial; k++) {
			T->ogr->name[k] = strdup (GMT->common.a.name[k]);
			T->ogr->type[k] = GMT->common.a.type[k];
		}
		for (seg = 0; seg < T->n_segments; seg++) {	/* For each segment in the table */
			S = T->segment[seg];
			gmt_alloc_ogr_seg (GMT, S, T->ogr->n_aspatial);	/* Copy over any feature-specific values */
			for (k = 0; k < T->ogr->n_aspatial; k++) {	/* For each column to turn into a constant aspatial value */
				col = GMT->common.a.col[k];
				if (col < 0) continue;	/* Multisegment header entry instead */
				for (row = 1, stop = false; !stop && row < S->n_rows; ++row) {
					if (!doubleAlmostEqualZero (S->coord[col][row], S->coord[col][row-1]))
						stop = true;
				}
				if (stop) {
					////GMT_Report (GMT->parent, GMT_MSG_NORMAL, "The -a option specified a constant column but its contents vary!\n");
					return (GMT_RUNTIME_ERROR);
				}
				else
					S->ogr->dvalue[k] = S->coord[col][0];
			}
		}
		/* OK, successfully passed the constant column tests, if any */
		for (seg = col = 0; seg < T->n_segments; seg++) {	/* Free up columns now stored as aspatial values */
			S = T->segment[seg];
			for (k = 0; k < T->ogr->n_aspatial; k++) if (GMT->common.a.col[k] > 0) GMT_free (GMT, S->coord[GMT->common.a.col[k]]);
			for (k = col = 0; k < T->n_columns; k++) {
				while (!S->coord[k]) k++;	/* Next available column */
				S->coord[col++] = S->coord[k];	/* Update pointers */
			}
			S->n_columns = col;	/* May have lost some columns now */
		}
		T->n_columns = D->n_columns = col;	/* May have lost some columns now */
	}
	if (T->ogr->geometry == GMT_IS_POLYGON || T->ogr->geometry == GMT_IS_MULTIPOLYGON) {	/* Must check consistency */
		for (seg = 0; seg < T->n_segments; seg++) {	/* For each segment in the table */
			if ((T->ogr->geometry == GMT_IS_POLYGON || T->ogr->geometry == GMT_IS_MULTIPOLYGON) && GMT_polygon_is_open (GMT, T->segment[seg]->coord[GMT_X], T->segment[seg]->coord[GMT_Y], T->segment[seg]->n_rows)) {
				////GMT_Report (GMT->parent, GMT_MSG_NORMAL, "The -a option specified [M]POLY but open segments were detected!\n");
				GMT_Destroy_Data (GMT->parent, &D[GMT_OUT]);
				return (GMT_RUNTIME_ERROR);
			}
			gmt_alloc_ogr_seg (GMT, T->segment[seg], T->ogr->n_aspatial);	/* Copy over any feature-specific values */
			T->segment[seg]->ogr->pol_mode = GMT_IS_PERIMETER;
			GMT_set_seg_minmax (GMT, T->segment[seg]);	/* Make sure min/max are set per polygon */

		}
		/* OK, they are all polygons.  Determine any polygon holes: if a point is fully inside another polygon (not on the edge) */
		for (seg1 = 0; seg1 < T->n_segments; seg1++) {	/* For each segment in the table */
			for (seg2 = seg1 + 1; seg2 < T->n_segments; seg2++) {	/* For each segment in the table */
				if (GMT_inonout (GMT, T->segment[seg1]->coord[GMT_X][0], T->segment[seg1]->coord[GMT_Y][0], T->segment[seg2]) == GMT_INSIDE) T->segment[seg1]->ogr->pol_mode = GMT_IS_HOLE;
				if (GMT_inonout (GMT, T->segment[seg2]->coord[GMT_X][0], T->segment[seg2]->coord[GMT_Y][0], T->segment[seg1]) == GMT_INSIDE) T->segment[seg2]->ogr->pol_mode = GMT_IS_HOLE;
			}
		}
	}
	if (T->ogr->geometry > GMT_IS_POINT && T->ogr->geometry != GMT_IS_MULTIPOINT) {	/* Must check for Dateline crossings */
		unsigned int n_split;
		uint64_t n_segs = T->n_segments;
		struct GMT_DATASEGMENT **L = NULL;

		for (seg = 0; seg < T->n_segments; seg++) {	/* For each segment in the table */
			if (!GMT_crossing_dateline (GMT, T->segment[seg])) continue;	/* GIS-safe feature! */
			////GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "Feature %" PRIu64 " crosses the Dateline\n", seg);
			if (!GMT->common.a.clip) continue;	/* Not asked to clip */
			/* Here we must split into east and west part(s) */
			if (T->ogr->geometry == GMT_IS_POLYGON || T->ogr->geometry == GMT_IS_MULTIPOLYGON) {	/* Clipping must add dateline segments */
				/* Clip into two closed polygons.  Eventually, perhaps return more (eliminate bridges) */
				n_split = GMT_split_poly_at_dateline (GMT, T->segment[seg], &L);
			}
			else {	/* Clipping just needs to add crossing points */
				/* Truncate into two or more line segments */
				n_split = GMT_split_line_at_dateline (GMT, T->segment[seg], &L);
			}
			T->segment = GMT_memory (GMT, T->segment, n_segs + n_split - 1, struct GMT_DATASEGMENT *);	/* Allow more space for new segments */
			GMT_free_segment (GMT, &(T->segment[seg]), D->alloc_mode);	/* Delete the old one */
			T->segment[seg] = L[0];			/* Hook in the first replacement */
			for (k = 1; k < n_split; k++) T->segment[n_segs++] = L[k];	/* Add the remaining segments to the end */
			GMT_free (GMT, L);
		}
		D->n_segments = T->n_segments = n_segs;	/* Update number of segments */

	}
	GMT->current.io.geo.range = GMT_IS_M180_TO_P180_RANGE;	/* Select the -180/180 output range format */
	return (0);
}

int GMT_write_dataset (struct GMT_CTRL *GMT, void *dest, unsigned int dest_type, struct GMT_DATASET *D, bool use_GMT_io, int table)
{	/* Writes an entire data set to file or stream */
	unsigned int tbl, u_table;
	bool close_file = false;
	int error, append = 0;
	int *fd = NULL;
	char file[GMT_BUFSIZ] = {""}, tmpfile[GMT_BUFSIZ] = {""}, open_mode[4] = {""}, *out_file = tmpfile;
	FILE *fp = NULL;

	if (dest_type == GMT_IS_FILE && dest && ((char *)dest)[0] == '>') append = 1;	/* Want to append to existing file */
	if (use_GMT_io)	/* Use GMT->current.io.info settings to determine if input is ascii/binary, else it defaults to ascii */
		strcpy (open_mode, (append) ? GMT->current.io.a_mode : GMT->current.io.w_mode);
	else			/* Force ASCII mode */
		strcpy (open_mode, (append) ? "a" : "w");

	/* Convert any destination type to stream */

	switch (dest_type) {
		case GMT_IS_FILE:	/* dest is a file name */
			strncpy (file, (const char *)dest, GMT_BUFSIZ);
			if (D->io_mode < GMT_WRITE_TABLE) {	/* Only need one destination */
				if ((fp = GMT_fopen (GMT, &file[append], open_mode)) == NULL) {
					////GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Cannot open file %s\n", &file[append]);
					return (EXIT_FAILURE);
				}
				close_file = true;	/* We only close files we have opened here */
				////GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "Write Data Table to file %s\n", &file[append]);
			}
			break;
		case GMT_IS_STREAM:	/* Open file pointer given, just copy */
			fp = (FILE *)dest;
			if (fp == NULL) fp = GMT->session.std[GMT_OUT];	/* Default destination */
			if (fp == GMT->session.std[GMT_OUT])
				strcpy (file, "<stdout>");
			else
				strcpy (file, "<output stream>");
			////GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "Write Data Table to %s\n", file);
			break;
		case GMT_IS_FDESC:		/* Open file descriptor given, just convert to file pointer */
			fd = dest;
			if (fd && (fp = fdopen (*fd, open_mode)) == NULL) {
				////GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Cannot convert file descriptor %d to stream in GMT_write_table\n", *fd);
				return (EXIT_FAILURE);
			}
			if (fd == NULL) fp = GMT->session.std[GMT_OUT];	/* Default destination */
			if (fp == GMT->session.std[GMT_OUT])
				strcpy (file, "<stdout>");
			else
				strcpy (file, "<output file descriptor>");
			////GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "Write Data Table to %s\n", file);
			break;
		default:
			////GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Unrecognized source type %d in GMT_write_table\n", dest_type);
			return (EXIT_FAILURE);
			break;
	}

	if (D->io_mode == GMT_WRITE_OGR && gmt_prep_ogr_output (GMT, D)) {	/* Must preprocess aspatial information and set metadata */
		////GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Failed to prepare for OGR output formatting\n");
		return (EXIT_FAILURE);
	}
	for (tbl = 0; tbl < D->n_tables; tbl++) {
		if (table != GMT_NOTSET && (u_table = table) != tbl) continue;	/* Selected a specific table */
		/*if (D->io_mode > GMT_WRITE_TABLE) {	 Write segments to separate files; must pass original file name in case a template
			if ((error = GMT_write_table (GMT, dest, GMT_IS_FILE, D->table[tbl], use_GMT_io, D->io_mode))) return (error);
		}*/
		else if (D->io_mode == GMT_WRITE_TABLE) {	/* Must write this table a its own file */
			if (D->table[tbl]->file[GMT_OUT])
				out_file = D->table[tbl]->file[GMT_OUT];
			else
				sprintf (tmpfile, file, D->table[tbl]->id);
			////GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "Write Data Table to %s\n", out_file);
			//if ((error = GMT_write_table (GMT, out_file, GMT_IS_FILE, D->table[tbl], use_GMT_io, D->io_mode))) return (error);
		}
		else {	/* Write to stream we set up earlier */
			//if ((error = GMT_write_table (GMT, fp, GMT_IS_STREAM, D->table[tbl], use_GMT_io, D->io_mode))) return (error);
		}
	}

	if (close_file) GMT_fclose (GMT, fp);

	return (0);	/* OK status */
}

static const char *GMT_type[GMT_N_TYPES] = {"byte", "byte", "integer", "integer", "integer", "integer", "integer", "integer", "double", "double", "string", "datetime"};

/* This version of fgets will check for input record truncation, that is,
 * the input record is longer than the given size.  Since calls to GMT_fgets
 * ASSUME they get a logical record, we will give a warning if truncation
 * occurs and read until we have consumed the linefeed, thus making the
 * i/o machinery ready for the next logical record.
 */

char *GMT_fgets (struct GMT_CTRL *GMT, char *str, int size, FILE *stream)
{
	
	str[size-2] = '\0'; /* Set last but one record to 0 */
	//if(stream!=NULL)
	//{
	//	printf("\n\n\n\nher eis the size %d len %d", size, strlen(str));
	//	fflush(stdout);
	//}
	if (!fgets (str, size, stream))
	{
	
		return (NULL); /* Got nothing */
	}

	/* fgets will always set str[size-1] = '\0' if more data than str can handle is found.
	 * Thus, we examine str[size-2].  If this is neither '\0' nor '\n' then we have only
	 * read a portion of a logical record that is longer than size.
	 */
	
	if (!(str[size-2] == '\0' || str[size-2] == '\n')) {
		/* Only got part of a record */
		int c, n = 0;
		/* Read char-by-char until newline is consumed */
		while ((c = fgetc (stream)) != '\n' && c != EOF)
			++n;
		if (c == '\n')
			/* We expect fgets to retain '\n', so add it */
			str[size-2] = '\n';
		else
			/* EOF without '\n' */
			--n;
		/* This will report wrong lengths if last line has no '\n' but we don't care */
		////GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Long input record (%d bytes) was truncated to first %d bytes!\n", size+n, size-2);
	}
	return (str);
}

bool gmt_ogr_parser (struct GMT_CTRL *GMT, char *record)
{	/* Parsing of the GMT/OGR vector specification (v 1.0). See Appendix R */
	return (GMT->current.io.ogr_parser (GMT, record));	/* We call either the header or data parser depending on pointer */
}

void GMT_set_segmentheader (struct GMT_CTRL *GMT, int direction, bool true_false)
{	/* Enable/Disable multi-segment headers for either input or output */
	GMT->current.io.multi_segments[direction] = true_false;
}

static inline void gmt_update_prev_rec (struct GMT_CTRL *GMT, uint64_t n_use) {
	/* Update previous record before reading the new record*/
	if (GMT->current.io.need_previous) GMT_memcpy (GMT->current.io.prev_rec, GMT->current.io.curr_rec, n_use, double);
}

unsigned int gmt_assign_aspatial_cols (struct GMT_CTRL *GMT)
{	/* This function will load input columns with aspatial data as requested by -a.
 	 * It will then handle any possible -i scalings/offsets as well for those columns.
 	 * This is how the @D values end up in the input data record we read. */

	unsigned int k, n;
	double value;
	if (GMT->current.io.ogr != GMT_OGR_TRUE) return (0);	/* No point checking further since file is not GMT/OGR */
	for (k = n = 0; k < GMT->common.a.n_aspatial; k++) {	/* For each item specified in -a */
		if (GMT->common.a.col[k] < 0) continue;	/* Not meant for data columns */
		value = GMT->current.io.OGR->dvalue[GMT->common.a.ogr[k]];
		gmt_convert_col (GMT->current.io.col[GMT_IN][GMT->common.a.col[k]], value);
		GMT->current.io.curr_rec[GMT->common.a.col[k]] = value;
		n++;
	}
	return (n);
}

/* Like strsep but ignores empty fields */
char *strsepz (char **stringp, const char *delim) {
	char *c;
	while ( (c = strsep(stringp, delim)) != NULL && *c == '\0' );
	return c;
}

static inline uint64_t gmt_n_cols_needed_for_gaps (struct GMT_CTRL *GMT, uint64_t n) {
	/* Return the actual items needed (which may be more than n if gap testing demands it) */
	if (GMT->common.g.active) return (MAX (n, GMT->common.g.n_col));	/* n or n_col (if larger) */
	return (n);	/* No gap checking, n it is */
}

bool gmt_gap_detected (struct GMT_CTRL *GMT)
{	/* Determine if two points are "far enough apart" to constitude a data gap and thus "pen up" */
	uint64_t i;

	if (!GMT->common.g.active || GMT->current.io.pt_no == 0) return (false);	/* Not active or on first point in a segment */
	/* Here we must determine if any or all of the selected gap criteria [see gmt_set_gap_param] are met */
	for (i = 0; i < GMT->common.g.n_methods; i++) {	/* Go through each criterion */
		if ((GMT->common.g.get_dist[i] (GMT, GMT->common.g.col[i]) > GMT->common.g.gap[i]) != GMT->common.g.match_all) return (!GMT->common.g.match_all);
	}
	return (GMT->common.g.match_all);
}

int gmt_set_gap (struct GMT_CTRL *GMT) {	/* Data gaps are special since there is no multiple-segment header flagging the gap; thus next time the record is already read */
	GMT->current.io.status = GMT_IO_GAP;
	GMT->current.io.seg_no++;
	//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "Data gap detected via -g; Segment header inserted near/at line # %" PRIu64 "\n", GMT->current.io.rec_no);
	sprintf (GMT->current.io.segment_header, "Data gap detected via -g; Segment header inserted");
	return (0);
}

void gmt_adjust_periodic (struct GMT_CTRL *GMT) {
	while (GMT->current.io.curr_rec[GMT_X] > GMT->common.R.wesn[XHI] && (GMT->current.io.curr_rec[GMT_X] - 360.0) >= GMT->common.R.wesn[XLO]) GMT->current.io.curr_rec[GMT_X] -= 360.0;
	while (GMT->current.io.curr_rec[GMT_X] < GMT->common.R.wesn[XLO] && (GMT->current.io.curr_rec[GMT_X] + 360.0) <= GMT->common.R.wesn[XLO]) GMT->current.io.curr_rec[GMT_X] += 360.0;
	/* If data is not inside the given range it will satisfy (lon > east) */
	/* Now it will be outside the region on the same side it started out at */
}

void gmt_adjust_projected (struct GMT_CTRL *GMT) {
	/* Case of incoming projected map coordinates that we wish to rever to lon/lat */
	if (GMT->current.proj.inv_coord_unit != GMT_IS_METER) {	/* Must first scale to meters */
		GMT->current.io.curr_rec[GMT_X] *= GMT->current.proj.m_per_unit[GMT->current.proj.inv_coord_unit];
		GMT->current.io.curr_rec[GMT_Y] *= GMT->current.proj.m_per_unit[GMT->current.proj.inv_coord_unit];
	}
	(*GMT->current.proj.inv) (GMT, &GMT->current.io.curr_rec[GMT_X], &GMT->current.io.curr_rec[GMT_Y], GMT->current.io.curr_rec[GMT_X], GMT->current.io.curr_rec[GMT_Y]);
}

void * gmt_ascii_input (struct GMT_CTRL *GMT, FILE *fp, uint64_t *n, int *status)
{
	uint64_t pos, col_no = 0, col_pos, n_convert, n_ok = 0, kind, add, n_use = 0;
	int64_t in_col;
	bool done = false, bad_record, set_nan_flag = false;
	char line[GMT_BUFSIZ] = {""}, *p = NULL, *token, *stringp;
	double val;

	/* gmt_ascii_input will skip blank lines and shell comment lines which start
	 * with #.  Fields may be separated by spaces, tabs, or commas.  The routine returns
	 * the actual number of items read [or 0 for segment header and -1 for EOF]
	 * If *n is passed as GMT_BUFSIZ it will be reset to the actual number of fields.
	 * If gap checking is in effect and one of the checks involves a column beyond
	 * the ones otherwise needed by the program we extend the reading so we may
	 * examin the column needed in the gap test.
	 * *status returns the number of fields read, 0 for header records, -1 for EOF.
	 * We return NULL (headers or errors) or pointer to GMT->current.io.curr_rec.
	 */
	//printf("file : %s line : %d func: %s %d \n",__FILE__,__LINE__,__func__, sizeof(line));
	while (!done) {	/* Done becomes true when we successfully have read a data record */

		/* First read until we get a non-blank, non-comment record, or reach EOF */
		//printf("-----------> came here %d\n", strlen(line));
		GMT->current.io.rec_no++;		/* Counts up, regardless of what this record is (data, junk, segment header, etc) */
		GMT->current.io.rec_in_tbl_no++;	/* Counts up, regardless of what this record is (data, junk, segment header, etc) */
		//printf("file : %s line : %d func: %s %d \n",__FILE__,__LINE__,__func__, sizeof(line));
		if (GMT->current.setting.io_header[GMT_IN] && GMT->current.io.rec_in_tbl_no <= GMT->current.setting.io_n_header_items) {	/* Must treat first io_n_header_items as headers */
			//printf("file : %s line : %d func: %s %d \n",__FILE__,__LINE__,__func__, sizeof(line));
			p = GMT_fgets (GMT, line, GMT_BUFSIZ, fp);	/* Get the line */
			if (GMT->common.h.mode == GMT_COMMENT_IS_RESET) continue;	/* Simplest way to replace headers on output is to ignore them on input */
			strncpy (GMT->current.io.current_record, line, GMT_BUFSIZ);
			GMT->current.io.status = GMT_IO_TABLE_HEADER;
			//GMT->current.setting.io_header[GMT_OUT] = true;	/* Turn on table headers on output PW: No! If we get here via -hi then no header output was requested */
			*status = 0;
			return (NULL);
		}
		/* Here we are done with any header records implied by -h */
		if (GMT->current.setting.io_blankline[GMT_IN]) {	/* Treat blank lines as segment markers, so only read a single line */
			p = GMT_fgets (GMT, line, GMT_BUFSIZ, fp);
			GMT->current.io.rec_no++, GMT->current.io.rec_in_tbl_no++;
		}
		else {	/* Default is to skip all blank lines until we get something else (or hit EOF) */
			//printf("file : %s line : %d func: %s %d %p \n",__FILE__,__LINE__,__func__, sizeof(line),fp);
			//if(fp)
				//printf("===========> not null \n");
			//printf("-----------> came here %d\n", strlen(line));
			while ((p = GMT_fgets (GMT, line, GMT_BUFSIZ, fp)) && GMT_is_a_blank_line (line)) GMT->current.io.rec_no++, GMT->current.io.rec_in_tbl_no++;
		}
		if (!p) {	/* Ran out of records, which can happen if file ends in a comment record */
			GMT->current.io.status = GMT_IO_EOF;
			if (GMT->current.io.give_report && GMT->current.io.n_bad_records) {	/* Report summary and reset counters */
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "This file had %" PRIu64 " data records with invalid x and/or y values\n", GMT->current.io.n_bad_records);
				GMT->current.io.n_bad_records = GMT->current.io.pt_no = GMT->current.io.n_clean_rec = 0;
				GMT->current.io.rec_no = GMT->current.io.rec_in_tbl_no = 0;
			}
			*status = -1;
			return (NULL);
		}
		if (gmt_ogr_parser (GMT, line)) continue;	/* If we parsed a GMT/OGR record we must go up to top of loop and get the next record */
		if (line[0] == '#') {	/* Got a file header, copy it and return */
			if (GMT->common.h.mode == GMT_COMMENT_IS_RESET) continue;	/* Simplest way to replace headers on output is to ignore them on input */
			strncpy (GMT->current.io.current_record, line, GMT_BUFSIZ);
			GMT->current.io.status = GMT_IO_TABLE_HEADER;
			*status = 0;
			return (NULL);
		}

		if ((kind = gmt_is_segment_header (GMT, line))) {	/* Got a segment header, take action and return */
			GMT->current.io.status = GMT_IO_SEGMENT_HEADER;
			GMT_set_segmentheader (GMT, GMT_OUT, true);	/* Turn on segment headers on output */
			GMT->current.io.seg_no++;
			GMT->current.io.segment_header[0] = '\0';
			if (kind == 1) {
				/* Just save the header content, not the marker and leading whitespace */
				strncpy (GMT->current.io.segment_header, GMT_trim_segheader (GMT, line), GMT_BUFSIZ);
			}
			/* else we got a segment break instead - and header was set to NULL */
			*status = 0;
			return (NULL);
		}

		/* Here we know we are processing a data record */

		if (GMT->common.a.active && GMT->current.io.ogr == GMT_OGR_FALSE) {	/* Cannot give -a and not be reading an OGR/GMT file */
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Aspatial associations set with -a but input file is not in OGR/GMT format!\n");
			GMT_exit (GMT, EXIT_FAILURE); return NULL;
		}

		n_use = gmt_n_cols_needed_for_gaps (GMT, *n);	/* Gives the actual columns we need (which may > *n if gap checking is active; if gap check we also update prev_rec) */
		gmt_update_prev_rec (GMT, n_use);

		/* First chop off trailing whitespace and commas */

		GMT_strstrip (line, false); /* Eliminate DOS endings and trailing white space, add linefeed */

		bad_record = set_nan_flag = false;		/* Initialize flags */
		strncpy (GMT->current.io.current_record, line, GMT_BUFSIZ);	/* Keep copy of current record around */
		col_no = pos = n_ok = 0;			/* Initialize counters */
		in_col = -1;					/* Since we will increment right away inside the loop */

		stringp = line;
		while (!bad_record && col_no < n_use && (token = strsepz (&stringp, " \t,")) != NULL) {	/* Get one field at the time until we run out or have issues */
			++in_col;	/* This is the actual column number in the input file */
			if (GMT->common.i.active) {	/* Must do special column-based processing since the -i option was set */
				if (GMT->current.io.col_skip[in_col]) continue;		/* Just skip and not even count this column */
				col_pos = GMT->current.io.col[GMT_IN][col_no].order;	/* Which data column will receive this value */
			}
			else				/* Default column order */
				col_pos = col_no;
			if ((n_convert = GMT_scanf (GMT, token, GMT->current.io.col_type[GMT_IN][col_pos], &val)) == GMT_IS_NAN) {	/* Got a NaN or it failed to decode the string */
				if (GMT->current.setting.io_nan_records || !GMT->current.io.skip_if_NaN[col_pos]) {	/* This field (or all fields) can be NaN so we pass it on */
					GMT->current.io.curr_rec[col_pos] = GMT->session.d_NaN;
					n_ok++;	/* Since NaN is considered an OK result */
				}
				else	/* Cannot have NaN in this column, flag record as bad */
					bad_record = true;
				if (GMT->current.io.skip_if_NaN[col_pos]) set_nan_flag = true;	/* Flag that we found NaN in a column that means we should skip */
			}
			else {					/* Successful decode, assign the value to the input array */
				gmt_convert_col (GMT->current.io.col[GMT_IN][col_no], val);
				GMT->current.io.curr_rec[col_pos] = val;
				n_ok++;
			}
			col_no++;		/* Count up number of columns found */
		}
		if ((add = gmt_assign_aspatial_cols (GMT))) {	/* We appended <add> columns given via aspatial OGR/GMT values */
			col_no += add;
			n_ok += add;
		}
		if (bad_record) {	/* This record failed our test and had NaNs */
			GMT->current.io.n_bad_records++;
			if (GMT->current.io.give_report && (GMT->current.io.n_bad_records == 1)) {	/* Report 1st occurrence of bad record */
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Encountered first invalid ASCII data record near/at line # %" PRIu64 "\n", GMT->current.io.rec_no);
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Likely causes:\n");
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "(1) Invalid x and/or y values, i.e. NaNs or garbage in text strings.\n");
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "(2) Incorrect data type assumed if -J, -f are not set or set incorrectly.\n");
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "(3) The -: switch is implied but not set.\n");
			}
		}
		else if (GMT->current.io.skip_duplicates && GMT->current.io.pt_no) {	/* Test to determine if we should skip repeated duplicate records with same x,y */
			done = !(GMT->current.io.curr_rec[GMT_X] == GMT->current.io.prev_rec[GMT_X] && GMT->current.io.curr_rec[GMT_Y] == GMT->current.io.prev_rec[GMT_Y]);	/* Yes, duplicate */
		}
		else
			done = true;	/* Success, we can get out of this loop and return what we got */
	}
	GMT->current.io.status = (GMT->current.io.read_mixed || n_ok == n_use || *n == GMT_MAX_COLUMNS) ? 0 : GMT_IO_MISMATCH;	/* Hopefully set status to 0 (OK) */
	if (*n == GMT_MAX_COLUMNS) *n = n_ok;							/* Update the number of expected fields */
	if (GMT_REC_IS_ERROR (GMT)) 
		{
		printf( "Mismatch between actual (%d) and expected (%d) fields near line %" PRIu64 "\n", col_no, *n, GMT->current.io.rec_no);
		
	}
	if (GMT->current.setting.io_lonlat_toggle[GMT_IN] && col_no >= 2) double_swap (GMT->current.io.curr_rec[GMT_X], GMT->current.io.curr_rec[GMT_Y]);	/* Got lat/lon instead of lon/lat */
	if (GMT->current.proj.inv_coordinates) gmt_adjust_projected (GMT);	/* Must apply inverse projection to get lon, lat */
	if (GMT->current.io.col_type[GMT_IN][GMT_X] & GMT_IS_GEO) gmt_adjust_periodic (GMT);	/* Must account for periodicity in 360 as per current rule*/

	if (gmt_gap_detected (GMT)) {*status = gmt_set_gap (GMT); return (GMT->current.io.curr_rec); }	/* A gap between this an previous record was detected (see -g) so we set status and return 0 */

	GMT->current.io.pt_no++;	/* Got a valid data record (which is true even if it was a gap) */
	*status = (int)n_ok;			/* Return the number of fields successfully read */
	if (set_nan_flag) {
		GMT->current.io.status |= GMT_IO_NAN;	/* Say we found NaNs */
		return (GMT->current.io.curr_rec);	/* Pass back pointer to data array */
	}
	return ((GMT->current.io.status) ? NULL : GMT->current.io.curr_rec);	/* Pass back pointer to data array */
}

bool gmt_skip_output (struct GMT_CTRL *GMT, double *cols, uint64_t n_cols)
{	/* Consult the -s[<cols>][a|r] setting and the cols values to determine if this record should be output */
	uint64_t c, n_nan;

	if (n_cols > GMT_MAX_COLUMNS) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Number of output data columns (%d) exceeds limit (GMT_MAX_COLUMS = %d)\n", n_cols, GMT_MAX_COLUMNS);
		return (true);	/* Skip record since we cannot access that many columns */
	}
	if (GMT->current.setting.io_nan_mode == GMT_IO_NAN_OK) return (false);				/* Normal case; output the record */
	if (GMT->current.setting.io_nan_mode == GMT_IO_NAN_ONE) {	/* -sa: Skip records if any NaNs are found */
		for (c = 0; c < n_cols; c++) if (GMT_is_dnan (cols[c]))  return (true);	/* Found a NaN so we skip */
		return (false);	/* No NaNs, output record */
	}
	for (c = n_nan = 0; c < GMT->current.io.io_nan_ncols; c++) {			/* Check each of the specified columns set via -s */
		if (GMT->current.io.io_nan_col[c] >= n_cols) continue;			/* Input record does not have this column */
		if (GMT_is_dnan (cols[GMT->current.io.io_nan_col[c]])) n_nan++;		/* Count the nan columns found */
	}
	if (n_nan < GMT->current.io.io_nan_ncols  && GMT->current.setting.io_nan_mode == GMT_IO_NAN_KEEP) return (true);	/* Skip records if -sr and not enough NaNs found */
	if (n_nan == GMT->current.io.io_nan_ncols && GMT->current.setting.io_nan_mode == GMT_IO_NAN_SKIP) return (true);	/* Skip records if -s and NaNs in specified columns */
	return (false);	/* No match, output record */
}

/* This is the lowest-most input function in GMT.  All ASCII table data are read via
 * gmt_ascii_input.  Changes here affect all programs that read such data. */

int GMT_ascii_output (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *ptr)
{
	uint64_t i, col, last, n_out;
	int e = 0, wn = 0;
	double val;

	//printf("%f   ------ \n ",*ptr);

	//printf("getting 11111111111111 1 in ascii: %f 2: %f, 3:%f, 4:%f \n", ptr[0], ptr[1], ptr[2], ptr[3]);
	
	if (gmt_skip_output (GMT, ptr, n)) return (0);	/* Record was skipped via -s[a|r] */
	n_out = (GMT->common.o.active) ? GMT->common.o.n_cols : n;
	//printf("nishita >>>>>>>>>>>>>>>>>>>>>>>>>>>>n %d n_out %d \n\n",n,n_out);
	last = n_out - 1;				/* Last filed, need to output linefeed instead of delimiter */

	for (i = 0; i < n_out && e >= 0; i++) {		/* Keep writing all fields unless there is a read error (e == -1) */
		if (GMT->common.o.active)	/* Which data column to pick */
			col = GMT->current.io.col[GMT_OUT][i].col;
		else if (GMT->current.setting.io_lonlat_toggle[GMT_OUT] && i < 2)
			col = 1 - i;	/* Write lat/lon instead of lon/lat */
		else
			col = i;	/* Just goto next column */
		//printf("col %d ",col);
		val = (col >= n) ? GMT->session.d_NaN : ptr[col];	/* If we request beyond length of array, return NaN */
		e = GMT_ascii_output_col (GMT, fp, val, col);	/* Write one item without any separator at the end */

		if (i == last)					/* This is the last field, must add newline */
			putc ('\n', fp);
		else if (GMT->current.setting.io_col_separator[0])		/* Not last field, and a separator is required */
			fprintf (fp, "%s", GMT->current.setting.io_col_separator);

		wn += e;
	}
	//printf("nishita 22>>>>>>>>>>>>>>>>>>>>>>>>>>>> \n\n");
	return ((e < 0) ? e : wn);
}



bool gmt_ogr_data_parser (struct GMT_CTRL *GMT, char *record)
{	/* Parsing of the GMT/OGR vector specification (v 1.0) for data feature records.
 	 * We KNOW GMT->current.io.ogr == GMT_OGR_TRUE, i.e., current file is a GMT/OGR file.
	 * We also KNOW that GMT->current.io.OGR has been allocated by gmt_ogr_header_parser.
	 * For GMT/OGR files we must parse and store the metadata in GMT->current.io.OGR,
	 * from where higher-level functions can access it.  GMT_End_IO will free the structure.
	 * This function returns true if we parsed a GMT/OGR record and false otherwise.
	 * If we encounter a parsing error we stop parsing any further by setting GMT->current.io.ogr = GMT_OGR_FALSE.
	 * We loop until all @<info> tags have been processed on this record.
	 */

	unsigned int n_aspatial;
	bool quote;
	char *p = NULL;
	struct GMT_OGR *S = NULL;

	if (record[0] != '#') return (false);			/* Not a comment record so no point looking further */
	if (!(p = strchr (record, '@'))) return (false);	/* Not an OGR/GMT record since @ was not found */

	/* Here we are reasonably sure that @? strings are OGR/GMT feature specifications */

	GMT_chop (record);	/* Get rid of linefeed etc */

	S = GMT->current.io.OGR;	/* Set S shorthand */
	quote = false;

	while (*p == '@') {
		++p;	/* Move to first char after @ */
		switch (p[0]) {	/* These are the feature tags only: @D, @P, @H */
			case 'D':	/* Aspatial data values, store in segment header  */
				//if (!S->geometry) { GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Bad OGR/GMT: @D given but no geometry set\n"); return (false);}
				//n_aspatial = gmt_ogr_decode_aspatial_values (GMT, &p[1], S);
				//if (S->n_aspatial != n_aspatial) {
					//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "OGR/GMT: Some @D items not specified (set to NULL)\n");
				//}
				break;

			case 'P':	/* Polygon perimeter, store in segment header  */
				//if (!(S->geometry == GMT_IS_POLYGON || S->geometry == GMT_IS_MULTIPOLYGON)) {
				//	GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Bad OGR/GMT: @P only valid for polygons\n");
					//GMT->current.io.ogr = GMT_OGR_FALSE;
					//return (false);
				//}
				//S->pol_mode = GMT_IS_PERIMETER;
				break;

			case 'H':	/* Polygon hole, store in segment header  */
				//if (!(S->geometry == GMT_IS_POLYGON || S->geometry == GMT_IS_MULTIPOLYGON)) {
					//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Bad OGR/GMT: @H only valid for polygons\n");
					//GMT->current.io.ogr = GMT_OGR_FALSE;
					//return (false);
				//}
				//S->pol_mode = GMT_IS_HOLE;
				//break;

			default:	/* Bad OGR record? */
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Bad OGR/GMT: Cannot have @%c after FEATURE_DATA\n", p[0]);
				//GMT->current.io.ogr = GMT_OGR_FALSE;
				break;
		}
		while (*p && (quote || *p != '@')) if (*p++ == '\"') quote = !quote;	/* Wind to next @ except skip if inside double quotes */
	}
	return (true);
}

int gmt_get_ogr_id (struct GMT_OGR *G, char *name)
{
	unsigned int k;
	for (k = 0; k < G->n_aspatial; k++) if (!strcmp (name, G->name[k])) return (k);
	return (GMT_NOTSET);
}

void gmt_align_ogr_values (struct GMT_CTRL *GMT)
{	/* Simplify aspatial data grabbing when -a is used */
	unsigned int k;
	int id;
	if (!GMT->common.a.active) return;	/* Nothing selected with -a */
	for (k = 0; k < GMT->common.a.n_aspatial; k++) {	/* Process the requested columns */
		id = gmt_get_ogr_id (GMT->current.io.OGR, GMT->common.a.name[k]);	/* See what order in the OGR struct this -a column appear */
		GMT->common.a.ogr[k] = id;
	}
}

void gmt_copy_and_truncate (char *out, char *in)
{	/* Duplicate in to out, then find the first space not inside quotes and truncate string there */
	bool quote = false;
	while (*in && (quote || *in != ' ')) {
		*out++ = *in;	/* Copy char */
		if (*in++ == ' ') quote = !quote;	/* Wind to next space except skip if inside double quotes */
	}
	*out = '\0';	/* Terminate string */
}

int gmt_ogr_get_type (char *item)
{
	if (!strcmp (item, "double") || !strcmp (item, "DOUBLE")) return (GMT_DOUBLE);
	if (!strcmp (item, "float") || !strcmp (item, "FLOAT")) return (GMT_FLOAT);
	if (!strcmp (item, "integer") || !strcmp (item, "INTEGER")) return (GMT_INT);
	if (!strcmp (item, "char") || !strcmp (item, "CHAR")) return (GMT_CHAR);
	if (!strcmp (item, "string") || !strcmp (item, "STRING")) return (GMT_TEXT);
	if (!strcmp (item, "datetime") || !strcmp (item, "DATETIME")) return (GMT_DATETIME);
	if (!strcmp (item, "logical") || !strcmp (item, "LOGICAL")) return (GMT_UCHAR);
	return (GMT_NOTSET);
}

unsigned int gmt_ogr_decode_aspatial_types (struct GMT_CTRL *GMT, char *record, struct GMT_OGR *S)
{	/* Parse @T aspatial types; this is done once per dataset and follows @N */
	unsigned int pos = 0, col = 0;
	size_t n_alloc;
	char buffer[GMT_BUFSIZ] = {""}, p[GMT_BUFSIZ];

	n_alloc = (S->type) ? GMT_BUFSIZ : 0;
	gmt_copy_and_truncate (buffer, record);
	while ((GMT_strtok (buffer, "|", &pos, p))) {
		if (col >= S->n_aspatial) {
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Bad OGR/GMT: @T record has more items than declared by @N\n");
			continue;
		}
		if (col == n_alloc) S->type = GMT_memory (GMT, S->type, n_alloc += GMT_TINY_CHUNK, unsigned int);
		S->type[col++] = gmt_ogr_get_type (p);
	}
	if (n_alloc < GMT_BUFSIZ && col < n_alloc) S->type = GMT_memory (GMT, S->type, col, unsigned int);
	return (col);
}

unsigned int gmt_ogr_decode_aspatial_names (struct GMT_CTRL *GMT, char *record, struct GMT_OGR *S)
{	/* Decode @N aspatial names; this is done once per dataset */
	unsigned int pos = 0, col = 0;
	size_t n_alloc;
	char buffer[GMT_BUFSIZ] = {""}, p[GMT_BUFSIZ] = {""};

	n_alloc = (S->type) ? GMT_BUFSIZ : 0;
	gmt_copy_and_truncate (buffer, record);
	while ((GMT_strtok (buffer, "|", &pos, p))) {
		if (col == n_alloc) S->name = GMT_memory (GMT, S->name, n_alloc += GMT_TINY_CHUNK, char *);
		S->name[col++] = strdup (p);
	}
	if (n_alloc < GMT_BUFSIZ && col < n_alloc) S->name = GMT_memory (GMT, S->name, col, char *);
	return (col);
}

bool gmt_ogr_header_parser (struct GMT_CTRL *GMT, char *record)
{	/* Parsing of the GMT/OGR vector specification (v 1.0).
 	 * GMT->current.io.ogr can have three states:
	 *	GMT_OGR_UNKNOWN (-1) if not yet set [this is how it is initialized in GMTAPI_Begin_IO].
	 *	GMT_OGR_FALSE    (0) if file has been determined NOT to be a GMT/OGR file.
	 *	GMT_OGR_TRUE    (+1) if it has met the criteria and is a GMT/OGR file.
	 * For GMT/OGR files we must parse and store the metadata in GMT->current.io.OGR,
	 * from where higher-level functions can access it.  GMT_End_IO will free the structure.
	 * This function returns true if we parsed a GMT/OGR record and false otherwise.
	 * If we encounter a parsing error we stop parsing any further by setting GMT->current.io.ogr = GMT_OGR_FALSE.
	 * We loop until all @<info> tags have been processed on this record.
	 * gmt_ogr_parser will point to this function until the header has been parsed, then it is
	 * set to point to gmt_ogr_data_parser instead, to speed up data record processing.
	 */

	unsigned int n_aspatial, k, geometry = 0;
	bool quote;
	char *p = NULL;
	struct GMT_OGR *S = NULL;

	if (GMT->current.io.ogr == GMT_OGR_FALSE) return (false);	/* No point parsing further if we KNOW it is not OGR */
	if (record[0] != '#') return (false);			/* Not a comment record so no point looking any further */
	if (GMT->current.io.ogr == GMT_OGR_TRUE && !strncmp (record, "# FEATURE_DATA", 14)) {	/* It IS an OGR file and we found end of OGR header section and start of feature data */
		GMT->current.io.ogr_parser = &gmt_ogr_data_parser;	/* From now on only parse for feature tags */
		gmt_align_ogr_values (GMT);	/* Simplify copy from aspatial values to input columns as per -a option */
		return (true);
	}
	if (!(p = strchr (record, '@'))) return (false);	/* Not an OGR/GMT record since @ was not found */

	if (GMT->current.io.ogr == GMT_OGR_UNKNOWN && !strncmp (p, "@VGMT", 5)) {	/* Found the OGR version identifier, look for @G if on the same record */
		if (GMT->common.a.output) {	/* Cannot read OGR files when -a is used to define output */
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Cannot read OGR/GMT files when -a is used to define output format\n");
			GMT_exit (GMT, EXIT_FAILURE); return false;
		}
		GMT->current.io.ogr = GMT_OGR_TRUE;		/* File is now known to be a GMT/OGR geospatial file */
		if (!(p = strchr (&p[5], '@'))) return (true);	/* No more @ codes; goto next record */
	}
	if (GMT->current.io.ogr != GMT_OGR_TRUE) return (false);	/* No point parsing further since file is not GMT/OGR (at least not yet) */

	/* Here we are reasonably sure that @? strings are OGR/GMT header specifications */

	GMT_chop (record);	/* Get rid of linefeed etc */

	/* Allocate S the first time we get here */

	if (!GMT->current.io.OGR) GMT->current.io.OGR = GMT_memory (GMT, NULL, 1, struct GMT_OGR);
	S = GMT->current.io.OGR;
	quote = false;

	while (*p == '@') {
		++p;	/* Move to first char after @ */

		switch (p[0]) {	/* These are the header tags */

			case 'G':	/* Geometry */
				if (!strncmp (&p[1], "LINESTRING", 10))
					geometry = GMT_IS_LINESTRING;
				else if (p[1] == 'P') {
					if (!strncmp (&p[2], "OLYGON", 6))
						geometry = GMT_IS_POLYGON;
					else if (!strncmp (&p[2], "OINT", 4))
						geometry = GMT_IS_POINT;
				}
				else if (!strncmp (&p[1], "MULTI", 5)) {
					if (!strncmp (&p[6], "POINT", 5))
						geometry = GMT_IS_MULTIPOINT;
					else if (!strncmp (&p[6], "LINESTRING", 10))
						geometry = GMT_IS_MULTILINESTRING;
					else if (!strncmp (&p[6], "POLYGON", 7))
						geometry = GMT_IS_MULTIPOLYGON;
					else {
						//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Bad OGR/GMT: @G unrecognized geometry\n");
						GMT->current.io.ogr = GMT_OGR_FALSE;
						return (false);
					}
				}
				if (!S->geometry)
					S->geometry = geometry;
				else if (S->geometry != geometry) {
					//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Bad OGR/GMT: @G cannot have different geometries\n");
					GMT->current.io.ogr = GMT_OGR_FALSE;
				}
				break;

			case 'N':	/* Aspatial name fields, store in table header */
				if (!S->geometry)
					{ 
					//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Bad OGR/GMT: @N given but no geometry set\n"); 
					return (false);
					}
				if (S->name) {	/* Already set */
					//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Bad OGR/GMT: @N Cannot have more than one per segment\n");
					GMT->current.io.ogr = GMT_OGR_FALSE;
					return (false);
				}
				n_aspatial = gmt_ogr_decode_aspatial_names (GMT, &p[1], S);
				if (S->n_aspatial == 0)
					S->n_aspatial = n_aspatial;
				else if (S->n_aspatial != n_aspatial) {
					//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Bad OGR/GMT: @N number of items vary\n");
					GMT->current.io.ogr = GMT_OGR_FALSE;
				}
				break;

			case 'J':	/* Dataset projection strings (one of 4 kinds) */
				switch (p[1]) {
					case 'e': k = 0;	break;	/* EPSG code */
					case 'g': k = 1;	break;	/* GMT proj code */
					case 'p': k = 2;	break;	/* Proj.4 code */
					case 'w': k = 3;	break;	/* OGR WKT representation */
					default:
						//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Bad OGR/GMT: @J given unknown format (%c)\n", (int)p[1]);
						GMT->current.io.ogr = GMT_OGR_FALSE;
						return (false);
				}
				S->proj[k] = strdup (&p[2]);
				break;

			case 'R':	/* Dataset region */
				if (S->region) { /* Already set */
					//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Bad OGR/GMT: @R can only appear once\n");
					GMT->current.io.ogr = GMT_OGR_FALSE;
					return (false);
				}
				S->region = strdup (&p[1]);
				break;

			case 'T':	/* Aspatial field types, store in table header  */
				if (!S->geometry) { 
					//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Bad OGR/GMT: @T given but no geometry set\n"); 
					return (false);}
				if (S->type) {	/* Already set */
					//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Bad OGR/GMT: @T Cannot have more than one per segment\n");
					GMT->current.io.ogr = GMT_OGR_FALSE;
					return (false);
				}
				n_aspatial = gmt_ogr_decode_aspatial_types (GMT, &p[1], S);
				if (S->n_aspatial == 0)
					S->n_aspatial = n_aspatial;
				else if (S->n_aspatial != n_aspatial) {
					//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Bad OGR/GMT: @T number of items vary\n");
					GMT->current.io.ogr = GMT_OGR_FALSE;
				}
				break;

			default:	/* Just record, probably means this is NOT a GMT/OGR file after all */
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Bad OGR/GMT: @%c not allowed before FEATURE_DATA\n", (int)p[0]);
				GMT->current.io.ogr = GMT_OGR_FALSE;
				break;
		}

		while (*p && (quote || *p != '@')) if (*p++ == '\"') quote = !quote;	/* Wind to next @ except skip if inside double quotes */
	}
	return (true);
}

void * GMT_ascii_textinput (struct GMT_CTRL *GMT, FILE *fp, uint64_t *n, int *status)
{
	bool more = true;
	char line[GMT_BUFSIZ] = {""}, *p = NULL;

	/* GMT_ascii_textinput will read one text line and return it, setting
	 * header or segment flags in the process.
	 */

	while (more) {
		/* First read until we get a non-blank, non-comment record, or reach EOF */

		GMT->current.io.rec_no++;		/* Counts up, regardless of what this record is (data, junk, segment header, etc) */
		GMT->current.io.rec_in_tbl_no++;	/* Counts up, regardless of what this record is (data, junk, segment header, etc) */
		while ((p = GMT_fgets (GMT, line, GMT_BUFSIZ, fp)) && gmt_ogr_parser (GMT, line)) {	/* Exits loop when we successfully have read a data record */
			GMT->current.io.rec_no++;		/* Counts up, regardless of what this record is (data, junk, segment header, etc) */
			GMT->current.io.rec_in_tbl_no++;	/* Counts up, regardless of what this record is (data, junk, segment header, etc) */
		}
		/* Here we come once any OGR headers have been parsed and we have a real (non-OGR header) record */
		if (GMT->current.setting.io_header[GMT_IN] && GMT->current.io.rec_in_tbl_no <= GMT->current.setting.io_n_header_items) {	/* Must treat first io_n_header_items as headers */
			if (GMT->common.h.mode == GMT_COMMENT_IS_RESET) continue;	/* Simplest way to replace headers on output is to ignore them on input */
			strncpy (GMT->current.io.current_record, line, GMT_BUFSIZ);
			GMT->current.io.status = GMT_IO_TABLE_HEADER;
			*status = 0;
			return (NULL);
		}
		if (!p) {	/* Ran out of records */
			GMT->current.io.status = GMT_IO_EOF;
			*n = 0ULL;
			*status = -1;
			return (NULL);
		}
		if (line[0] == '#') {	/* Got a file header, take action and return */
			if (GMT->common.h.mode == GMT_COMMENT_IS_RESET) continue;	/* Simplest way to replace headers on output is to ignore them on input */
			strncpy (GMT->current.io.current_record, line, GMT_BUFSIZ);
			GMT->current.io.status = GMT_IO_TABLE_HEADER;
			*n = 1ULL;
			*status = 0;
			return (NULL);
		}

		if (line[0] == GMT->current.setting.io_seg_marker[GMT_IN]) {	/* Got a segment header, take action and return */
			GMT->current.io.status = GMT_IO_SEGMENT_HEADER;
			GMT_set_segmentheader (GMT, GMT_OUT, true);	/* Turn on segment headers on output */
			GMT->current.io.seg_no++;
			/* Just save the header content, not the marker and leading whitespace */
			strncpy (GMT->current.io.segment_header, GMT_trim_segheader (GMT, line), GMT_BUFSIZ);
			*n = 1ULL;
			*status = 0;
			return (NULL);
		}
		more = false;	/* Got a valid record */
	}

	/* Normal data record */

	/* First chop off trailing whitespace and commas */

	GMT_strstrip (line, false); /* Eliminate DOS endings and trailing white space */

	strncpy (GMT->current.io.current_record, line, GMT_BUFSIZ);

	GMT->current.io.status = 0;
	GMT->current.io.pt_no++;	/* Got a valid text record */
	*n = 1ULL;			/* We always return 1 item as there are no columns */
	*status = 1;
	return (GMT->current.io.current_record);
}

bool GMT_parse_segment_item (struct GMT_CTRL *GMT, char *in_string, char *pattern, char *out_string)
{
	/* Scans the in_string for the occurrence of an option switch (e.g, -L) and
	 * if found, extracts the argument and returns it via out_string.  Function
	 * return true if the pattern was found and false otherwise.
	 * out_string must be allocated and have space for the copying */
	char *t = NULL;
	size_t k;
	if (!in_string || !pattern) return (false);	/* No string or pattern passed */
	if (!(t = strstr (in_string, pattern))) return (false);	/* Option not present */
	if (!out_string) return (true);	/* If NULL is passed as out_string then we just return true if we find the option */
	out_string[0] = '\0';	/* Reset string to empty before we try to set it below */
	k = (size_t)t - (size_t)in_string; /* Position of pattern in in_string */
	if (k && !(in_string[k-1] == ' ' || in_string[k-1] == '\t')) return (false);	/* Option not first or preceeded by whitespace */
	t += 2;	/* Position of the argument */
	if (t[0] == '\"')	/* Double quoted argument, must scan from next character until terminal quote */
		sscanf (++t, "%[^\"]", out_string);
	else if (t[0] == '\'')	/* Single quoted argument, must scan from next character until terminal quote */
		sscanf (++t, "%[^\']", out_string);
	else	/* Scan until next white space; stop also when there is leading white space, indicating no argument at all! */
		sscanf (t, "%[^ \t]", out_string);
	return (true);
}

struct GMT_TEXTTABLE * GMT_read_texttable (struct GMT_CTRL *GMT, void *source, unsigned int source_type)
{
	/* Reads an entire segment text data set into memory */

	bool close_file = false, header = true, no_segments, first_seg = true;
	int status;
	size_t n_row_alloc = GMT_CHUNK, n_seg_alloc = GMT_CHUNK, n_head_alloc = GMT_TINY_CHUNK;
	uint64_t row = 0, n_read = 0, seg = 0, ncol = 0;
	char file[GMT_BUFSIZ] = {""}, *in = NULL;
	FILE *fp = NULL;
	struct GMT_TEXTTABLE *T = NULL;

	/* Determine input source */
	printf("source_type %s ",source_type);
	if (source_type == GMT_IS_FILE) {	/* source is a file name */
		strncpy (file, (const char *)source, GMT_BUFSIZ);
		if ((fp = GMT_fopen (GMT, file, "r")) == NULL) {
			////GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Cannot open file %s\n", file);
			return (NULL);
		}
		close_file = true;	/* We only close files we have opened here */
	}
	else if (source_type == GMT_IS_STREAM) {	/* Open file pointer given, just copy */
		fp = (FILE *)source;
		if (fp == NULL) fp = GMT->session.std[GMT_IN];	/* Default input */
		if (fp == GMT->session.std[GMT_IN])
			strcpy (file, "<stdin>");
		else
			strcpy (file, "<input stream>");
	}
	else if (source_type == GMT_IS_FDESC) {		/* Open file descriptor given, just convert to file pointer */
		int *fd = source;
		if (fd && (fp = fdopen (*fd, "r")) == NULL) {
			////GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Cannot convert file descriptor %d to stream in GMT_read_texttable\n", *fd);
			return (NULL);
		}
		if (fd == NULL) fp = GMT->session.std[GMT_IN];	/* Default input */
		if (fp == GMT->session.std[GMT_IN])
			strcpy (file, "<stdin>");
		else
			strcpy (file, "<input file descriptor>");
	}
	else {
		////GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Unrecognized source type %d in GMT_read_texttable\n", source_type);
		return (NULL);
	}

	in = GMT_ascii_textinput (GMT, fp, &ncol, &status);	/* Get first record */
	n_read++;
	if (GMT_REC_IS_EOF (GMT)) {
		//GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "File %s is empty!\n", file);
		return (NULL);
	}

	/* Allocate the Table structure */

	T = GMT_memory (GMT, NULL, 1, struct GMT_TEXTTABLE);
	T->file[GMT_IN] = strdup (file);
	T->segment = GMT_memory (GMT, NULL, n_seg_alloc, struct GMT_TEXTSEGMENT *);
	T->header  = GMT_memory (GMT, NULL, n_head_alloc, char *);

	while (status >= 0 && !GMT_REC_IS_EOF (GMT)) {	/* Not yet EOF */
		if (header) {
			while ((GMT->current.setting.io_header[GMT_IN] && n_read <= GMT->current.setting.io_n_header_items) || GMT_REC_IS_TABLE_HEADER (GMT)) { /* Process headers */
				T->header[T->n_headers] = strdup (GMT->current.io.current_record);
				T->n_headers++;
				if (T->n_headers == n_head_alloc) {
					n_head_alloc <<= 1;
					T->header = GMT_memory (GMT, T->header, n_head_alloc, char *);
				}
				in = GMT_ascii_textinput (GMT, fp, &ncol, &status);
				n_read++;
			}
			if (T->n_headers)
				T->header = GMT_memory (GMT, T->header, T->n_headers, char *);
			else {	/* No header records found */
				GMT_free (GMT, T->header);
				T->header = NULL;
			}
			header = false;	/* Done processing header block; other comments are GIS/OGR encoded comments */
		}

		no_segments = !GMT_REC_IS_SEGMENT_HEADER (GMT);	/* Not a multi-segment file.  We then assume file has only one segment */

		while (no_segments || (GMT_REC_IS_SEGMENT_HEADER (GMT) && !GMT_REC_IS_EOF (GMT))) {
			/* PW: This will need to change to allow OGR comments to follow segment header */
			/* To use different line-distances for each segment, place the distance in the segment header */
			if (first_seg || T->segment[seg]->n_rows > 0) {
				if (!first_seg) seg++;	/* Only advance segment if last had any points */
				T->segment[seg] = GMT_memory (GMT, NULL, 1, struct GMT_TEXTSEGMENT);
				first_seg = false;
			}
			n_read++;
			/* Segment initialization */
			n_row_alloc = GMT_CHUNK;
			row = 0;
			if (!no_segments) {
				in = GMT_ascii_textinput (GMT, fp, &ncol, &status);	/* Don't read if we didnt read a segment header up front */
				n_read++;
			}
			no_segments = false;	/* This has now served its purpose */
		}
		if (GMT_REC_IS_EOF (GMT)) continue;	/* At EOF; get out of this loop */
		if (!no_segments) {			/* Handle info stored in multi-seg header record */
			char buffer[GMT_BUFSIZ] = {""};
			if (GMT_parse_segment_item (GMT, GMT->current.io.segment_header, "-L", buffer)) T->segment[seg]->label = strdup (buffer);
			if (strlen (GMT->current.io.segment_header)) T->segment[seg]->header = strdup (GMT->current.io.segment_header);
		}

		T->segment[seg]->record = GMT_memory (GMT, NULL, n_row_alloc, char *);

		while (!(GMT->current.io.status & (GMT_IO_SEGMENT_HEADER | GMT_IO_EOF))) {	/* Keep going until false or find a new segment header */

			if (in) T->segment[seg]->record[row++] = strdup (in);	/* in might be NULL if comment record is found - these are skipped */

			if (row == n_row_alloc) {
				n_row_alloc <<= 1;
				T->segment[seg]->record = GMT_memory (GMT, T->segment[seg]->record, n_row_alloc, char *);
			}
			in = GMT_ascii_textinput (GMT, fp, &ncol, &status);
			n_read++;
		}
		T->segment[seg]->n_rows = row;	/* Number of records in this segment */
		T->n_records += row;		/* Total number of records so far */
		T->segment[seg]->id = seg;	/* Internal segment number */

		/* Reallocate to free up some memory */

		T->segment[seg]->record = GMT_memory (GMT, T->segment[seg]->record, T->segment[seg]->n_rows, char *);

		if (T->segment[seg]->n_rows == 0) {	/* Empty segment; we delete to avoid problems downstream in applications */
			GMT_free (GMT, T->segment[seg]);
			seg--;	/* Go back to where we were */
		}

		if (seg == (n_seg_alloc-1)) {
			size_t n_old_alloc = n_seg_alloc;
			n_seg_alloc <<= 1;
			T->segment = GMT_memory (GMT, T->segment, n_seg_alloc, struct GMT_TEXTSEGMENT *);
			GMT_memset (&(T->segment[n_old_alloc]), n_seg_alloc - n_old_alloc, struct GMT_TEXTSEGMENT *);	/* Set to NULL */
		}
	}
	if (close_file) GMT_fclose (GMT, fp);

	if (T->segment[seg]->n_rows == 0)	/* Last segment was empty; we delete to avoid problems downstream in applications */
		GMT_free (GMT, T->segment[seg]);
	else
		seg++;
	if (seg < n_seg_alloc) T->segment = GMT_memory (GMT, T->segment, seg, struct GMT_TEXTSEGMENT *);
	T->n_segments = seg;

	return (T);
}

int GMT_adjust_loose_wesn (struct GMT_CTRL *GMT, double wesn[], struct GMT_GRID_HEADER *header)
{
	/* Used to ensure that sloppy w,e,s,n values are rounded to the gridlines or pixels in the referenced grid.
	 * Upon entry, the boundaries w,e,s,n are given as a rough approximation of the actual subset needed.
	 * The routine will limit the boundaries to the grids region and round w,e,s,n to the nearest gridline or
	 * pixel boundaries (depending on the grid orientation).
	 * Warnings are produced when the w,e,s,n boundaries are adjusted, so this routine is currently not
	 * intended to throw just any values at it (although one could).
	 */

	bool global, error = false;
	double val, dx, small;
	char format[GMT_LEN64] = {""};

	switch (GMT_minmaxinc_verify (GMT, wesn[XLO], wesn[XHI], header->inc[GMT_X], GMT_SMALL)) {	/* Check if range is compatible with x_inc */
		case 3:
			return (GMT_GRDIO_BAD_XINC);
			break;
		case 2:
			return (GMT_GRDIO_BAD_XRANGE);
			break;
		default:
			/* Everything is seemingly OK */
			break;
	}
	switch (GMT_minmaxinc_verify (GMT, wesn[YLO], wesn[YHI], header->inc[GMT_Y], GMT_SMALL)) {	/* Check if range is compatible with y_inc */
		case 3:
			return (GMT_GRDIO_BAD_YINC);
			break;
		case 2:
			return (GMT_GRDIO_BAD_YRANGE);
			break;
		default:
			/* Everything is OK */
			break;
	}
	global = GMT_grd_is_global (GMT, header);

	if (!global) {
		if (GMT_x_is_lon (GMT, GMT_IN)) {
			/* If longitudes are all west of range or all east of range, try moving them by 360 degrees east or west */
			if (wesn[XHI] < header->wesn[XLO])
				wesn[XLO] += 360.0, wesn[XHI] += 360.0;
			else if (wesn[XLO] > header->wesn[XHI])
				wesn[XLO] -= 360.0, wesn[XHI] -= 360.0;
		}
		if (header->wesn[XLO] - wesn[XLO] > GMT_SMALL) wesn[XLO] = header->wesn[XLO], error = true;
		if (wesn[XHI] - header->wesn[XHI] > GMT_SMALL) wesn[XHI] = header->wesn[XHI], error = true;
	}
	if (header->wesn[YLO] - wesn[YLO] > GMT_SMALL) wesn[YLO] = header->wesn[YLO], error = true;
	if (wesn[YHI] - header->wesn[YHI] > GMT_SMALL) wesn[YHI] = header->wesn[YHI], error = true;
	//if (error)
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: Region exceeds grid domain. Region reduced to grid domain.\n");

	if (!(GMT_x_is_lon (GMT, GMT_IN) && GMT_360_RANGE (wesn[XLO], wesn[XHI]) && global)) {    /* Do this unless a 360 longitude wrap */
		small = GMT_SMALL * header->inc[GMT_X];

		val = header->wesn[XLO] + lrint ((wesn[XLO] - header->wesn[XLO]) * header->r_inc[GMT_X]) * header->inc[GMT_X];
		dx = fabs (wesn[XLO] - val);
		if (GMT_x_is_lon (GMT, GMT_IN)) dx = fmod (dx, 360.0);
		if (dx > small) {
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: (w - x_min) must equal (NX + eps) * x_inc), where NX is an integer and |eps| <= %g.\n", GMT_SMALL);
			sprintf (format, "Warning: w reset from %s to %s\n", GMT->current.setting.format_float_out, GMT->current.setting.format_float_out);
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, format, wesn[XLO], val);
			wesn[XLO] = val;
		}

		val = header->wesn[XLO] + lrint ((wesn[XHI] - header->wesn[XLO]) * header->r_inc[GMT_X]) * header->inc[GMT_X];
		dx = fabs (wesn[XHI] - val);
		if (GMT_x_is_lon (GMT, GMT_IN)) dx = fmod (dx, 360.0);
		if (dx > small) {
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: (e - x_min) must equal (NX + eps) * x_inc), where NX is an integer and |eps| <= %g.\n", GMT_SMALL);
			sprintf (format, "Warning: e reset from %s to %s\n", GMT->current.setting.format_float_out, GMT->current.setting.format_float_out);
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, format, wesn[XHI], val);
			wesn[XHI] = val;
		}
	}

	/* Check if s,n are a multiple of y_inc offset from y_min - if not adjust s, n */
	small = GMT_SMALL * header->inc[GMT_Y];

	val = header->wesn[YLO] + lrint ((wesn[YLO] - header->wesn[YLO]) * header->r_inc[GMT_Y]) * header->inc[GMT_Y];
	if (fabs (wesn[YLO] - val) > small) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: (s - y_min) must equal (NY + eps) * y_inc), where NY is an integer and |eps| <= %g.\n", GMT_SMALL);
		sprintf (format, "Warning: s reset from %s to %s\n", GMT->current.setting.format_float_out, GMT->current.setting.format_float_out);
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, format, wesn[YLO], val);
		wesn[YLO] = val;
	}

	val = header->wesn[YLO] + lrint ((wesn[YHI] - header->wesn[YLO]) * header->r_inc[GMT_Y]) * header->inc[GMT_Y];
	if (fabs (wesn[YHI] - val) > small) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: (n - y_min) must equal (NY + eps) * y_inc), where NY is an integer and |eps| <= %g.\n", GMT_SMALL);
		sprintf (format, "Warning: n reset from %s to %s\n", GMT->current.setting.format_float_out, GMT->current.setting.format_float_out);
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, format, wesn[YHI], val);
		wesn[YHI] = val;
	}
	return (GMT_NOERROR);
}

void gmt_copy_ogr_seg (struct GMT_CTRL *GMT, struct GMT_DATASEGMENT *S, struct GMT_OGR *G)
{	/* Allocates the OGR structure for a given segment and copies current values from table OGR segment */
	unsigned int col;

	gmt_alloc_ogr_seg (GMT, S, G->n_aspatial);
	for (col = 0; col < G->n_aspatial; col++) {
		if (G->tvalue != NULL && G->tvalue[col])
			S->ogr->tvalue[col] = strdup (G->tvalue[col]);
		if (G->dvalue != NULL && G->dvalue[col])
			S->ogr->dvalue[col] = G->dvalue[col];
	}
	S->ogr->pol_mode = G->pol_mode;
}

void GMT_prep_tmp_arrays (struct GMT_CTRL *GMT, size_t row, size_t n_cols)
{
	size_t col;

	/* Check if this is the very first time, if so we initialize the arrays */
	if (GMT->hidden.mem_cols == 0)
		gmt_init_tmp_arrays (GMT, n_cols);	/* First time we get here */

	/* Check if we are exceeding our column count so far, if so we must allocate more columns */
	else if (n_cols > GMT->hidden.mem_cols) {	/* Must allocate more columns, this is expected to happen rarely */
		GMT->hidden.mem_coord = GMT_memory (GMT, GMT->hidden.mem_coord, n_cols, double *);	/* New ones are NOT NULL */
		for (col = GMT->hidden.mem_cols; col < n_cols; col++)	/* Explicitly allocate the new additions */
			GMT->hidden.mem_coord[col] = GMT_memory (GMT, NULL, GMT->hidden.mem_rows, double);
		GMT->hidden.mem_cols = n_cols;		/* Updated column count */
	}

	/* Check if we are exceeding our allocated count for this column.  If so allocate more rows */

	if (row < GMT->hidden.mem_rows) return;	/* Nothing to do */

	/* Here we must allocate more rows, this is expected to happen rarely given the large initial allocation */

	while (row >= GMT->hidden.mem_rows) GMT->hidden.mem_rows <<= 1;	/* Double up until enough */
	for (col = 0; col < GMT->hidden.mem_cols; col++)	/* Add more memory via realloc */
		GMT->hidden.mem_coord[col] = GMT_memory (GMT, GMT->hidden.mem_coord[col], GMT->hidden.mem_rows, double);

	//GMT_Report (GMT->parent, GMT_MSG_DEBUG, "GMT memory: Increase %" PRIuS " temporary column arrays to new length : %" PRIuS "\n", GMT->hidden.mem_cols, GMT->hidden.mem_rows);
	/* Note: Any additions to these arrays are not guaranteed to be set to zero */
}

/* To avoid lots of alloc and realloc calls we prefer to allocate a sizeable array
 * per coordinate axes once, then use that temporary space for reading and
 * calculations, and then alloc permanent space elsewhere and call memcpy to
 * place the final memory there.  We assume that for most purposes we will
 * need GMT_INITIAL_MEM_COL_ALLOC columns [2] and allocate GMT_INITIAL_MEM_ROW_ALLOC
 * [2097152U] rows for each column.  This is 32 Mb for double precision data.
 * These arrays are expected to hardly ever beeing reallocated as that would
 * only happen for very long segments, a rare occurance. For most typical data
 * we may have lots of smaller segments but rarely do any segment exceed the
 * 1048576U length initialized above.  Thus, reallocs are generally avoided.
 * Note: (1) All columns share a signle n_alloc counter and the code belows will
 *           check whenever arrays need to be extended.
 *	 (2) We chose to maintain a small set of column vectors rather than a single
 *	     item since GMT tends to use columns vectors and thus the book-keeping is
 *	     simpler and the number of columns is typically very small (2-3).
 */

void gmt_init_tmp_arrays (struct GMT_CTRL *GMT, size_t n_cols)
{
	/* Initialization of GMT coordinate temp arrays - this is called at most once per GMT session  */

	size_t col;

	if (n_cols == 0) n_cols = GMT_INITIAL_MEM_COL_ALLOC;	/* Allocate at least this many */
	GMT->hidden.mem_coord  = GMT_memory (GMT, GMT->hidden.mem_coord, n_cols, double *);	/* These are all NULL */
	GMT->hidden.mem_cols = n_cols;	/* How many columns we have initialized */
	for (col = 0; col < n_cols; col++)	/* For each column, reallocate space for n_rows */
		GMT->hidden.mem_coord[col] = GMT_memory (GMT, NULL, GMT_INITIAL_MEM_ROW_ALLOC, double);
	GMT->hidden.mem_rows = GMT_INITIAL_MEM_ROW_ALLOC;
	//GMT_Report (GMT->parent, GMT_MSG_DEBUG, "GMT memory: Initialize %" PRIuS " temporary column arrays, each of length : %" PRIuS "\n", GMT->hidden.mem_cols, GMT->hidden.mem_rows);
}

void GMT_set_seg_polar (struct GMT_CTRL *GMT, struct GMT_DATASEGMENT *S)
{	/* Must check if polygon is a polar cap.  We use the idea that the sum of
 	 * the change in angle along the polygon is either -/+360 (if pole is inside)
	 * or 0 if outside.  The sign (not used here) gives the handedness of the polygon.
	 * Which pole (S or N) is determined by computng the average latitude and
	 * assuming the pole is in the heimsphere most visited.  This may not be
	 * true of course. */
	uint64_t row;
	int n_360;
	double dlon, lon_sum = 0.0, lat_sum = 0.0;
	static char *pole[3] = {"south", "no", "north"};

	if (GMT_polygon_is_open (GMT, S->coord[GMT_X], S->coord[GMT_Y], S->n_rows)) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Error: Cannot call GMT_set_seg_polar on an open polygon\n");
		return;
	}
	for (row = 0; row < S->n_rows - 1; row++) {
		GMT_set_delta_lon (S->coord[GMT_X][row], S->coord[GMT_X][row+1], dlon);
		lon_sum += dlon;
		lat_sum += S->coord[GMT_Y][row];
	}
	n_360 = irint (lon_sum / 360.0);	/* This is either -1, 0, or +1 since lon_sum is either -360, 0, +360 plus some noise */
	if (n_360) {	/* true if contains a pole; adjust rectangular bounds and set pole flag */
		S->pole = irint (copysign (1.0, lat_sum));	/* So, 0 means not polar */
		S->min[GMT_X] = 0.0;	S->max[GMT_X] = 360.0;
		if (S->pole == -1) S->lat_limit = S->min[GMT_Y], S->min[GMT_Y] = -90.0;
		if (S->pole == +1) S->lat_limit = S->max[GMT_Y], S->max[GMT_Y] = +90.0;
	}
	else
		S->pole = 0;
	//GMT_Report (GMT->parent, GMT_MSG_DEBUG, "GMT_set_seg_polar: N = %" PRIu64 " Multiples of 360: %d  Residual: %g Polygon contains %s pole.\n", S->n_rows, n_360, lon_sum - n_360 * 360.0, pole[S->pole+1]);
}

struct GMT_DATATABLE * GMT_create_table (struct GMT_CTRL *GMT, uint64_t n_segments, uint64_t n_rows, uint64_t n_columns, bool alloc_only)
{
	/* Allocate the new Table structure given the specified dimensions.
	 * If n_columns == 0 it means we don't know that dimension yet.
	 * If alloc_only is true then we do NOT set the corresponding counters (i.e., n_segments).  */
	uint64_t seg;
	struct GMT_DATATABLE *T = NULL;

	T = GMT_memory (GMT, NULL, 1, struct GMT_DATATABLE);
	if (!alloc_only) T->n_segments = n_segments;
	if (!alloc_only) T->n_records = n_segments * n_rows;
	T->n_alloc = n_segments;
	if (n_columns) {
		T->min = GMT_memory (GMT, NULL, n_columns, double);
		T->max = GMT_memory (GMT, NULL, n_columns, double);
	}
	T->n_columns = n_columns;
	if (n_segments) {
		T->segment = GMT_memory (GMT, NULL, n_segments, struct GMT_DATASEGMENT *);
		for (seg = 0; n_columns && seg < n_segments; seg++) {
			T->segment[seg] = GMT_memory (GMT, NULL, 1, struct GMT_DATASEGMENT);
			if (GMT_alloc_segment (GMT, T->segment[seg], n_rows, n_columns, true)) return (NULL);
		}
	}

	return (T);
}

struct GMT_DATATABLE * GMT_read_table (struct GMT_CTRL *GMT, void *source, unsigned int source_type, bool greenwich, unsigned int *geometry, bool use_GMT_io)
{
	/* Reads an entire data set into a single table in memory with any number of segments */

	bool ascii, close_file = false, header = true, no_segments, first_seg = true, poly, check_geometry;
	int status;
	uint64_t n_expected_fields;
	uint64_t n_read = 0, row = 0, seg = 0, col;
	size_t n_head_alloc = GMT_TINY_CHUNK;
	char open_mode[4] = {""}, file[GMT_BUFSIZ] = {""}, line[GMT_LEN64] = {""};
	double d, *in = NULL;
	FILE *fp = NULL;
	struct GMT_DATATABLE *T = NULL;
	void * (*psave) (struct GMT_CTRL *, FILE *, uint64_t *, int *) = NULL;	/* Pointer to function reading tables */

	if (use_GMT_io) {	/* Use GMT->current.io.info settings to determine if input is ascii/binary, else it defaults to ascii */
		n_expected_fields = GMT->common.b.active[GMT_IN] ? GMT->common.b.ncol[GMT_IN] : GMT_MAX_COLUMNS;
		strcpy (open_mode, GMT->current.io.r_mode);
		ascii = !GMT->common.b.active[GMT_IN];
	}
	else {			/* Force ASCII mode */
		n_expected_fields = GMT_MAX_COLUMNS;	/* GMT->current.io.input will return the number of columns */
		strcpy (open_mode, "r");
		ascii = true;
		psave = GMT->current.io.input;			/* Save the previous pointer since we need to change it back at the end */
		GMT->current.io.input = GMT->session.input_ascii;	/* Override and use ascii mode */
	}

#ifdef SET_IO_MODE
	if (!ascii) GMT_setmode (GMT, GMT_IN);
#endif

	check_geometry = ((*geometry & GMT_IS_POLY) && (*geometry & GMT_IS_LINE));	/* Have to determine if these are closed polygons or not */
	poly = (((*geometry & GMT_IS_POLY) || *geometry == GMT_IS_MULTIPOLYGON) && (*geometry & GMT_IS_LINE) == 0);	/* To enable polar cap assessment in i/o */

	/* Determine input source */

	if (source_type == GMT_IS_FILE) {	/* source is a file name */
		strncpy (file, (const char *)source, GMT_BUFSIZ);
		if ((fp = GMT_fopen (GMT, file, open_mode)) == NULL) {
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Cannot open file %s\n", file);
			if (!use_GMT_io) GMT->current.io.input = psave;	/* Restore previous setting */
			return (NULL);
		}
		close_file = true;	/* We only close files we have opened here */
	}
	else if (source_type == GMT_IS_STREAM) {	/* Open file pointer given, just copy */
		fp = (FILE *)source;
		if (fp == NULL) fp = GMT->session.std[GMT_IN];	/* Default input */
		if (fp == GMT->session.std[GMT_IN])
			strcpy (file, "<stdin>");
		else
			strcpy (file, "<input stream>");
	}
	else if (source_type == GMT_IS_FDESC) {		/* Open file descriptor given, just convert to file pointer */
		int *fd = source;
		if (fd && (fp = fdopen (*fd, open_mode)) == NULL) {
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Cannot convert file descriptor %d to stream in GMT_read_table\n", *fd);
			if (!use_GMT_io) GMT->current.io.input = psave;	/* Restore previous setting */
			return (NULL);
		}
		if (fd == NULL) fp = GMT->session.std[GMT_IN];	/* Default input */
		if (fp == GMT->session.std[GMT_IN])
			strcpy (file, "<stdin>");
		else
			strcpy (file, "<input file descriptor>");
	}
	else {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Unrecognized source type %d in GMT_read_table\n", source_type);
		if (!use_GMT_io) GMT->current.io.input = psave;	/* Restore previous setting */
		return (NULL);
	}

	in = GMT->current.io.input (GMT, fp, &n_expected_fields, &status);	/* Get first record */
	n_read++;
	if (GMT_REC_IS_EOF(GMT)) {
		//GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "File %s is empty!\n", file);
		if (!use_GMT_io) GMT->current.io.input = psave;	/* Restore previous setting */
		return (NULL);
	}
	/* Allocate the Table structure with GMT_CHUNK segments, but none has any rows or columns */

	T = GMT_create_table (GMT, GMT_CHUNK, 0U, 0U, false);

	T->file[GMT_IN] = strdup (file);
	T->header = GMT_memory (GMT, NULL, n_head_alloc, char *);

	while (status >= 0 && !GMT_REC_IS_EOF (GMT)) {	/* Not yet EOF */
		if (header) {
			while ((GMT->current.setting.io_header[GMT_IN] && n_read <= GMT->current.setting.io_n_header_items) || GMT_REC_IS_TABLE_HEADER (GMT)) { /* Process headers */
				T->header[T->n_headers] = strdup (GMT->current.io.current_record);
				T->n_headers++;
				if (T->n_headers == n_head_alloc) {
					n_head_alloc <<= 1;
					T->header = GMT_memory (GMT, T->header, n_head_alloc, char *);
				}
				in = GMT->current.io.input (GMT, fp, &n_expected_fields, &status);
				n_read++;
			}
			if (T->n_headers)
				T->header = GMT_memory (GMT, T->header, T->n_headers, char *);
			else {	/* No header records found */
				GMT_free (GMT, T->header);
				T->header = NULL;
			}
			header = false;	/* Done processing header block; other comments are GIS/OGR encoded comments */
		}

		if (GMT_REC_IS_EOF (GMT)) continue;	/* Got EOF after headers */

		no_segments = !GMT_REC_IS_SEGMENT_HEADER (GMT);	/* Not a multi-segment file.  We then assume file has only one segment */

		while (no_segments || (GMT_REC_IS_SEGMENT_HEADER (GMT) && !GMT_REC_IS_EOF (GMT))) {
			/* To use different line-distances for each segment, place the distance in the segment header */
			if (first_seg || T->segment[seg]->n_rows > 0) {
				if (!first_seg) seg++;	/* Only advance segment if last had any points */
				T->segment[seg] = GMT_memory (GMT, NULL, 1, struct GMT_DATASEGMENT);
				first_seg = false;
			}
			n_read++;
			if (ascii && !no_segments) {	/* Only ascii files can have info stored in multi-seg header records */
				if (GMT_parse_segment_item (GMT, GMT->current.io.segment_header, "-D", line)) {	/* Found a potential -D<dist> option in the header */
					if (sscanf (line, "%lg", &d) == 1) T->segment[seg]->dist = d;	/* If readable, assign it to dist, else leave as zero */
				}
			}
			/* Segment initialization */
			row = 0;
			if (!no_segments) {	/* Read data if we read a segment header up front, but guard against headers which sets in = NULL */
				while (!GMT_REC_IS_EOF (GMT) && (in = GMT->current.io.input (GMT, fp, &n_expected_fields, &status)) == NULL) n_read++;
			}
			T->segment[seg]->n_columns = n_expected_fields;	/* This is where number of columns are determined */
			no_segments = false;	/* This has now served its purpose */
		}
		if (GMT_REC_IS_EOF (GMT)) continue;	/* At EOF; get out of this loop */
		if (ascii && !no_segments) {	/* Only ascii files can have info stored in multi-seg header record */
			char buffer[GMT_BUFSIZ] = {""};
			if (strlen (GMT->current.io.segment_header)) {
				T->segment[seg]->header = strdup (GMT->current.io.segment_header);
				if (GMT_parse_segment_item (GMT, GMT->current.io.segment_header, "-L", buffer)) T->segment[seg]->label = strdup (buffer);
			}
			if (GMT->current.io.ogr == GMT_OGR_TRUE) gmt_copy_ogr_seg (GMT, T->segment[seg], GMT->current.io.OGR);	/* Copy over any feature-specific values */
		}

		if (poly && T->segment[seg]->n_columns < 2) {
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "File %s does not have at least 2 columns required for polygons (found %d)\n", file, T->segment[seg]->n_columns);
			if (!use_GMT_io) GMT->current.io.input = psave;	/* Restore previous setting */
			return (NULL);
		}

		while (! (GMT->current.io.status & (GMT_IO_SEGMENT_HEADER | GMT_IO_GAP | GMT_IO_EOF))) {	/* Keep going until false or find a new segment header */
			if (GMT->current.io.status & GMT_IO_MISMATCH) {
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Mismatch between actual (%d) and expected (%d) fields near line %" PRIu64 "\n", status, n_expected_fields, n_read);
				if (!use_GMT_io) GMT->current.io.input = psave;	/* Restore previous setting */
				return (NULL);
			}

			GMT_prep_tmp_arrays (GMT, row, T->segment[seg]->n_columns);	/* Init or reallocate tmp read vectors */
			for (col = 0; col < T->segment[seg]->n_columns; col++) {
				GMT->hidden.mem_coord[col][row] = in[col];
				if (GMT->current.io.col_type[GMT_IN][col] & GMT_IS_LON) {	/* Must handle greenwich/dateline alignments */
					if (greenwich && GMT->hidden.mem_coord[col][row] > 180.0) GMT->hidden.mem_coord[col][row] -= 360.0;
					if (!greenwich && GMT->hidden.mem_coord[col][row] < 0.0)  GMT->hidden.mem_coord[col][row] += 360.0;
				}
			}

			row++;
			in = GMT->current.io.input (GMT, fp, &n_expected_fields, &status);
			while (GMT_REC_IS_TABLE_HEADER (GMT)) in = GMT->current.io.input (GMT, fp, &n_expected_fields, &status);	/* Just wind past other comments */
			n_read++;
		}

		if (check_geometry) {	/* Determine if dealing with closed polygons or lines based on first segment only */
			if (!GMT_polygon_is_open (GMT, GMT->hidden.mem_coord[GMT_X], GMT->hidden.mem_coord[GMT_Y], row)) poly = true;
			check_geometry = false;	/* Done with one-time checking */
			*geometry = (poly) ? GMT_IS_POLY : GMT_IS_LINE;	/* Update the geometry setting */
		}
		if (poly) {	/* If file contains a polygon then we must close it if needed */
			if (GMT->current.io.col_type[GMT_IN][GMT_X] & GMT_IS_GEO) {	/* Must check for polar cap */
				double dlon = GMT->hidden.mem_coord[GMT_X][0] - GMT->hidden.mem_coord[GMT_X][row-1];
				if (!((fabs (dlon) == 0.0 || fabs (dlon) == 360.0) && GMT->hidden.mem_coord[GMT_Y][0] == GMT->hidden.mem_coord[GMT_Y][row-1])) {
					GMT_prep_tmp_arrays (GMT, row, T->segment[seg]->n_columns);	/* Maybe reallocate tmp read vectors */
					for (col = 0; col < T->segment[seg]->n_columns; col++) GMT->hidden.mem_coord[col][row] = GMT->hidden.mem_coord[col][0];
					row++;	/* Explicitly close polygon */
				}
			}
			else if (GMT_polygon_is_open (GMT, GMT->hidden.mem_coord[GMT_X], GMT->hidden.mem_coord[GMT_Y], row)) {	/* Cartesian closure */
				GMT_prep_tmp_arrays (GMT, row, T->segment[seg]->n_columns);	/* Init or update tmp read vectors */
				for (col = 0; col < T->segment[seg]->n_columns; col++) GMT->hidden.mem_coord[col][row] = GMT->hidden.mem_coord[col][0];
				row++;	/* Explicitly close polygon */
			}
			if (GMT_parse_segment_item (GMT, T->segment[seg]->header, "-Ph", NULL)) T->segment[seg]->pol_mode = GMT_IS_HOLE;
			/* If this is a hole then set link from previous segment to this one */
			if (seg && GMT_polygon_is_hole (T->segment[seg])) T->segment[seg-1]->next = T->segment[seg];
		}

		if (row == 0) {	/* Empty segment; we delete to avoid problems downstream in applications */
			GMT_free (GMT, T->segment[seg]);
			seg--;	/* Go back to where we were */
		}
		else {	/* OK to populate segment and increment counters */
			GMT_assign_segment (GMT, T->segment[seg], row, T->segment[seg]->n_columns);	/* Allocate and place arrays into segment */
			GMT_set_seg_minmax (GMT, T->segment[seg]);	/* Set min/max */
			if (poly && (GMT->current.io.col_type[GMT_IN][GMT_X] & GMT_IS_GEO)) GMT_set_seg_polar (GMT, T->segment[seg]);
			T->n_records += row;		/* Total number of records so far */
			T->segment[seg]->id = seg;	/* Internal segment number */
		}
		/* Reallocate to free up some memory */

		if (seg == (T->n_alloc-1)) {	/* Need to allocate more segments */
			size_t n_old_alloc = T->n_alloc;
			T->n_alloc <<= 1;
			T->segment = GMT_memory (GMT, T->segment, T->n_alloc, struct GMT_DATASEGMENT *);
			GMT_memset (&(T->segment[n_old_alloc]), T->n_alloc - n_old_alloc, struct GMT_DATASEGMENT *);	/* Set to NULL */
		}

		/* If a gap was detected, forget about it now, so we can use the data for the next segment */

		GMT->current.io.status -= (GMT->current.io.status & GMT_IO_GAP);
	}
	if (close_file) GMT_fclose (GMT, fp);
	if (!use_GMT_io) GMT->current.io.input = psave;	/* Restore previous setting */

	if (first_seg) {	/* Never saw any segment or data records */
		GMT_free_table (GMT, T, GMT_ALLOCATED_BY_GMT);
		return (NULL);
	}
	if (T->segment[seg]->n_rows == 0) {	/* Last segment was empty; we delete to avoid problems downstream in applications */
		GMT_free (GMT, T->segment[seg]);
		if (seg == 0) {	/* Happens when we just read 1 segment header with no data */
			GMT_free_table (GMT, T, GMT_ALLOCATED_BY_GMT);
			return (NULL);
		}
	}
	else
		seg++;
	T->segment = GMT_memory (GMT, T->segment, seg, struct GMT_DATASEGMENT *);
	T->n_segments = seg;
	T->n_columns = T->segment[0]->n_columns;
	/* Determine table min,max values */
	T->min = GMT_memory (GMT, NULL, T->n_columns, double);
	T->max = GMT_memory (GMT, NULL, T->n_columns, double);
	for (col = 0; col < T->n_columns; col++) {T->min[col] = DBL_MAX; T->max[col] = -DBL_MAX;}
	for (seg = 0; seg < T->n_segments; seg++) {
		for (col = 0; col < T->n_columns; col++) {
			T->min[col] = MIN (T->min[col], T->segment[seg]->min[col]);
			T->max[col] = MAX (T->max[col], T->segment[seg]->max[col]);
		}
		if (T->segment[seg]->pole) {T->min[GMT_X] = 0.0; T->max[GMT_X] = 360.0;}
	}

	return (T);
}
int GMT_scanf_arg (struct GMT_CTRL *GMT, char *s, unsigned int expectation, double *val)
{
	/* Version of GMT_scanf used for cpt & command line arguments only (not data records).
	 * It differs from GMT_scanf in that if the expectation is GMT_IS_UNKNOWN it will
	 * check to see if the argument is (1) an absolute time string, (2) a geographical
	 * location string, or if not (3) a floating point string.  To ensure backward
	 * compatibility: if we encounter geographic data it will also set the GMT->current.io.type[]
	 * variable accordingly so that data i/o will work as in 3.4
	 */

	char c;

	if (expectation == GMT_IS_UNKNOWN) {		/* Expectation for this column not set - must be determined if possible */
		c = s[strlen(s)-1];
		if (strchr (s, (int)'T'))		/* Found a T in the argument - assume Absolute time */
			expectation = GMT_IS_ARGTIME;
		else if (c == 't')			/* Found trailing t - assume Relative time */
			expectation = GMT_IS_ARGTIME;
		else if (strchr ("WwEe", (int)c))	/* Found trailing W or E - assume Geographic longitudes */
			expectation = GMT_IS_LON;
		else if (strchr ("SsNn", (int)c))	/* Found trailing S or N - assume Geographic latitudes */
			expectation = GMT_IS_LAT;
		else if (strchr ("DdGg", (int)c))	/* Found trailing G or D - assume Geographic coordinate */
			expectation = GMT_IS_GEO;
		else if (strchr (s, (int)':'))		/* Found a : in the argument - assume Geographic coordinates */
			expectation = GMT_IS_GEO;
		else 					/* Found nothing - assume floating point */
			expectation = GMT_IS_FLOAT;
	}

	/* OK, here we have an expectation, now call GMT_scanf */

	return (GMT_scanf (GMT, s, expectation, val));
}

bool GMT_is_a_blank_line (char *line) {
	/* Returns true if we should skip this line (because it is blank) */
	unsigned int i = 0;
	while (line[i] && (line[i] == ' ' || line[i] == '\t')) i++;	/* Wind past leading whitespace or tabs */
	if (line[i] == '\n' || line[i] == '\r' || line[i] == '\0') return (true);
	return (false);
}

bool GMT_is_a_NaN_line (struct GMT_CTRL *GMT, char *line)
{	/* Returns true if record is NaN NaN [NaN NaN] etc */
	unsigned int pos = 0;
	char p[GMT_LEN256] = {""};

	while ((GMT_strtok (line, " \t,", &pos, p))) {
		GMT_str_tolower (p);
		if (strncmp (p, "nan", 3U)) return (false);
	}
	return (true);
}

unsigned int gmt_is_segment_header (struct GMT_CTRL *GMT, char *line)
{	/* Returns 1 if this record is a GMT segment header;
	 * Returns 2 if this record is a segment breaker;
	 * Otherwise returns 0 */
	if (GMT->current.setting.io_blankline[GMT_IN] && GMT_is_a_blank_line (line)) return (2);	/* Treat blank line as segment break */
	if (GMT->current.setting.io_nanline[GMT_IN] && GMT_is_a_NaN_line (GMT, line)) return (2);	/* Treat NaN-records as segment break */
	if (line[0] == GMT->current.setting.io_seg_marker[GMT_IN]) return (1);	/* Got a regular GMT segment header */
	return (0);	/* Not a segment header */
}

unsigned int GMT_verify_expectations (struct GMT_CTRL *GMT, unsigned int wanted, unsigned int got, char *item)
{	/* Compare what we wanted with what we got and see if it is OK */
	unsigned int error = 0;

	if (wanted == GMT_IS_UNKNOWN) {	/* No expectations set */
		switch (got) {
			case GMT_IS_ABSTIME:	/* Found a T in the string - ABSTIME ? */
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Error: %s appears to be an Absolute Time String: ", item);
				//if (GMT_is_geographic (GMT, GMT_IN))
					//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "This is not allowed for a map projection\n");
				//else
					//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "You must specify time data type with option -f.\n");
				error++;
				break;

			case GMT_IS_GEO:	/* Found a : in the string - GEO ? */
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: %s appears to be a Geographical Location String: ", item);
				//if (GMT->current.proj.projection == GMT_LINEAR)
					//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "You should append d to the -Jx or -JX projection for geographical data.\n");
				//else
					//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "You should specify geographical data type with option -f.\n");
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Will proceed assuming geographical input data.\n");
				break;

			case GMT_IS_LON:	/* Found a : in the string and then W or E - LON ? */
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: %s appears to be a Geographical Longitude String: ", item);
				//if (GMT->current.proj.projection == GMT_LINEAR)
					//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "You should append d to the -Jx or -JX projection for geographical data.\n");
				//else
					//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "You should specify geographical data type with option -f.\n");
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Will proceed assuming geographical input data.\n");
				printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
				break;

			case GMT_IS_LAT:	/* Found a : in the string and then S or N - LAT ? */
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: %s appears to be a Geographical Latitude String: ", item);
				//if (GMT->current.proj.projection == GMT_LINEAR)
					//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "You should append d to the -Jx or -JX projection for geographical data.\n");
				//else
					//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "You should specify geographical data type with option -f.\n");
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Will proceed assuming geographical input data.\n");
				printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
				break;

			case GMT_IS_FLOAT:
				break;
			default:
				break;
		}
	}
	else {
		switch (got) {
			case GMT_IS_NAN:
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Error: Could not decode %s, return NaN.\n", item);
				printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
				error++;
				break;

			case GMT_IS_LAT:
				if (wanted == GMT_IS_LON) {
					//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Error: Expected longitude, but %s is a latitude!\n", item);
					//printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
					error++;
				}
				break;

			case GMT_IS_LON:
				if (wanted == GMT_IS_LAT) {
					//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Error: Expected latitude, but %s is a longitude!\n", item);
					//printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
					error++;
				}
				break;
			default:
				break;
		}
	}
	return (error);
}

struct GMT_OGR * GMT_duplicate_ogr (struct GMT_CTRL *GMT, struct GMT_OGR *G)
{	/* Duplicate GMT/OGR structure, if used */
	unsigned int k;
	struct GMT_OGR *G_dup = NULL;
	if (!G) return (NULL);	/* Nothing to do */
	G_dup = GMT_memory (GMT, NULL, 1, struct GMT_OGR);
	if (G->region) G_dup->region = strdup (G->region);
	for (k = 0; k < 4; k++) if (G->proj[k]) G_dup->proj[k] = strdup (G->proj[k]);
	G_dup->geometry = G->geometry;
	if (G->n_aspatial) {
		G_dup->n_aspatial = G->n_aspatial;
		G_dup->name = GMT_memory (GMT, NULL, G->n_aspatial, char *);
		for (k = 0; k < G->n_aspatial; k++) if (G->name[k]) G_dup->name[k] = strdup (G->name[k]);
		G_dup->type = GMT_memory (GMT, NULL, G->n_aspatial, unsigned int);
		GMT_memcpy (G_dup->type, G->type, G->n_aspatial, int);
	}
	return (G_dup);
}

int GMT_nc_get_att_text (struct GMT_CTRL *GMT, int ncid, int varid, char *name, char *text, size_t textlen)
{	/* This function is a replacement for nc_get_att_text that avoids overflow of text
	 * ncid, varid, name, text	: as in nc_get_att_text
	 * textlen			: maximum number of characters to copy to string text
	 */
	int status;
	size_t attlen;
	char *att = NULL;

	status = nc_inq_attlen (ncid, varid, name, &attlen);
	if (status != NC_NOERR) {
		*text = '\0';
		return status;
	}
	att = GMT_memory (GMT, NULL, attlen, char);
	status = nc_get_att_text (ncid, varid, name, att);
	if (status == NC_NOERR) {
		attlen = MIN (attlen, textlen-1); /* attlen does not include terminating '\0') */
		strncpy (text, att, attlen); /* Copy att to text */
		text[attlen] = '\0'; /* Terminate string */
	}
	else
		*text = '\0';
	GMT_free (GMT, att);
	return status;
}
FILE *gmt_nc_fopen (struct GMT_CTRL *GMT, const char *filename, const char *mode)
/* Open a netCDF file for column I/O. Append ?var1/var2/... to indicate the requested columns.
 * Currently only reading is supported.
 * The routine returns a fake file pointer (in fact the netCDF file ID), but stores
 * all the relevant information in the GMT->current.io struct (ncid, ndim, nrec, varid, add_offset,
 * scale_factor, missing_value). Some of these are allocated here, and have to be
 * deallocated upon GMT_fclose.
 * Also asigns GMT->current.io.col_type[GMT_IN] based on the variable attributes.
 */
{
	char file[GMT_BUFSIZ] = {""}, path[GMT_BUFSIZ] = {""};
	int i, j, nvars, dimids[5] = {-1, -1, -1, -1, -1}, ndims, in, id;
	size_t n, item[2];
	size_t tmp_pointer; /* To avoid "cast from pointer to integer of different size" */
	double t_value[5], dummy[2];
	char varnm[20][GMT_LEN64], long_name[GMT_LEN256] = {""}, units[GMT_LEN256] = {""};
	char varname[GMT_LEN64] = {""}, dimname[GMT_LEN64] = {""};
	struct GMT_TIME_SYSTEM time_system;
	bool by_value;

	if (mode[0] != 'r') {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "GMT_fopen does not support netCDF writing mode\n");
		printf("GMT_fopen does not support netCDF writing mode\n");
		GMT_exit (GMT, EXIT_FAILURE); return NULL;
	}

	GMT_memset (varnm, 20 * GMT_LEN64, char);

	nvars = sscanf (filename,
		"%[^?]?%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]",
		file, varnm[0], varnm[1], varnm[2], varnm[3], varnm[4], varnm[5], varnm[6], varnm[7], varnm[8], varnm[9], varnm[10],
		varnm[11], varnm[12], varnm[13], varnm[14], varnm[15], varnm[16], varnm[17], varnm[18], varnm[19]) - 1;
	//char *a = GMT_getdatapath (GMT, file, path, R_OK);
	//a[10] = '\0';
	//fflush(stdout);
	//return;
	if (nc_open (GMT_getdatapath (GMT, file, path, R_OK), NC_NOWRITE, &GMT->current.io.ncid)) return (NULL);
	if (GMT_compat_check (GMT, 4)) {
		if (nvars <= 0) nvars = sscanf (GMT->common.b.varnames,
			"%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]",
			varnm[0], varnm[1], varnm[2], varnm[3], varnm[4], varnm[5], varnm[6], varnm[7], varnm[8], varnm[9], varnm[10],
			varnm[11], varnm[12], varnm[13], varnm[14], varnm[15], varnm[16], varnm[17], varnm[18], varnm[19]);
	}
	if (nvars <= 0)
		nc_inq_nvars (GMT->current.io.ncid, &GMT->current.io.nvars);
	else
		GMT->current.io.nvars = nvars;
	GMT->current.io.varid = GMT_memory (GMT, NULL, GMT->current.io.nvars, int);
	GMT->current.io.scale_factor = GMT_memory (GMT, NULL, GMT->current.io.nvars, double);
	GMT->current.io.add_offset = GMT_memory (GMT, NULL, GMT->current.io.nvars, double);
	GMT->current.io.missing_value = GMT_memory (GMT, NULL, GMT->current.io.nvars, double);
	GMT->current.io.ndim = GMT->current.io.nrec = 0;

	for (i = 0; i < GMT->current.io.nvars; i++) {

		/* Check for indices */
		for (j = 0; j < 5; j++) GMT->current.io.t_index[i][j] = 0, GMT->current.io.count[i][j] = 1;
		j = in = 0, by_value = false;
		while (varnm[i][j] && varnm[i][j] != '(' && varnm[i][j] != '[') j++;
		if (varnm[i][j] == '(') {
			in = sscanf (&varnm[i][j+1], "%lf,%lf,%lf,%lf", &t_value[1], &t_value[2], &t_value[3], &t_value[4]);
			varnm[i][j] = '\0';
			by_value = true;
		}
		else if (varnm[i][j] == '[') {
			in = sscanf (&varnm[i][j+1], "%u,%u,%u,%u", &GMT->current.io.t_index[i][1], &GMT->current.io.t_index[i][2], &GMT->current.io.t_index[i][3], &GMT->current.io.t_index[i][4]);
			varnm[i][j] = '\0';
		}

		/* Get variable ID and variable name */
		if (nvars <= 0)
			GMT->current.io.varid[i] = i;
		else
		{
			//GMT_err_fail (GMT, nc_inq_varid (GMT->current.io.ncid, varnm[i], &GMT->current.io.varid[i]), file);
			printf("\n\n**Error**\n\n");
			exit(0);
		}
		nc_inq_varname (GMT->current.io.ncid, GMT->current.io.varid[i], varname);

		/* Check number of dimensions */
		nc_inq_varndims (GMT->current.io.ncid, GMT->current.io.varid[i], &ndims);
		if (ndims > 5) {
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "NetCDF variable %s has too many dimensions (%d)\n", varname, j);
			printf("NetCDF variable %s has too many dimensions (%d)\n", varname, j);
			GMT_exit (GMT, EXIT_FAILURE); return NULL;
		}
		if (ndims - in < 1) {
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "NetCDF variable %s has %" PRIuS " dimensions, cannot specify more than %d indices; ignoring remainder\n", varname, ndims, ndims-1);
			printf("NetCDF variable %s has %d dimensions, cannot specify more than %d indices; ignoring remainder\n", varname, ndims, ndims-1);
			for (j = in; j < ndims; j++) GMT->current.io.t_index[i][j] = 1;
		}
		if (ndims - in > 2)
		{
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "NetCDF variable %s has %" PRIuS " dimensions, showing only 2\n", varname, ndims);
			printf("NetCDF variable %s has %d dimensions, showing only 2\n", varname, ndims);
		}
		/* Get information of the first two dimensions */
		nc_inq_vardimid(GMT->current.io.ncid, GMT->current.io.varid[i], dimids);
		nc_inq_dimlen(GMT->current.io.ncid, dimids[0], &n);
		if (GMT->current.io.ndim != 0 && GMT->current.io.ndim != n) {
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "NetCDF variable %s has different dimension (%" PRIuS ") from others (%" PRIuS ")\n", varname, n, GMT->current.io.ndim);
			printf("NetCDF variable %s has different dimension (%d) from others (%d)\n", varname, n, GMT->current.io.ndim);
			GMT_exit (GMT, EXIT_FAILURE); return NULL;
		}
		GMT->current.io.ndim = n;
		if (dimids[1] >= 0 && ndims - in > 1) {
			nc_inq_dimlen(GMT->current.io.ncid, dimids[1], &n);
		}
		else
			n = 1;
		GMT->current.io.count[i][1] = (int)n;
		GMT->current.io.ncols += (int)n;

		/* If selected by value instead of index */
		for (j = 1; by_value && j <= in; j++) {
			nc_inq_dim (GMT->current.io.ncid, dimids[j], dimname, &n);
			nc_inq_varid (GMT->current.io.ncid, dimname, &id);
			item[0] = 0, item[1] = n-1;
			if (nc_get_att_double (GMT->current.io.ncid, id, "actual_range", dummy)) {
				nc_get_var1_double (GMT->current.io.ncid, id, &item[0], &dummy[0]);
				nc_get_var1_double (GMT->current.io.ncid, id, &item[1], &dummy[1]);
			}
			GMT->current.io.t_index[i][j] = lrint((t_value[j] - dummy[0]) / (dummy[1] - dummy[0]));
		}

		/* Get scales, offsets and missing values */
		if (nc_get_att_double(GMT->current.io.ncid, GMT->current.io.varid[i], "scale_factor", &GMT->current.io.scale_factor[i])) GMT->current.io.scale_factor[i] = 1.0;
		if (nc_get_att_double(GMT->current.io.ncid, GMT->current.io.varid[i], "add_offset", &GMT->current.io.add_offset[i])) GMT->current.io.add_offset[i] = 0.0;
		if (nc_get_att_double (GMT->current.io.ncid, GMT->current.io.varid[i], "_FillValue", &GMT->current.io.missing_value[i]) &&
		    nc_get_att_double (GMT->current.io.ncid, GMT->current.io.varid[i], "missing_value", &GMT->current.io.missing_value[i])) GMT->current.io.missing_value[i] = GMT->session.d_NaN;

		/* Scan for geographical or time units */
		if (GMT_nc_get_att_text (GMT, GMT->current.io.ncid, GMT->current.io.varid[i], "long_name", long_name, GMT_LEN256)) long_name[0] = 0;
		if (GMT_nc_get_att_text (GMT, GMT->current.io.ncid, GMT->current.io.varid[i], "units", units, GMT_LEN256)) units[0] = 0;
		GMT_str_tolower (long_name); GMT_str_tolower (units);

		if (GMT->current.io.col_type[GMT_IN][i] == GMT_IS_FLOAT)
			{ /* Float type is preset, do not alter */ }
		else if (!strcmp (long_name, "longitude") || strstr (units, "degrees_e"))
			GMT->current.io.col_type[GMT_IN][i] = GMT_IS_LON;
		else if (!strcmp (long_name, "latitude") || strstr (units, "degrees_n"))
			GMT->current.io.col_type[GMT_IN][i] = GMT_IS_LAT;
		else if (!strcmp (long_name, "time") || !strcmp (varname, "time")) {
			GMT->current.io.col_type[GMT_IN][i] = GMT_IS_RELTIME;
			GMT_memcpy (&time_system, &GMT->current.setting.time_system, 1, struct GMT_TIME_SYSTEM);
			if (GMT_get_time_system (GMT, units, &time_system) || GMT_init_time_system_structure (GMT, &time_system))
			{
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: Time units [%s] in NetCDF file not recognised, defaulting to gmt.conf.\n", units);
				printf("Warning: Time units [%s] in NetCDF file not recognised, defaulting to gmt.conf.\n", units);
			}
			/* Determine scale between data and internal time system, as well as the offset (in internal units) */
			GMT->current.io.scale_factor[i] = GMT->current.io.scale_factor[i] * time_system.scale * GMT->current.setting.time_system.i_scale;
			GMT->current.io.add_offset[i] *= time_system.scale;	/* Offset in seconds */
			GMT->current.io.add_offset[i] += GMT_DAY2SEC_F * ((time_system.rata_die - GMT->current.setting.time_system.rata_die) + (time_system.epoch_t0 - GMT->current.setting.time_system.epoch_t0));
			GMT->current.io.add_offset[i] *= GMT->current.setting.time_system.i_scale;	/* Offset in internal time units */
		}
		else if (GMT->current.io.col_type[GMT_IN][i] == GMT_IS_UNKNOWN)
			GMT->current.io.col_type[GMT_IN][i] = GMT_IS_FLOAT;
	}

	//GMT->current.io.input = gmt_nc_input;
	tmp_pointer = (size_t)(-GMT->current.io.ncid);
	return ((FILE *)tmp_pointer);
}

int GMT_z_output (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *data)
{
	int err;
	if (gmt_skip_output (GMT, data, n)) return (0);	/* Record was skipped via -s[a|r] */
	err = GMT->current.io.write_item (GMT, fp, n, data);
	/* Cast below since the output functions are declared with uint64_t but cannot really exceed 4096... SHould change uint64_t to uint32_t */
	return (err ? -1 : (int)n);	/* Return -1 if failed, else n items written */
}



int GMT_init_z_io (struct GMT_CTRL *GMT, char format[], bool repeat[], enum GMT_swap_direction swab, off_t skip, char type, struct GMT_Z_IO *r)
{
	bool first = true;
	unsigned int k;

	GMT_memset (r, 1, struct GMT_Z_IO);

	for (k = 0; k < 2; k++) {	/* Loop over the two format flags */
		switch (format[k]) {
			/* These 4 cases will set the format orientation for input */
			case 'T':
				if (first) r->format = GMT_IS_ROW_FORMAT;
				r->y_step = 1;	first = false;	break;
			case 'B':
				if (first) r->format = GMT_IS_ROW_FORMAT;
				r->y_step = -1;	first = false;	break;
			case 'L':
				if (first) r->format = GMT_IS_COL_FORMAT;
				r->x_step = 1;	first = false;	break;
			case 'R':
				if (first) r->format = GMT_IS_COL_FORMAT;
				r->x_step = -1;	first = false;	break;
			default:
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Syntax error -Z: %c not a valid format specifier!\n", format[k]);
				printf("Syntax error -Z: %c not a valid format specifier!\n", format[k]);
				GMT_exit (GMT, EXIT_FAILURE); return EXIT_FAILURE;
				break;
		}
	}

	if (!strchr ("AacuhHiIlLfd", type)) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Syntax error -Z: %c not a valid data type!\n", type);
		printf("Syntax error -Z: %c not a valid data type!\n", type);
		GMT_exit (GMT, EXIT_FAILURE); return EXIT_FAILURE;

	}

	r->x_missing = (repeat[GMT_X]) ? 1 : 0;	r->y_missing = (repeat[GMT_Y]) ? 1 : 0;
	r->skip = skip;			r->swab = swab;
	r->binary = (strchr ("Aa", type)) ? false : true;
	GMT->current.io.read_item  = GMT_get_io_ptr (GMT, GMT_IN,  swab, type);	/* Set read pointer depending on data format */
	GMT->current.io.write_item = GMT_get_io_ptr (GMT, GMT_OUT, swab, type);	/* Set write pointer depending on data format */
	GMT->common.b.type[GMT_IN] = GMT->common.b.type[GMT_OUT] = type;		/* Since -b is not setting this */
	if (r->binary) {	/* Use the binary modes (which only matters under Windoze)  */
		strcpy (GMT->current.io.r_mode, "rb");
		strcpy (GMT->current.io.w_mode, "wb");
		strcpy (GMT->current.io.a_mode, "ab+");
	}
	return (GMT_OK);
}

int GMT_parse_z_io (struct GMT_CTRL *GMT, char *txt, struct GMT_PARSE_Z_IO *z)
{
	int value;
	unsigned int i, k = 0, start;

	if (!txt) return (EXIT_FAILURE);	/* Must give a non-NULL argument */
	if (!txt[0]) return (0);		/* Default -ZTLa */

	for (start = 0; !z->not_grid && txt[start] && start < 2; start++) {	/* Loop over the first 2 flags unless dataset is not a grid */

		switch (txt[start]) {

			/* These 4 cases will set the format orientation for input */

			case 'T':
			case 'B':
			case 'L':
			case 'R':
				if (k > 2) {
					//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Syntax error -Z: Choose format from [TBLR][TBLR]!\n");
					printf("Syntax error -Z: Choose format from [TBLR][TBLR]!\n");
					return (EXIT_FAILURE);
				}
				z->format[k++] = txt[start];
				break;
			default:
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Syntax error -Z: Must begin with [TBLR][TBLR]!\n");
				printf("Syntax error -Z: Must begin with [TBLR][TBLR]!\n");
				return (EXIT_FAILURE);
				break;
		}
	}

	for (i = start; txt[i]; i++) {	/* Loop over remaining flags */

		switch (txt[i]) {

			/* Set this if file is periodic, is grid registered, but repeating column or row is missing from input */

			case 'x':
				z->repeat[GMT_X] = true;	break;
			case 'y':
				z->repeat[GMT_Y] = true;	break;

			/* Optionally skip the given number of bytes before reading data */

			case 's':
				i++;
				if (txt[i]) {	/* Read the byte count for skipping */
					value = atoi (&txt[i]);
					if (value < 0) {
						//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Syntax error -Z: Skip must be positive\n");
						printf("Syntax error -Z: Skip must be positive\n");
						return (EXIT_FAILURE);
					}
					z->skip = value;
					while (txt[i] && isdigit ((int)txt[i])) i++;
					i--;
				}
				break;

			case 'w':
				z->swab = (k_swap_in | k_swap_out); 	break;	/* Default is swap both input and output when selected */

			/* Set read pointer depending on data format */

			case 'A': /* ASCII (next regular float (%lg) from the stream) */
			case 'a': /* ASCII (1 per record) */
			case 'c': /* Binary int8_t */
			case 'u': /* Binary uint8_t */
			case 'h': /* Binary int16_t */
			case 'H': /* Binary uint16_t */
			case 'i': /* Binary int32_t */
			case 'I': /* Binary uint32_t */
			case 'l': /* Binary int64_t */
			case 'L': /* Binary uint64_t */
			case 'f': /* Binary 4-byte float */
			case 'd': /* Binary 8-byte double */
				z->type = txt[i];
				break;

			default:
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Syntax error -Z: %c not a valid modifier!\n", txt[i]);
				printf("Syntax error -Z: %c not a valid modifier!\n", txt[i]);
				return (EXIT_FAILURE);
				break;
		}
	}

	return (0);
}
uint64_t gmt_col_ij (struct GMT_Z_IO *r, struct GMT_GRID *G, uint64_t ij)
{
	/* Translates incoming ij (no padding) to gmt_ij (includes padding) for column-structured data */

	r->gmt_j = (unsigned int)(r->start_row + r->y_step * (ij % r->y_period));
	r->gmt_i = (unsigned int)(r->start_col + r->x_step * (ij / r->y_period));

	return (GMT_IJP (G->header, r->gmt_j, r->gmt_i));
}

uint64_t gmt_row_ij (struct GMT_Z_IO *r, struct GMT_GRID *G, uint64_t ij)
{
	/* Translates incoming ij (no padding) to gmt_ij (includes padding) for row-structured data */

	r->gmt_j = (unsigned int)(r->start_row + r->y_step * (ij / r->x_period));
	r->gmt_i = (unsigned int)(r->start_col + r->x_step * (ij % r->x_period));

	return (GMT_IJP (G->header, r->gmt_j, r->gmt_i));
}

int GMT_set_z_io (struct GMT_CTRL *GMT, struct GMT_Z_IO *r, struct GMT_GRID *G)
{
	/* THIS SHOULD NOT BE FATAL!
	if ((r->x_missing || r->y_missing) && G->header->registration == GMT_GRID_PIXEL_REG) return (GMT_GRDIO_RI_NOREPEAT);
	*/
	r->start_col = ((r->x_step == 1) ? 0 : G->header->nx - 1 - r->x_missing);
	r->start_row = ((r->y_step == 1) ? r->y_missing : G->header->ny - 1);
	r->get_gmt_ij = (r->format == GMT_IS_COL_FORMAT) ? gmt_col_ij : gmt_row_ij;
	r->x_period = G->header->nx - r->x_missing;
	r->y_period = G->header->ny - r->y_missing;
	r->n_expected = r->x_period * r->y_period;
	return (GMT_NOERROR);
}
void GMT_set_cartesian (struct GMT_CTRL *GMT, unsigned int dir)
{
	/* Eliminate lots of repeated statements to do this: */
	GMT->current.io.col_type[dir][GMT_X] = GMT_IS_FLOAT;
	GMT->current.io.col_type[dir][GMT_Y] = GMT_IS_FLOAT;
}

void GMT_set_xycolnames (struct GMT_CTRL *GMT, char *string)
{
	char *xy[2][2] = {{"x", "y"}, {"lon", "lat"}};
	unsigned int mode = (GMT_is_geographic (GMT, GMT_OUT)) ? 1 : 0;
	unsigned int ix = (GMT->current.setting.io_lonlat_toggle[GMT_OUT]) ? 1 : 0, iy;
	iy = 1 - ix;
	sprintf (string, "%s[0]\t%s[1]", xy[mode][ix], xy[mode][iy]);
}

void gmt_write_multilines (struct GMT_CTRL *GMT, FILE *fp, char *text, char *prefix) {
	/* Optional title(s) or remarks provided; could be several lines separated by \n */
	char p[GMT_BUFSIZ] = {""}, line[GMT_BUFSIZ] = {""};
	unsigned int pos = 0, k = 0;

	while (GMT_strtok (text, "\\", &pos, p)) {
		sprintf (line, "# %7s : %s", prefix, &p[k]);
		GMT_write_tableheader (GMT, fp, line);
		k = 1;	/* Need k to skip the n in \n */
	}
}

void GMT_write_tableheader (struct GMT_CTRL *GMT, FILE *fp, char *txt)
{
	/* Output ASCII segment header; skip if mode is binary.
	 * We append a newline (\n) if not is present */

	if (!GMT->current.setting.io_header[GMT_OUT]) return;	/* No output headers requested */
	if (GMT_binary_header (GMT, GMT_OUT))		/* Must write a binary header */
		GMT_io_binary_header (GMT, fp, GMT_OUT);
	else if (!txt || !txt[0])				/* Blank header */
		fprintf (fp, "#\n");
	else {
		if (txt[0] != '#') fputc ('#', fp);	/* Make sure we have # at start */
		fprintf (fp, "%s", txt);
		if (txt[strlen(txt)-1] != '\n') fputc ('\n', fp);	/* Make sure we have \n at end */
	}
}

void GMT_write_newheaders (struct GMT_CTRL *GMT, FILE *fp, uint64_t n_cols)
{	/* Common ascii header records added on output */
	if (GMT->common.b.active[GMT_OUT]) return;		/* No output headers for binary files */
	if (!GMT->current.setting.io_header[GMT_OUT]) return;	/* No output headers requested, so don't bother */
	if (GMT->common.h.title) {	/* Optional title(s) provided; could be several lines separated by \n */
		gmt_write_multilines (GMT, fp, GMT->common.h.title, "Title");
	}
	/* Always write command line */
	GMT_write_tableheader (GMT, fp, GMT_create_header_item (GMT->parent, GMT_COMMENT_IS_COMMAND | GMT_COMMENT_IS_OPTION, GMT->current.options));
	if (GMT->common.h.remark) {	/* Optional remark(s) provided; could be several lines separated by \n */
		gmt_write_multilines (GMT, fp, GMT->common.h.remark, "Remark");
	}
	if (GMT->common.h.add_colnames) {	/* Want output comment with column names */
		if (GMT->common.h.colnames)	/* Optional column names already provided */
			GMT_write_tableheader (GMT, fp, GMT->common.h.colnames);
		else if (n_cols) {	/* Generate names col1[0], col2[1] etc */
			uint64_t col, first = 1;
			char record[GMT_BUFSIZ] = {""}, txt[GMT_LEN64] = {""};
			if (n_cols >= 2) {	/* Place x and y first */
				GMT_set_xycolnames (GMT, record);
				first++;
			}
			else
				sprintf (record, "col1[0]");
			for (col = first; col < n_cols; col++) {
				sprintf (txt, "\tcol%" PRIu64 "[%" PRIu64 "]", col+1, col);
				strcat (record, txt);
			}
			GMT_write_tableheader (GMT, fp, record);
		}
	}
}


