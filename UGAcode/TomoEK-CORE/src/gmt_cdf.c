/*--------------------------------------------------------------------
 *	$Id: gmt_cdf.c 12822 2014-01-31 23:39:56Z remko $
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
 *
 *	G M T _ C D F . C   R O U T I N E S
 *
 * Takes care of all grd input/output built on NCAR's netCDF routines (which is
 * an XDR implementation)
 * Most functions will return with error message if an internal error is returned.
 * There functions are only called indirectly via the GMT_* grdio functions.
 *
 * Author:	Paul Wessel
 * Date:	1-JAN-2010
 * Version:	5
 *
 * Public functions:
 *
 *	GMT_cdf_read_grd_info :		Read header from file
 *	GMT_cdf_read_grd :		Read header and data set from file
 *	GMT_cdf_update_grd_info :	Update header in existing file
 *	GMT_cdf_write_grd_info :	Write header to new file
 *	GMT_cdf_write_grd :		Write header and data set to new file
 *
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
//#include "gmt_io.h"
#include <netcdf.h>
#include <limits.h>
#include "gmt_grdio.h"
#include "gmt_type.h"
#include "memory.h"
#include "gmt_macros.h"
#include "gmt_resources.h"
#include "gmt_grd.h"
#include "gmt_common.h"


#define PRIo64       "I64o"
#define PRIu64       "I64u"
#define PRIx64       "I64x"
#define PRIX64       "I64X"


int gmt_cdf_grd_info (struct GMT_CTRL *GMT, int ncid, struct GMT_GRID_HEADER *header, char job)
{
	int err;	/* Implicity by GMT_err_trap */
	int nm[2];
	double dummy[2];
	char text[GMT_GRID_COMMAND_LEN320+GMT_GRID_REMARK_LEN160];
	nc_type z_type;

	/* Dimension ids, varibale ids, etc. */
	int side_dim, xysize_dim, x_range_id, y_range_id, z_range_id, inc_id, nm_id, z_id, dims[1];

	/* Define and get dimensions and variables */

	//GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Enter gmt_cdf_grd_info with argument %c\n", (int)job);

	if (job == 'w') {
		GMT_err_trap (nc_def_dim (ncid, "side", 2U, &side_dim));
		GMT_err_trap (nc_def_dim (ncid, "xysize", header->nm, &xysize_dim));

		dims[0]	= side_dim;
		GMT_err_trap (nc_def_var (ncid, "x_range", NC_DOUBLE, 1, dims, &x_range_id));
		GMT_err_trap (nc_def_var (ncid, "y_range", NC_DOUBLE, 1, dims, &y_range_id));
		GMT_err_trap (nc_def_var (ncid, "z_range", NC_DOUBLE, 1, dims, &z_range_id));
		GMT_err_trap (nc_def_var (ncid, "spacing", NC_DOUBLE, 1, dims, &inc_id));
		GMT_err_trap (nc_def_var (ncid, "dimension", NC_LONG, 1, dims, &nm_id));

		switch (header->type) {
			case GMT_GRID_IS_CB: z_type = NC_BYTE; break;
			case GMT_GRID_IS_CS: z_type = NC_SHORT; break;
			case GMT_GRID_IS_CI: z_type = NC_INT; break;
			case GMT_GRID_IS_CF: z_type = NC_FLOAT; break;
			case GMT_GRID_IS_CD: z_type = NC_DOUBLE; break;
			default:			z_type = NC_NAT;
		}

		dims[0]	= xysize_dim;
		GMT_err_trap (nc_def_var (ncid, "z", z_type, 1, dims, &z_id));
	}
	else {
		GMT_err_trap (nc_inq_varid (ncid, "x_range", &x_range_id));
		GMT_err_trap (nc_inq_varid (ncid, "y_range", &y_range_id));
		GMT_err_trap (nc_inq_varid (ncid, "z_range", &z_range_id));
		GMT_err_trap (nc_inq_varid (ncid, "spacing", &inc_id));
		GMT_err_trap (nc_inq_varid (ncid, "dimension", &nm_id));
		GMT_err_trap (nc_inq_varid (ncid, "z", &z_id));
		GMT_err_trap (nc_inq_vartype (ncid, z_id, &z_type));
		switch (z_type) {
			case NC_BYTE:   header->type = GMT_GRID_IS_CB; break;
			case NC_SHORT:  header->type = GMT_GRID_IS_CS; break;
			case NC_INT:    header->type = GMT_GRID_IS_CI; break;
			case NC_FLOAT:  header->type = GMT_GRID_IS_CF; break;
			case NC_DOUBLE: header->type = GMT_GRID_IS_CD; break;
			default:        header->type = k_grd_unknown_fmt; break;
		}
	}
	header->z_id = z_id;

	/* Get or assign attributes */

	GMT_memset (text, GMT_GRID_COMMAND_LEN320+GMT_GRID_REMARK_LEN160, char);

	if (job == 'u') GMT_err_trap (nc_redef (ncid));

	if (job == 'r') {
		int reg;
		GMT_err_trap (nc_get_att_text (ncid, x_range_id, "units", header->x_units));
		GMT_err_trap (nc_get_att_text (ncid, y_range_id, "units", header->y_units));
		GMT_err_trap (nc_get_att_text (ncid, z_range_id, "units", header->z_units));
		GMT_err_trap (nc_get_att_double (ncid, z_id, "scale_factor", &header->z_scale_factor));
		GMT_err_trap (nc_get_att_double (ncid, z_id, "add_offset", &header->z_add_offset));
		GMT_err_trap (nc_get_att_int (ncid, z_id, "node_offset", &reg));
		header->registration = reg;
		nc_get_att_float (ncid, z_id, "_FillValue", &header->nan_value);
		GMT_err_trap (nc_get_att_text (ncid, NC_GLOBAL, "title", header->title));
		GMT_err_trap (nc_get_att_text (ncid, NC_GLOBAL, "source", text));
		strncpy (header->command, text, GMT_GRID_COMMAND_LEN320);
		strncpy (header->remark, &text[GMT_GRID_COMMAND_LEN320], GMT_GRID_REMARK_LEN160);

		GMT_err_trap (nc_get_var_double (ncid, x_range_id, dummy));
		header->wesn[XLO] = dummy[0];
		header->wesn[XHI] = dummy[1];
		GMT_err_trap (nc_get_var_double (ncid, y_range_id, dummy));
		header->wesn[YLO] = dummy[0];
		header->wesn[YHI] = dummy[1];
		GMT_err_trap (nc_get_var_double (ncid, inc_id, dummy));
		header->inc[GMT_X] = dummy[0];
		header->inc[GMT_Y] = dummy[1];
		GMT_err_trap (nc_get_var_int (ncid, nm_id, nm));
		header->nx = nm[0];
		header->ny = nm[1];
		GMT_err_trap (nc_get_var_double (ncid, z_range_id, dummy));
		header->z_min = dummy[0];
		header->z_max = dummy[1];
		header->row_order = k_nc_start_north; /* N->S */
	}
	else {
		int reg;
		strncpy (text, header->command, GMT_GRID_COMMAND_LEN320);
		strncpy (&text[GMT_GRID_COMMAND_LEN320], header->remark, GMT_GRID_REMARK_LEN160);
		GMT_err_trap (nc_put_att_text (ncid, x_range_id, "units", GMT_GRID_UNIT_LEN80, header->x_units));
		GMT_err_trap (nc_put_att_text (ncid, y_range_id, "units", GMT_GRID_UNIT_LEN80, header->y_units));
		GMT_err_trap (nc_put_att_text (ncid, z_range_id, "units", GMT_GRID_UNIT_LEN80, header->z_units));
		GMT_err_trap (nc_put_att_double (ncid, z_id, "scale_factor", NC_DOUBLE, 1U, &header->z_scale_factor));
		GMT_err_trap (nc_put_att_double (ncid, z_id, "add_offset", NC_DOUBLE, 1U, &header->z_add_offset));
		if (z_type != NC_FLOAT && z_type != NC_DOUBLE)
			header->nan_value = rintf (header->nan_value); /* round to integer */
		GMT_err_trap (nc_put_att_float (ncid, z_id, "_FillValue", z_type, 1U, &header->nan_value));
		reg = header->registration;
		GMT_err_trap (nc_put_att_int (ncid, z_id, "node_offset", NC_LONG, 1U, &reg));
		GMT_err_trap (nc_put_att_text (ncid, NC_GLOBAL, "title", GMT_GRID_TITLE_LEN80, header->title));
		GMT_err_trap (nc_put_att_text (ncid, NC_GLOBAL, "source", GMT_GRID_COMMAND_LEN320+GMT_GRID_REMARK_LEN160, text));

		GMT_err_trap (nc_enddef (ncid));

		dummy[0] = header->wesn[XLO];	dummy[1] = header->wesn[XHI];
		GMT_err_trap (nc_put_var_double (ncid, x_range_id, dummy));
		dummy[0] = header->wesn[YLO];	dummy[1] = header->wesn[YHI];
		GMT_err_trap (nc_put_var_double (ncid, y_range_id, dummy));
		dummy[0] = header->inc[GMT_X];	dummy[1] = header->inc[GMT_Y];
		GMT_err_trap (nc_put_var_double (ncid, inc_id, dummy));
		nm[0] = header->nx;	nm[1] = header->ny;
		GMT_err_trap (nc_put_var_int (ncid, nm_id, nm));
		if (header->z_min <= header->z_max) {
			dummy[0] = header->z_min; dummy[1] = header->z_max;
		}
		else {
			dummy[0] = 0.0; dummy[1] = 0.0;
		}
		GMT_err_trap (nc_put_var_double (ncid, z_range_id, dummy));
	}
	return (GMT_NOERROR);
}

int GMT_cdf_read_grd_info (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *header)
{
	int ncid;
	int err;
	if (!strcmp (header->name,"=")) return (GMT_GRDIO_NC_NO_PIPE);
	GMT_err_trap (nc_open (header->name, NC_NOWRITE, &ncid));
	GMT_err_trap (gmt_cdf_grd_info (GMT, ncid, header, 'r'));
	GMT_err_trap (nc_close (ncid));
	return (GMT_NOERROR);
}

int GMT_cdf_update_grd_info (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *header)
{
	int ncid, old_fill_mode;
	int err;
	if (!strcmp (header->name,"=")) return (GMT_GRDIO_NC_NO_PIPE);
	GMT_err_trap (nc_open (header->name, NC_WRITE, &ncid));
	GMT_err_trap (nc_set_fill (ncid, NC_NOFILL, &old_fill_mode)); 
	GMT_err_trap (gmt_cdf_grd_info (GMT, ncid, header, 'u'));
	GMT_err_trap (nc_close (ncid));
	return (GMT_NOERROR);
}

int GMT_cdf_write_grd_info (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *header)
{
	int ncid, old_fill_mode;
	int err;
	if (!strcmp (header->name,"=")) return (GMT_GRDIO_NC_NO_PIPE);
	GMT_err_trap (nc_create (header->name, NC_CLOBBER, &ncid));
	GMT_err_trap (nc_set_fill (ncid, NC_NOFILL, &old_fill_mode)); 
	GMT_err_trap (gmt_cdf_grd_info (GMT, ncid, header, 'w'));
	GMT_err_trap (nc_close (ncid));
	return (GMT_NOERROR);
}

bool GMT_init_complex (struct GMT_GRID_HEADER *header, unsigned int complex_mode, uint64_t *imag_offset)
{	/* Sets complex-related parameters based on the input complex_mode variable:
	 * If complex_mode & GMT_GRID_NO_HEADER then we do NOT want to write a header [output only; only some formats]
	 * If grid is the imaginary components of a complex grid then we compute the offset
	 * from the start of the complex array where the first imaginary value goes, using the serial arrangement.
	 */

	bool do_header = !(complex_mode & GMT_GRID_NO_HEADER);	/* Want no header if this bit is set */
	/* Imaginary components are stored after the real components if complex */
	*imag_offset = (complex_mode & GMT_GRID_IS_COMPLEX_IMAG) ? header->size / 2ULL : 0ULL;
	
	return (do_header);
}


int GMT_cdf_read_grd (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *header, float *grid, double wesn[], unsigned int *pad, unsigned int complex_mode)
{	/* header:	grid structure header
	 * grid:	array with final grid
	 * wesn:	Sub-region to extract  [Use entire file if 0,0,0,0]
	 * padding:	# of empty rows/columns to add on w, e, s, n of grid, respectively
	 * complex_mode:	&4 | &8 if complex array is to hold real (4) and imaginary (8) parts (otherwise read as real only)
	 *		Note: The file has only real values, we simply allow space in the complex array
	 *		for real and imaginary parts when processed by grdfft etc.
	 *
	 * Reads a subset of a grid file and optionally pads the array with extra rows and columns
	 * header values for nx and ny are reset to reflect the dimensions of the logical array,
	 * not the physical size (i.e., the padding is not counted in nx and ny)
	 */

	int  ncid;
	size_t start[1], edge[1];
	bool check;
	int j, err;
	int first_col, last_col, first_row, last_row;
	unsigned int i, width_in, height_in;
	unsigned int width_out, *actual_col = NULL;
	uint64_t ij, kk, imag_offset;
	float *tmp = NULL;

	//GMT_err_pass (GMT, GMT_grd_prep_io (GMT, header, wesn, &width_in, &height_in, &first_col, &last_col, &first_row, &last_row, &actual_col), header->name);
	(void)GMT_init_complex (header, complex_mode, &imag_offset);	/* Set offset for imaginary complex component */

	width_out = width_in;		/* Width of output array */
	if (pad[XLO] > 0) width_out += pad[XLO];
	if (pad[XHI] > 0) width_out += pad[XHI];

	/* Open the NetCDF file */

	if (!strcmp (header->name,"=")) return (GMT_GRDIO_NC_NO_PIPE);
	GMT_err_trap (nc_open (header->name, NC_NOWRITE, &ncid));
	check = !isnan (header->nan_value);

	/* Load data row by row. The data in the file is stored in the same
	 * "upside down" fashion as within GMT. The first row is the top row */

	tmp = GMT_memory (GMT, NULL, header->nx, float);

	edge[0] = header->nx;
	ij = imag_offset + pad[YHI] * width_out + pad[XLO];
	header->z_min =  DBL_MAX;
	header->z_max = -DBL_MAX;

	for (j = first_row; j <= last_row; j++, ij += width_out) {
		start[0] = j * header->nx;
		GMT_err_trap (nc_get_vara_float (ncid, header->z_id, start, edge, tmp));	/* Get one row */
		for (i = 0; i < width_in; i++) {	/* Check for and handle NaN proxies */
			kk = ij+i;
			grid[kk] = tmp[actual_col[i]];
			if (check && grid[kk] == header->nan_value)
				grid[kk] = GMT->session.f_NaN;
			if (GMT_is_fnan (grid[kk])) continue;
			header->z_min = MIN (header->z_min, (double)grid[kk]);
			header->z_max = MAX (header->z_max, (double)grid[kk]);
		}
	}

	header->nx = width_in;
	header->ny = height_in;
	GMT_memcpy (header->wesn, wesn, 4, double);

	GMT_err_trap (nc_close (ncid));

	GMT_free (GMT, actual_col);
	GMT_free (GMT, tmp);
	return (GMT_NOERROR);
}

int GMT_cdf_write_grd (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *header, float *grid, double wesn[], unsigned int *pad, unsigned int complex_mode)
{	/* header:	grid structure header
	 * grid:	array with final grid
	 * wesn:	Sub-region to write out  [Use entire file if 0,0,0,0]
	 * padding:	# of empty rows/columns to add on w, e, s, n of grid, respectively
	 * complex_mode:	&1 | &2 if complex array is to hold real (1) and imaginary (2) parts (otherwise read as real only)
	 *		Note: The file has only real values, we simply allow space in the complex array
	 *		for real and imaginary parts when processed by grdfft etc.
	 */

	size_t start[1], edge[1];
	int ncid, old_fill_mode;
	long *tmp_i = NULL;
	int err;
	unsigned int i, *actual_col = NULL;
	unsigned int j, width_out, height_out, width_in;
	int first_col, last_col, first_row, last_row;
	uint64_t ij, nr_oor = 0, imag_offset;
	float *tmp_f = NULL;
	double limit[2] = {-FLT_MAX, FLT_MAX};
	float value;
	nc_type z_type;

	/* Determine the value to be assigned to missing data, if not already done so */

	switch (header->type) {
		case GMT_GRID_IS_CB:
			if (isnan (header->nan_value)) header->nan_value = CHAR_MIN;
			limit[0] = CHAR_MIN - 0.5; limit[1] = CHAR_MAX + 0.5;
			z_type = NC_BYTE; break;
		case GMT_GRID_IS_CS:
			if (isnan (header->nan_value)) header->nan_value = SHRT_MIN;
			limit[0] = SHRT_MIN - 0.5; limit[1] = SHRT_MAX + 0.5;
			z_type = NC_SHORT; break;
		case GMT_GRID_IS_CI:
			if (isnan (header->nan_value)) header->nan_value = INT_MIN;
			limit[0] = INT_MIN - 0.5; limit[1] = INT_MAX + 0.5;
			z_type = NC_INT; break;
		case GMT_GRID_IS_CF:
			z_type = NC_FLOAT; break;
		case GMT_GRID_IS_CD:
			z_type = NC_DOUBLE; break;
		default:
			z_type = NC_NAT;
	}

	//GMT_err_pass (GMT, GMT_grd_prep_io (GMT, header, wesn, &width_out, &height_out, &first_col, &last_col, &first_row, &last_row, &actual_col), header->name);
	(void)GMT_init_complex (header, complex_mode, &imag_offset);	/* Set offset for imaginary complex component */

	width_in = width_out;		/* Physical width of input array */
	if (pad[XLO] > 0) width_in += pad[XLO];
	if (pad[XHI] > 0) width_in += pad[XHI];

	GMT_memcpy (header->wesn, wesn, 4, double);
	header->nx = width_out;
	header->ny = height_out;

	/* Write grid header */

	if (!strcmp (header->name,"=")) return (GMT_GRDIO_NC_NO_PIPE);
	GMT_err_trap (nc_create (header->name, NC_CLOBBER, &ncid));
	GMT_err_trap (nc_set_fill (ncid, NC_NOFILL, &old_fill_mode)); 
	GMT_err_trap (gmt_cdf_grd_info (GMT, ncid, header, 'w'));

	/* Set start position for writing grid */

	edge[0] = width_out;
	ij = first_col + pad[XLO] + (first_row + pad[YHI]) * width_in;
	header->z_min =  DBL_MAX;
	header->z_max = -DBL_MAX;

	/* Store z-variable */

	if (z_type == NC_FLOAT || z_type == NC_DOUBLE) {
		tmp_f = GMT_memory (GMT, NULL, width_in, float);
		for (j = 0; j < height_out; j++, ij += width_in) {
			start[0] = j * width_out;
			for (i = 0; i < width_out; i++) {
				value = grid[ij+actual_col[i]+imag_offset];
				if (!isfinite (value)) {
					if (isinf(value))
						nr_oor++; /* out of float range */
					tmp_f[i] = header->nan_value;
				}
				else {
					tmp_f[i] = value;
					header->z_min = MIN (header->z_min, (double)tmp_f[i]);
					header->z_max = MAX (header->z_max, (double)tmp_f[i]);
				}
			}
			GMT_err_trap (nc_put_vara_float (ncid, header->z_id, start, edge, tmp_f));
		}
		GMT_free (GMT, tmp_f);
	}
	else {
		tmp_i = GMT_memory (GMT, NULL, width_in, long);
		for (j = 0; j < height_out; j++, ij += width_in) {
			start[0] = j * width_out;
			for (i = 0; i < width_out; i++) {
				value = grid[ij+actual_col[i]+imag_offset];
				if (!isfinite (value)) {
					if (isinf(value))
						nr_oor++; /* out of float range */
					tmp_i[i] = lrintf (header->nan_value);
				}
				else if (value <= limit[0] || value >= limit[1]) {
					tmp_i[i] = lrintf (header->nan_value);
					nr_oor++;
				}
				else {
					tmp_i[i] = lrintf (value);
					header->z_min = MIN (header->z_min, (double)tmp_i[i]);
					header->z_max = MAX (header->z_max, (double)tmp_i[i]);
				}
			}
			GMT_err_trap (nc_put_vara_long (ncid, header->z_id, start, edge, tmp_i));
		}
		GMT_free (GMT, tmp_i);
	}

	if (nr_oor > 0) 
		printf ("Warning: %" PRIu64 " out-of-range grid values converted to _FillValue [%s]\n", nr_oor, header->name);

	GMT_free (GMT, actual_col);

	if (header->z_min <= header->z_max) {
		limit[0] = header->z_min; limit[1] = header->z_max;
	}
	else {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: No valid values in grid [%s]\n", header->name);
		limit[0] = 0.0; limit[1] = 0.0;
	}
	GMT_err_trap (nc_put_var_double (ncid, header->z_id - 3, limit));

	/* Close grid */

	GMT_err_trap (nc_close (ncid));

	return (GMT_NOERROR);
}
