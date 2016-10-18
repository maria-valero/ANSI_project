/*
 * grd_io.h
 *
 *  Created on: Feb 22, 2015
 *      Author: nishita
 */

#ifndef GRD_IO_H_
#define GRD_IO_H_
#include <stdbool.h>
#include <stdio.h>
#include "gmt_resources.h"

#define GMT_to_inch(GMT,value) GMT_convert_units (GMT, value, GMT->current.setting.proj_length_unit, GMT_INCH)

static inline void scale_and_offset_f (float *data, size_t length, float scale, float offset) {
	/* Routine that scales and offsets the data in a vector
	 *  data:   Single-precision real input vector
	 *  length: The number of elements to process
	 * This function uses the vDSP portion of the Accelerate framework if possible */
	size_t n;
	if (scale == 1.0) /* offset only */
		for (n = 0; n < length; ++n)
			data[n] += offset;
	else if (offset == 0.0) /* scale only */
		for (n = 0; n < length; ++n)
			data[n] *= scale;
	else /* scale + offset */
		for (n = 0; n < length; ++n)
			data[n] = data[n] * scale + offset;
}
int GMT_access (struct GMT_CTRL *GMT, const char* filename, int mode);
int GMT_scanf (struct GMT_CTRL *GMT, char *s, unsigned int expectation, double *val);
int GMT_write_grd_info (struct GMT_CTRL *GMT, char *file, struct GMT_GRID_HEADER *header);
void gmt_grd_get_units (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *header);
void GMT_grd_zminmax (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *h, float *z);
void GMT_io_binary_header (struct GMT_CTRL *GMT, FILE *fp, unsigned int dir);
int GMT_adjust_loose_wesn (struct GMT_CTRL *GMT, double wesn[], struct GMT_GRID_HEADER *header);
struct GMT_DATASET * GMT_duplicate_dataset (struct GMT_CTRL *GMT, struct GMT_DATASET *Din, unsigned int mode, unsigned int *geometry);
struct GMT_GRID *GMT_duplicate_grid (struct GMT_CTRL *GMT, struct GMT_GRID *G, unsigned int mode);
int GMT_read_grd (struct GMT_CTRL *GMT, char *file, struct GMT_GRID_HEADER *header, float *grid, double *wesn, unsigned int *pad, int complex_mode);

struct GMT_MATRIX * GMT_duplicate_matrix (struct GMT_CTRL *GMT, struct GMT_MATRIX *M_in, bool duplicate_data);
int GMT_fclose (struct GMT_CTRL *GMT, FILE *stream);
FILE *GMT_fopen (struct GMT_CTRL *GMT, const char *filename, const char *mode);
void GMT_free_ogr (struct GMT_CTRL *GMT, struct GMT_OGR **G, unsigned int mode);
void GMT_free_table (struct GMT_CTRL *GMT, struct GMT_DATATABLE *table, enum GMT_enum_alloc alloc_mode);
void GMT_grd_pad_on (struct GMT_CTRL *GMT, struct GMT_GRID *G, unsigned int *pad);
void GMT_set_dataset_minmax (struct GMT_CTRL *GMT, struct GMT_DATASET *D);
char *GMT_trim_segheader (struct GMT_CTRL *GMT, char *line);
int GMT_update_grd_info (struct GMT_CTRL *GMT, char *file, struct GMT_GRID_HEADER *header);
int GMT_alloc_segment (struct GMT_CTRL *GMT, struct GMT_DATASEGMENT *S, uint64_t n_rows, uint64_t n_columns, bool first);
int GMT_alloc_univector (struct GMT_CTRL *GMT, union GMT_UNIVECTOR *u, unsigned int type, uint64_t n_rows);
int gmt_alloc_vectors (struct GMT_CTRL *GMT, struct GMT_VECTOR *V);
void GMT_assign_segment (struct GMT_CTRL *GMT, struct GMT_DATASEGMENT *S, uint64_t n_rows, uint64_t n_columns);
char *GMT_getdatapath (struct GMT_CTRL *GMT, const char *stem, char *path, int mode);
bool GMT_check_url_name (char *fname);
struct GMT_GRID * GMT_create_grid (struct GMT_CTRL *GMT);
#endif /* GRD_IO_H_ */
