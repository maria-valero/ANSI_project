/*
 * gmt_io.h
 *
 *  Created on: Feb 23, 2015
 *      Author: nishita
 */

#ifndef GMT_IO_H_
#define GMT_IO_H_
#include "common_math.h"
#include "gmt_resources.h"
#include "gmt_type.h"
#include "gmt_constants.h"


/* Macros to simplify check for return status */
#define GMT_REC_IS_TABLE_HEADER(C)	(C->current.io.status & GMT_IO_TABLE_HEADER)
#define GMT_REC_IS_SEGMENT_HEADER(C)	(C->current.io.status & GMT_IO_SEGMENT_HEADER)
#define GMT_REC_IS_ANY_HEADER(C)	(C->current.io.status & GMT_IO_ANY_HEADER)
#define GMT_REC_IS_ERROR(C)		(C->current.io.status & GMT_IO_MISMATCH)
#define GMT_REC_IS_EOF(C)		(C->current.io.status & GMT_IO_EOF)
#define GMT_REC_IS_NAN(C)		(C->current.io.status & GMT_IO_NAN)
#define GMT_REC_IS_GAP(C)		(C->current.io.status & GMT_IO_GAP)
#define GMT_REC_IS_NEW_SEGMENT(C)	(C->current.io.status & GMT_IO_NEW_SEGMENT)
#define GMT_REC_IS_LINE_BREAK(C)	(C->current.io.status & GMT_IO_LINE_BREAK)
#define GMT_REC_IS_FILE_BREAK(C)	(C->current.io.status & GMT_IO_NEXT_FILE)
#define GMT_REC_IS_DATA(C)		(C->current.io.status == 0 || C->current.io.status == GMT_IO_NAN)

#define GMT_polygon_is_hole(S) (S->pol_mode == GMT_IS_HOLE || (S->ogr && S->ogr->pol_mode == GMT_IS_HOLE))
/* Must add M, m, E, Z, and/or S to the common option processing list */
#define GMT_OPT(opt) opt
enum GMT_enum_ogr {
	GMT_IS_LINESTRING = 2,
	GMT_IS_POLYGON,
	GMT_IS_MULTIPOINT,
	GMT_IS_MULTILINESTRING,
	GMT_IS_MULTIPOLYGON};

/* Return codes from GMT_inonout */
enum GMT_enum_inside {
	GMT_OUTSIDE = 0,
	GMT_ONEDGE,
	GMT_INSIDE};

	/* How to handle NaNs in records */

enum GMT_io_nan_enum {
	GMT_IO_NAN_OK = 0,	/* NaNs are fine; just ouput the record as is */
	GMT_IO_NAN_SKIP,	/* -s[cols]	: Skip records with z == NaN in selected cols [z-col only] */
	GMT_IO_NAN_KEEP,	/* -sr		: Skip records with z != NaN */
	GMT_IO_NAN_ONE};	/* -sa		: Skip records with at least one NaN */

struct GMT_PARSE_Z_IO {	/* -Z[<flags>] */
	bool active;		/* true if selected */
	bool not_grid;		/* false if binary data file is a grid so organization matters */
	bool repeat[2];		/* true if periodic in x|y and repeating row/col is missing */
	enum GMT_swap_direction swab;	/* k_swap_none = no byte swapping, k_swap_inswaps input, k_swap_out swaps output, combine to swap both */
	off_t skip;		/* Initial bytes to skip before reading */
	char type;		/* Data type flag A|a|c|u|h|H|i|I|l|L|f|d */
	char format[2];		/* 2-char code describing row/col organization for grids */
};

struct GMT_Z_IO {		/* Used when processing z(x,y) table input when (x,y) is implicit */
	bool swab;		/* true if we must swap byte-order */
	bool binary;		/* true if we are reading/writing binary data */
	bool input;		/* true if we are reading, false if we are writing */
	int x_step;	/* +1 if logical x values increase to right, else -1 */
	int y_step;	/* +1 if logical y values increase upwards, else -1 */
	unsigned int x_missing;	/* 1 if a periodic (right) column is implicit (i.e., not stored) */
	unsigned int y_missing;	/* 1 if a periodic (top) row is implicit (i.e., not stored) */
	unsigned int format;	/* Either GMT_IS_COL_FORMAT or GMT_IS_ROW_FORMAT */
	unsigned int x_period;	/* length of a row in the input data ( <= nx, see x_missing) */
	unsigned int y_period;	/* length of a col in the input data ( <= ny, see y_missing) */
	unsigned int start_col;	/* First logical column in file */
	unsigned int start_row;	/* First logical row in file */
	unsigned int gmt_i;		/* Current column number in the GMT registered grid */
	unsigned int gmt_j;		/* Current row number in the GMT registered grid */
	uint64_t n_expected;	/* Number of data element expected to be read */
	off_t skip;		/* Number of bytes to skip before reading data */
	uint64_t (*get_gmt_ij) (struct GMT_Z_IO *, struct GMT_GRID *, uint64_t);	/* Pointer to function that converts running number to GMT ij */
};

#define GMT_SAME_LATITUDE(A,B)  (doubleAlmostEqualZero (A,B))			/* A and B are the same latitude */
#define gmt_convert_col(S,x) {if (S.convert) x = ((S.convert == 2) ? log10 (x) : x) * S.scale + S.offset;}

int GMT_scanf_arg (struct GMT_CTRL *GMT, char *s, unsigned int expectation, double *val);
unsigned int GMT_verify_expectations (struct GMT_CTRL *GMT, unsigned int wanted, unsigned int got, char *item);
struct GMT_DATATABLE * GMT_read_table (struct GMT_CTRL *GMT, void *source, unsigned int source_type, bool greenwich, unsigned int *geometry, bool use_GMT_io);
int GMT_write_dataset (struct GMT_CTRL *GMT, void *dest, unsigned int dest_type, struct GMT_DATASET *D, bool use_GMT_io, int table);
void * GMT_ascii_textinput (struct GMT_CTRL *GMT, FILE *fp, uint64_t *n, int *status);
bool GMT_is_a_blank_line (char *line);
unsigned int gmt_is_segment_header (struct GMT_CTRL *GMT, char *line);

unsigned int gmt_assign_aspatial_cols (struct GMT_CTRL *GMT);
void * gmt_ascii_input (struct GMT_CTRL *GMT, FILE *fp, uint64_t *n, int *status);
bool gmt_ogr_header_parser (struct GMT_CTRL *GMT, char *record);
int GMT_ascii_output (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *ptr);

int GMT_nc_get_att_text (struct GMT_CTRL *GMT, int ncid, int varid, char *name, char *text, size_t textlen);
int GMT_z_output (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *data);

#endif /* GMT_IO_H_ */
