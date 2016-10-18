#include <stdarg.h>
#include <netcdf.h>
#include <stdio.h>
#include <netcdf.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <unistd.h>
#include <time.h>
#include "grd_io.h"
#include "gmt_io.h"
#include "gmt_init.h"
#include "gmt_private.h"
#include "surface.h"
#include "gmt_type.h"
#include "memory.h"
#include "gmt_macros.h"
#include "gmt_defaults.h"
#include "gmt_nan.h"
#include "gmt_api.h"

#define assert(e) ((void)0)
/* Two different i/o mode: GMT_Put|Get_Data vs GMT_Put|Get_Record */
enum GMT_enum_iomode {
	GMT_BY_SET 	= 0,	/* Default is to read the entire set */
	GMT_BY_REC	= 1};	/* Means we will access the registere files on a record-by-record basis */

extern const char* g_api_error_string[];
extern int GMT_surface (void *V_API, int mode, void *args);
static inline struct GMT_GRID    *gmt_get_grid_data (struct GMT_GRID *ptr) {return (ptr);}
static inline struct GMT_GRID_ROWBYROW * gmt_get_rbr_ptr (struct GMT_GRID_ROWBYROW *ptr) {return (ptr);}

static inline struct GMTAPI_CTRL * gmt_get_api_ptr (struct GMTAPI_CTRL *ptr) {return (ptr);} //nishita
//static inline struct GMT_PALETTE * gmt_get_cpt_ptr (struct GMT_PALETTE **ptr) {return (*ptr);}
static inline struct GMT_DATASET * gmt_get_dataset_ptr (struct GMT_DATASET **ptr) {return (*ptr);}
static inline struct GMT_TEXTSET * gmt_get_textset_ptr (struct GMT_TEXTSET **ptr) {return (*ptr);}
static inline struct GMT_GRID    * gmt_get_grid_ptr (struct GMT_GRID **ptr) {return (*ptr);}
//static inline struct GMT_MATRIX  * gmt_get_matrix_ptr (struct GMT_MATRIX **ptr) {return (*ptr);}
//static inline struct GMT_VECTOR  * gmt_get_vector_ptr (struct GMT_VECTOR **ptr) {return (*ptr);}
static inline double      	 * gmt_get_coord_ptr (double **ptr) {return (*ptr);}
/* Misc. local text strings needed in this file only, used when debug verbose is on (-Vd) */

static const char *GMT_method[] = {"File", "Stream", "File Descriptor", "Memory Copy", "Memory Reference"};
static const char *GMT_family[] = {"Data Table", "Text Table", "GMT Grid", "CPT Table", "GMT Image", "GMT Vector", "GMT Matrix", "GMT Coord"};
static const char *GMT_via[] = {"User Vector", "User Matrix"};
static const char *GMT_direction[] = {"Input", "Output"};
static const char *GMT_stream[] = {"Standard", "User-supplied"};
//static const char *GMT_status[] = {"Unused", "In-use", "Used"};
static const char *GMT_geometry[] = {"Not Set", "Point", "Line", "Polygon", "Point|Line|Poly", "Line|Poly", "Surface", "Non-Geographical"};

/* p_func_size_t is used as a pointer to functions that returns a size_t dimension */
typedef size_t (*p_func_size_t) (uint64_t row, uint64_t col, size_t dim);
#define GMT_REC_IS_DATA(C)		(C->current.io.status == 0 || C->current.io.status == GMT_IO_NAN)

/* We asked for subset of grid if the wesn pointer is not NULL but indicates a region */
#define full_region(wesn) (!wesn || (wesn[XLO] == wesn[XHI] && wesn[YLO] == wesn[YHI]))

/* Get current setting for in/out columns */


char *ptrvoid (char ** p)	/* Handle as char ** just to determine if address is of a NULL pointer */
	{ return *p; }

int GMTAPI_Decode_ID (char *filename)
{	/* Checking if filename contains a name with embedded GMTAPI Object ID.
	 * If found we return the ID, otherwise we return GMT_NOTSET.
 	*/
	int object_ID = GMT_NOTSET;

	if (GMT_File_Is_Memory (filename)) {	/* Passing ID of an already registered object */
		if (sscanf (&filename[9], "%d", &object_ID) != 1) return (GMT_NOTSET);	/* Get the object ID unless we fail scanning */
	}
	return (object_ID);	/* Returns GMT_NOTSET if no embedded ID was found */
}

int GMT_Message (void *V_API, unsigned int mode, char *format, ...)
{	/* Message independent of verbosity, optionally with timestamp.
	 * mode = 0:	No time stamp
	 * mode = 1:	Abs time stamp formatted via GMT_TIME_STAMP
	 * mode = 2:	Report elapsed time since last reset.
	 * mode = 4:	Reset elapsed time to 0, no time stamp.
	 * mode = 6:	Reset elapsed time and report it as well.
	 */
	time_t toc_abs;
	clock_t toc, S;
	size_t source_info_len;
	unsigned int H, M, milli;
	char message[4*GMT_BUFSIZ] = {""}, stamp[GMT_LEN256] = {""};
	struct GMTAPI_CTRL *API = NULL;
	va_list args;

	if (V_API == NULL) return_error (V_API, GMT_NOT_A_SESSION);
	if (format == NULL) return GMT_PTR_IS_NULL;	/* Format cannot be NULL */
	API = gmt_get_api_ptr (V_API);
	if (mode) toc = clock ();
	if (mode & 4) API->GMT->current.time.tic = toc;

	switch (mode) {
		case 1:
			strftime (stamp, sizeof(stamp), API->GMT->current.setting.format_time_stamp, localtime (&toc_abs));
			break;
		case 2:
		case 6:
			S = (toc - (clock_t)API->GMT->current.time.tic);	/* Elapsed time in ticks */
			milli = (unsigned)(((float)S / CLOCKS_PER_SEC - (int)(S / CLOCKS_PER_SEC)) * CLOCKS_PER_SEC);	/* millisec */
			S /= CLOCKS_PER_SEC;	/* Elapsed time in seconds */
			H = urint (floor (S * GMT_SEC2HR));
			S -= H * GMT_HR2SEC_I;
			M = urint (floor (S * GMT_SEC2MIN));
			S -= M * GMT_MIN2SEC_I;
			sprintf (stamp, "Elapsed time %2.2d:%2.2d:%2.2d.%2.2d", H, M, (int)S, milli);
			break;
		default: break;
	}
	va_start (args, format);
	if (mode % 4) sprintf (message, "%s | ", stamp);	/* Lead with the time stamp */
	source_info_len = strlen (message);
	vsnprintf (message + source_info_len, 4*GMT_BUFSIZ - source_info_len, format, args);
	va_end (args);
	assert (strlen (message) < 4*GMT_BUFSIZ);
	API->print_func (API->GMT->session.std[GMT_ERR], message);	/* Do the printing */
	return (GMT_NOERROR);
}

int GMTAPI_report_error (void *V_API, int error)
{	/* Write error message to log or stderr, then return error code back.
 	 * All functions can call this, even if API has not been initialized. */
	FILE *fp = NULL;
	bool report;
	char message[GMT_BUFSIZ];
	struct GMTAPI_CTRL *API = gmt_get_api_ptr (V_API);

	report = (API) ? API->error != API->last_error : true;
	if (report && error != GMT_OK) {	/* Report error */
		if (!API || !API->GMT || (fp = API->GMT->session.std[GMT_ERR]) == NULL) fp = stderr;
		if (API && API->session_tag) {
			sprintf (message, "[Session %s (%d)]: Error returned from GMT API: %s (%d)\n",
				API->session_tag, API->session_ID, g_api_error_string[error], error);
			GMT_Message (API, GMT_TIME_NONE, message);
		}
		else
			fprintf (fp, "Error returned from GMT API: %s (%d)\n", g_api_error_string[error], error);
	}
	if (API) API->last_error = API->error, API->error = error;	/* Update API error value if API exists */
	return (error);
}

unsigned int gmtry (unsigned int geometry)
{
	/* Return index to text representation in GMT_geometry[] */
	if (geometry == GMT_IS_POINT)   return 1;
	if (geometry == GMT_IS_LINE)    return 2;
	if (geometry == GMT_IS_POLY)    return 3;
	if (geometry == GMT_IS_PLP)     return 4;
	if ((geometry & GMT_IS_LINE) && (geometry & GMT_IS_POLY)) return 5;
	if (geometry == GMT_IS_SURFACE) return 6;
	if (geometry == GMT_IS_NONE)    return 7;
	return 0;
}

bool GMTAPI_Validate_Geometry (struct GMTAPI_CTRL *API, int family, int geometry)
{	/* Sanity check that geometry and family are compatible; note they may be -1 hence int */
	bool problem = false;
	if (geometry == GMT_NOTSET || family == GMT_NOTSET) return false;	/* No errors if nothing to check */
	switch (family) {
		case GMT_IS_TEXTSET: if (!(geometry == GMT_IS_NONE || (geometry & GMT_IS_PLP))) problem = true; break;	/* Textsets can hold many things... */
		case GMT_IS_DATASET: if (!(geometry == GMT_IS_NONE || (geometry & GMT_IS_PLP))) problem = true; break;	/* Datasets can hold many things... */
		case GMT_IS_GRID:    if (geometry != GMT_IS_SURFACE) problem = true; break;	/* Only surface is valid */
		case GMT_IS_IMAGE:   if (geometry != GMT_IS_SURFACE) problem = true; break;	/* Only surface is valid */
		case GMT_IS_CPT:     if (geometry != GMT_IS_NONE) problem = true;    break;	/* Only text is valid */
		case GMT_IS_VECTOR:  if (geometry != GMT_IS_POINT) problem = true;   break;	/* Only point is valid */
		case GMT_IS_MATRIX:  if (geometry == GMT_IS_NONE) problem = true;    break;	/* Matrix can hold surfaces or TEXTSETs */
		case GMT_IS_COORD:   if (geometry != GMT_IS_NONE) problem = true;    break;	/* Only text is valid */
	}
	return (problem);
}

int GMTAPI_Validate_ID (struct GMTAPI_CTRL *API, int family, int object_ID, int direction)
{	/* Checks to see if the given object_ID is listed and of the right direction.  If so
 	 * we return the item number; otherwise return GMT_NOTSET and set API->error to the error code.
	 * Note: int arguments MAY be GMT_NOTSET, hence signed ints.  If object_ID == GMT_NOTSET
	 * then we only look for TEXTSETS or DATASETS .*/
	unsigned int i;
	int item, s_value;
	//printf("GMTAPI_Validate_ID: %s %s %d   object_ID :%d \n",__FILE__ ,__func__,__LINE__,object_ID);
	 /* Search for the object in the active list.  However, if object_ID == GMT_NOTSET we pick the first in that direction */

	for (i = 0, item = GMT_NOTSET; item == GMT_NOTSET && i < API->n_objects; i++) {
		if (!API->object[i]) continue;									/* Empty object */
		if (direction != GMT_NOTSET && API->object[i]->status != GMT_IS_UNUSED) continue;		/* Already used this object */
		if (!(family == GMT_NOTSET || (s_value = API->object[i]->family) == family)) continue;		/* Not the required data type */
		if (object_ID == GMT_NOTSET && (s_value = API->object[i]->direction) == direction) item = i;	/* Pick the first object with the specified direction */
		if (object_ID == GMT_NOTSET && !(API->object[i]->family == GMT_IS_DATASET || API->object[i]->family == GMT_IS_TEXTSET)) continue;	/* Must be data/text-set */
		else if (direction == GMT_NOTSET && (s_value = API->object[i]->ID) == object_ID) item = i;	/* Pick the requested object regardless of direction */
		else if ((s_value = API->object[i]->ID) == object_ID) item = i;					/* Pick the requested object */
	}
	if (item == GMT_NOTSET) { 
		//printf("GMTAPI_Validate_ID: %s %s %d\n",__FILE__ ,__func__,__LINE__);
		API->error = GMT_NOT_A_VALID_ID; return (GMT_NOTSET); }			/* No such object found */

	/* OK, we found the object; is it the right kind (input or output)? */
	if (direction != GMT_NOTSET && (s_value = API->object[item]->direction) != direction) {
		/* Passing an input object but it is listed as output, or vice versa */
		if (direction == GMT_IN)  { API->error = GMT_NOT_INPUT_OBJECT;  
			//printf("GMTAPI_Validate_ID: %s %s %d\n",__FILE__ ,__func__,__LINE__);			
			return (GMT_NOTSET); }
		if (direction == GMT_OUT) { API->error = GMT_NOT_OUTPUT_OBJECT; 
			//printf("GMTAPI_Validate_ID: %s %s %d\n",__FILE__ ,__func__,__LINE__);
			return (GMT_NOTSET); }
	}
	/* Here we have been successful in finding the right object */
	//printf("GMTAPI_Validate_ID: %s %s %d\n",__FILE__ ,__func__,__LINE__);
	return (item);
}
int GMTAPI_Memory_Registered (struct GMTAPI_CTRL *API, enum GMT_enum_family family, unsigned int direction, void *resource) {
	/* Determine if resource is a filename and that it has already been registered */
	int object_ID = 0, item;

	if (family == GMT_IS_COORD) return (GMT_NOTSET);	/* Coordinate arrays are never a registered memory resource */
	if ((object_ID = GMTAPI_Decode_ID ((char *)resource)) == GMT_NOTSET) return (GMT_NOTSET);	/* Not a registered resource */
	if ((item = GMTAPI_Validate_ID (API, family, object_ID, direction)) == GMT_NOTSET) return (GMT_NOTSET);	/* Not the right attributes */
	return (object_ID);	/* resource is a registered and valid item */
}

int GMTAPI_is_registered (struct GMTAPI_CTRL *API, enum GMT_enum_family family, unsigned int geometry, unsigned int direction, unsigned int mode, char *filename, void *resource)
{	/* Checks to see if the given data pointer has already been registered.
 	 * This can happen for grids which first gets registered reading the header
 	 * and then is registered again when reading the whole grid.  In those cases
	 * we dont want to register them twice.
	 */
	unsigned int i;
	int item;

	if (API->n_objects == 0) return (GMT_NOTSET);	/* There are no known resources yet */

	 /* Search for the object in the active list.  However, if object_ID == GMT_NOTSET we pick the first in that direction */

	for (i = 0, item = GMT_NOTSET; item == GMT_NOTSET && i < API->n_objects; i++) {
		if (!API->object[i]) continue;					/* Empty object */
		if (API->object[i]->status != GMT_IS_UNUSED) {	/* Has already been read - do we wish to reset this count ? */
			if (family == GMT_IS_GRID && (mode & GMT_GRID_DATA_ONLY)) {
				if (mode & GMT_GRID_IS_COMPLEX_MASK) {
					/* Check if complex grid already has one layer and we are reading the next one */
					struct GMT_GRID *G = gmt_get_grid_data (resource);	/* Get pointer to grid */
					unsigned int cmplx = mode & GMT_GRID_IS_COMPLEX_MASK;
					if (G->header->complex_mode & GMT_GRID_IS_COMPLEX_MASK && G->header->complex_mode != cmplx && filename) {
						/* Apparently so, either had real and now getting image or vice versa. */
						free ((void *)API->object[i]->filename);	/* Free previous grid name and replace with current name */
						API->object[i]->filename = strdup (filename);
						mode |= GMT_IO_RESET;	/* Reset so we may read in the 2nd component grid */
					}
				}
				else {	/* Just read the header earlier, do the reset */
					mode |= GMT_IO_RESET;	/* Reset so we may read in the grid data */
				}
			}

			if (!(mode & GMT_IO_RESET)) continue;	/* No reset so we refuse */
			API->object[i]->status = GMT_IS_UNUSED;	/* Reset so we may continue to read it */
		}
		if (API->object[i]->direction != direction) continue;		/* Wrong direction */
		if (API->object[i]->family != family) continue;			/* Wrong family */
		if (API->object[i]->geometry != geometry) continue;		/* Wrong geometry */
		if (resource && API->object[i]->resource == resource) item = API->object[i]->ID;	/* Yes: already registered */
		if (resource && API->object[i]->data == resource) item = API->object[i]->ID;		/* Yes: already registered */
	}
	return (item);		/* The ID of the object (or -1) */
}

struct GMTAPI_DATA_OBJECT * GMTAPI_Make_DataObject (struct GMTAPI_CTRL *API, enum GMT_enum_family family, unsigned int method, unsigned int geometry, void *resource, unsigned int direction)
{	/* Simply the creation and initialization of this DATA_OBJECT structure */
	struct GMTAPI_DATA_OBJECT *S_obj = (struct GMTAPI_DATA_OBJECT *) GMT_memory (API->GMT, NULL, 1, struct GMTAPI_DATA_OBJECT);

	S_obj->family = family;
	S_obj->method = method;
	S_obj->geometry = (enum GMT_enum_geometry)geometry;
	S_obj->resource = resource;
	S_obj->direction = (enum GMT_io_enum) direction;

	return (S_obj);
}

int GMTAPI_Add_Data_Object (struct GMTAPI_CTRL *API, struct GMTAPI_DATA_OBJECT *object)
{	/* Hook object to end of linked list and assign unique id (> 0) which is returned */

	/* Find the first entry in the API->object array which is unoccupied, and if
	 * they are all occupied then reallocate the array to make more space.
	 * We thus find and return the lowest available ID. */
	int object_ID;

	API->n_objects++;		/* Add one more entry to the tally */
	if (API->n_objects == API->n_objects_alloc) {	/* Must allocate more space for more data descriptors */
		size_t old_n_alloc = API->n_objects_alloc;
		API->n_objects_alloc += GMT_SMALL_CHUNK;
		API->object = (struct GMTAPI_DATA_OBJECT **) GMT_memory (API->GMT, API->object, API->n_objects_alloc, struct GMTAPI_DATA_OBJECT *);
		GMT_memset (&(API->object[old_n_alloc]), API->n_objects_alloc - old_n_alloc, struct GMTAPI_DATA_OBJECT *);	/* Set to NULL */
		if (!(API->object)) {		/* Failed to allocate more memory */
			API->n_objects--;	/* Undo our premature increment */
			return_value (API, GMT_MEMORY_ERROR, GMT_NOTSET);
			//printf("GMT_MEMORY_ERROR, GMT_NOTSET");
		}
	}
	object_ID = object->ID = API->unique_ID++;	/* Assign a unique object ID */
	API->object[API->n_objects-1] = object;		/* Hook the current object onto the end of the list */

	return (object_ID);
}

void GMT_str_toupper (char *value)
{
	/* Convert entire string to upper case */
	int i, c;
	for (i = 0; value[i]; i++) {
		c = (int)value[i];
		value[i] = (char) toupper (c);
	}
}

unsigned int GMTAPI_Add_Existing (struct GMTAPI_CTRL *API, enum GMT_enum_family family, unsigned int geometry, unsigned int direction, int *first_ID)
{	/* In this mode, we find all registrered resources of matching family,geometry,direction that are unused and turn select to true. */
	unsigned int i, n;

	*first_ID = GMT_NOTSET;	/* Not found yet */
	for (i = n = 0; i < API->n_objects; i++) {
		if (!API->object[i]) continue;				/* A freed object, skip */
		if (API->object[i]->direction != direction) continue;	/* Wrong direction */
		if (API->object[i]->geometry != geometry) continue;	/* Wrong geometry */
		if (API->object[i]->status != GMT_IS_UNUSED) continue;	/* Already used */
		if (family != API->object[i]->family) continue;		/* Wrong data type */
		n++;
		if (*first_ID == GMT_NOTSET) *first_ID = API->object[i]->ID;
		API->object[i]->selected = true;
	}
	return (n);
}

bool GMTAPI_Not_Used (struct GMTAPI_CTRL *API, char *name)
{
	/* See if this file has already been registered and used.  If so, do not add it again */
	unsigned int item = 0;
	bool not_used = true;
	while (item < API->n_objects && not_used) {
		if (API->object[item] && API->object[item]->direction == GMT_IN && API->object[item]->status != GMT_IS_UNUSED && API->object[item]->filename && !strcmp (API->object[item]->filename, name))
			/* Used resource with same name */
			not_used = false;	/* Got item with same name, but used */
		else
			item++;	/* No, keep looking */
	}
	return (not_used);
}

int GMTAPI_Init_Import (struct GMTAPI_CTRL *API, enum GMT_enum_family family, unsigned int geometry, unsigned int mode, struct GMT_OPTION *head)
{	/* Handle registration of data files given with option arguments and/or stdin as input sources.
	 * These are the possible actions taken:
	 * 1. If (mode | GMT_ADD_FILES_IF_NONE) is true and NO resources have previously been registered, then we scan the option list for files (option == '<' (input)).
	 *    For each file found we register the item as a resource.
	 * 2. If (mode | GMT_ADD_FILES_ALWAYS) is true then we always scan the option list for files (option == '<' (input)).
	 *    For each file found we register the item as a resource.
	 * 3. If (mode & GMT_ADD_STDIO_IF_NONE) is true we will register stdin as an input source only if there are NO input items registered.
	 * 4. If (mode & GMT_ADD_STDIO_ALWAYS) is true we will register stdin as an input source, regardless of other items already registered.
	 */

	int object_ID, first_ID = GMT_NOTSET, item;//nishita didn't do that
 	unsigned int n_reg = 0;
	struct GMT_OPTION *current = NULL;
	double *wesn = NULL;
	//printf("*******************> %s %s %d\n",__FILE__ ,__func__,__LINE__);
	//printf ( "GMTAPI_Init_Import: Passed family = %s and geometry = %s\n", GMT_family[family], GMT_geometry[gmtry(geometry)]);

	if (mode & GMT_ADD_EXISTING) {
		n_reg = GMTAPI_Add_Existing (API, family, geometry, GMT_IN, &first_ID);
	}
	//printf("*******************> %s %s %d first_ID::: %d \n",__FILE__ ,__func__,__LINE__,first_ID);
	if ((mode & GMT_ADD_FILES_ALWAYS) || ((mode & GMT_ADD_FILES_IF_NONE))) {	/* Wish to register all input file args as sources */
		//printf("*******************> %s %s %d\n",__FILE__ ,__func__,__LINE__);
		current = head;
		while (current) {		/* Loop over the list and look for input files */
			//printf("*******************> %s %s %d\n",__FILE__ ,__func__,__LINE__);
			if (current->option == GMT_OPT_INFILE && GMTAPI_Not_Used (API, current->arg)) {	/* File given, register it if not already used */
				//printf("*******************> %s %s %d\n",__FILE__ ,__func__,__LINE__);
				if (geometry == GMT_IS_SURFACE) {	/* Grids may require a subset */
					if (API->GMT->common.R.active) {	/* Global subset may have been specified (it might also match the grids domain) */
						//printf("*******************> %s %s %d\n",__FILE__ ,__func__,__LINE__);
						wesn = (double *)GMT_memory (API->GMT, NULL, 4, double);
						GMT_memcpy (wesn, API->GMT->common.R.wesn, 4, double);
					}
				}
				if ((object_ID = GMT_Register_IO (API, family, GMT_IS_FILE, geometry, GMT_IN, wesn, current->arg)) == GMT_NOTSET)
		{
					//return_value (API, API->error, GMT_NOTSET);	/* Failure to register */
					//printf("*******************> %s %s %d\n",__FILE__ ,__func__,__LINE__);
					return (GMT_NOTSET);// nishita
		}
				//printf("*******************> %s %s %d\n",__FILE__ ,__func__,__LINE__);
				n_reg++;	/* Count of new items registered */
				if (wesn) GMT_free (API->GMT, wesn);
				if (first_ID == GMT_NOTSET) first_ID = object_ID;	/* Found our first ID */
				if ((item = GMTAPI_Validate_ID (API, family, object_ID, GMT_IN)) == GMT_NOTSET)
				{
					//return_value (API, API->error, GMT_NOTSET);	/* Some internal error... */
					//printf("*******************> %s %s %d object_ID %d \n",__FILE__ ,__func__,__LINE__,object_ID);
					return (GMT_NOTSET); //nishita
				}
				API->object[item]->selected = true;
			}
			current = current->next;	/* Go to next option */
		}
		//printf ( "GMTAPI_Init_Import: Added %d new sources\n", n_reg);
	}

	/* Note that n_reg can have changed if we added file args above */

	if ((mode & GMT_ADD_STDIO_ALWAYS) || ((mode & GMT_ADD_STDIO_IF_NONE) && n_reg == 0)) {	/* Wish to register stdin pointer as a source */
		if ((object_ID = GMT_Register_IO (API, family, GMT_IS_STREAM, geometry, GMT_IN, NULL, API->GMT->session.std[GMT_IN])) == GMT_NOTSET){
			//return_value (API, API->error, GMT_NOTSET);	/* Failure to register stdin */
			//printf("*******************> %s %s %d\n",__FILE__ ,__func__,__LINE__);
			return (GMT_NOTSET);//nishita
		}
		n_reg++;		/* Add the single item */
		if (first_ID == GMT_NOTSET) first_ID = object_ID;	/* Found our first ID */
		if ((item = GMTAPI_Validate_ID (API, family, object_ID, GMT_IN)) == GMT_NOTSET)
		{
			//return_value (API, API->error, GMT_NOTSET);	/* Some internal error... */
			//printf("*******************> %s %s %d\n",__FILE__ ,__func__,__LINE__);
			return (GMT_NOTSET); //nishita
		}
		API->object[item]->selected = true;
		printf ( "GMTAPI_Init_Import: Added stdin to registered sources\n");
	}
	//printf("*******************> %s %s %d\n",__FILE__ ,__func__,__LINE__);
	return (first_ID);
}

int GMTAPI_Init_Export (struct GMTAPI_CTRL *API, enum GMT_enum_family family, unsigned int geometry, unsigned int mode, struct GMT_OPTION *head)
{	/* Handle registration of output file given with option arguments and/or stdout as output destinations.
	 * Only a single output may be considered.  These are the possible actions taken:
	 * 1. If (mode | GMT_ADD_FILES_IF_NONE) is true and NO destinations have previously been registered,
	 *    then we scan the option list for files (option == '>' (output)).
	 *    Only one file can be registered as a destination; finding more than one results in an error.
	 * 2. If (mode | GMT_ADD_FILES_ALWAYS) is true then we always scan the option list for files (option == '>' (output)).
	 *    Only one file can be registered as a destination; finding more than one results in an error.
	 * 3. If (mode & GMT_ADD_STDIO_IF_NONE) is true we will register stdout as the only destination if there is NO output item registered.
	 * 4. If (mode & GMT_ADD_STDIO_ALWAYS) is true we will register stdout as an destination,
	 *    and give error if other output items have already been registered.
	 */

	unsigned int n_reg = 0;
	int object_ID, item;
	struct GMT_OPTION *current = NULL;

	//printf ("GMTAPI_Init_Export: Passed family = %s and geometry = %s\n", GMT_family[family], GMT_geometry[gmtry(geometry)]);

	if (mode & GMT_ADD_EXISTING) {
		n_reg = GMTAPI_Add_Existing (API, family, geometry, GMT_OUT, &object_ID);
	}
	if (n_reg > 1)
		//return_value (API, GMT_ONLY_ONE_ALLOWED, GMT_NOTSET);	/* Only one output allowed */
		return (GMT_NOTSET); //nishita

	if ((mode & GMT_ADD_FILES_ALWAYS) || (mode & GMT_ADD_FILES_IF_NONE)) {	/* Wish to register a single output file arg as destination */
		current = head;
		while (current) {		/* Loop over the list and look for input files */
			if (current->option == GMT_OPT_OUTFILE) n_reg++;	/* File given, count it */
			current = current->next;				/* Go to next option */
		}
		if (n_reg > 1)
			//return_value (API, GMT_ONLY_ONE_ALLOWED, GMT_NOTSET);	/* Only one output allowed */
			return (GMT_NOTSET); //nishita

		if (n_reg == 1) {	/* Register the single output file found above */
			current = head;
			while (current) {		/* Loop over the list and look for output files (we know there is only one) */
				if (current->option == GMT_OPT_OUTFILE) {	/* File given, register it */
					if ((object_ID = GMT_Register_IO (API, family, GMT_IS_FILE, geometry, GMT_OUT, NULL, current->arg)) == GMT_NOTSET)
						//return_value (API, API->error, GMT_NOTSET);	/* Failure to register */
						return (GMT_NOTSET);//nishita
					if ((item = GMTAPI_Validate_ID (API, family, object_ID, GMT_OUT)) == GMT_NOTSET)
						//return_value (API, API->error, GMT_NOTSET);	/* Some internal error... */
						return(GMT_NOTSET); //nishita
					API->object[item]->selected = true;
					printf ("GMTAPI_Init_Export: Added 1 new destination\n");
				}
				current = current->next;	/* Go to next option */
			}
		}
	}
	/* Note that n_reg can have changed if we added file arg */

	if ((mode & GMT_ADD_STDIO_ALWAYS) && n_reg == 1)
		//return_value (API, GMT_ONLY_ONE_ALLOWED, GMT_NOTSET);	/* Only one output destination allowed at once */
		return (GMT_NOTSET); //nishita
	if (n_reg == 0 && ((mode & GMT_ADD_STDIO_ALWAYS) || (mode & GMT_ADD_STDIO_IF_NONE))) {	/* Wish to register stdout pointer as a destination */
		if ((object_ID = GMT_Register_IO (API, family, GMT_IS_STREAM, geometry, GMT_OUT, NULL, API->GMT->session.std[GMT_OUT])) == GMT_NOTSET)
			//return_value (API, API->error, GMT_NOTSET);	/* Failure to register stdout?*/
			return (GMT_NOTSET); //nishita
		if ((item = GMTAPI_Validate_ID (API, family, object_ID, GMT_OUT)) == GMT_NOTSET)
			//return_value (API, API->error, GMT_NOTSET);	/* Some internal error... */
			return (GMT_NOTSET); //nishita
		API->object[item]->selected = true;
		printf ( "GMTAPI_Init_Export: Added stdout to registered destinations\n");
		n_reg = 1;	/* Only have one item */
	}
	if (n_reg == 0)
		//return_value (API, GMT_OUTPUT_NOT_SET, GMT_NOTSET);	/* No output set */
		return (GMT_NOTSET); //nishita
	return (object_ID);
}

int GMTAPI_Next_IO_Source (struct GMTAPI_CTRL *API, unsigned int direction)
{	/* Get ready for the next source/destination (open file, initialize counters, etc.).
	 * Note this is only a mechanism for dataset and textset files where it is common
	 * to give many files on the command line (e.g., *.txt) and we do rec-by-rec processing. */
	int *fd = NULL;	/* !!! Must be int* due to nature of Unix system function */
	int error = 0, kind, via = 0;
	static const char *dir[2] = {"from", "to"};
	static const char *operation[2] = {"Reading", "Writing"};
	char *mode = NULL;
	struct GMT_MATRIX *M_obj = NULL;
	struct GMT_VECTOR *V_obj = NULL;
	struct GMTAPI_DATA_OBJECT *S_obj = NULL;

	S_obj = API->object[API->current_item[direction]];		/* For shorthand purposes only */
	//printf(/*API, GMT_MSG_DEBUG,*/ "GMTAPI_Next_IO_Source: Selected object %d\n", S_obj->ID);
	mode = (direction == GMT_IN) ? API->GMT->current.io.r_mode : API->GMT->current.io.w_mode;	/* Reading or writing */
	S_obj->close_file = false;		/* Do not want to close file pointers passed to us unless WE open them below */
	/* Either use binary n_columns settings or initialize to unknown, i.e., GMT_MAX_COLUMNS */
	S_obj->n_expected_fields = (API->GMT->common.b.ncol[direction]) ? API->GMT->common.b.ncol[direction] : GMT_MAX_COLUMNS;
	GMT_memset (API->GMT->current.io.curr_pos[direction], 3U, uint64_t);	/* Reset file, seg, point counters */
	if (S_obj->method >= GMT_VIA_VECTOR) via = (S_obj->method / GMT_VIA_VECTOR) - 1;
	
	switch (S_obj->method) {	/* File, array, stream etc ? */
		case GMT_IS_FILE:	/* Filename given; we must open file ourselves */
			if (S_obj->family == GMT_IS_GRID) return (/*GMTAPI_report_error (API,*/ GMT_NOT_A_VALID_TYPE/*)*/);	/* Grids not allowed here */
			if ((S_obj->fp = GMT_fopen (API->GMT, S_obj->filename, mode)) == NULL) {	/* Trouble opening file */
				printf (/*API, GMT_MSG_NORMAL, */ "Unable to open file %s for %s\n", S_obj->filename, GMT_direction[direction]);
				return (GMT_ERROR_ON_FOPEN);
			}
			S_obj->close_file = true;	/* We do want to close files we are opening, but later */
			strncpy (API->GMT->current.io.current_filename[direction], S_obj->filename, GMT_BUFSIZ);
			printf (/*API, GMT_MSG_LONG_VERBOSE,*/ "%s %s %s file %s\n",
				operation[direction], GMT_family[S_obj->family], dir[direction], S_obj->filename);
			if (GMT_binary_header (API->GMT, direction)) {
				GMT_io_binary_header (API->GMT, S_obj->fp, direction);
				printf (/*API, GMT_MSG_NORMAL, */ "%s %d bytes of header %s binary file %s\n",
					operation[direction], API->GMT->current.setting.io_n_header_items, dir[direction], S_obj->filename);
			}
			break;

		case GMT_IS_STREAM:	/* Given a stream; no need to open (or close) anything */
//#ifdef SET_IO_MODE
	//		if (S_obj->family == GMT_IS_DATASET && S_obj->fp == API->GMT->session.std[direction])
		//		GMT_setmode (API->GMT, (int)direction);	/* Windows may need to have its read mode changed from text to binary */
//#endif
			kind = (S_obj->fp == API->GMT->session.std[direction]) ? 0 : 1;	/* 0 if stdin/out, 1 otherwise for user pointer */
			sprintf (API->GMT->current.io.current_filename[direction], "<%s %s>", GMT_stream[kind], GMT_direction[direction]);
			printf (/*API, GMT_MSG_LONG_VERBOSE,*/ "%s %s %s %s %s stream\n",
				operation[direction], GMT_family[S_obj->family], dir[direction], GMT_stream[kind], GMT_direction[direction]);
			if (GMT_binary_header (API->GMT, direction)) {
				GMT_io_binary_header (API->GMT, S_obj->fp, direction);
				printf (/*API, GMT_MSG_NORMAL, */ "%s %d bytes of header %s binary %s stream\n",
					operation[direction], API->GMT->current.setting.io_n_header_items, dir[direction], GMT_stream[kind]);
			}
			break;

		case GMT_IS_FDESC:	/* Given a pointer to a file handle; otherwise same as stream */
			fd = (int *)S_obj->fp;
			if ((S_obj->fp = GMT_fdopen (*fd, mode)) == NULL) {	/* Reopen handle as stream */
				printf (/*API, GMT_MSG_NORMAL, */ "Unable to open file descriptor %d for %s\n", *fd, GMT_direction[direction]);
				return (GMT_ERROR_ON_FDOPEN);
			}
			kind = (S_obj->fp == API->GMT->session.std[direction]) ? 0 : 1;	/* 0 if stdin/out, 1 otherwise for user pointer */
			sprintf (API->GMT->current.io.current_filename[direction], "<%s %s>", GMT_stream[kind], GMT_direction[direction]);
			printf (/*API, GMT_MSG_LONG_VERBOSE,*/ "%s %s %s %s %s stream via supplied file descriptor\n",
				operation[direction], GMT_family[S_obj->family], dir[direction], GMT_stream[kind], GMT_direction[direction]);
			if (GMT_binary_header (API->GMT, direction)) {
				GMT_io_binary_header (API->GMT, S_obj->fp, direction);
				printf (/*API, GMT_MSG_NORMAL, */ "%s %d bytes of header %s binary %s stream via supplied file descriptor\n",
					operation[direction], API->GMT->current.setting.io_n_header_items, dir[direction], GMT_stream[kind]);
			}
			break;

	 	case GMT_IS_DUPLICATE:	/* Copy, nothing to do [PW: not tested] */
			printf (/*API, GMT_MSG_LONG_VERBOSE,*/ "%s %s %s memory copy supplied by pointer\n",
				operation[direction], GMT_family[S_obj->family], dir[direction]);
			break;

	 	case GMT_IS_REFERENCE:	/* Reference, nothing to do [PW: not tested] */
			printf (/*API, GMT_MSG_LONG_VERBOSE,*/ "%s %s %s memory reference supplied by pointer\n",
				operation[direction], GMT_family[S_obj->family], dir[direction]);
			break;

	 	case GMT_IS_DUPLICATE + GMT_VIA_MATRIX:	/* This means reading or writing a dataset record-by-record via a user matrix [PW: not tested] */
		case GMT_IS_REFERENCE + GMT_VIA_MATRIX:
			if (S_obj->family != GMT_IS_DATASET) return (/*GMTAPI_report_error (API,*/ GMT_NOT_A_VALID_TYPE/*)*/);
			printf (/*API, GMT_MSG_LONG_VERBOSE,*/ "%s %s %s %s memory location via %s\n",
				operation[direction], GMT_family[S_obj->family], dir[direction], GMT_direction[direction], GMT_via[via]);
			if (direction == GMT_IN) {	/* Hard-wired limit as pass in from calling program */
				M_obj = (struct GMT_MATRIX *) S_obj->resource;
				S_obj->n_rows = M_obj->n_rows;
				S_obj->n_columns = M_obj->n_columns;
				API->GMT->common.b.ncol[direction] = M_obj->n_columns;	/* Basically, we are doing what GMT calls binary i/o */
			}
			API->GMT->common.b.active[direction] = true;
			strcpy (API->GMT->current.io.current_filename[direction], "<memory>");
			break;

		 case GMT_IS_DUPLICATE + GMT_VIA_VECTOR:	/* These 2 means reading a dataset record-by-record via user vector arrays [PW: not tested] */
		 case GMT_IS_REFERENCE + GMT_VIA_VECTOR:
			if (S_obj->family != GMT_IS_DATASET) return (/*GMTAPI_report_error (API,*/  GMT_NOT_A_VALID_TYPE/*)*/);
			printf (/*API, GMT_MSG_LONG_VERBOSE,*/ "%s %s %s %s memory location via %s\n",
					operation[direction], GMT_family[S_obj->family], dir[direction], GMT_direction[direction], GMT_via[via]);
			V_obj = S_obj->resource;
			if (direction == GMT_OUT && V_obj->alloc_mode == GMT_ALLOCATED_BY_GMT) {	/* Must allocate output space */
				S_obj->n_alloc = GMT_CHUNK;
				/* S_obj->n_rows is 0 which means we are allocating more space as we need it */
				V_obj->n_rows = S_obj->n_alloc;
				if ((error = gmt_alloc_vectors (API->GMT, V_obj)) != GMT_OK) return (/*GMTAPI_report_error (API,*/ error/*)*/);
			}
			else
				S_obj->n_rows = V_obj->n_rows;	/* Hard-wired limit as passed in from calling program */
			S_obj->n_columns = V_obj->n_columns;
			API->GMT->common.b.ncol[direction] = V_obj->n_columns;	/* Basically, we are doing what GMT calls binary i/o */
			API->GMT->common.b.active[direction] = true;
			strcpy (API->GMT->current.io.current_filename[direction], "<memory>");
			break;

		default:
			printf (/*API, GMT_MSG_NORMAL, */ "GMTAPI: Internal error: GMTAPI_Next_IO_Source called with illegal method\n");
			break;
	}

	/* A few things pertaining only to data/text tables */
	API->GMT->current.io.rec_in_tbl_no = 0;	/* Start on new table */
	S_obj->import = (S_obj->family == GMT_IS_TEXTSET) ? &GMT_ascii_textinput : API->GMT->current.io.input;	/* The latter may point to ascii or binary input functions */

	return (GMT_OK);
}

int GMTAPI_Next_Data_Object (struct GMTAPI_CTRL *API, enum GMT_enum_family family, unsigned int direction)
{	/* Sets up current_item to be the next unused item of the required direction; or return EOF.
	 * When EOF is returned, API->current_item[direction] holds the last object ID used. */
	bool found = false;
	unsigned int item;

	item = API->current_item[direction] + 1;	/* Advance to next item, if possible */
	while (item < API->n_objects && !found) {
		if (API->object[item] && API->object[item]->selected && API->object[item]->status == GMT_IS_UNUSED && API->object[item]->direction == direction && family == API->object[item]->family)
			found = true;	/* Got item that is selected and unused, has correct direction and family */
		else
			item++;	/* No, keep looking */
	}
	if (found) {	/* Update to use next item */
		API->current_item[direction] = item;	/* The next item */
		return (GMTAPI_Next_IO_Source (API, direction));	/* Initialize the next source/destination */
	}
	else
		return (EOF);	/* No more objects available for this direction; return EOF */
}


int GMT_Put_Record (void *V_API, unsigned int mode, void *record)
{	/* Writes a single data record to destimation.
	 * We use mode to signal the kind of record:
	 *   GMT_WRITE_TABLE_HEADER: Write an ASCII table header
	 *   GMT_WRITE_SEGMENT_HEADER: Write an ASCII or binary segment header
	 *   GMT_WRITE_DOUBLE:    Write an ASCII or binary data record
	 *   GMT_WRITE_TEXT:      Write an ASCII data record
	 * For text: If record == NULL use internal current record or header.
	 * Returns 0 if a record was written successfully (See what -s[r] can do).
	 * If an error occurs we return -1 and set API->error.
	 */
	int error = 0;
	uint64_t *p = NULL, col, ij;
	char *s = NULL;
	double *d = NULL;
	p_func_size_t GMT_2D_to_index = NULL;
	struct GMTAPI_DATA_OBJECT *S_obj = NULL;
	struct GMT_MATRIX *M_obj = NULL;
	struct GMT_VECTOR *V_obj = NULL;
	struct GMTAPI_CTRL *API = NULL;
	//myfunc3();

	//double * out = (double *)record;
	//printf("getting 11111111111111 1: %f 2: %f, 3:%f, 4:%f \n", out[0], out[1], out[2], out[3]);

	
	if (V_API == NULL) return_error (V_API, GMT_NOT_A_SESSION);
	API = gmt_get_api_ptr (V_API);
	if (!API->io_enabled[GMT_OUT]) return_error (API, GMT_ACCESS_NOT_ENABLED);

	S_obj = API->object[API->current_item[GMT_OUT]];	/* Shorthand for the data source we are working on */
	if (S_obj->status == GMT_IS_USED) return_error (API, GMT_WRITTEN_ONCE);	/* Only allow writing of a data set once [unless we reset status] */
	//printf("nishita __________ method %d mode%d \n ",S_obj->method,mode);
	switch (S_obj->method) {	/* File, array, stream etc ? */
		case GMT_IS_FILE:
	 	case GMT_IS_STREAM:
	 	case GMT_IS_FDESC:
			switch (mode) {
				case GMT_WRITE_DOUBLE:		/* Export either a formatted ASCII data record or a binary record */
					if (API->GMT->common.b.ncol[GMT_OUT] == UINT_MAX) API->GMT->common.b.ncol[GMT_OUT] = API->GMT->common.b.ncol[GMT_IN];
					error = API->GMT->current.io.output (API->GMT, S_obj->fp, API->GMT->common.b.ncol[GMT_OUT], record);
					break;
				case GMT_WRITE_TABLE_START:	/* Write title and command to start of file; skip if binary */
					GMT_write_newheaders (API->GMT, S_obj->fp, S_obj->n_columns);	error = 1;	/* Write one item */
					break;
				default:
					//GMT_Report (API, GMT_MSG_NORMAL, "GMTAPI: Internal error: GMT_Put_Record called with illegal mode %u\n", mode);
					return_error (API, GMT_NOT_A_VALID_IO_MODE);
					break;
			}
			break;
		default:
//			GMT_Report (API, GMT_MSG_NORMAL, "GMTAPI: Internal error: GMT_Put_Record called with illegal method\n");
			return_error (API, GMT_NOT_A_VALID_METHOD);
			break;
	}

	if (!error && (mode == GMT_WRITE_DOUBLE || mode == GMT_WRITE_TEXT)) API->current_rec[GMT_OUT]++;	/* Only increment if we placed a data record on the output */

	if (S_obj->n_alloc && API->current_rec[GMT_OUT] == S_obj->n_alloc) {	/* Must allocate more memory */
		size_t size;
		S_obj->n_alloc += GMT_CHUNK;
		size = S_obj->n_alloc;
		if (S_obj->method == (GMT_IS_DUPLICATE + GMT_VIA_MATRIX)) {
			size *= API->GMT->common.b.ncol[GMT_OUT];
			if ((error = GMT_alloc_univector (API->GMT, &(M_obj->data), M_obj->type, size)) != GMT_OK) return (error);
		}
		else {
			V_obj->n_rows = size;
			if ((error = gmt_alloc_vectors (API->GMT, V_obj)) != GMT_OK) return (error);
		}
	}
	S_obj->status = GMT_IS_USING;	/* Have started writing to this destination */

	return ((error) ? -1 : 0);
}


//check end


int GMT_Begin_IO (void *V_API, unsigned int family, unsigned int direction, unsigned int header)
{
	/* Initializes the rec-by-rec i/o mechanism for either input or output (given by direction).
	 * GMT_Begin_IO must be called before any data i/o is allowed.
	 * family:	The kind of data must be GMT_IS_DATASET or TEXTSET.
	 * direction:	Either GMT_IN or GMT_OUT.
	 * header:	Either GMT_HEADER_ON|OFF, controls the writing of the table start header info block
	 * Returns:	false if successfull, true if error.
	 */
	int error;
	struct GMTAPI_CTRL *API = NULL;

	if (V_API == NULL)
		//return_error (V_API, GMT_NOT_A_SESSION);
		return (GMT_NOT_A_SESSION); //nishita
	if (!(direction == GMT_IN || direction == GMT_OUT))
		//return_error (V_API, GMT_NOT_A_VALID_DIRECTION);
		return (GMT_NOT_A_VALID_DIRECTION); //nishita
	if (!(family == GMT_IS_DATASET || family == GMT_IS_TEXTSET))
		//return_error (V_API, GMT_NOT_A_VALID_IO_ACCESS);
		return (GMT_NOT_A_VALID_IO_ACCESS); //nishita
	API = gmt_get_api_ptr ((struct GMTAPI_CTRL *)V_API);
	if (!API->registered[direction]) printf(/*API, GMT_MSG_DEBUG,*/ "GMT_Begin_IO: Warning: No %s resources registered\n", GMT_direction[direction]);

	/* Must initialize record-by-record machinery for dataset or textset */
	//printf(/*API, GMT_MSG_DEBUG,*/ "GMT_Begin_IO: Initialize record-by-record access for %s\n", GMT_direction[direction]);
	API->current_item[direction] = -1;	/* GMTAPI_Next_Data_Object (below) will wind it to the first item >= 0 */
	if ((error = GMTAPI_Next_Data_Object (API, (enum GMT_enum_family)family, direction)))
		//return_error (API, GMT_NO_RESOURCES);	/* Something went bad */
		return (GMT_NO_RESOURCES);	//nshita
	API->io_mode[direction] = GMT_BY_REC;
	API->io_enabled[direction] = true;	/* OK to access resources */
	API->GMT->current.io.need_previous = (API->GMT->common.g.active || API->GMT->current.io.skip_duplicates);
	API->GMT->current.io.ogr = GMT_OGR_UNKNOWN;
	API->GMT->current.io.segment_header[0] = API->GMT->current.io.current_record[0] = 0;
	//printf(/*API, GMT_MSG_DEBUG,*/ "GMT_Begin_IO: %s resource access is now enabled [record-by-record]\n", GMT_direction[direction]);
	if (direction == GMT_OUT && header == GMT_HEADER_ON && !API->GMT->common.b.active[GMT_OUT])
		GMT_Put_Record (API, GMT_WRITE_TABLE_START, NULL);	/* Write standard ascii header block */

	return (GMT_OK);	/* No error encountered */
}

int GMT_Init_IO (void *V_API, unsigned int family, unsigned int geometry, unsigned int direction, unsigned int mode, unsigned int n_args, void *args)
{
	/* Registers program option file arguments as sources/destinations for the current module.
	 * All modules planning to use std* and/or command-line file args must call GMT_Init_IO to register these resources.
	 * family:	The kind of data (GMT_IS_DATASET|TEXTSET|CPT|GRID)
	 * geometry:	Either GMT_IS_NONE|POINT|LINE|POLYGON|SURFACE
	 * direction:	Either GMT_IN or GMT_OUT
	 * mode:	Bitflags composed of 1 = add command line (option) files, 2 = add std* if no other input/output,
	 *		4 = add std* regardless.  mode must be > 0.
	 * n_args:	Either 0 if we pass linked option structs or argc if we pass argv[]
	 * args:	Either linked list of program option arguments (n_args == 0) or char *argv[].
	 *
	 * Returns:	false if successfull, true if error.
	 */
	int object_ID;	/* ID of first object [only for debug purposes - not used in this function; ignore -Wunused-but-set-variable warning */
	struct GMT_OPTION *head = NULL;
	struct GMTAPI_CTRL *API = NULL;

	if (V_API == NULL)
		//return_error (V_API, GMT_NOT_A_SESSION);
		return (GMT_NOT_A_SESSION); //nishita
	API = gmt_get_api_ptr ((struct GMTAPI_CTRL *)(struct GMTAPI_CTRL *)V_API);
	if (GMTAPI_Validate_Geometry (API, family, geometry))
		//return_error (API, GMT_BAD_GEOMETRY);
		return (GMT_BAD_GEOMETRY); //nishita
	if (!(direction == GMT_IN || direction == GMT_OUT))
		//return_error (API, GMT_NOT_A_VALID_DIRECTION);
		return (GMT_NOT_A_VALID_DIRECTION);//nishita
	if (!((mode & GMT_ADD_FILES_IF_NONE) || (mode & GMT_ADD_FILES_ALWAYS) || (mode & GMT_ADD_STDIO_IF_NONE) || (mode & GMT_ADD_STDIO_ALWAYS) || (mode & GMT_ADD_EXISTING)))
		//return_error (API, GMT_NOT_A_VALID_MODE);
		return (GMT_NOT_A_VALID_MODE); //nishita

	if (n_args == 0) /* Passed the head of linked option structures */
		{
		
		head =  args;
		}
	else		/* Passed argc, argv, likely from Fortran */
		{
		
		head = GMT_Create_Options (API, n_args, args);
		}
	//GMT_io_banner (API->GMT, direction);	/* Message for binary i/o */
	if (direction == GMT_IN)
		object_ID = GMTAPI_Init_Import (API, (enum GMT_enum_family)family, geometry, mode, head);
	else
		object_ID = GMTAPI_Init_Export (API, (enum GMT_enum_family) family, geometry, mode, head);
	//printf ("GMT_Init_IO: Returned first %s object ID = %d\n", GMT_direction[direction], object_ID);
	return (API->error);
}

/* Mapping of internal [row][col] indices to a single 1-D index.
 * Internally, row and col starts at 0.  These will be accessed
 * via pointers to these functions, hence they are not macros.
 */

size_t GMTAPI_2D_to_index_C_normal (uint64_t row, uint64_t col, size_t dim)
{	/* Maps (row,col) to 1-D index for C normal grid */
	return (((size_t)row * dim) + (size_t)col);	/* Normal grid */
}

size_t GMTAPI_2D_to_index_C_cplx_real (uint64_t row, uint64_t col, size_t dim)
{	/* Maps (row,col) to 1-D index for C complex grid, real component */
	return (2*((size_t)row * dim) + (size_t)col);	/* Complex grid, real(1) component */
}

size_t GMTAPI_2D_to_index_C_cplx_imag (uint64_t row, uint64_t col, size_t dim)
{	/* Maps (row,col) to 1-D index for C complex grid, imaginary component */
	return (2*((size_t)row * dim) + (size_t)col + 1ULL);	/* Complex grid, imag(2) component */
}

size_t GMTAPI_2D_to_index_F_normal (uint64_t row, uint64_t col, size_t dim)
{	/* Maps (row,col) to 1-D index for Fortran */
	return (((size_t)col * dim) + (size_t)row);
}

size_t GMTAPI_2D_to_index_F_cplx_real (uint64_t row, uint64_t col, size_t dim)
{	/* Maps (row,col) to 1-D index for Fortran complex grid, real component */
	return (2*((size_t)col * dim) + (size_t)row);	/* Complex grid, real(1) */
}

size_t GMTAPI_2D_to_index_F_cplx_imag (uint64_t row, uint64_t col, size_t dim)
{	/* Maps (row,col) to 1-D index for Fortran complex grid, imaginary component  */
	return (2*((size_t)col * dim) + (size_t)row + 1ULL);	/* Complex grid, imag(2) component */
}



p_func_size_t GMTAPI_get_2D_to_index (struct GMTAPI_CTRL *API, enum GMT_enum_fmt shape, unsigned int mode)
{
	/* Return pointer to the required 2D-index function above.  Here
	 * shape is either GMT_IS_ROW_FORMAT (C) or GMT_IS_COL_FORMAT (FORTRAN);
	 * mode is either 0 (regular grid), GMT_GRID_IS_COMPLEX_REAL (complex real) or GMT_GRID_IS_COMPLEX_IMAG (complex imag)
	 */
	p_func_size_t p = NULL;

	switch (mode & GMT_GRID_IS_COMPLEX_MASK) {
		case GMT_GRID_IS_REAL:
			p = (shape == GMT_IS_ROW_FORMAT) ? GMTAPI_2D_to_index_C_normal : GMTAPI_2D_to_index_F_normal;
			break;
		case GMT_GRID_IS_COMPLEX_REAL:
			p = (shape == GMT_IS_ROW_FORMAT) ? GMTAPI_2D_to_index_C_cplx_real : GMTAPI_2D_to_index_F_cplx_real;
			break;
		case GMT_GRID_IS_COMPLEX_IMAG:
			p = (shape == GMT_IS_ROW_FORMAT) ? GMTAPI_2D_to_index_C_cplx_imag : GMTAPI_2D_to_index_F_cplx_imag;
			break;
		default:
			printf (/*API, GMT_MSG_NORMAL, */ "GMTAPI_get_2D_to_index: Illegal mode passed - aborting\n");
			//API_exit (API, EXIT_FAILURE);
			return (NULL); //nishita
	}
	return (p);
}

/* Note: Many/all of these do not need to check if API == NULL since they are called from functions that do. */
/* Private functions used by this library only.  These are not accessed outside this file. */

double GMTAPI_get_val (struct GMTAPI_CTRL *API, union GMT_UNIVECTOR *u, uint64_t row, unsigned int type)
{	/* Returns a double value from the <type> column array pointed to by the union pointer *u, at row position row.
 	 * Used in GMTAPI_Import_Dataset and GMTAPI_Import_Grid. */
	double val;

	switch (type) {	/* Use type to select the correct array from which to extract a value */
		case GMT_UCHAR:		val = u->uc1[row];		break;
		case GMT_CHAR:		val = u->sc1[row];		break;
		case GMT_USHORT:	val = u->ui2[row];		break;
		case GMT_SHORT:		val = u->si2[row];		break;
		case GMT_UINT:		val = u->ui4[row];		break;
		case GMT_INT:		val = u->si4[row];		break;
		case GMT_ULONG:		val = (double)u->ui8[row];	break;
		case GMT_LONG:		val = (double)u->si8[row];	break;
		case GMT_FLOAT:		val = u->f4[row];		break;
		case GMT_DOUBLE:	val = u->f8[row];		break;
		default:
			printf (/*API, GMT_MSG_NORMAL, */ "Internal error in GMTAPI_get_val: Passed bad type (%d), returning NaN\n", type);
			val = API->GMT->session.d_NaN;
			API->error = GMT_NOT_A_VALID_TYPE;
			break;
	}
	return (val);
}

void * GMT_Get_Record (void *V_API, unsigned int mode, int *retval)
{
	/* Retrieves the next data record from the virtual input source and
	 * returns the number of columns found via *retval (unless retval == NULL).
	 * If current record is a segment header then we return 0.
	 * If we reach EOF then we return EOF.
	 * mode is either GMT_READ_DOUBLE (data columns), GMT_READ_TEXT (text string) or
	 *	GMT_READ_MIXED (expect data but tolerate read errors).
	 * Also, if (mode | GMT_FILE_BREAK) is true then we will return empty-handed
	 *	when we get to the end of a file except the final file (which is EOF).
	 *	The calling module can then take actions appropriate between data files.
	 * The double array OR text string is returned via the pointer *record.
	 * If not a data record we return NULL, and pass status via API->GMT->current.io.status.
	 */

	int status;
	int64_t n_fields = 0;
	uint64_t *p = NULL, col, ij, n_nan;
	bool get_next_record;
	char *t_record = NULL;
	void *record = NULL;
	p_func_size_t GMT_2D_to_index = NULL;
	struct GMTAPI_DATA_OBJECT *S_obj = NULL;
	struct GMT_TEXTSET *DT_obj = NULL;
	struct GMT_DATASET *DS_obj = NULL;
	struct GMT_MATRIX *M_obj = NULL;
	struct GMT_VECTOR *V_obj = NULL;
	struct GMTAPI_CTRL *API = NULL;
	
	if (V_API == NULL)
		return_null (V_API, GMT_NOT_A_SESSION);
		//return (GMT_NOT_A_SESSION); //nishita
	
	API = gmt_get_api_ptr ((struct GMTAPI_CTRL *)(struct GMTAPI_CTRL *)V_API);
	if (retval) *retval = 0;
	
	if (!API->io_enabled[GMT_IN])
	{
		return_null (API, GMT_ACCESS_NOT_ENABLED);
	}
		//return (GMT_ACCESS_NOT_ENABLED); //nishita
	S_obj = API->object[API->current_item[GMT_IN]];	/* Shorthand for the current data source we are working on */
	API->GMT->current.io.read_mixed = (mode == GMT_READ_MIXED);	/* Cannot worry about constant # of columns if text is present */

	do {	/* We do this until we can secure the next record or we run out of records (and return EOF) */
		get_next_record = false;	/* We expect to read one data record and return */
		API->GMT->current.io.status = 0;	/* Initialize status to OK */
		if (S_obj->status == GMT_IS_USED) {		/* Finished reading from this resource, go to next resource */
			if (API->GMT->current.io.ogr == GMT_OGR_TRUE)
				return_null (API, GMT_OGR_ONE_TABLE_ONLY);	/* Only allow single tables if GMT/OGR */
				//return (GMT_OGR_ONE_TABLE_ONLY); //nishita
			if (GMTAPI_Next_Data_Object (API, S_obj->family, GMT_IN) == EOF)	/* That was the last source, return */
				n_fields = EOF;
			else {
				S_obj = API->object[API->current_item[GMT_IN]];	/* Shorthand for the next data source to work on */
				get_next_record = true;				/* Since we haven't read the next record yet */
			}
			continue;
		}
		switch (S_obj->method) {
			case GMT_IS_FILE:	/* File, stream, and fd are all the same for us, regardless of data or text input */
		 	case GMT_IS_STREAM:
		 	case GMT_IS_FDESC:
				

				{
					char filename[0xFFF];
					int fno;
    					ssize_t r;
					int MAXSIZE = 0xFFF;
    					char proclnk[0xFFF];
					
					fno = fileno(S_obj->fp);
					sprintf(proclnk, "/proc/self/fd/%d", fno);
					r = readlink(proclnk, filename, MAXSIZE);
					if (r < 0)
					{
					    printf("failed to readlink\n");
					    exit(1);
					}
					filename[r] = '\0';
				}


		 		record = S_obj->import (API->GMT, S_obj->fp, &(S_obj->n_expected_fields), &status);	/* Get that next record */
				
				n_fields = S_obj->n_columns = status;	/* Get that next record */
				if (API->GMT->current.io.status & GMT_IO_EOF) {			/* End-of-file in current file (but there may be many files) */
					S_obj->status = GMT_IS_USED;	/* Mark as read */
					if (S_obj->close_file) {	/* Close if it was a file that we opened earlier */
						GMT_fclose (API->GMT, S_obj->fp);
						S_obj->close_file = false;
					}
					if (GMTAPI_Next_Data_Object (API, S_obj->family, GMT_IN) == EOF)	/* That was the last source, return */
						n_fields = EOF;					/* EOF is ONLY returned when we reach the end of the LAST data file */
					else if (mode & GMT_FILE_BREAK) {			/* Return empty handed to indicate a break between files */
						n_fields = GMT_IO_NEXT_FILE;			/* We flag this situation with a special return value */
						API->GMT->current.io.status = GMT_IO_NEXT_FILE;
					}
					else {	/* Get ready to read the next data file */
						S_obj = API->object[API->current_item[GMT_IN]];	/* Shorthand for the next data source to work on */
						get_next_record = true;				/* Since we haven't read the next record yet */
					}
					API->GMT->current.io.tbl_no++;				/* Update number of tables we have processed */
				}
				else
					S_obj->status = GMT_IS_USING;				/* Mark as being read */

				if (GMT_REC_IS_DATA (API->GMT) && S_obj->n_expected_fields != GMT_MAX_COLUMNS) API->GMT->common.b.ncol[GMT_IN] = S_obj->n_expected_fields;	/* Set the actual column count */
				break;

			case GMT_IS_DUPLICATE + GMT_VIA_MATRIX:	/* Here we copy/read from a user memory location */
			case GMT_IS_REFERENCE + GMT_VIA_MATRIX:
				if (API->current_rec[GMT_IN] >= S_obj->n_rows) {	/* Our only way of knowing we are done is to quit when we reach the number of rows that was registered */
					API->GMT->current.io.status = GMT_IO_EOF;
					S_obj->status = GMT_IS_USED;	/* Mark as read */
					if (retval) *retval = EOF;
					return (NULL);	/* Done with this array */
				}
				else
					S_obj->status = GMT_IS_USING;				/* Mark as being read */
				M_obj = S_obj->resource;
				GMT_2D_to_index = GMTAPI_get_2D_to_index (API, M_obj->shape, GMT_GRID_IS_REAL);
				for (col = n_nan = 0; col < S_obj->n_columns; col++) {	/* We know the number of columns from registration */
					ij = GMT_2D_to_index (API->current_rec[GMT_IN], col, M_obj->dim);
					API->GMT->current.io.curr_rec[col] = GMTAPI_get_val (API, &(M_obj->data), ij, M_obj->type);
					if (GMT_is_dnan (API->GMT->current.io.curr_rec[col])) n_nan++;
				}
				if (n_nan == S_obj->n_columns) {
					API->GMT->current.io.status = GMT_IO_SEGMENT_HEADER;	/* Flag as segment header */
					record = NULL;
				}
				else
					record = API->GMT->current.io.curr_rec;
				n_fields = S_obj->n_columns;
				break;

			 case GMT_IS_DUPLICATE + GMT_VIA_VECTOR:	/* Here we copy from a user memory location that points to an array of column vectors */
			 case GMT_IS_REFERENCE + GMT_VIA_VECTOR:
				if (API->current_rec[GMT_IN] >= S_obj->n_rows) {	/* Our only way of knowing we are done is to quit when we reach the number or rows that was registered */
					API->GMT->current.io.status = GMT_IO_EOF;
					S_obj->status = GMT_IS_USED;	/* Mark as read */
					if (retval) *retval = EOF;
					return (NULL);	/* Done with this array */
				}
				else
					S_obj->status = GMT_IS_USING;				/* Mark as being read */
				V_obj =  S_obj->resource;
				for (col = n_nan = 0; col < S_obj->n_columns; col++) {	/* We know the number of columns from registration */
					API->GMT->current.io.curr_rec[col] = GMTAPI_get_val (API, &(V_obj->data[col]), API->current_rec[GMT_IN], V_obj->type[col]);
					if (GMT_is_dnan (API->GMT->current.io.curr_rec[col])) n_nan++;
				}
				if (n_nan == S_obj->n_columns) {
					API->GMT->current.io.status = GMT_IO_SEGMENT_HEADER;	/* Flag as segment header */
					record = NULL;
				}
				else
					record = API->GMT->current.io.curr_rec;
				n_fields = S_obj->n_columns;
				break;

			case GMT_IS_REFERENCE:	/* Only for textsets and datasets */
				p = API->GMT->current.io.curr_pos[GMT_IN];	/* Shorthand used below */
				if (S_obj->family == GMT_IS_DATASET) {
					DS_obj =  S_obj->resource;

					status = 0;
					if (p[2] == DS_obj->table[p[0]]->segment[p[1]]->n_rows) {	/* Reached end of current segment */
						p[1]++, p[2] = 0;				/* Advance to next segments 1st row */
						status = GMT_IO_SEGMENT_HEADER;			/* Indicates a segment boundary */
					}
					if (p[1] == DS_obj->table[p[0]]->n_segments) {		/* Also the end of a table ("file") */
						p[0]++, p[1] = 0;
						if (mode & GMT_FILE_BREAK) {			/* Return empty handed to indicate a break between files */
							status = GMT_IO_NEXT_FILE;
							n_fields = GMT_IO_NEXT_FILE;
							record = NULL;
						}
					}
					if (p[0] == (uint64_t)DS_obj->n_tables) {	/* End of entire data set */
						status = GMT_IO_EOF;
						n_fields = EOF;
						record = NULL;
						S_obj->status = GMT_IS_USED;	/* Mark as read */
					}
					if (!status) {	/* OK get the record */
						for (col = n_nan = 0; col < DS_obj->n_columns; col++) {
							API->GMT->current.io.curr_rec[col] = DS_obj->table[p[0]]->segment[p[1]]->coord[col][p[2]];
							if (GMT_is_dnan (API->GMT->current.io.curr_rec[col])) n_nan++;
						}
						p[2]++;
						n_fields = API->GMT->common.b.ncol[GMT_IN] = DS_obj->n_columns;
						if (n_nan == DS_obj->n_columns) {
							API->GMT->current.io.status = GMT_IO_SEGMENT_HEADER;	/* Flag as segment header */
							record = NULL;
						}
						else
							record = API->GMT->current.io.curr_rec;
						S_obj->status = GMT_IS_USING;	/* Mark as read */
					}
					API->GMT->current.io.status = status;
				}
				if (S_obj->family == GMT_IS_TEXTSET) {
					DT_obj = S_obj->resource;
					if (p[2] == DT_obj->table[p[0]]->segment[p[1]]->n_rows) {p[1]++, p[2] = 0;}
					if (p[1] == DT_obj->table[p[0]]->n_segments) {p[0]++, p[1] = 0;}
					if (p[0] == DT_obj->n_tables) {
						n_fields = EOF;
						API->GMT->current.io.status = GMT_IO_EOF;
						record = NULL;
						S_obj->status = GMT_IS_USED;	/* Mark as read */
					}
					else {
						t_record = DT_obj->table[p[0]]->segment[p[1]]->record[p[2]++];
						API->GMT->current.io.status = 0;
						if (t_record[0] == API->GMT->current.setting.io_seg_marker[GMT_IN]) {
							/* Segment header: Just save the header content, not the
							 *                 marker and leading whitespace
							 */
							strncpy (API->GMT->current.io.segment_header,
								GMT_trim_segheader (API->GMT, t_record), GMT_BUFSIZ);
							API->GMT->current.io.status = GMT_IO_SEGMENT_HEADER;
							record = NULL;
						}
						else {	/* Regular record */
							strncpy (API->GMT->current.io.current_record, t_record, GMT_BUFSIZ);
							record = t_record;
						}
						n_fields = 1;
						S_obj->status = GMT_IS_USING;	/* Mark as read */
					}
				}
				break;
			default:
				printf (/*API, GMT_MSG_NORMAL, */ "GMTAPI: Internal error: GMT_Get_Record called with illegal method\n");
				break;
		}
	} while (get_next_record);

	if (!(n_fields == EOF || n_fields == GMT_IO_NEXT_FILE)) API->current_rec[GMT_IN]++;	/* Increase record count, unless EOF */

	if (retval) *retval = (int)n_fields;	/* Requested we return the number of fields found */
	return (record);	/* Return pointer to current record */
}

int GMT_End_IO (void *V_API, unsigned int direction, unsigned int mode)
{
	/* Terminates the i/o mechanism for either input or output (given by direction).
	 * GMT_End_IO must be called after all data i/o is completed.
	 * direction:	Either GMT_IN or GMT_OUT
	 * mode:	Either GMT_IO_DONE (nothing), GMT_IO_RESET (let all resources be accessible again), or GMT_IO_UNREG (unreg all accessed resources).
	 * NOTE: 	Mode not yet implemented until we see a use.
	 * Returns:	false if successfull, true if error.
	 */
	int error = 0;
	unsigned int item, method = 0, via = 0;
	struct GMTAPI_DATA_OBJECT *S_obj = NULL;
	struct GMTAPI_CTRL *API = NULL;

	if (V_API == NULL)
		//return_error (V_API, GMT_NOT_A_SESSION);
		return (GMT_NOT_A_SESSION); //nishita
	if (!(direction == GMT_IN || direction == GMT_OUT))
		//return_error (V_API, GMT_NOT_A_VALID_DIRECTION);
		return (GMT_NOT_A_VALID_DIRECTION); //nishita
	if (mode > GMT_IO_UNREG)
		//return_error (V_API, GMT_NOT_A_VALID_IO_MODE);
		return (GMT_NOT_A_VALID_IO_MODE); //nishita
	API = gmt_get_api_ptr ((struct GMTAPI_CTRL *)(struct GMTAPI_CTRL *)V_API);
	GMT_free_ogr (API->GMT, &(API->GMT->current.io.OGR), 0);	/* Free segment-related array */
	if (direction == GMT_OUT) {		/* Finalize output issues */
		S_obj = API->object[API->current_item[GMT_OUT]];	/* Shorthand for the data source we are working on */
		if (S_obj) {	/* Dealt with file i/o */
			S_obj->status = GMT_IS_USED;	/* Done writing to this destination */
			if (S_obj->method >= GMT_VIA_VECTOR) {
				via = (S_obj->method / GMT_VIA_VECTOR) - 1;
				method = S_obj->method - (via + 1) * GMT_VIA_VECTOR;	/* Array index that have any GMT_VIA_* removed */
			}
			else
				method = S_obj->method;
			if ((method == GMT_IS_DUPLICATE || method == GMT_IS_REFERENCE) && API->io_mode[GMT_OUT] == GMT_BY_REC) {	/* GMT_Put_Record: Must realloc last segment and the tables segment array */
				if (S_obj->actual_family == GMT_IS_DATASET) {	/* Dataset type */
					struct GMT_DATASET *D_obj =  S_obj->resource;
					if (D_obj && D_obj->table && D_obj->table[0]) {
						struct GMT_DATATABLE *T_obj = D_obj->table[0];
						uint64_t *p = API->GMT->current.io.curr_pos[GMT_OUT];	/* Short-hand for counts of tbl, seg, rows */
						if (!T_obj->segment[p[GMT_SEG]]) T_obj->segment[p[GMT_SEG]] = GMT_memory (API->GMT, NULL, 1, struct GMT_DATASEGMENT);
						GMT_assign_segment (API->GMT, T_obj->segment[p[GMT_SEG]], p[GMT_ROW], T_obj->n_columns);	/* Allocate and place arrays into segment */
						if (API->GMT->current.io.current_record[0]) T_obj->segment[p[GMT_SEG]]->header = strdup (API->GMT->current.io.current_record);
						p[GMT_SEG]++;	/* Total number of segments */
						T_obj->n_segments++;
						/* Realloc final number of segments */
						if (p[GMT_SEG] < T_obj->n_alloc) T_obj->segment = GMT_memory (API->GMT, T_obj->segment, T_obj->n_segments, struct GMT_DATASEGMENT *);
						D_obj->n_segments = T_obj->n_segments;
						GMT_set_tbl_minmax (API->GMT, T_obj);
						GMT_set_dataset_minmax (API->GMT, D_obj);
					}
				}
				else if (S_obj->actual_family == GMT_IS_MATRIX) {	/* Matrix type */
					if (S_obj->n_alloc != API->current_rec[GMT_OUT]) {	/* Must finalize memory */
						struct GMT_MATRIX *M_obj = S_obj->resource;
						size_t size = S_obj->n_alloc = API->current_rec[GMT_OUT];
						size *= API->GMT->common.b.ncol[GMT_OUT];
						if ((error = GMT_alloc_univector (API->GMT, &(M_obj->data), M_obj->type, size)) != GMT_OK)
							//return_error (V_API, error);
							return ( error);  //nishita
					}
				}
				else if (S_obj->actual_family == GMT_IS_VECTOR) {	/* Vector type */
					if (S_obj->n_alloc != API->current_rec[GMT_OUT]) {	/* Must finalize memory */
						struct GMT_VECTOR *V_obj = S_obj->resource;
						size_t size = S_obj->n_alloc = API->current_rec[GMT_OUT];
						V_obj->n_rows = size;
						if ((error = gmt_alloc_vectors (API->GMT, V_obj)) != GMT_OK)
							//return_error (V_API, error);
							return (error);
					}
				}
				else if (S_obj->actual_family == GMT_IS_TEXTSET) {	/* Textset type */
					struct GMT_TEXTSET *D_obj = S_obj->resource;
					if (D_obj && D_obj->table && D_obj->table[0]) {
						struct GMT_TEXTTABLE *T_obj = D_obj->table[0];
						uint64_t *p = API->GMT->current.io.curr_pos[GMT_OUT];
						if (p[GMT_SEG] > 0) T_obj->segment[p[GMT_SEG]]->record = GMT_memory (API->GMT, T_obj->segment[p[GMT_SEG]]->record, T_obj->segment[p[GMT_SEG]]->n_rows, char *);	/* Last segment */
						T_obj->segment = GMT_memory (API->GMT, T_obj->segment, T_obj->n_segments, struct GMT_TEXTSEGMENT *);
						D_obj->n_segments = T_obj->n_segments;
					}
				}
			}
			if (S_obj->close_file) {	/* Close file that we opened earlier */
				GMT_fclose (API->GMT, S_obj->fp);
				S_obj->close_file = false;
			}
		}
	}
	API->io_enabled[direction] = false;	/* No longer OK to access resources */
	API->current_rec[direction] = 0;	/* Reset for next use */
	for (item = 0; item < API->n_objects; item++) {
		if (!API->object[item]) continue;	/* Skip empty object */
		if (API->object[item]->direction != direction) continue;	/* Not the required direction */
		if (API->object[item]->selected) API->object[item]->selected = false;	/* No longer a selected resource */
	}

	//printf(/*API, GMT_MSG_DEBUG,*/ "GMT_End_IO: %s resource access is now disabled\n", GMT_direction[direction]);

	return (GMT_OK);	/* No error encountered */
}
void * GMT_Duplicate_Data (void *V_API, unsigned int family, unsigned int mode, void *data)
{
	/* Create an duplicate container of the requested kind and optionally allocate space
	 * or duplicate content.
	 * The known families are GMT_IS_{DATASET,TEXTSET,GRID,CPT,IMAGE}.
 	 * Pass mode as one of GMT_DUPLICATE_{NONE|ALLOC|DATA} to just duplicate the
	 * container and header structures, allocate space of same dimensions as original,
	 * or allocate space and duplicate contents.  For GMT_IS_{DATA|TEXT}SET you may add
	 * the modifiers GMT_ALLOC_VERTICAL or GMT_ALLOC_HORIZONTAL. Also, for GMT_IS_DATASET
	 * you can manipulate the incoming data->dim to overwrite the number of items allocated.
	 * [By default we follow the dimensions of hte incoming data].
	 *
	 * Return: Pointer to new resource, or NULL if an error (set via API->error).
	 */

	int object_ID, item;
	unsigned int geometry = 0U, pmode = 0U;
	void *new_obj = NULL;
	struct GMTAPI_CTRL *API = NULL;

	if (V_API == NULL)
		//return_null (V_API, GMT_NOT_A_SESSION);
		return (GMT_NOT_A_SESSION);
	API = gmt_get_api_ptr ((struct GMTAPI_CTRL *)V_API);

	switch (family) {	/* dataset, cpt, text, grid , image, vector, matrix */
		case GMT_IS_GRID:	/* GMT grid, allocate header but not data array */
			new_obj = GMT_duplicate_grid (API->GMT, (struct GMT_GRID *)data, mode);
			geometry = GMT_IS_SURFACE;
			break;
//#ifdef HAVE_GDAL
	//	case GMT_IS_IMAGE:	/* GMT image, allocate header but not data array */
		//	new_obj = GMT_duplicate_image (API->GMT, data, mode);
			//geometry = GMT_IS_SURFACE;
			//break;
//#endif
		/*case GMT_IS_DATASET:	 GMT dataset, allocate the requested tables, segments, rows, and columns
			pmode = (mode & (GMT_ALLOC_VERTICAL + GMT_ALLOC_HORIZONTAL));	 Just isolate any special allocation modes
			mode -= pmode;	 Remove the hor/ver flags from the rest of mode
			if (mode == GMT_DUPLICATE_DATA)
				new_obj = GMT_duplicate_dataset (API->GMT, data, pmode, &geometry);
			else if (mode == GMT_DUPLICATE_ALLOC) {	 Allocate data set of same size, possibly modulated by Din->dim (of > 0) and pmode
				struct GMT_DATASET *Din = data;	 We know this is a GMT_DATASET pointer
				new_obj = GMT_alloc_dataset (API->GMT, data, Din->dim[GMT_ROW], Din->dim[GMT_COL], pmode);
				geometry = Din->geometry;
				GMT_memset (Din->dim, 4, uint64_t);	 Reset alloc dimensions
			}
			else {	 Just want a dataset structure
				struct GMT_DATASET *Din = data;	 We know this is a GMT_DATASET pointer
				new_obj = GMT_memory (API->GMT, NULL, 1, struct GMT_DATASET);
				geometry = Din->geometry;
			}
			break;*/
		/*case GMT_IS_TEXTSET:	 GMT text dataset, allocate the requested tables, segments, and rows
			pmode = (mode & (GMT_ALLOC_VERTICAL + GMT_ALLOC_HORIZONTAL));	 Just isolate any special allocation modes
			mode -= pmode;	 Remove the hor/ver flags from the rest of mode
			if (mode == GMT_DUPLICATE_DATA)
				new_obj = GMT_duplicate_textset (API->GMT, data, pmode);
			else if (mode == GMT_DUPLICATE_ALLOC)	 Allocate text set of same size, possibly modulated by pmode
				new_obj =  GMT_alloc_textset (API->GMT, data, pmode);
			else	 Just want a dataset structure
				new_obj = GMT_memory (API->GMT, NULL, 1, struct GMT_TEXTSET);
			geometry = GMT_IS_NONE;
			break;
		case GMT_IS_CPT:	 GMT CPT table, allocate one with space for dim[0] color entries
			//new_obj = GMT_duplicate_palette (API->GMT, data, 0);
			//geometry = GMT_IS_NONE;
			//break;*/
		default:
			API->error = GMT_NOT_A_VALID_FAMILY;
			break;
	}
	if (API->error)
		//return_null (API, API->error);
		return(API->error); //nishita

	/* Now register this dataset so it can be deleted by GMT_Destroy_Data */
	if ((object_ID = GMT_Register_IO (API, family, GMT_IS_REFERENCE, geometry, GMT_IN, NULL, new_obj)) == GMT_NOTSET)
		//return_null (API, API->error);	/* Failure to register */
		return(API->error); //nishita
	if ((item = GMTAPI_Validate_ID (API, family, object_ID, GMT_IN)) == GMT_NOTSET)
		//return_null (API, API->error);
		return(API->error); //nishita
	API->object[item]->geometry = geometry;	/* Ensure same geometry */
	API->object[item]->data = new_obj;		/* Retain pointer to the allocated data so we use garbage collection later */

	printf (/*API, GMT_MSG_LONG_VERBOSE,*/ "Successfully duplicated a %s\n", GMT_family[family]);

	return (new_obj);
}

int GMTAPI_Begin_IO (struct GMTAPI_CTRL *API, unsigned int direction)
{
	/* Initializes the i/o mechanism for either input or output (given by direction).
	 * GMTAPI_Begin_IO must be called before any bulk data i/o is allowed.
	 * direction:	Either GMT_IN or GMT_OUT.
	 * Returns:	false if successfull, true if error.
	 */

	if (API == NULL) return_error (API, GMT_NOT_A_SESSION);
	if (!(direction == GMT_IN || direction == GMT_OUT)) return_error (API, GMT_NOT_A_VALID_DIRECTION);
	if (!API->registered[direction]) printf(/*API, GMT_MSG_DEBUG,*/ "GMTAPI_Begin_IO: Warning: No %s resources registered\n", GMT_direction[direction]);

	API->io_mode[direction] = GMT_BY_SET;
	API->io_enabled[direction] = true;	/* OK to access resources */
	API->GMT->current.io.ogr = GMT_OGR_UNKNOWN;
	API->GMT->current.io.read_mixed = false;
	API->GMT->current.io.need_previous = (API->GMT->common.g.active || API->GMT->current.io.skip_duplicates);
	API->GMT->current.io.segment_header[0] = API->GMT->current.io.current_record[0] = 0;
	printf(/*API, GMT_MSG_DEBUG,*/ "GMTAPI_Begin_IO: %s resource access is now enabled [container]\n", GMT_direction[direction]);

	return (GMT_OK);	/* No error encountered */
}

size_t GMTAPI_set_grdarray_size (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *h, unsigned int mode, double *wesn)
{	/* Determines size of grid given grid spacing and grid domain in h.
 	 * However, if wesn is given and not empty we use that sub-region instead.
 	 * Finally, the current pad is used when calculating the grid size.
	 * NOTE: This function leaves h unchanged by testing on a temporary header. */
	struct GMT_GRID_HEADER *h_tmp = NULL;
	size_t size;

	/* Must duplicate header and possibly reset wesn, then set pad and recalculate the dims */
	h_tmp = GMT_memory (GMT, NULL, 1, struct GMT_GRID_HEADER);
	GMT_memcpy (h_tmp, h, 1, struct GMT_GRID_HEADER);
	h_tmp->complex_mode |= mode;	/* Set the mode-to-be so that if complex the size is doubled */

	if (!full_region (wesn)) {
		GMT_memcpy (h_tmp->wesn, wesn, 4, double);	/* Use wesn instead of header info */
		GMT_adjust_loose_wesn (GMT, wesn, h);		/* Subset requested; make sure wesn matches header spacing */
	}
	GMT_grd_setpad (GMT, h_tmp, GMT->current.io.pad);	/* Use the system pad setting by default */
	GMT_set_grddim (GMT, h_tmp);			/* Computes all integer parameters */
	size = h_tmp->size;				/* This is the size needed to hold grid + padding */
	GMT_free (GMT, h_tmp);
	return (size);
}

bool GMTAPI_adjust_grdpadding (struct GMT_GRID_HEADER *h, unsigned int *pad)
{	/* Compares current grid pad status to output pad requested.  If we need
	 * to adjust a pad we return true here, otherwise false. */
	unsigned int side;

	for (side = 0; side < 4; side++) if (h->pad[side] != pad[side]) return (true);
	return (false);
}

int gmt_open_grd (struct GMT_CTRL *GMT, char *file, struct GMT_GRID *G, char mode, unsigned int access_mode)
{
	/* Read or write the header structure and initialize row-by-row machinery.
	 * We fill the GMT_GRID_ROWBYROW structure with all the required information.
	 * mode can be w or r.  Upper case W or R refers to headerless
	 * grdraster-type files.  The access_mode dictates if we automatically advance
	 * row counter to next row after read/write or if we use the rec_no to seek
	 * first.
	 */

	int r_w, err;
	bool header = true, magic = true, alloc = false;
	int cdf_mode[3] = { NC_NOWRITE, NC_WRITE, NC_WRITE};	/* MUST be ints */
	char *bin_mode[3] = { "rb", "rb+", "wb"};
	char *fmt = NULL;
	struct GMT_GRID_ROWBYROW *R = gmt_get_rbr_ptr ((struct GMT_GRID_ROWBYROW *)G->extra);	/* Shorthand to row-by-row book-keeping structure */

	if (mode == 'r' || mode == 'R') {	/* Open file for reading */
		if (mode == 'R') {	/* File has no header; can only work if G->header has been set already, somehow */
			header = false;
			if (G->header->nx == 0 || G->header->ny == 0) {
				printf(/*GMT->parent, GMT_MSG_NORMAL,*/ "Unable to read header-less grid file %s without a preset header structure\n", file);
				return (GMT_GRDIO_OPEN_FAILED);
			}
		}
		r_w = 0;	mode = 'r';
	}
	else if (mode == 'W') {	/* Write headerless grid */
		r_w = 2;	mode = 'w';
		header = magic = false;
	}
	else {	/* Regular writing of grid with header */
		r_w = 1;
		magic = false;
	}
	if (header) {
		if (mode == 'r' && !R->open)	/* First time reading the info */
			GMT_read_grd_info (GMT, file, G->header);
		else if (R->open)		/* Coming back to update the header */
			GMT_update_grd_info (GMT, file, G->header);
		else				/* First time writing the header */
			GMT_write_grd_info (GMT, file, G->header);
	}
	else /* Fallback to existing header */
		GMT_err_trap (GMT_grd_get_format (GMT, file, G->header, magic));
	if (R->open) return (GMT_NOERROR);	/* Already set the first time */
	fmt = GMT->session.grdformat[G->header->type];
	if (fmt[0] == 'c') {		/* Open netCDF file, old format */
		GMT_err_trap (nc_open (G->header->name, cdf_mode[r_w], &R->fid));
		R->edge[0] = G->header->nx;
		R->start[0] = 0;
		R->start[1] = 0;
	}
	else if (fmt[0] == 'n') {	/* Open netCDF file, COARDS-compliant format */
		GMT_err_trap (nc_open (G->header->name, cdf_mode[r_w], &R->fid));
		R->edge[0] = 1;
		R->edge[1] = G->header->nx;
		R->start[0] = G->header->ny-1;
		R->start[1] = 0;
	}
	else {		/* Regular binary file with/w.o standard GMT header, or Sun rasterfile */
	//(1) //nishita
	//
	if (r_w == 0) {	/* Open for plain reading */
			if ((R->fp = GMT_fopen (GMT, G->header->name, bin_mode[0])) == NULL)
				return (GMT_GRDIO_OPEN_FAILED);
		}
		else if ((R->fp = GMT_fopen (GMT, G->header->name, bin_mode[r_w])) == NULL)
			return (GMT_GRDIO_CREATE_FAILED);
		/* Seek past the grid header, unless there is none */
	
		if (header && fseek (R->fp, (off_t)GMT_GRID_HEADER_SIZE, SEEK_SET)) return (GMT_GRDIO_SEEK_FAILED);
		alloc = (fmt[1] != 'f');	/* Only need to allocate the v_row array if grid is not float */
	}
	R->size = GMT_grd_data_size (GMT, G->header->type, &G->header->nan_value);
	R->check = !isnan (G->header->nan_value);
	R->open = true;

	if (fmt[1] == 'm')	/* Bit mask */
		R->n_byte = lrint (ceil (G->header->nx / 32.0)) * R->size;
	else if (fmt[0] == 'r' && fmt[1] == 'b')	/* Sun Raster uses multiple of 2 bytes */
		R->n_byte = lrint (ceil (G->header->nx / 2.0)) * 2 * R->size;
	else	/* All other */
		R->n_byte = G->header->nx * R->size;

	if (alloc) R->v_row = GMT_memory (GMT, NULL, R->n_byte, char);

	R->row = 0;
	R->auto_advance = (access_mode & GMT_GRID_ROW_BY_ROW_MANUAL) ? false : true;	/* Read sequentially or random-access rows */
	return (GMT_NOERROR);
}

void GMTAPI_info_to_grdheader (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *h, struct GMT_MATRIX *M_obj)
{	/* Unpacks the necessary items into the grid header from the matrix parameters */
	h->nx = (unsigned int)M_obj->n_columns;
	h->ny = (unsigned int)M_obj->n_rows;
	h->registration = M_obj->registration;
	GMT_memcpy (h->wesn, M_obj->range, 4, double);
	/* Compute xy_off and increments */
	h->xy_off = (h->registration == GMT_GRID_NODE_REG) ? 0.0 : 0.5;
	h->inc[GMT_X] = GMT_get_inc (GMT, h->wesn[XLO], h->wesn[XHI], h->nx, h->registration);
	h->inc[GMT_Y] = GMT_get_inc (GMT, h->wesn[YLO], h->wesn[YHI], h->ny, h->registration);
}

struct GMT_GRID * GMTAPI_Import_Grid (struct GMTAPI_CTRL *API, int object_ID, unsigned int mode, struct GMT_GRID *grid)
{	/* Handles the reading of a 2-D grid given in one of several ways.
	 * Get the entire grid:
 	 * 	mode = GMT_GRID_ALL reads both header and grid;
	 * Get a subset of the grid:  Call GMTAPI_Import_Grid twice:
	 * 	1. first with mode = GMT_GRID_HEADER_ONLY which reads header only.  Then, pass
	 *	   the new S_obj-> wesn to match your desired subregion
	 *	2. 2nd with mode = GMT_GRID_DATA_ONLY, which reads grid based on header's settings
	 * If the grid->data array is NULL it will be allocated for you.
	 */

	int item, new_item, new_ID;
	bool done = true, new = false, row_by_row, via = false;
 	uint64_t row, col, i0, i1, j0, j1, ij, ij_orig;
	size_t size;
	enum GMT_enum_gridio both_set = (GMT_GRID_HEADER_ONLY | GMT_GRID_DATA_ONLY);
	double dx, dy;
	p_func_size_t GMT_2D_to_index = NULL;
	struct GMT_GRID *G_obj = NULL, *G_orig = NULL;
	struct GMT_MATRIX *M_obj = NULL;
	struct GMTAPI_DATA_OBJECT *S_obj = NULL;

	//GMT_Report (API, GMT_MSG_DEBUG, "GMTAPI_Import_Grid: Passed ID = %d and mode = %d\n", object_ID, mode);
	
	if ((item = GMTAPI_Validate_ID (API, GMT_IS_GRID, object_ID, GMT_IN)) == GMT_NOTSET) return_null (API, API->error);

	S_obj = API->object[item];		/* Current data object */
	//if (S_obj->status != GMT_IS_UNUSED && !(mode & GMT_IO_RESET)) return_null (API, GMT_READ_ONCE);	/* Already read this resources before, so fail unless overridden by mode */
	if (S_obj->status != GMT_IS_UNUSED && S_obj->method == GMT_IS_FILE && !(mode & GMT_IO_RESET)) return_null (API, GMT_READ_ONCE);	/* Already read this file before, so fail unless overridden by mode */
	if ((mode & both_set) == both_set) mode -= both_set;	/* Allow users to have set GMT_GRID_HEADER_ONLY | GMT_GRID_DATA_ONLY; reset to GMT_GRID_ALL */
	row_by_row = ((mode & GMT_GRID_ROW_BY_ROW) || (mode & GMT_GRID_ROW_BY_ROW_MANUAL));
	if (row_by_row && S_obj->method != GMT_IS_FILE)
	{
		//GMT_Report (API, GMT_MSG_NORMAL, "Can only use method GMT_IS_FILE when row-by-row reading of grid is selected\n");
		return_null (API, GMT_NOT_A_VALID_METHOD);
	}

	if (S_obj->region && grid) {	/* See if this is really a subset or just the same region as the grid */
		if (grid->header->wesn[XLO] == S_obj->wesn[XLO] && grid->header->wesn[XHI] == S_obj->wesn[XHI] && grid->header->wesn[YLO] == S_obj->wesn[YLO] && grid->header->wesn[YHI] == S_obj->wesn[YHI]) S_obj->region = false;
	}
	
	switch (S_obj->method) {
		case GMT_IS_FILE:	/* Name of a grid file on disk */
			if (grid == NULL) {	/* Only allocate grid struct when not already allocated */
				if (mode & GMT_GRID_DATA_ONLY) return_null (API, GMT_NO_GRDHEADER);		/* For mode & GMT_GRID_DATA_ONLY grid must already be allocated */
				G_obj = GMT_create_grid (API->GMT);
				new = true;
			}
			else
				G_obj = grid;	/* We are working on a grid already allocated */
			done = (mode & GMT_GRID_HEADER_ONLY) ? false : true;	/* Not done until we read grid */
			
			if (! (mode & GMT_GRID_DATA_ONLY)) {		/* Must init header and read the header information from file */
				//printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
				if (row_by_row) {	/* Special row-by-row processing mode */
					char r_mode = (mode & GMT_GRID_NO_HEADER) ? 'R' : 'r';
					/* If we get here more than once we only allocate extra once */
					if (G_obj->extra == NULL) G_obj->extra = GMT_memory (API->GMT, NULL, 1, struct GMT_GRID_ROWBYROW);
					if (gmt_open_grd (API->GMT, S_obj->filename, G_obj, r_mode, mode)) {	/* Open the grid for incremental row reading */
						if (new) GMT_free_grid (API->GMT, &G_obj, false);
						return_null (API, GMT_GRID_READ_ERROR);
					}
				}
				else if (/*GMT_err_pass (API->GMT,*/GMT_read_grd_info (API->GMT, S_obj->filename, G_obj->header) !=0) { //nishita
					if (new) GMT_free_grid (API->GMT, &G_obj, false);
					return_null (API, GMT_GRID_READ_ERROR);
				}
				if (mode & GMT_GRID_HEADER_ONLY) break;	/* Just needed the header, get out of here */
			}
			//printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
			/* Here we will read the grid data themselves. */
			/* To get a subset we use wesn that is not NULL or contain 0/0/0/0.
			 * Otherwise we extract the entire file domain */
			size = GMTAPI_set_grdarray_size (API->GMT, G_obj->header, mode, S_obj->wesn);	/* Get array dimension only, which includes padding */
			if (!G_obj->data) {	/* Array is not allocated yet, do so now. We only expect header (and possibly w/e/s/n subset) to have been set correctly */
				G_obj->header->size = size;
				G_obj->data = GMT_memory_aligned (API->GMT, NULL, G_obj->header->size, float);
			}
			else {	/* Already have allocated space; check that it is enough */
				if (size > G_obj->header->size) return_null (API, GMT_GRID_READ_ERROR);
			}
			//GMT_Report (API, GMT_MSG_LONG_VERBOSE, "Reading grid from file %s\n", S_obj->filename);
			
			//if (GMT_err_pass (API->GMT, GMT_read_grd (API->GMT, S_obj->filename, G_obj->header, G_obj->data, S_obj->wesn,
				//			API->GMT->current.io.pad, mode), S_obj->filename))
			if(GMT_read_grd (API->GMT, S_obj->filename, G_obj->header, G_obj->data, S_obj->wesn,
							API->GMT->current.io.pad, mode) != 0) //added by nishita
				return_null (API, GMT_GRID_READ_ERROR);
			//if (GMT_err_pass (API->GMT, GMT_grd_BC_set (API->GMT, G_obj, GMT_IN), S_obj->filename))
			if(GMT_grd_BC_set (API->GMT, G_obj, GMT_IN)  != 0) //added by nishita
				return_null (API, GMT_GRID_BC_ERROR);	/* Set boundary conditions */
			G_obj->alloc_mode = GMT_ALLOCATED_BY_GMT;
			S_obj->resource = G_obj;	/* Set resource pointer to the grid */
			break;

	 	case GMT_IS_DUPLICATE:	/* GMT grid and header in a GMT_GRID container object. */
			if ((G_orig = S_obj->resource) == NULL) return_null (API, GMT_PTR_IS_NULL);
			if (grid == NULL) {	/* Only allocate when not already allocated */
				if (mode & GMT_GRID_DATA_ONLY) return_null (API, GMT_NO_GRDHEADER);		/* For mode & GMT_GRID_DATA_ONLY grid must already be allocated */
				G_obj = GMT_create_grid (API->GMT);
			}
			else
				G_obj = grid;	/* We are passing in a grid already */
			done = (mode & GMT_GRID_HEADER_ONLY) ? false : true;	/* Not done until we read grid */
			if (! (mode & GMT_GRID_DATA_ONLY)) {	/* Must init header and copy the header information from the existing grid */
				GMT_memcpy (G_obj->header, G_orig->header, 1, struct GMT_GRID_HEADER);
				if (mode & GMT_GRID_HEADER_ONLY) break;	/* Just needed the header, get out of here */
			}
			/* Here we will read grid data. */
			/* To get a subset we use wesn that is not NULL or contain 0/0/0/0.
			 * Otherwise we use everything passed in */
			 
			//GMT_Report (API, GMT_MSG_LONG_VERBOSE, "Duplicating grid data from GMT_GRID memory location\n");
			
			if (!G_obj->data) {	/* Array is not allocated, do so now. We only expect header (and possibly subset w/e/s/n) to have been set correctly */
				G_obj->header->size = GMTAPI_set_grdarray_size (API->GMT, G_obj->header, mode, S_obj->wesn);	/* Get array dimension only, which may include padding */
				G_obj->data = GMT_memory_aligned (API->GMT, NULL, G_obj->header->size, float);
			}
			G_obj->alloc_mode = GMT_ALLOCATED_BY_GMT;
			if (!S_obj->region && GMT_grd_pad_status (API->GMT, G_obj->header, API->GMT->current.io.pad)) {	/* Want an exact copy with no subset and same padding */
				GMT_memcpy (G_obj->data, G_orig->data, G_orig->header->size, float);
				break;		/* Done with this grid */
			}
			/* Here we need to do more work: Either extract subset or add/change padding, or both. */
			/* Get start/stop row/cols for subset (or the entire domain) */
			dx = G_obj->header->inc[GMT_X] * G_obj->header->xy_off;	dy = G_obj->header->inc[GMT_Y] * G_obj->header->xy_off;
			j1 = (unsigned int)GMT_grd_y_to_row (API->GMT, G_obj->header->wesn[YLO]+dy, G_orig->header);
			j0 = (unsigned int)GMT_grd_y_to_row (API->GMT, G_obj->header->wesn[YHI]-dy, G_orig->header);
			i0 = (unsigned int)GMT_grd_x_to_col (API->GMT, G_obj->header->wesn[XLO]+dx, G_orig->header);
			i1 = (unsigned int)GMT_grd_x_to_col (API->GMT, G_obj->header->wesn[XHI]-dx, G_orig->header);
			GMT_memcpy (G_obj->header->pad, API->GMT->current.io.pad, 4, int);	/* Set desired padding */
			for (row = j0; row <= j1; row++) {
				for (col = i0; col <= i1; col++, ij++) {
					ij_orig = GMT_IJP (G_orig->header, row, col);	/* Position of this (row,col) in original grid organization */
					ij = GMT_IJP (G_obj->header, row, col);		/* Position of this (row,col) in output grid organization */
					G_obj->data[ij] = G_orig->data[ij_orig];
				}
			}
			GMT_BC_init (API->GMT, G_obj->header);	/* Initialize grid interpolation and boundary condition parameters */
			//if (GMT_err_pass (API->GMT, GMT_grd_BC_set (API->GMT, G_obj, GMT_IN), "Grid memory")) 
			if(GMT_grd_BC_set (API->GMT, G_obj, GMT_IN) != 0) //nishita
				return_null (API, GMT_GRID_BC_ERROR);	/* Set boundary conditions */
			break;

	 	case GMT_IS_REFERENCE:	/* GMT grid and header in a GMT_GRID container object by reference */
			if (S_obj->region) return_null (API, GMT_SUBSET_NOT_ALLOWED);
			//GMT_Report (API, GMT_MSG_LONG_VERBOSE, "Referencing grid data from GMT_GRID memory location\n");
			if ((G_obj = S_obj->resource) == NULL) return_null (API, GMT_PTR_IS_NULL);
			done = (mode & GMT_GRID_HEADER_ONLY) ? false : true;	/* Not done until we read grid */
			//GMT_Report (API, GMT_MSG_DEBUG, "GMTAPI_Import_Grid: Change alloc mode\n");
			G_obj->alloc_mode = G_obj->alloc_mode;
			//GMT_Report (API, GMT_MSG_DEBUG, "GMTAPI_Import_Grid: Check pad\n");
			GMT_BC_init (API->GMT, G_obj->header);	/* Initialize grid interpolation and boundary condition parameters */
			//if (GMT_err_pass (API->GMT, GMT_grd_BC_set (API->GMT, G_obj, GMT_IN), "Grid memory")) 
				if( GMT_grd_BC_set (API->GMT, G_obj, GMT_IN) !=0) //nishita
				return_null (API, GMT_GRID_BC_ERROR);	/* Set boundary conditions */
			if (!GMTAPI_adjust_grdpadding (G_obj->header, API->GMT->current.io.pad)) break;	/* Pad is correct so we are done */
			/* Here we extend G_obj->data to allow for padding, then rearrange rows */
			if (G_obj->alloc_mode == GMT_ALLOCATED_EXTERNALLY) return_null (API, GMT_PADDING_NOT_ALLOWED);
			//GMT_Report (API, GMT_MSG_DEBUG, "GMTAPI_Import_Grid: Add pad\n");
			GMT_grd_pad_on (API->GMT, G_obj, API->GMT->current.io.pad);
			//GMT_Report (API, GMT_MSG_DEBUG, "GMTAPI_Import_Grid: Return from GMT_IS_REFERENCE\n");
			break;

	 	case GMT_IS_DUPLICATE + GMT_VIA_MATRIX:	/* The user's 2-D grid array of some sort, + info in the args [NOT YET FULLY TESTED] */
			if ((M_obj = S_obj->resource) == NULL) return_null (API, GMT_PTR_IS_NULL);
			if (S_obj->region) return_null (API, GMT_SUBSET_NOT_ALLOWED);
			G_obj = (grid == NULL) ? GMT_create_grid (API->GMT) : grid;	/* Only allocate when not already allocated */
			G_obj->header->complex_mode = mode;	/* Set the complex mode */
			if (! (mode & GMT_GRID_DATA_ONLY)) {
				GMTAPI_info_to_grdheader (API->GMT, G_obj->header, M_obj);	/* Populate a GRD header structure */
				if (mode & GMT_GRID_HEADER_ONLY) break;	/* Just needed the header */
			}
			G_obj->alloc_mode = GMT_ALLOCATED_BY_GMT;
			/* Must convert to new array */
			//GMT_Report (API, GMT_MSG_LONG_VERBOSE, "Importing grid data from user memory location\n");
			GMT_set_grddim (API->GMT, G_obj->header);	/* Set all dimensions */
			G_obj->data = GMT_memory_aligned (API->GMT, NULL, G_obj->header->size, float);
			GMT_2D_to_index = GMTAPI_get_2D_to_index (API, M_obj->shape, GMT_GRID_IS_REAL);
			GMT_grd_loop (API->GMT, G_obj, row, col, ij) {
				ij_orig = GMT_2D_to_index (row, col, M_obj->dim);
				G_obj->data[ij] = (float)GMTAPI_get_val (API, &(M_obj->data), ij_orig, M_obj->type);
			}
			GMT_BC_init (API->GMT, G_obj->header);	/* Initialize grid interpolation and boundary condition parameters */
			
			//if (GMT_err_pass (API->GMT, GMT_grd_BC_set (API->GMT, G_obj, GMT_IN), "Grid memory")) 
			if( GMT_grd_BC_set (API->GMT, G_obj, GMT_IN) != 0) //nishita
				return_null (API, GMT_GRID_BC_ERROR);	/* Set boundary conditions */

			
			new_ID = GMT_Register_IO (API, GMT_IS_GRID, GMT_IS_DUPLICATE, S_obj->geometry, GMT_IN, NULL, G_obj);	/* Register a new resource to hold G_obj */
			if ((new_item = GMTAPI_Validate_ID (API, GMT_IS_GRID, new_ID, GMT_IN)) == GMT_NOTSET) return_null (API, GMT_NOTSET);	/* Some internal error... */
			API->object[new_item]->data = G_obj;
			API->object[new_item]->status = GMT_IS_USED;	/* Mark as read */
			G_obj->alloc_level = API->object[new_item]->alloc_level;	/* Since allocated here */
			via = true;
			break;

	 	case GMT_IS_REFERENCE + GMT_VIA_MATRIX:	/* The user's 2-D grid array of some sort, + info in the args [NOT YET FULLY TESTED] */
			if ((M_obj = S_obj->resource) == NULL) return_null (API, GMT_PTR_IS_NULL);
			if (S_obj->region) return_null (API, GMT_SUBSET_NOT_ALLOWED);
			G_obj = (grid == NULL) ? GMT_create_grid (API->GMT) : grid;	/* Only allocate when not already allocated */
			G_obj->header->complex_mode = mode;	/* Set the complex mode */
			if (! (mode & GMT_GRID_DATA_ONLY)) {
				GMTAPI_info_to_grdheader (API->GMT, G_obj->header, M_obj);	/* Populate a GRD header structure */
				if (mode & GMT_GRID_HEADER_ONLY) break;	/* Just needed the header */
			}
			if (!(M_obj->shape == GMT_IS_ROW_FORMAT && M_obj->type == GMT_FLOAT && M_obj->alloc_mode == GMT_ALLOCATED_EXTERNALLY && (mode & GMT_GRID_IS_COMPLEX_MASK)))
				 return_null (API, GMT_NOT_A_VALID_IO_ACCESS);
			//GMT_Report (API, GMT_MSG_LONG_VERBOSE, "Referencing grid data from user memory location\n");
			G_obj->data = M_obj->data.f4;
			S_obj->alloc_mode = M_obj->alloc_mode;	/* Pass on alloc_mode of matrix */
			G_obj->alloc_mode = M_obj->alloc_mode;
			GMT_BC_init (API->GMT, G_obj->header);	/* Initialize grid interpolation and boundary condition parameters */
			//if (GMT_err_pass (API->GMT, GMT_grd_BC_set (API->GMT, G_obj, GMT_IN), "Grid memory")) 
			if(GMT_grd_BC_set (API->GMT, G_obj, GMT_IN) !=0) //nishita
				return_null (API, GMT_GRID_BC_ERROR);	/* Set boundary conditions */
			if (!GMTAPI_adjust_grdpadding (G_obj->header, API->GMT->current.io.pad)) break;	/* Pad is correct so we are done */
			if (G_obj->alloc_mode == GMT_ALLOCATED_EXTERNALLY) return_null (API, GMT_PADDING_NOT_ALLOWED);
			/* Here we extend G_obj->data to allow for padding, then rearrange rows */
			GMT_grd_pad_on (API->GMT, G_obj, API->GMT->current.io.pad);
			new_ID = GMT_Register_IO (API, GMT_IS_GRID, GMT_IS_REFERENCE, S_obj->geometry, GMT_IN, NULL, G_obj);	/* Register a new resource to hold G_obj */
			if ((new_item = GMTAPI_Validate_ID (API, GMT_IS_GRID, new_ID, GMT_IN)) == GMT_NOTSET) return_null (API, GMT_NOTSET);	/* Some internal error... */
			API->object[new_item]->data = G_obj;
			API->object[new_item]->status = GMT_IS_USED;	/* Mark as read */
			G_obj->alloc_level = API->object[new_item]->alloc_level;	/* Since allocated here */
			via = true;
			break;

		default:
			//GMT_Report (API, GMT_MSG_NORMAL, "Wrong method used to import grid\n");
			return_null (API, GMT_NOT_A_VALID_METHOD);
			break;
	}

	if (done) S_obj->status = GMT_IS_USED;	/* Mark as read (unless we just got the header) */
	if (!via) S_obj->data = G_obj;		/* Retain pointer to the allocated data so we use garbage collection later */

	return (G_obj);	/* Pass back out what we have so far */
}


void * GMTAPI_Import_Data (struct GMTAPI_CTRL *API, enum GMT_enum_family family, int object_ID, unsigned int mode, void *data)
{
	/* Function that will import the data object referred to by the object_ID (or all registered inputs if object_ID == GMT_NOTSET).
	 * This is a wrapper functions for CPT, Dataset, Textset, Grid and Image imports; see the specific functions
	 * for details on the arguments, in particular the mode setting (or see the GMT API documentation).
	 */
	int item;
	void *new_obj = NULL;

	if (API == NULL) return_null (API, GMT_NOT_A_SESSION);			/* GMT_Create_Session has not been called */
	if (!API->registered[GMT_IN]) return_null (API, GMT_NO_INPUT);		/* No sources registered yet */

	/* Get information about this resource first */
	if ((item = GMTAPI_Validate_ID (API, family, object_ID, GMT_IN)) == GMT_NOTSET) return_null (API, API->error);

	/* The case where object_ID is not set but a virtual (memory) file is found is a special case: we must supply the correct object_ID */
	if (object_ID == GMT_NOTSET && item && API->object[item]->method != GMT_IS_FILE) object_ID = API->object[item]->ID;	/* Found virtual file; set actual object_ID */

	switch (family) {	/* CPT, Dataset, or Grid */
		case GMT_IS_CPT:
			//new_obj = GMTAPI_Import_CPT (API, object_ID, mode);			/* Try to import a CPT */
			//break;
		case GMT_IS_DATASET:
			//new_obj = GMTAPI_Import_Dataset (API, object_ID, mode);		/* Try to import data tables */
			break;
		case GMT_IS_TEXTSET:
			//new_obj = GMTAPI_Import_Textset (API, object_ID, mode);		/* Try to import text tables */
			break;
		case GMT_IS_GRID:
			new_obj = GMTAPI_Import_Grid (API, object_ID, mode, data);		/* Try to import a grid */
			break;
//#ifdef HAVE_GDAL
	//	case GMT_IS_IMAGE:
		//	new_obj = GMTAPI_Import_Image (API, object_ID, mode, data);		/* Try to import a image */
			//break;
//#endif
		default:
			API->error = GMT_NOT_A_VALID_FAMILY;
			break;
	}
	if (new_obj == NULL) return_null (API, API->error);	/* Return NULL as something went wrong */
	return (new_obj);	/* Successful, return pointer */
}

void * GMT_Get_Data (void *V_API, int object_ID, unsigned int mode, void *data)
{
	/* Function to import registered data sources directly into program memory as a set (not record-by-record).
	 * data is pointer to an existing grid container when we read a grid in two steps, otherwise use NULL.
	 * ID is the registered resource from which to import.
	 * Return: Pointer to data container, or NULL if there were errors (passed back via API->error).
	 */
	int item, family;
	bool was_enabled;
	void *new_obj = NULL;
	struct GMTAPI_CTRL *API = NULL;

	if (V_API == NULL) return_null (V_API, GMT_NOT_A_SESSION);

	/* Determine the item in the object list that matches this ID and direction */
	API = gmt_get_api_ptr ((struct GMTAPI_CTRL *)V_API);
	if (object_ID == GMT_NOTSET) {	/* Must pick up the family from the shelf */
		family = API->shelf;
		API->shelf = GMT_NOTSET;
	}
	else
		family = GMT_NOTSET;
	if ((item = GMTAPI_Validate_ID (API, family, object_ID, GMT_IN)) == GMT_NOTSET) {
		return_null (API, API->error);
	}

	was_enabled = API->io_enabled[GMT_IN];
	if (!was_enabled && GMTAPI_Begin_IO (API, GMT_IN) != GMT_OK) {	/* Enables data input if not already set and sets access mode */
		return_null (API, API->error);
	}
	API->object[item]->selected = true;	/* Make sure it the requested data set is selected */

	/* OK, try to do the importing */
	if ((new_obj = GMTAPI_Import_Data (API, API->object[item]->family, object_ID, mode, data)) == NULL) {
		return_null (API, API->error);
	}

	if (!was_enabled && GMT_End_IO (API, GMT_IN, 0) != GMT_OK) {	/* Disables data input if we had to set it in this function */
		return_null (API, API->error);
	}
//#ifdef DEBUG
	//GMTAPI_Set_Object (API, API->object[item]);
	//GMT_list_API (API, "GMT_Get_Data");
//#endif
	return (new_obj);		/* Return pointer to the data container */
}

void * GMT_Read_Data (void *V_API, unsigned int family, unsigned int method, unsigned int geometry, unsigned int mode, double wesn[], char *input, void *data)
{
	/* Function to read data files directly into program memory as a set (not record-by-record).
	 * We can combine the <register resource - import resource > sequence in
	 * one combined function.  See GMT_Register_IO for details on arguments.
	 * data is pointer to an existing grid container when we read a grid in two steps, otherwise use NULL.
	 * Case 1: input != NULL: Register input as the source and import data.
	 * Case 2: input == NULL: Register stdin as the source and import data.
	 * Case 3: geometry == 0: Loop over all previously registered AND unread sources and combine as virtual dataset [DATASET|TEXTSET only]
	 * Case 4: family is GRID|IMAGE and method = GMT_GRID_DATA_ONLY: Just find already registered resource
	 * Return: Pointer to data container, or NULL if there were errors (passed back via API->error).
	 */
	int in_ID = GMT_NOTSET, item;
	bool just_get_data;
	void *new_obj = NULL;
	struct GMTAPI_CTRL *API = NULL;

	if (V_API == NULL) return_null (V_API, GMT_NOT_A_SESSION);

	API = gmt_get_api_ptr ((struct GMTAPI_CTRL *)V_API);
	just_get_data = GMT_File_Is_Memory (input);

	if ((family == GMT_IS_GRID || family == GMT_IS_IMAGE) && (mode & GMT_GRID_DATA_ONLY)) {	/* Case 4: Already registered when we obtained header, find object ID */
		if ((in_ID = GMTAPI_is_registered (API, (enum GMT_enum_family)family, geometry, GMT_IN, mode, input, data)) == GMT_NOTSET) 
			{
			
			return_null (API, GMT_OBJECT_NOT_FOUND);	/* Could not find it */
			}
		if (!full_region (wesn)) {	/* Must update subset selection */
			int item;
			if ((item = GMTAPI_Validate_ID (API, family, in_ID, GMT_IN)) == GMT_NOTSET) 
				{
				
				return_null (API, API->error);
				}
			GMT_memcpy (API->object[item]->wesn, wesn, 4, double);
		}
	}
	else if (input) {	/* Case 1: Load from a single, given source. Register it first. */
		if ((in_ID = GMT_Register_IO (API, family, method, geometry, GMT_IN, wesn, input)) == GMT_NOTSET) 
			{
			
			return_null (API, API->error);
			}
	}
	else if (input == NULL && geometry) {	/* Case 2: Load from stdin.  Register stdin first */
		if ((in_ID = GMT_Register_IO (API, family, GMT_IS_STREAM, geometry, GMT_IN, wesn, API->GMT->session.std[GMT_IN])) == GMT_NOTSET){
			

			return_null (API, API->error);	/* Failure to register std??? */
			}
	}
	else {	/* Case 3: input == NULL && geometry == 0, so use all previously registered sources (unless already used). */
		if (!(family == GMT_IS_DATASET || family == GMT_IS_TEXTSET)) return_null (API, GMT_ONLY_ONE_ALLOWED);	/* Virtual source only applies to data and text tables */
		API->shelf = family;	/* Save which one it is so we know in GMT_Get_Data */
	}
	if (just_get_data) {
		if ((item = GMTAPI_Validate_ID (API, GMT_NOTSET, in_ID, GMT_NOTSET)) == GMT_NOTSET) {
			return_null (API, API->error);
		}
#ifdef DEBUG
		//GMTAPI_Set_Object (API, API->object[item]);
#endif
		return ((API->object[item]->data) ? API->object[item]->data : API->object[item]->resource);	/* Return pointer to the data */
	}

	/* OK, try to do the importing */
	if (in_ID != GMT_NOTSET) {	/* Make sure we select the item we just registered */
		if ((item = GMTAPI_Validate_ID (API, GMT_NOTSET, in_ID, GMT_NOTSET)) == GMT_NOTSET) {
				return_null (API, API->error);
		}
		API->object[item]->selected = true;	/* Make sure the item we want is now selected */
	}
	if ((new_obj = GMT_Get_Data (API, in_ID, mode, data)) == NULL) return_null (API, API->error);

#ifdef DEBUG
	//GMT_list_API (API, "GMT_Read_Data");
#endif

	return (new_obj);		/* Return pointer to the data container */
}

void GMTAPI_increment_D (struct GMT_DATASET *D_obj, uint64_t n_rows, uint64_t n_columns)
{	/* Increment dimensions for this single dataset/segment */
	D_obj->table[D_obj->n_tables]->segment[0]->n_rows = n_rows;
	D_obj->table[D_obj->n_tables]->segment[0]->n_columns = D_obj->table[D_obj->n_tables]->n_columns = n_columns;
	D_obj->table[D_obj->n_tables]->n_records += n_rows;
	D_obj->table[D_obj->n_tables]->n_segments = 1;
	D_obj->n_tables++;	/* Since we just read one table */
}

struct GMT_DATASET * GMTAPI_Import_Dataset (struct GMTAPI_CTRL *API, int object_ID, unsigned int mode)
{	/* Does the actual work of loading in the entire virtual data set (possibly via many sources)
	 * If object_ID == GMT_NOTSET we get all registered input tables, otherwise we just get the one requested.
	 * Note: Memory is allocated for the Dataset except for method GMT_IS_REFERENCE.
	 */

	int item, first_item = 0, this_item = GMT_NOTSET, last_item, new_item, new_ID;
	unsigned int geometry;
	bool allocate = false, update = false, all_D, use_GMT_io, greenwich = true, via;
	size_t n_alloc;
	uint64_t row, seg, col, ij;
	p_func_size_t GMT_2D_to_index = NULL;

	struct GMT_DATASET *D_obj = NULL, *Din_obj = NULL;
	struct GMT_MATRIX *M_obj = NULL;
	struct GMT_VECTOR *V_obj = NULL;
	struct GMTAPI_DATA_OBJECT *S_obj = NULL;

	

	if (object_ID == GMT_NOTSET) {	/* Means there is more than one source: Merge all registered data tables into a single virtual data set */
		last_item = API->n_objects - 1;	/* Must check all objects */
		allocate = true;
		n_alloc = GMT_TINY_CHUNK;
	}
	else {		/* Requested a single, specific data table */
		if ((first_item = GMTAPI_Validate_ID (API, GMT_IS_DATASET, object_ID, GMT_IN)) == GMT_NOTSET) return_null (API, API->error);
		last_item = first_item;
		n_alloc = 1;
	}

	/* Allocate data set and an initial list of tables */
	D_obj = GMT_memory (API->GMT, NULL, 1, struct GMT_DATASET);
	D_obj->table = GMT_memory (API->GMT, NULL, n_alloc, struct GMT_DATATABLE *);
	D_obj->alloc_mode = GMT_ALLOCATED_BY_GMT;	/* So GMT_* modules can free this memory (may override below) */
	D_obj->alloc_level = API->GMT->hidden.func_level;	/* So GMT_* modules can free this memory (may override below) */
	use_GMT_io = !(mode & GMT_IO_ASCII);		/* false if we insist on ASCII reading */
	API->GMT->current.io.seg_no = API->GMT->current.io.rec_no = API->GMT->current.io.rec_in_tbl_no = 0;	/* Reset for each new dataset */
	if (API->GMT->common.R.active && API->GMT->common.R.wesn[XLO] < -180.0 && API->GMT->common.R.wesn[XHI] > -180.0) greenwich = false;

	for (item = first_item; item <= last_item; item++) {	/* Look through all sources for registered inputs (or just one) */
		S_obj = API->object[item];	/* S_obj is the current data object */
		if (!S_obj) {	/* Probably not a good sign */
			printf(/*API, GMT_MSG_DEBUG,*/ "GMTAPI_Import_Dataset: Skipped empty object (item = %d)\n", item);
			continue;
		}
		if (!S_obj->selected) continue;			/* Registered, but not selected */
		if (S_obj->direction == GMT_OUT) continue;	/* We're doing reading here, so skip output objects */
		if (S_obj->family != GMT_IS_DATASET) continue;	/* We're doing datasets here, so skip other data types */
		if (S_obj->status != GMT_IS_UNUSED) { 	/* Already read this resource before; are we allowed to re-read? */
			if (S_obj->method == GMT_IS_STREAM || S_obj->method == GMT_IS_FDESC)
				//return_null (API, GMT_READ_ONCE);	/* Not allowed to re-read streams */
				return (GMT_READ_ONCE);//nishita
			if (!(mode & GMT_IO_RESET))
				//return_null (API, GMT_READ_ONCE);	/* Not authorized to re-read */
				return (GMT_READ_ONCE);//nishita
		}
		if (this_item == GMT_NOTSET) this_item = item;	/* First item that worked */
		via = false;
		geometry = (API->GMT->common.a.output) ? API->GMT->common.a.geometry : S_obj->geometry;	/* When reading GMT and writing OGR/GMT we must make sure we set this first */
		switch (S_obj->method) {	/* File, array, stream etc ? */
	 		case GMT_IS_FILE:	/* Import all the segments, then count total number of records */
//#ifdef SET_IO_MODE
	//			if (item == first_item) GMT_setmode (API->GMT, GMT_IN);	/* Windows may need to switch read mode from text to binary */
//#endif
				/* GMT_read_table will report where it is reading from if level is GMT_MSG_LONG_VERBOSE */
				if (API->GMT->current.io.ogr == GMT_OGR_TRUE && D_obj->n_tables > 0)	/* Only single tables if GMT/OGR */
					return_null (API, GMT_OGR_ONE_TABLE_ONLY);
				printf (/*API, GMT_MSG_LONG_VERBOSE,*/ "Reading %s from %s %s\n", GMT_family[S_obj->family], GMT_method[S_obj->method], S_obj->filename);
				if ((D_obj->table[D_obj->n_tables] = GMT_read_table (API->GMT, S_obj->filename, S_obj->method, greenwich, &geometry, use_GMT_io)) == NULL) continue;		/* Ran into an empty file (e.g., /dev/null or equivalent). Skip to next item, */
				D_obj->table[D_obj->n_tables]->id = D_obj->n_tables;	/* Give sequential internal object_ID numbers to tables */
				D_obj->n_tables++;	/* Since we just read one */
				update = true;
				break;

			case GMT_IS_STREAM:	/* Import all the segments, then count total number of records */
	 		case GMT_IS_FDESC:
				/* GMT_read_table will report where it is reading from if level is GMT_MSG_LONG_VERBOSE */
//#ifdef SET_IO_MODE
	//			if (item == first_item) GMT_setmode (API->GMT, GMT_IN);	/* Windows may need to switch read mode from text to binary */
//#endif
				if (API->GMT->current.io.ogr == GMT_OGR_TRUE && D_obj->n_tables > 0)	/* Only single tables if GMT/OGR */
					return_null (API, GMT_OGR_ONE_TABLE_ONLY);
				printf (/*API, GMT_MSG_LONG_VERBOSE,*/ "Reading %s from %s "/*%" PRIxS*/ "\n", GMT_family[S_obj->family], GMT_method[S_obj->method], (size_t)S_obj->fp);
				if ((D_obj->table[D_obj->n_tables] = GMT_read_table (API->GMT, S_obj->fp, S_obj->method, greenwich, &geometry, use_GMT_io)) == NULL) continue;		/* Ran into an empty file (e.g., /dev/null or equivalent). Skip to next item, */
				D_obj->table[D_obj->n_tables]->id = D_obj->n_tables;	/* Give sequential internal object_ID numbers to tables */
				D_obj->n_tables++;	/* Since we just read one */
				update = true;
				break;

			case GMT_IS_DUPLICATE:	/* Duplicate the input dataset */
				if (n_alloc > 1) return_null (API, GMT_ONLY_ONE_ALLOWED);
				if ((Din_obj = S_obj->resource) == NULL) return_null (API, GMT_PTR_IS_NULL);
				printf (/*API, GMT_MSG_LONG_VERBOSE,*/ "Duplicating data table from GMT_DATASET memory location\n");
				GMT_free (API->GMT, D_obj->table);	/* Free up what we allocated earlier since GMT_alloc_dataset does it all */
				GMT_free (API->GMT, D_obj);
				D_obj = GMT_duplicate_dataset (API->GMT, Din_obj, GMT_ALLOC_NORMAL, NULL);
				break;

			case GMT_IS_REFERENCE:	/* Just pass memory location, so free up what we allocated first */
				if (n_alloc > 1)
					//return_null (API, GMT_ONLY_ONE_ALLOWED);
					return (GMT_ONLY_ONE_ALLOWED);//nishita
				printf (/*API, GMT_MSG_LONG_VERBOSE,*/ "Referencing data table from GMT_DATASET memory location\n");
				GMT_free (API->GMT, D_obj->table);	/* Free up what we allocated up front since we just wish to pass the pointer */
				GMT_free (API->GMT, D_obj);
				if ((D_obj = S_obj->resource) == NULL)
					//return_null (API, GMT_PTR_IS_NULL);
					return (GMT_PTR_IS_NULL);//nishita
				break;

	 		case GMT_IS_DUPLICATE + GMT_VIA_MATRIX:
				/* Each array source becomes a separate table with a single segment */
				if ((M_obj = S_obj->resource) == NULL) return_null (API, GMT_PTR_IS_NULL);
				printf (/*API, GMT_MSG_LONG_VERBOSE,*/ "Duplicating data table from user array location\n");
				D_obj->table[D_obj->n_tables] = GMT_memory (API->GMT, NULL, 1, struct GMT_DATATABLE);
				D_obj->table[D_obj->n_tables]->segment = GMT_memory (API->GMT, NULL, 1, struct GMT_DATASEGMENT *);
				D_obj->table[D_obj->n_tables]->segment[0] = GMT_memory (API->GMT, NULL, 1, struct GMT_DATASEGMENT);
				GMT_alloc_segment (API->GMT, D_obj->table[D_obj->n_tables]->segment[0], M_obj->n_rows, M_obj->n_columns, true);
				GMT_2D_to_index = GMTAPI_get_2D_to_index (API, M_obj->shape, GMT_GRID_IS_REAL);
				for (row = 0; row < M_obj->n_rows; row++) {
					for (col = 0; col < M_obj->n_columns; col++) {
						ij = GMT_2D_to_index (row, col, M_obj->dim);
						D_obj->table[D_obj->n_tables]->segment[0]->coord[col][row] = GMTAPI_get_val (API, &(M_obj->data), ij, M_obj->type);
					}
				}
				GMTAPI_increment_D (D_obj, M_obj->n_rows, M_obj->n_columns);	/* Update counters for D_obj */
				new_ID = GMT_Register_IO (API, GMT_IS_DATASET, GMT_IS_DUPLICATE, geometry, GMT_IN, NULL, D_obj);	/* Register a new resource to hold D_obj */
				if ((new_item = GMTAPI_Validate_ID (API, GMT_IS_DATASET, new_ID, GMT_IN)) == GMT_NOTSET) return_null (API, GMT_NOTSET);	/* Some internal error... */
				API->object[new_item]->data = D_obj;
				API->object[new_item]->status = GMT_IS_USED;	/* Mark as read */
				D_obj->alloc_level = API->object[new_item]->alloc_level;	/* Since allocated here */
				update = via = true;
				break;

	 		case GMT_IS_DUPLICATE + GMT_VIA_VECTOR:
				/* Each column array source becomes column arrays in a separate table with a single segment */
				if ((V_obj = S_obj->resource) == NULL) return_null (API, GMT_PTR_IS_NULL);
				printf (/*API, GMT_MSG_LONG_VERBOSE,*/ "Duplicating data table from user "/*%" PRIu64*/ " column arrays of length  "/*%" PRIu64*/ "\n", V_obj->n_columns, V_obj->n_rows);
				D_obj->table[D_obj->n_tables] = GMT_memory (API->GMT, NULL, 1, struct GMT_DATATABLE);
				D_obj->table[D_obj->n_tables]->segment = GMT_memory (API->GMT, NULL, 1, struct GMT_DATASEGMENT *);
				D_obj->table[D_obj->n_tables]->segment[0] = GMT_memory (API->GMT, NULL, 1, struct GMT_DATASEGMENT);
				GMT_alloc_segment (API->GMT, D_obj->table[D_obj->n_tables]->segment[0], V_obj->n_rows, V_obj->n_columns, true);
				for (col = 0, all_D = true; all_D && col < V_obj->n_columns; col++) if (V_obj->type[col] != GMT_DOUBLE) all_D = false;
				if (all_D) {	/* Can use fast memcpy */
					for (col = 0; col < V_obj->n_columns; col++)
						GMT_memcpy (D_obj->table[D_obj->n_tables]->segment[0]->coord[col], V_obj->data[col].f8, V_obj->n_rows, double);
				}
				else {	/* Must copy items individually */
					for (row = 0; row < V_obj->n_rows; row++) {
						for (col = 0; col < V_obj->n_columns; col++)
							D_obj->table[D_obj->n_tables]->segment[0]->coord[col][row] = GMTAPI_get_val (API, &(V_obj->data[col]), row, V_obj->type[col]);
					}
				}
				GMTAPI_increment_D (D_obj, V_obj->n_rows, V_obj->n_columns);	/* Update counters for D_obj */
				new_ID = GMT_Register_IO (API, GMT_IS_DATASET, GMT_IS_DUPLICATE, geometry, GMT_IN, NULL, D_obj);	/* Register a new resource to hold D_obj */
				if ((new_item = GMTAPI_Validate_ID (API, GMT_IS_DATASET, new_ID, GMT_IN)) == GMT_NOTSET) return_null (API, GMT_NOTSET);	/* Some internal error... */
				API->object[new_item]->data = D_obj;
				API->object[new_item]->status = GMT_IS_USED;	/* Mark as read */
				D_obj->alloc_level = API->object[new_item]->alloc_level;	/* Since allocated here */
				update = via = true;
				break;

		 	case GMT_IS_REFERENCE + GMT_VIA_VECTOR:
				if ((V_obj = S_obj->resource) == NULL) return_null (API, GMT_PTR_IS_NULL);
				if (V_obj->type[0] != GMT_DOUBLE) return_null (API, GMT_NOT_A_VALID_TYPE);
				/* Each column array source becomes preallocated column arrays in a separate table with a single segment */
				printf (/*API, GMT_MSG_LONG_VERBOSE,*/ "Referencing data table from user "/*%" PRIu64*/ " column arrays of length"/* %" PRIu64*/ "\n", V_obj->n_columns, V_obj->n_rows);
				D_obj->table[D_obj->n_tables] = GMT_memory (API->GMT, NULL, 1, struct GMT_DATATABLE);
				D_obj->table[D_obj->n_tables]->segment = GMT_memory (API->GMT, NULL, 1, struct GMT_DATASEGMENT *);
				D_obj->table[D_obj->n_tables]->segment[0] = GMT_memory (API->GMT, NULL, 1, struct GMT_DATASEGMENT);
				GMT_alloc_segment (API->GMT, D_obj->table[D_obj->n_tables]->segment[0], 0, V_obj->n_columns, true);
				for (col = 0; col < V_obj->n_columns; col++)
					D_obj->table[D_obj->n_tables]->segment[0]->coord[col] = V_obj->data[col].f8;
				GMTAPI_increment_D (D_obj, V_obj->n_rows, V_obj->n_columns);	/* Update counters for D_obj */
				D_obj->alloc_mode = GMT_ALLOCATED_EXTERNALLY;	/* Since we just hooked on the arrays */
				new_ID = GMT_Register_IO (API, GMT_IS_DATASET, GMT_IS_REFERENCE, geometry, GMT_IN, NULL, D_obj);	/* Register a new resource to hold D_obj */
				if ((new_item = GMTAPI_Validate_ID (API, GMT_IS_DATASET, new_ID, GMT_IN)) == GMT_NOTSET) return_null (API, GMT_NOTSET);	/* Some internal error... */
				API->object[new_item]->data = D_obj;
				API->object[new_item]->status = GMT_IS_USED;	/* Mark as read */
				D_obj->alloc_level = API->object[new_item]->alloc_level;	/* Since allocated here */
				S_obj->family = GMT_IS_VECTOR;	/* Done with the via business now */
				update = via = true;
				break;

			default:	/* Barking up the wrong tree here... */
				printf (/*API, GMT_MSG_NORMAL, */ "Wrong method used to import data tables\n");
				GMT_free (API->GMT, D_obj->table);
				GMT_free (API->GMT, D_obj);
				return_null (API, GMT_NOT_A_VALID_METHOD);
				break;
		}
		if (update) {	/* Total up statistics */
			D_obj->n_segments += D_obj->table[D_obj->n_tables-1]->n_segments;	/* Sum up total number of segments across the data set */
			D_obj->n_records += D_obj->table[D_obj->n_tables-1]->n_records;	/* Sum up total number of records across the data set */
			/* Update segment IDs so they are sequential across many tables (GMT_read_table sets the ids relative to current table). */
			if (D_obj->n_tables > 1) {
				for (seg = 0; seg < D_obj->table[D_obj->n_tables-1]->n_segments; seg++)
					D_obj->table[D_obj->n_tables-1]->segment[seg]->id += D_obj->table[D_obj->n_tables-2]->n_segments;
			}
			if (allocate && D_obj->n_tables == n_alloc) {	/* Must allocate space for more tables */
				size_t old_n_alloc = n_alloc;
				n_alloc += GMT_TINY_CHUNK;
				D_obj->table = GMT_memory (API->GMT, D_obj->table, n_alloc, struct GMT_DATATABLE *);
				GMT_memset (&(D_obj->table[old_n_alloc]), n_alloc - old_n_alloc, struct GMT_DATATABLE *);	/* Set to NULL */
			}
		}
		S_obj->alloc_mode = D_obj->alloc_mode;	/* Clarify allocation mode for this entity */
#if 0
		if (col_check (D_obj->table[D_obj->n_tables-1], &n_cols)) {	/* Different tables have different number of columns, which is not good */
			return_null (API, GMT_N_COLS_VARY);
		}
#endif
		S_obj->status = GMT_IS_USED;	/* Mark as read */
	}
	if (D_obj->n_tables == 0) {	/* Only found empty files (e.g., /dev/null) and we have nothing to show for our efforts.  Return an single empty table with no segments. */
		D_obj->table = GMT_memory (API->GMT, D_obj->table, 1, struct GMT_DATATABLE *);
		D_obj->table[0] = GMT_memory (API->GMT, NULL, 1, struct GMT_DATATABLE);
		D_obj->n_tables = 1;	/* But we must indicate we found one (empty) table */
	}
	else {	/* Found one or more tables */
		if (allocate && D_obj->n_tables < n_alloc) D_obj->table = GMT_memory (API->GMT, D_obj->table, D_obj->n_tables, struct GMT_DATATABLE *);
		D_obj->n_columns = D_obj->table[0]->n_columns;
		if (!D_obj->min) D_obj->min = GMT_memory (API->GMT, NULL, D_obj->n_columns, double);
		if (!D_obj->max) D_obj->max = GMT_memory (API->GMT, NULL, D_obj->n_columns, double);
	}
	GMT_set_dataset_minmax (API->GMT, D_obj);	/* Set the min/max values for the entire dataset */
	D_obj->geometry = geometry;	/* Since GMT_read_table may have updated it */
	if (!via) API->object[this_item]->data = D_obj;	/* Retain pointer to the allocated data so we use garbage collection later */
	return (D_obj);
}

/*struct GMT_TEXTSET * GMTAPI_Import_Textset (struct GMTAPI_CTRL *API, int object_ID, unsigned int mode)
{	 Does the actual work of loading in the entire virtual text set (possibly via many sources)
	 * If object_ID == GMT_NOTSET we get all registered input tables, otherwise we just get the one requested.
	 * Note: Memory is allocated for the Dataset except for GMT_IS_REFERENCE.


	int item, first_item = 0, last_item, this_item = GMT_NOTSET, new_item, new_ID;
	bool update = false, allocate = false, via;
	size_t n_alloc;
	uint64_t row, seg;
	char *t_ptr = NULL;
	struct GMT_TEXTSET *T_obj = NULL;
	struct GMT_MATRIX *M_obj = NULL;
	struct GMTAPI_DATA_OBJECT *S_obj = NULL;

	printf(API, GMT_MSG_DEBUG, "GMTAPI_Import_Textset: Passed ID = %d and mode = %d\n", object_ID, mode);

	T_obj = GMT_memory (API->GMT, NULL, 1, struct GMT_TEXTSET);

	if (object_ID == GMT_NOTSET) {	 More than one source: Merge all registered data tables into a single virtual text set
		last_item = API->n_objects - 1;	 Must check all objects
		allocate  = true;
		n_alloc = GMT_TINY_CHUNK;
	}
	else {		 Requested a single, specific data table
		if ((first_item = GMTAPI_Validate_ID (API, GMT_IS_TEXTSET, object_ID, GMT_IN)) == GMT_NOTSET) return_null (API, API->error);
		last_item  = first_item;
		n_alloc = 1;
	}
	T_obj->table = GMT_memory (API->GMT, NULL, n_alloc, struct GMT_TEXTTABLE *);
	T_obj->alloc_mode = GMT_ALLOCATED_BY_GMT;	 So GMT_* modules can free this memory (may override below)
	T_obj->alloc_level = API->GMT->hidden.func_level;	 So GMT_* modules can free this memory (may override below)

	for (item = first_item; item <= last_item; item++) {	 Look through all sources for registered inputs (or just one)
		S_obj = API->object[item];	 Current object
		if (!S_obj) {	 Probably not a good sign
			printf(API, GMT_MSG_DEBUG, "GMTAPI_Import_Textset: Skipped empty object (item = %d)\n", item);
			continue;
		}
		if (!S_obj->selected) continue;			 Registered, but not selected
		if (S_obj->direction == GMT_OUT) continue;	 We're doing reading here, so bugger off!
		if (S_obj->family != GMT_IS_TEXTSET) continue;	 We're doing textsets here, so skip other things
		if (S_obj->status != GMT_IS_UNUSED) {	 Already read this resource before; are we allowed to re-read?
			if (S_obj->method == GMT_IS_STREAM || S_obj->method == GMT_IS_FDESC) return_null (API, GMT_READ_ONCE);	 Not allowed to re-read streams
			if (!(mode & GMT_IO_RESET)) return_null (API, GMT_READ_ONCE);	 Not authorized to re-read
		}
		if (this_item == GMT_NOTSET) this_item = item;	 First item that worked
		via = false;
		switch (S_obj->method) {	 File, array, stream etc ?
			case GMT_IS_FILE:	 Import all the segments, then count total number of records
//#ifdef SET_IO_MODE
	//			if (item == first_item) GMT_setmode (API->GMT, GMT_IN);
//#endif
				 GMT_read_texttable will report where it is reading from if level is GMT_MSG_LONG_VERBOSE
				printf (API, GMT_MSG_LONG_VERBOSE, "Reading %s from %s %s\n", GMT_family[S_obj->family], GMT_method[S_obj->method], S_obj->filename);
				if ((T_obj->table[T_obj->n_tables] = GMT_read_texttable (API->GMT, S_obj->filename, S_obj->method)) == NULL) continue;	 Ran into an empty file (e.g., /dev/null or equivalent). Skip to next item,
				T_obj->table[T_obj->n_tables]->id = T_obj->n_tables;	 Give internal object_ID numbers to tables
				update = true;
				break;
	 		case GMT_IS_STREAM:
	 		case GMT_IS_FDESC:
//#ifdef SET_IO_MODE
	//			if (item == first_item) GMT_setmode (API->GMT, GMT_IN);
//#endif
				 GMT_read_texttable will report where it is reading from if level is GMT_MSG_LONG_VERBOSE
				printf (API, GMT_MSG_LONG_VERBOSE, "Reading %s from %s %" PRIxS "\n", GMT_family[S_obj->family], GMT_method[S_obj->method], (size_t)S_obj->fp);
				if ((T_obj->table[T_obj->n_tables] = GMT_read_texttable (API->GMT, S_obj->fp, S_obj->method)) == NULL) continue;	 Ran into an empty file (e.g., /dev/null or equivalent). Skip to next item,
				T_obj->table[T_obj->n_tables]->id = T_obj->n_tables;	 Give internal object_ID numbers to tables
				update = true;
				break;
			case GMT_IS_DUPLICATE:	 Duplicate the input dataset
				if (n_alloc > 1)
					//return_null (API, GMT_ONLY_ONE_ALLOWED);
					return (GMT_ONLY_ONE_ALLOWED); //nishita
				printf (API, GMT_MSG_LONG_VERBOSE, "Duplicating text table from GMT_TEXTSET memory location\n");
				GMT_free (API->GMT, T_obj->table);	 Free up what we allocated since GMT_alloc_dataset does it all
				GMT_free (API->GMT, T_obj);
				if (S_obj->resource == NULL) return_null (API, GMT_PTR_IS_NULL);
				T_obj = GMT_duplicate_textset (API->GMT, S_obj->resource, GMT_ALLOC_NORMAL);
				break;
			case GMT_IS_REFERENCE:	 Just pass memory location, so free up what we allocated first
				if (n_alloc > 1) return_null (API, GMT_ONLY_ONE_ALLOWED);
				printf (API, GMT_MSG_LONG_VERBOSE,  "Referencing data table from GMT_TEXTSET memory location\n");
				GMT_free (API->GMT, T_obj->table);	 Free up what we allocated since GMT_alloc_textset does it all
				GMT_free (API->GMT, T_obj);
				if ((T_obj = S_obj->resource) == NULL) return_null (API, GMT_PTR_IS_NULL);
				break;
	 		case GMT_IS_DUPLICATE + GMT_VIA_MATRIX:
				 Each matrix source becomes a separate table with one segment
			 	if ((M_obj = S_obj->resource) == NULL)
			 		//return_null (API, GMT_PTR_IS_NULL);
			 		return(GMT_PTR_IS_NULL); //nishita
				printf (API, GMT_MSG_LONG_VERBOSE, "Duplicating text table from user matrix location\n");
				T_obj->table[T_obj->n_tables] = GMT_memory (API->GMT, NULL, 1, struct GMT_TEXTTABLE);
				T_obj->table[T_obj->n_tables]->segment = GMT_memory (API->GMT, NULL, 1, struct GMT_TEXTSEGMENT *);
				T_obj->table[T_obj->n_tables]->segment[0] = GMT_memory (API->GMT, NULL, 1, struct GMT_TEXTSEGMENT);
				T_obj->table[T_obj->n_tables]->segment[0]->record = GMT_memory (API->GMT, NULL, M_obj->n_rows, char *);
				t_ptr = (char *)M_obj->data.sc1;
				for (row = 0; row < (uint64_t)M_obj->n_rows; row++) {
					T_obj->table[T_obj->n_tables]->segment[0]->record[row] = strdup (&t_ptr[row*M_obj->dim]);
				}
				T_obj->table[T_obj->n_tables]->segment[0]->n_rows = M_obj->n_rows;
				T_obj->table[T_obj->n_tables]->n_records += M_obj->n_rows;
				T_obj->table[T_obj->n_tables]->n_segments = 1;
				new_ID = GMT_Register_IO (API, GMT_IS_TEXTSET, GMT_IS_DUPLICATE, S_obj->geometry, GMT_IN, NULL, T_obj);	 Register a new resource to hold T_obj
				if ((new_item = GMTAPI_Validate_ID (API, GMT_IS_DATASET, new_ID, GMT_IN)) == GMT_NOTSET) return_null (API, GMT_NOTSET);	 Some internal error...
				API->object[new_item]->data = T_obj;
				API->object[new_item]->status = GMT_IS_USED;	 Mark as read
				T_obj->alloc_level = API->object[new_item]->alloc_level;	 Since allocated here
				update = via = true;
				break;
			default:	 Barking up the wrong tree here...
				printf (API, GMT_MSG_NORMAL,  "Wrong method used to import data tables\n");
				GMT_free (API->GMT, T_obj->table);
				GMT_free (API->GMT, T_obj);
				return_null (API, GMT_NOT_A_VALID_METHOD);
				break;
		}
		if (update) {
			T_obj->n_segments += T_obj->table[T_obj->n_tables]->n_segments;	 Sum up total number of segments across the data set
			T_obj->n_records += T_obj->table[T_obj->n_tables]->n_records;	 Sum up total number of records across the data set
			 Update segment object_IDs so they are sequential across many tables (GMT_read_table sets the ids relative to current table).
			if (T_obj->n_tables > 0)
				for (seg = 0; seg < T_obj->table[T_obj->n_tables]->n_segments; seg++)
					T_obj->table[T_obj->n_tables]->segment[seg]->id += T_obj->table[T_obj->n_tables-1]->n_segments;
			T_obj->n_tables++;
		}
		if (allocate && T_obj->n_tables == n_alloc) {	 Must allocate space for more tables
			size_t old_n_alloc = n_alloc;
			n_alloc += GMT_TINY_CHUNK;
			T_obj->table = GMT_memory (API->GMT, T_obj->table, n_alloc, struct GMT_TEXTTABLE *);
			GMT_memset (&(T_obj->table[old_n_alloc]), n_alloc - old_n_alloc, struct GMT_TEXTTABLE *);	 Set to NULL
		}
		S_obj->alloc_mode = T_obj->alloc_mode;	 Clarify allocation mode for this entity
		S_obj->status = GMT_IS_USED;	 Mark as read
	}

	if (T_obj->n_tables == 0) {	 Only found empty files (e.g., /dev/null) and we have nothing to show for our efforts.  Return an single empty table with no segments.
		T_obj->table = GMT_memory (API->GMT, T_obj->table, 1, struct GMT_TEXTTABLE *);
		T_obj->table[0] = GMT_memory (API->GMT, NULL, 1, struct GMT_TEXTTABLE);
		T_obj->n_tables = 1;	 But we must indicate we found one (empty) table
	}
	else {	 Found one or more tables
		if (allocate && T_obj->n_tables < n_alloc) T_obj->table = GMT_memory (API->GMT, T_obj->table, T_obj->n_tables, struct GMT_TEXTTABLE *);
	}
	if (!via) T_obj->geometry = API->object[this_item]->geometry;
	API->object[this_item]->data = T_obj;		 Retain pointer to the allocated data so we use garbage collection later

	return (T_obj);
}*/


unsigned int GMTAPI_count_objects (struct GMTAPI_CTRL *API, enum GMT_enum_family family, unsigned int geometry, unsigned int direction, int *first_ID)
{	/* Count how many data sets of the given family are currently registered and unused for the given direction (GMT_IN|GMT_OUT).
 	 * Also return the ID of the first unused data object for the given direction, geometry, and family (GMT_NOTSET if not found).
	 */
	unsigned int i, n;

	*first_ID = GMT_NOTSET;	/* Not found yet */
	for (i = n = 0; i < API->n_objects; i++) {
		if (!API->object[i]) continue;				/* A freed object, skip */
		if (API->object[i]->direction != direction) continue;	/* Wrong direction */
		if (API->object[i]->geometry != geometry) continue;	/* Wrong geometry */
		if (API->object[i]->status != GMT_IS_UNUSED) continue;	/* Already used */
		if (family != API->object[i]->family) continue;		/* Wrong data type */
		n++;
		if (*first_ID == GMT_NOTSET) *first_ID = API->object[i]->ID;
	}
	return (n);
}

void gmt_close_grd (struct GMT_CTRL *GMT, struct GMT_GRID *G)
{
	struct GMT_GRID_ROWBYROW *R = gmt_get_rbr_ptr ((struct GMT_GRID_ROWBYROW *)G->extra);	/* Shorthand to row-by-row book-keeping structure */
	if (R->v_row) GMT_free (GMT, R->v_row);
	if (GMT->session.grdformat[G->header->type][0] == 'c' || GMT->session.grdformat[G->header->type][0] == 'n')
		nc_close (R->fid);
	else
		GMT_fclose (GMT, R->fp);
	GMT_free (GMT, G->extra);
}

unsigned int GMT_free_grid_ptr (struct GMT_CTRL *GMT, struct GMT_GRID *G, bool free_grid)
{	/* By taking a reference to the grid pointer we can set it to NULL when done */
	if (!G) return 0;	/* Nothing to deallocate */
	/* Only free G->data if allocated by GMT AND free_grid is true */
	if (G->data && free_grid) {
		if (G->alloc_mode == GMT_ALLOCATED_BY_GMT) GMT_free_aligned (GMT, G->data);
		G->data = NULL;	/* This will remove reference to external memory since GMT_free_aligned would not have been called */
	}
	if (G->extra) gmt_close_grd (GMT, G);	/* Close input file used for row-by-row i/o */
	//if (G->header && G->alloc_mode == GMT_ALLOCATED_BY_GMT) GMT_free (GMT, G->header);
	if (G->header) {	/* Free the header structure and anything allocated by it */
		if (G->header->pocket) free (G->header->pocket);
		GMT_free (GMT, G->header);
	}
	return (G->alloc_mode);
}

void GMT_free_dataset_ptr (struct GMT_CTRL *GMT, struct GMT_DATASET *data)
{	/* This takes pointer to data array and thus can return it as NULL */
	unsigned int tbl, k;
	if (!data) return;	/* Do not try to free NULL pointer */
	if (!data->table) return;	/* Do not try to free NULL pointer of tables */
	for (tbl = 0; tbl < data->n_tables; tbl++) {
		GMT_free_table (GMT, data->table[tbl], data->alloc_mode);
	}
	if (data->min) GMT_free (GMT, data->min);
	if (data->max) GMT_free (GMT, data->max);
	GMT_free (GMT, data->table);
	for (k = 0; k < 2; k++) if (data->file[k]) free (data->file[k]);
}

void gmt_free_textsegment (struct GMT_CTRL *GMT, struct GMT_TEXTSEGMENT *segment)
{
	/* Free memory allocated by GMT_read_texttable */

	uint64_t row;
	unsigned int k;
	if (!segment) return;	/* Do not try to free NULL pointer */
	for (row = 0; row < segment->n_rows; row++) if (segment->record[row]) free (segment->record[row]);
	GMT_free (GMT, segment->record);
	if (segment->label) free ( segment->label);
	if (segment->header) free ( segment->header);
	for (k = 0; k < 2; k++) if (segment->file[k]) free (segment->file[k]);
	GMT_free (GMT, segment);
}

void gmt_free_texttable (struct GMT_CTRL *GMT, struct GMT_TEXTTABLE *table)
{
	unsigned int k;
	uint64_t seg;
	if (!table) return;	/* Do not try to free NULL pointer */
	for (seg = 0; seg < table->n_segments; seg++) gmt_free_textsegment (GMT, table->segment[seg]);
	for (k = 0; k < table->n_headers; k++) free (table->header[k]);
	if (table->n_headers) GMT_free (GMT, table->header);
	if (table->segment) GMT_free (GMT, table->segment);
	for (k = 0; k < 2; k++) if (table->file[k]) free (table->file[k]);
	GMT_free (GMT, table);
}

void GMT_free_textset_ptr (struct GMT_CTRL *GMT, struct GMT_TEXTSET *data)
{	/* This takes pointer to data array and thus can return it as NULL */

	unsigned int tbl, k;
	for (tbl = 0; tbl < data->n_tables; tbl++) gmt_free_texttable (GMT, data->table[tbl]);
	GMT_free (GMT, data->table);
	for (k = 0; k < 2; k++) if (data->file[k]) free (data->file[k]);
}

int GMTAPI_destroy_data_ptr (struct GMTAPI_CTRL *API, enum GMT_enum_family family, void *ptr)
{
	/* Like GMT_Destroy_Data but takes pointer to data rather than address of pointer.
	 * We pass true to make sure we free the memory.  Some objects (grid, matrix, vector) may
	 * point to externally allocated memory so we return the alloc_mode for those items.
	 * This is mostly for information since the pointers to such external memory have now
	 * been set to NULL instead of being freed.
	 * The containers are always allocated by GMT so those are freed at the end.
	 */

	if (API == NULL) return (GMT_NOT_A_SESSION);
	if (!ptr) return (GMT_OK);	/* Null pointer */

	switch (family) {	/* dataset, cpt, text table or grid */
		case GMT_IS_GRID:	/* GMT grid; return alloc mode of data array in case it was allocated externally */
			GMT_free_grid_ptr (API->GMT, (struct GMT_GRID *)ptr, true);
			break;
		case GMT_IS_DATASET:
			GMT_free_dataset_ptr (API->GMT, (struct GMT_DATASET *)ptr);
			break;
		case GMT_IS_TEXTSET:
			GMT_free_textset_ptr (API->GMT, (struct GMT_TEXTSET *)ptr);
			break;
		//case GMT_IS_CPT:
			//GMT_free_cpt_ptr (API->GMT, ptr);
			//break;
//#ifdef HAVE_GDAL
	//	case GMT_IS_IMAGE:
		//	GMT_free_image_ptr (API->GMT, ptr, true);
			//break;
//#endif
		case GMT_IS_COORD:
			/* Nothing to do as GMT_free below will do it */
			break;

		/* Also allow destoying of intermediate vector and matrix containers */
		case GMT_IS_MATRIX:	/* GMT matrix; return alloc mode of data array in case it was allocated externally */
			//GMT_free_matrix_ptr (API->GMT, ptr, true);
			break;
		case GMT_IS_VECTOR:	/* GMT vector; return alloc mode of data array in case it was allocated externally */
			//GMT_free_vector_ptr (API->GMT, ptr, true);
			break;
		default:
			return (GMTAPI_report_error (API, GMT_NOT_A_VALID_FAMILY));
			break;
	}
	GMT_free (API->GMT, ptr);	/* OK to free container */
	return (GMT_OK);	/* Null pointer */
}

void GMTAPI_put_val (struct GMTAPI_CTRL *API, union GMT_UNIVECTOR *u, double val, uint64_t row, unsigned int type)
{ /* Places a double value in the <type> column array[i] pointed to by the union pointer *u, at row position row.
	 * No check to see if the type can hold the value is performed, so truncation may result.
	 * Used in GMTAPI_Export_Dataset and GMTAPI_Export_Grid. */

	switch (type) {	/* Use type to select the correct array to which we will put a value */
		case GMT_UCHAR:  u->uc1[row] = (uint8_t)val;  break;
		case GMT_CHAR:   u->sc1[row] = (int8_t)val;   break;
		case GMT_USHORT: u->ui2[row] = (uint16_t)val; break;
		case GMT_SHORT:  u->si2[row] = (int16_t)val;  break;
		case GMT_UINT:   u->ui4[row] = (uint32_t)val; break;
		case GMT_INT:    u->si4[row] = (int32_t)val;  break;
		case GMT_ULONG:  u->ui8[row] = (uint64_t)val; break;
		case GMT_LONG:   u->si8[row] = (int64_t)val;  break;
		case GMT_FLOAT:  u->f4[row]  = (float)val;    break;
		case GMT_DOUBLE: u->f8[row]  = val;           break;
		default:
			printf (/*API, GMT_MSG_NORMAL, */ "Internal error in GMTAPI_get_val: Passed bad type (%d)\n", type);
			API->error = GMT_NOT_A_VALID_TYPE;
			break;
	}
}

int GMTAPI_Export_Dataset (struct GMTAPI_CTRL *API, int object_ID, unsigned int mode, struct GMT_DATASET *D_obj)
{	/* Does the actual work of writing out the specified data set to one destination.
	 * If object_ID == GMT_NOTSET we use the first registered output destination, otherwise we just use the one requested.
	 * See the GMTAPI documentation for how mode is used to create multiple files from segments, etc.
	 */
	int item, error, default_method;
	uint64_t tbl, col, offset;
	uint64_t row, seg, ij;
	p_func_size_t GMT_2D_to_index = NULL;
	struct GMTAPI_DATA_OBJECT *S_obj = NULL;
	struct GMT_DATASET *D_copy = NULL;
	struct GMT_MATRIX *M_obj = NULL;
	struct GMT_VECTOR *V_obj = NULL;
	void *ptr = NULL;

	//printf(/*API, GMT_MSG_DEBUG,*/ "GMTAPI_Export_Dataset: Passed ID = %d and mode = %d\n", object_ID, mode);

	if (object_ID == GMT_NOTSET) return (GMTAPI_report_error (API, GMT_OUTPUT_NOT_SET));
	if ((item = GMTAPI_Validate_ID (API, GMT_IS_DATASET, object_ID, GMT_OUT)) == GMT_NOTSET) return (GMTAPI_report_error (API, API->error));

	S_obj = API->object[item];	/* S is the object whose data we will export */
	if (S_obj->family != GMT_IS_DATASET) return (GMTAPI_report_error (API, GMT_NOT_A_VALID_FAMILY));	/* Called with wrong data type */
	if (mode >= GMT_WRITE_TABLE && !S_obj->filename) return (GMTAPI_report_error (API, GMT_OUTPUT_NOT_SET));	/* Must have filename when segments are to be written */
	if (S_obj->status != GMT_IS_UNUSED && !(mode & GMT_IO_RESET))	/* Only allow writing of a data set once unless overridden by mode */
		return (GMTAPI_report_error (API, GMT_WRITTEN_ONCE));
	default_method = GMT_IS_FILE;
	if (S_obj->filename)	/* Write to this file */
		ptr = S_obj->filename;
	else {			/* No filename so we switch to writing to the stream or fdesc */
		default_method = (S_obj->method == GMT_IS_FILE) ? GMT_IS_STREAM : S_obj->method;
		ptr = S_obj->fp;
//#ifdef SET_IO_MODE
	//	GMT_setmode (API->GMT, GMT_OUT);	/* Windows may need to switch write mode from text to binary */
//#endif
	}
	D_obj->io_mode = mode;	/* Handles if tables or segments should be written to separate files */
	switch (S_obj->method) {	/* File, array, stream etc ? */
	 	case GMT_IS_STREAM:
//#ifdef SET_IO_MODE
	//		GMT_setmode (API->GMT, GMT_OUT);	/* Windows may need to switch write mode from text to binary */
//#endif
		case GMT_IS_FILE:
	 	case GMT_IS_FDESC:
			/* GMT_write_dataset (or lower) will report where it is reading from if level is GMT_MSG_LONG_VERBOSE */
			if ((error = GMT_write_dataset (API->GMT, ptr, default_method, D_obj, true, GMT_NOTSET))) return (GMTAPI_report_error (API, error));
			break;

		case GMT_IS_DUPLICATE:		/* Duplicate the input dataset */
			if (S_obj->resource) return (GMTAPI_report_error (API, GMT_PTR_NOT_NULL));	/* The output resource must be NULL */
			printf(/*API,GMT_MSG_LONG_VERBOSE, */ "Duplicating data table to GMT_DATASET memory location\n");
			D_copy = GMT_duplicate_dataset (API->GMT, D_obj, GMT_ALLOC_NORMAL, NULL);
			S_obj->resource = D_copy;	/* Set resource pointer from object to this dataset */
			break;

		case GMT_IS_REFERENCE:	/* Just pass memory location */
			if (S_obj->resource) return (GMTAPI_report_error (API, GMT_PTR_NOT_NULL));	/* The output resource must be NULL */
			printf(/*API, GMT_MSG_LONG_VERBOSE, */"Referencing data table to GMT_DATASET memory location\n");
			D_obj->alloc_level = S_obj->alloc_level;	/* Since we are passing it up to the caller */
			S_obj->alloc_mode = D_obj->alloc_mode;
			S_obj->resource = D_obj;			/* Set resource pointer from object to this dataset */
			break;

	 	case GMT_IS_DUPLICATE + GMT_VIA_MATRIX:
			if (S_obj->resource == NULL) return (GMTAPI_report_error (API, GMT_PTR_IS_NULL));	/* The output resource must initially have info needed to do the output */
			printf (/*API, GMT_MSG_LONG_VERBOSE, */"Duplicating data table to user array location\n");
			M_obj = GMT_duplicate_matrix (API->GMT, (struct GMT_MATRIX *)S_obj->resource, false);
			/* Must allocate output space */
			M_obj->n_rows = D_obj->n_records;	/* Number of rows needed to hold the data records */
			M_obj->n_columns = D_obj->n_columns;	/* Number of columns needed to hold the data records */
			if (API->GMT->current.io.multi_segments[GMT_OUT]) M_obj->n_rows += D_obj->n_segments;	/* Add one row for each segment header */
			if (M_obj->shape == GMT_IS_ROW_FORMAT) {	/* C-style matrix layout */
				if (M_obj->dim == 0) M_obj->dim = D_obj->n_columns;
				if (M_obj->dim < D_obj->n_columns) return (GMTAPI_report_error (API, GMT_DIM_TOO_SMALL));
				S_obj->n_alloc = M_obj->n_rows * M_obj->dim;	/* Get total number of elements as n_rows * dim */
			}
			else {	/* Fortran style */
				if (M_obj->dim == 0) M_obj->dim = D_obj->n_records;
				if (M_obj->dim < D_obj->n_records) return (GMTAPI_report_error (API, GMT_DIM_TOO_SMALL));
				S_obj->n_alloc = M_obj->n_columns * M_obj->dim;	/* Get total number of elements as n_columns * dim */
			}
			if ((error = GMT_alloc_univector (API->GMT, &(M_obj->data), M_obj->type, S_obj->n_alloc)) != GMT_OK) return (GMTAPI_report_error (API, error));
			GMT_2D_to_index = GMTAPI_get_2D_to_index (API, M_obj->shape, GMT_GRID_IS_REAL);

			for (tbl = offset = 0; tbl < D_obj->n_tables; tbl++) {
				for (seg = 0; seg < D_obj->table[tbl]->n_segments; seg++) {
					for (row = 0; row < D_obj->table[tbl]->segment[seg]->n_rows; row++) {
						for (col = 0; col < D_obj->table[tbl]->segment[seg]->n_columns; col++) {
							ij = GMT_2D_to_index (row + offset, col, M_obj->dim);
							GMTAPI_put_val (API, &(M_obj->data), D_obj->table[tbl]->segment[seg]->coord[col][row], ij, M_obj->type);
						}
					}
					offset += D_obj->table[tbl]->segment[seg]->n_rows;	/* Since row starts at 0 for each segment */
				}
			}
			S_obj->resource = M_obj;	/* Set resource pointer from object to this matrix */
			break;

		case GMT_IS_DUPLICATE + GMT_VIA_VECTOR:
		/*case GMT_IS_REFERENCE + GMT_VIA_VECTOR:
			if ((V_obj = S_obj->resource) == NULL) return (GMTAPI_report_error (API, GMT_PTR_IS_NULL));	 The output resource must initially have info needed to do the output
			printf (API, GMT_MSG_LONG_VERBOSE, "Duplicating data table to user column arrays location\n");
			V_obj = GMT_duplicate_vector (API->GMT, S_obj->resource, false);
			V_obj->n_rows = D_obj->n_records;
			if (API->GMT->current.io.multi_segments[GMT_OUT]) V_obj->n_rows += D_obj->n_segments;
			if ((error = gmt_alloc_vectors (API->GMT, V_obj)) != GMT_OK) return (GMTAPI_report_error (API, error));
			for (tbl = ij = 0; tbl < D_obj->n_tables; tbl++) {
				for (seg = 0; seg < D_obj->table[tbl]->n_segments; seg++) {
					for (row = 0; row < D_obj->table[tbl]->segment[seg]->n_rows; row++, ij++) {
						for (col = 0; col < D_obj->table[tbl]->segment[seg]->n_columns; col++) {
							GMTAPI_put_val (API, &(V_obj->data[col]), D_obj->table[tbl]->segment[seg]->coord[col][row], ij, V_obj->type[col]);
						}
					}
				}
			}
			S_obj->resource = V_obj;
			break;*/

		default:
			printf (/*API, GMT_MSG_NORMAL,*/ "Wrong method used to export data tables\n");
			return (GMTAPI_report_error (API, GMT_NOT_A_VALID_METHOD));
			break;
	}
	S_obj->alloc_mode = D_obj->alloc_mode;	/* Clarify allocation mode for this entity */
	S_obj->status = GMT_IS_USED;	/* Mark as written */
	S_obj->data = D_obj;		/* Retain pointer to the allocated data so we can find its object via the data pointer later */
	S_obj->data = NULL;

	return GMT_OK;
}

void gmt_grd_xy_scale (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *h, unsigned int direction)
{
	unsigned int k;
	/* Apply the scaling of wesn,inc as given by the header's xy_* settings.
	 * After reading a grid it will have wesn/inc in meters.
	 * Before writing a grid, it may have units changed back to original units
	 * or scaled to anoter set of units */
	
	if (direction == GMT_IN) {
		if (h->xy_adjust[direction] == 0) return;	/* Nothing to do */
		if (h->xy_adjust[GMT_IN] & 2) return;		/* Already scaled them */
		for (k = 0; k < 4; k++) h->wesn[k] *= h->xy_unit_to_meter[GMT_IN];
		for (k = 0; k < 2; k++) h->inc[k]  *= h->xy_unit_to_meter[GMT_IN];
		h->xy_adjust[GMT_IN] = 2;	/* Now the grid is ready for use and in meters */
		if (h->xy_mode[direction])
			printf (/*GMT->parent, GMT_MSG_LONG_VERBOSE,*/ "Input grid file x/y unit was converted from meters to %s after reading.\n", GMT->current.proj.unit_name[h->xy_unit[GMT_IN]]);
		else
			printf (/*GMT->parent, GMT_MSG_LONG_VERBOSE,*/ "Input grid file x/y unit was converted from %s to meters after reading.\n", GMT->current.proj.unit_name[h->xy_unit[GMT_IN]]);
	}
	else if (direction == GMT_OUT) {	/* grid x/y are assumed to be in meters */
		if (h->xy_adjust[GMT_OUT] & 1) {	/* Was given a new unit for output */
			for (k = 0; k < 4; k++) h->wesn[k] /= h->xy_unit_to_meter[GMT_OUT];
			for (k = 0; k < 2; k++) h->inc[k]  /= h->xy_unit_to_meter[GMT_OUT];
			h->xy_adjust[GMT_OUT] = 2;	/* Now we are ready for writing */
			if (h->xy_mode[GMT_OUT])
				printf (/*GMT->parent, GMT_MSG_LONG_VERBOSE,*/ "Output grid file x/y unit was converted from %s to meters before writing.\n", GMT->current.proj.unit_name[h->xy_unit[GMT_OUT]]);
			else
				printf (/*GMT->parent, GMT_MSG_LONG_VERBOSE,*/ "Output grid file x/y unit was converted from meters to %s before writing.\n", GMT->current.proj.unit_name[h->xy_unit[GMT_OUT]]);
		}
		else if (h->xy_adjust[GMT_IN] & 2) {	/* Just undo old scaling */
			for (k = 0; k < 4; k++) h->wesn[k] /= h->xy_unit_to_meter[GMT_IN];
			for (k = 0; k < 2; k++) h->inc[k]  /= h->xy_unit_to_meter[GMT_IN];
			h->xy_adjust[GMT_IN] -= 2;	/* Now it is back to where we started */
			if (h->xy_mode[GMT_OUT])
				printf (/*GMT->parent, GMT_MSG_LONG_VERBOSE,*/ "Output grid file x/y unit was reverted back to %s from meters before writing.\n", GMT->current.proj.unit_name[h->xy_unit[GMT_IN]]);
			else
				printf (/*GMT->parent, GMT_MSG_LONG_VERBOSE,*/ "Output grid file x/y unit was reverted back from meters to %s before writing.\n", GMT->current.proj.unit_name[h->xy_unit[GMT_IN]]);
		}
	}
}


int GMT_read_grd_info (struct GMT_CTRL *GMT, char *file, struct GMT_GRID_HEADER *header)
{	/* file:	File name
	 * header:	grid structure header
	 * Note: The header reflects what is actually in the file, and all the dimensions
	 * reflect the number of rows, cols, size, pads etc.  However, if GMT_read_grd is
	 * called requesting a subset then these will be reset accordingly.
	 */

	int err;	/* Implied by GMT_err_trap */
	double scale, offset;
	float invalid;

	/* Save parameters on file name suffix before issuing GMT->session.readinfo */
	//GMT_err_trap (GMT_grd_get_format (GMT, file, header, true));
	GMT_grd_get_format (GMT, file, header, true);
	
	/* remember scale, offset, and invalid: */
	scale = header->z_scale_factor;
	offset = header->z_add_offset;
	invalid = header->nan_value;
	

	//GMT_err_trap ((*GMT->session.readinfo[header->type]) (GMT, header));
	(*GMT->session.readinfo[header->type]) (GMT, header);

	gmt_grd_xy_scale (GMT, header, GMT_IN);	/* Possibly scale wesn,inc */

	/* restore non-default scale, offset, and invalid: */
	if (scale != 1.0)
		header->z_scale_factor = scale;
	if (offset != 0.0)
		header->z_add_offset = offset;
	if (isfinite(invalid))
		header->nan_value = invalid;

	gmt_grd_get_units (GMT, header);
	header->grdtype = gmt_get_grdtype (GMT, header);

	//GMT_err_pass (GMT, GMT_grd_RI_verify (GMT, header, 0), file);
	GMT_set_grddim (GMT, header);	/* Set all integer dimensions and xy_off */
	
	/* unpack z-range: */
	header->z_min = header->z_min * header->z_scale_factor + header->z_add_offset;
	header->z_max = header->z_max * header->z_scale_factor + header->z_add_offset;

	return (GMT_NOERROR);
}

void GMTAPI_grdheader_to_info (struct GMT_GRID_HEADER *h, struct GMT_MATRIX *M_obj)
{	/* Packs the necessary items of the grid header into the matrix parameters */
	M_obj->n_columns = h->nx;
	M_obj->n_rows = h->ny;
	M_obj->registration = h->registration;
	GMT_memcpy (M_obj->range, h->wesn, 4, double);
}

void GMT_free_grid (struct GMT_CTRL *GMT, struct GMT_GRID **G, bool free_grid)
{	/* By taking a reference to the grid pointer we can set it to NULL when done */
	(void)GMT_free_grid_ptr (GMT, *G, free_grid);
	if (*G) GMT_free (GMT, *G);
}



int GMTAPI_Export_Grid (struct GMTAPI_CTRL *API, int object_ID, unsigned int mode, struct GMT_GRID *G_obj)
{	/* Writes out a single grid to destination */
	int item, error;
	bool done = true, row_by_row;
	uint64_t row, col, i0, i1, j0, j1, ij, ijp, ij_orig;
	size_t size;
	double dx, dy;
	p_func_size_t GMT_2D_to_index = NULL;
	struct GMTAPI_DATA_OBJECT *S_obj = NULL;
	struct GMT_GRID *G_copy = NULL;
	struct GMT_MATRIX *M_obj = NULL;

	printf(/*API, GMT_MSG_DEBUG,*/ "GMTAPI_Export_Grid: Passed ID = %d and mode = %d\n", object_ID, mode);

	if (object_ID == GMT_NOTSET) return (GMTAPI_report_error (API, GMT_OUTPUT_NOT_SET));
	if ((item = GMTAPI_Validate_ID (API, GMT_IS_GRID, object_ID, GMT_OUT)) == GMT_NOTSET) return (GMTAPI_report_error (API, API->error));

	S_obj = API->object[item];	/* The current object whose data we will export */
	if (S_obj->status != GMT_IS_UNUSED && !(mode & GMT_IO_RESET)) return (GMTAPI_report_error (API, GMT_WRITTEN_ONCE));	/* Only allow writing of a data set once, unless overridden by mode */
	row_by_row = ((mode & GMT_GRID_ROW_BY_ROW) || (mode & GMT_GRID_ROW_BY_ROW_MANUAL));
	if (row_by_row && S_obj->method != GMT_IS_FILE) {
		printf (/*API, GMT_MSG_NORMAL, */"Can only use method GMT_IS_FILE when row-by-row writing of grid is selected\n");
		return (GMTAPI_report_error (API, GMT_NOT_A_VALID_METHOD));
	}
	if (S_obj->region && G_obj) {	/* See if this is really a subset or just the same region as the grid */
		if (G_obj->header->wesn[XLO] == S_obj->wesn[XLO] && G_obj->header->wesn[XHI] == S_obj->wesn[XHI] && G_obj->header->wesn[YLO] == S_obj->wesn[YLO] && G_obj->header->wesn[YHI] == S_obj->wesn[YHI]) S_obj->region = false;
	}
	switch (S_obj->method) {
		case GMT_IS_FILE:	/* Name of a grid file on disk */
			if (mode & GMT_GRID_HEADER_ONLY) {	/* Update header structure only */
				printf (/*API, GMT_MSG_LONG_VERBOSE,*/ "Updating grid header for file %s\n", S_obj->filename);
				if (row_by_row) {	/* Special row-by-row processing mode */
					char w_mode = (mode & GMT_GRID_NO_HEADER) ? 'W' : 'w';
					/* Since we may get here twice (initial write; later update) we only allocate extra if NULL */
					if (G_obj->extra == NULL) G_obj->extra = GMT_memory (API->GMT, NULL, 1, struct GMT_GRID_ROWBYROW);
					if (gmt_open_grd (API->GMT, S_obj->filename, G_obj, w_mode, mode))	/* Open the grid for incremental row writing */
						return (GMTAPI_report_error (API, GMT_GRID_WRITE_ERROR));
				}
				else if (GMT_update_grd_info (API->GMT, NULL, G_obj->header))
					return (GMTAPI_report_error (API, GMT_GRID_WRITE_ERROR));
				done = false;	/* Since we are not done with writing */
			}
			else {
				printf (/*API, GMT_MSG_LONG_VERBOSE, */"Writing grid to file %s\n", S_obj->filename);
				GMT_write_grd (API->GMT, S_obj->filename, G_obj->header, G_obj->data, S_obj->wesn, G_obj->header->pad, mode);
				//if (GMT_err_pass (API->GMT, GMT_write_grd (API->GMT, S_obj->filename, G_obj->header, G_obj->data, S_obj->wesn, G_obj->header->pad, mode), S_obj->filename))
					return (GMTAPI_report_error (API, GMT_GRID_WRITE_ERROR));
			}
			break;

	 	case GMT_IS_DUPLICATE:	/* Duplicate GMT grid and header to a GMT_GRID container object. Subset allowed */
			if (S_obj->resource) return (GMTAPI_report_error (API, GMT_PTR_NOT_NULL));	/* The ouput resource pointer must be NULL */
			if (mode & GMT_GRID_HEADER_ONLY) return (GMTAPI_report_error (API, GMT_NOT_A_VALID_MODE));
			printf (/*API, GMT_MSG_LONG_VERBOSE,*/ "Duplicating grid data to GMT_GRID memory location\n");
			if (!S_obj->region) {	/* No subset, possibly same padding */
				G_copy = GMT_Duplicate_Data (API, GMT_IS_GRID, GMT_DUPLICATE_DATA, G_obj);
				if (GMTAPI_adjust_grdpadding (G_copy->header, API->GMT->current.io.pad))
					GMT_grd_pad_on (API->GMT, G_copy, API->GMT->current.io.pad);
				GMT_BC_init (API->GMT, G_copy->header);	/* Initialize grid interpolation and boundary condition parameters */
				//if (GMT_err_pass (API->GMT, GMT_grd_BC_set (API->GMT, G_copy, GMT_OUT), "Grid memory")) return (GMTAPI_report_error (API, GMT_GRID_BC_ERROR));	/* Set boundary conditions */
				//S_obj->resource = G_copy;	/* Set resource pointer to the grid */
				//break;		/* Done with this grid */
			}
			/* Here we need to extract subset, and possibly change padding. */
			/* Get start/stop row/cols for subset (or the entire domain) */
			G_copy = GMT_create_grid (API->GMT);
			GMT_memcpy (G_copy->header, G_obj->header, 1, struct GMT_GRID_HEADER);
			GMT_memcpy (G_copy->header->wesn, S_obj->wesn, 4, double);
			dx = G_obj->header->inc[GMT_X] * G_obj->header->xy_off;	dy = G_obj->header->inc[GMT_Y] * G_obj->header->xy_off;
			j1 = (unsigned int) GMT_grd_y_to_row (API->GMT, G_obj->header->wesn[YLO]+dy, G_obj->header);
			j0 = (unsigned int) GMT_grd_y_to_row (API->GMT, G_obj->header->wesn[YHI]-dy, G_obj->header);
			i0 = (unsigned int) GMT_grd_x_to_col (API->GMT, G_obj->header->wesn[XLO]+dx, G_obj->header);
			i1 = (unsigned int) GMT_grd_x_to_col (API->GMT, G_obj->header->wesn[XHI]-dx, G_obj->header);
			GMT_memcpy (G_obj->header->pad, API->GMT->current.io.pad, 4, int);		/* Set desired padding */
			G_copy->header->size = GMTAPI_set_grdarray_size (API->GMT, G_obj->header, mode, S_obj->wesn);	/* Get array dimension only, which may include padding */
			G_copy->data = GMT_memory_aligned (API->GMT, NULL, G_copy->header->size, float);
			G_copy->header->z_min = DBL_MAX;	G_copy->header->z_max = -DBL_MAX;	/* Must set zmin/zmax since we are not writing */
			for (row = j0; row <= j1; row++) {
				for (col = i0; col <= i1; col++, ij++) {
					ij_orig = GMT_IJP (G_obj->header, row, col);	/* Position of this (row,col) in original grid organization */
					ij = GMT_IJP (G_copy->header, row, col);	/* Position of this (row,col) in output grid organization */
					G_copy->data[ij] = G_obj->data[ij_orig];
					if (GMT_is_fnan (G_copy->data[ij])) continue;
					/* Update z_min, z_max */
					G_copy->header->z_min = MIN (G_copy->header->z_min, (double)G_copy->data[ij]);
					G_copy->header->z_max = MAX (G_copy->header->z_max, (double)G_copy->data[ij]);
				}
			}
			GMT_BC_init (API->GMT, G_copy->header);	/* Initialize grid interpolation and boundary condition parameters */
			
			
			//if (GMT_err_pass (API->GMT, GMT_grd_BC_set (API->GMT, G_copy, GMT_OUT), "Grid memory")) return (GMTAPI_report_error (API, GMT_GRID_BC_ERROR));	/* Set boundary conditions */


			S_obj->resource = G_copy;	/* Set resource pointer to the grid */
			break;

	 	case GMT_IS_REFERENCE:	/* GMT grid and header in a GMT_GRID container object - just pass the reference */
			if (S_obj->region) return (GMTAPI_report_error (API, GMT_SUBSET_NOT_ALLOWED));
			if (mode & GMT_GRID_HEADER_ONLY) return (GMTAPI_report_error (API, GMT_NOT_A_VALID_MODE));
			printf (/*API, GMT_MSG_LONG_VERBOSE,*/ "Referencing grid data to GMT_GRID memory location\n");
			if (GMTAPI_adjust_grdpadding (G_obj->header, API->GMT->current.io.pad))
				GMT_grd_pad_on (API->GMT, G_obj, API->GMT->current.io.pad);	/* Adjust pad */
			GMT_grd_zminmax (API->GMT, G_obj->header, G_obj->data);	/* Must set zmin/zmax since we are not writing */
			GMT_BC_init (API->GMT, G_obj->header);	/* Initialize grid interpolation and boundary condition parameters */
			//if (GMT_err_pass (API->GMT, GMT_grd_BC_set (API->GMT, G_obj, GMT_OUT), "Grid memory")) return (GMTAPI_report_error (API, GMT_GRID_BC_ERROR));	/* Set boundary conditions */
			S_obj->resource = G_obj;	/* Set resource pointer to the grid */
			G_obj->alloc_level = S_obj->alloc_level;	/* Since we are passing it up to the caller */
			break;

	 	/*case GMT_IS_DUPLICATE + GMT_VIA_MATRIX:	 The user's 2-D grid array of some sort, + info in the args [NOT FULLY TESTED]
			if (S_obj->resource == NULL) return (GMTAPI_report_error (API, GMT_PTR_IS_NULL));	 The output resource pointer cannot be NULL for matrix
			if (mode & GMT_GRID_HEADER_ONLY) return (GMTAPI_report_error (API, GMT_NOT_A_VALID_MODE));
			M_obj = GMT_duplicate_matrix (API->GMT, S_obj->resource, false);
			GMTAPI_grdheader_to_info (G_obj->header, M_obj);	 Populate an array with GRD header information
			printf (API, GMT_MSG_LONG_VERBOSE, "Exporting grid data to user memory location\n");
			size = GMT_get_nm (API->GMT, G_obj->header->nx, G_obj->header->ny);
			if ((error = GMT_alloc_univector (API->GMT, &(M_obj->data), M_obj->type, size)) != GMT_OK) return (error);
			GMT_2D_to_index = GMTAPI_get_2D_to_index (API, M_obj->shape, GMT_GRID_IS_REAL);
			GMT_grd_loop (API->GMT, G_obj, row, col, ijp) {
				ij = GMT_2D_to_index (row, col, M_obj->dim);
				GMTAPI_put_val (API, &(M_obj->data), (double)G_obj->data[ijp], ij, M_obj->type);
			}
			S_obj->resource = M_obj;	 Set resource pointer to the matrix
			break;*/
		#if 0
	 	case GMT_IS_REFERENCE + GMT_VIA_MATRIX:	/* The user's 2-D grid array of some sort, + info in the args [NOT FULLY TESTED] */
			if (S_obj->resource == NULL) return (GMTAPI_report_error (API, GMT_PTR_IS_NULL));	/* The output resource pointer cannot be NULL for matrix */
			if (mode & GMT_GRID_HEADER_ONLY) return (GMTAPI_report_error (API, GMT_NOT_A_VALID_MODE));
			if (GMTAPI_adjust_grdpadding (G_obj->header, API->GMT->current.io.pad))
				GMT_grd_pad_on (API->GMT, G_obj, API->GMT->current.io.pad);	/* Adjust pad */
			M_obj = GMT_duplicate_matrix (API->GMT, S_obj->resource, false);
			if (!(M_obj->shape == GMT_IS_ROW_FORMAT && M_obj->type == GMT_FLOAT && M_obj->alloc_mode == GMT_ALLOCATED_EXTERNALLY && (mode & GMT_GRID_IS_COMPLEX_MASK)))
				return (GMTAPI_report_error (API, GMT_NOT_A_VALID_IO_ACCESS));
			GMTAPI_grdheader_to_info (G_obj->header, M_obj);	/* Populate an array with GRD header information */
			printf (/*API, GMT_MSG_LONG_VERBOSE,*/ "Referencing grid data to user memory location\n");
			M_obj->data.f4 = G_obj->data;
			S_obj->resource = M_obj;	/* Set resource pointer to the matrix */
			break;
		#endif
		default:
			printf (/*API, GMT_MSG_NORMAL, */"Wrong method used to export grids\n");
			return (GMTAPI_report_error (API, GMT_NOT_A_VALID_METHOD));
			break;
	}

	if (done) S_obj->status = GMT_IS_USED;	/* Mark as written (unless we only updated header) */
	S_obj->data = G_obj;		/* Retain pointer to the allocated data so we can find the object via its pointer later */
	S_obj->data = NULL;

	return (GMT_OK);
}

int GMTAPI_Export_Data (struct GMTAPI_CTRL *API, enum GMT_enum_family family, int object_ID, unsigned int mode, void *data)
{
	/* Function that will export the single data object referred to by the object_ID as registered by GMT_Register_IO.
	 * Note: While there is no GMTAPI_Export_Image, these are handles as grids via GMTAPI_Export_Grid.
	 */
	int error, item;

	if (API == NULL) return (GMT_NOT_A_SESSION);			/* GMT_Create_Session has not been called */
	if (!API->registered[GMT_OUT]) return (/*GMTAPI_report_error (API,*/ GMT_NO_OUTPUT/*)*/);		/* No destination registered yet */

	/* Get information about this resource first */
	if ((item = GMTAPI_Validate_ID (API, family, object_ID, GMT_OUT)) == GMT_NOTSET) return (/*GMTAPI_report_error (API,*/ API->error/*)*/);

	/* The case where object_ID is not set but a virtual (memory) file is found is a special case: we must supply the correct object_ID */
	if (object_ID == GMT_NOTSET && item && API->object[item]->method != GMT_IS_FILE) object_ID = API->object[item]->ID;	/* Found virtual file; set actual object_ID */

	/* Check if this is a container passed from the outside to capture output */
	if (API->object[item]->messenger && API->object[item]->data) {	/* Need to destroy the dummy container before passing data out */
		error = GMTAPI_destroy_data_ptr (API, API->object[item]->family, API->object[item]->data);	/* Do the dirty deed */
		API->object[item]->messenger = false;	/* OK, now clean for output */
	}

	switch (family) {	/* CPT, Dataset, Textfile, or Grid */
		//case GMT_IS_CPT:	/* Export a CPT */
			//error = GMTAPI_Export_CPT (API, object_ID, mode, data);
			//break;
		case GMT_IS_DATASET:	/* Export a Data set */
			error = GMTAPI_Export_Dataset (API, object_ID, mode, (struct GMT_DATASET *)data);
			break;
		case GMT_IS_TEXTSET:	/* Export a Text set */
			//error = GMTAPI_Export_Textset (API, object_ID, mode, data);
			break;
		case GMT_IS_GRID:	/* Export a GMT grid */
			//printf("File: gmt_api.c Line : 2931\n");
			error = GMTAPI_Export_Grid (API, object_ID, mode, (struct GMT_GRID *)data);
			break;
		default:
			error = GMT_NOT_A_VALID_FAMILY;
			break;
	}
	return (/*GMTAPI_report_error (API,*/ error/*)*/);	/* Return status */
}

int GMT_Put_Data (void *V_API, int object_ID, unsigned int mode, void *data)
{
	/* Function to write data directly from program memory as a set (not record-by-record).
	 * We can combine the <register resource - export resource > sequence in
	 * one combined function.  See GMT_Register_IO for details on arguments.
	 * Here, *data is the pointer to the data object to save (CPT, dataset, textset, Grid)
	 * ID is the registered destination.
	 * While only one output destination is allowed, for DATA|TEXTSETS one can
	 * have the tables and even segments be written to individual files (see the mode
	 * description in the documentation for how to enable this feature.)
	 * Return: false if all is well, true if there was an error (and set API->error).
	 */
	int item;
	bool was_enabled;
	struct GMTAPI_CTRL *API = NULL;
	//printf("kawegruy\n");
	if (V_API == NULL)
		//return_error (V_API, GMT_NOT_A_SESSION);
		return (GMT_NOT_A_SESSION); //nishita
	API = gmt_get_api_ptr ((struct GMTAPI_CTRL *)V_API);

	/* Determine the item in the object list that matches this ID and direction */
	if ((item = GMTAPI_Validate_ID (API, GMT_NOTSET, object_ID, GMT_OUT)) == GMT_NOTSET)
		return_error (API, API->error);
		//return (API->error); //nishita
	was_enabled = API->io_enabled[GMT_OUT];
	if (!was_enabled && GMTAPI_Begin_IO (API, GMT_OUT) != GMT_OK) {	/* Enables data output if not already set and sets access mode */
		return_error (API, API->error);
		//return(API->error); // nishita
	}
	if (GMTAPI_Export_Data (API, API->object[item]->family, object_ID, mode, data) != GMT_OK) return_error (API, API->error);

	if (!was_enabled && GMT_End_IO (API, GMT_OUT, 0) != GMT_OK) {	/* Disables data output if we had to set it in this function */
		return_error (API, API->error);
		//return(API->error); // nishita
	}
//#ifdef DEBUG
	//GMTAPI_Set_Object (API, API->object[item]);
	//GMT_list_API (API, "GMT_Put_Data");
//#endif
	return (GMT_OK);	/* No error encountered */
}

int GMTAPI_Get_Object (struct GMTAPI_CTRL *API, int sfamily, void *ptr)
{	/* Returns the ID of the first object whose data pointer matches ptr.
	 * Unless family is GMT_NOTSET the object must be of the specified family.
	 */
	unsigned int i;
	enum GMT_enum_family family;
	int object_ID = GMT_NOTSET;	/* Not found yet */

	if (sfamily != GMT_NOTSET) family = sfamily;
	for (i = 0; object_ID == GMT_NOTSET && i < API->n_objects; i++) {	/* Loop over all objects */
		if (!API->object[i]) continue;	/* Skip freed objects */
		if (API->object[i]->data == NULL) continue;	/* No data pointer */
		if (sfamily != GMT_NOTSET && API->object[i]->family != family) continue;	/* Not the right family */
		if (API->object[i]->data == ptr && object_ID == GMT_NOTSET) object_ID = API->object[i]->ID;	/* Found a matching data pointer */
	}
	return (object_ID);	/* Return ID or -1 if not found */
}

int GMT_Write_Data (void *V_API, unsigned int family, unsigned int method, unsigned int geometry, unsigned int mode, double wesn[], char *output, void *data)
{
	/* Function to write data directly from program memory as a set (not record-by-record).
	 * We can combine the <register resource - export resource > sequence in
	 * one combined function.  See GMT_Register_IO for details on arguments.
	 * Here, *data is the pointer to the data object to save (CPT, dataset, textset, Grid)
	 * Case 1: output != NULL: Register this as the destination and export data.
	 * Case 2: output == NULL: Register stdout as the destination and export data.
	 * Case 3: geometry == 0: Use a previously registered single destination.
	 * While only one output destination is allowed, for DATA|TEXTSETS one can
	 * have the tables and even segments be written to individual files (see the mode
	 * description in the documentation for how to enable this feature.)
	 * Return: false if all is well, true if there was an error (and set API->error).
	 */
	unsigned int n_reg;
	int out_ID;
	struct GMTAPI_CTRL *API = NULL;

	if (V_API == NULL)
		//return_error (V_API, GMT_NOT_A_SESSION);
		return (GMT_NOT_A_SESSION); //nishita
	API = gmt_get_api_ptr((struct GMTAPI_CTRL *)V_API);

	if (output) {	/* Case 1: Save to a single specified destination (file or memory).  Register it first. */
		if ((out_ID = GMTAPI_Memory_Registered (API, (enum GMT_enum_family)family, GMT_OUT, output)) != GMT_NOTSET) {
			/* Output is a memory resource, passed via a @GMTAPI@-###### file name, and ###### is the out_ID.
			   In this case we must make some further checks.  We need to find the API object that holds data.
			   We do this below and get in_ID (the id of the data to write), whie out_ID is the id of where
			   things go (the output "memory").  Having the in_ID we get the array index in_item that matches
			   this ID and of the correct family.  We set direction to GMT_NOTSET since otherwise we may be
			   denied a hit as we dont really know what the direction is for in_ID.  Once in_item has been
			   secured we transfer ownership of this data from the in_ID object to the out_ID object.  That
			   way we avoid accidental premature freeing of the data object via the in_ID object since it now
			   will live on via out_ID and outlive the current module.
			    */
			int in_ID = GMT_NOTSET,  in_item = GMT_NOTSET;
			in_ID = GMTAPI_Get_Object (API, family, data);	/* Get the object ID of the input source */
			if (in_ID != GMT_NOTSET) in_item = GMTAPI_Validate_ID (API, family, in_ID, GMT_NOTSET);	/* Get the item in the API array; pass dir = GMT_NOTSET to bypass status check */
			if (in_item != GMT_NOTSET) {
				printf(/*API, GMT_MSG_DEBUG,*/ "GMT_Write_Data: Writing %s to memory object %d from object %d which transfers ownership\n", GMT_family[family], out_ID, in_ID);
				API->object[in_item]->no_longer_owner = true;	/* Since we have passed the content onto an output object */
			}
		}	/* else it is a regular file and we just register it and get the new out_ID needed below */
		else if ((out_ID = GMT_Register_IO (API, family, method, geometry, GMT_OUT, wesn, output)) == GMT_NOTSET) return_error (API, API->error);
	}
	else if (output == NULL && geometry) {	/* Case 2: Save to stdout.  Register stdout first. */
		if (family == GMT_IS_GRID)
			//return_error (API, GMT_STREAM_NOT_ALLOWED);	/* Cannot write grids to stream */
			return (GMT_STREAM_NOT_ALLOWED); //nishita
		if ((out_ID = GMT_Register_IO (API, family, GMT_IS_STREAM, geometry, GMT_OUT, wesn, API->GMT->session.std[GMT_OUT])) == GMT_NOTSET) return_error (API, API->error);	/* Failure to register std??? */
	}
	else {	/* Case 3: output == NULL && geometry == 0, so use the previously registered destination */
		if ((n_reg = GMTAPI_count_objects (API, (enum GMT_enum_family)family, geometry, GMT_OUT, &out_ID)) != 1)
			return_error (API, GMT_NO_OUTPUT);	/* There is no registered output */
			//return(GMT_NO_OUTPUT); // nishita
	}
	/* With out_ID in hand we can now put the data where it should go */
	if (GMT_Put_Data (API, out_ID, mode, data) != GMT_OK)
		return_error (API, API->error);
		//return (API->error); //nishita

//#ifdef DEBUG
	//GMT_list_API (API, "GMT_Write_Data");
//#endif
	return (GMT_OK);	/* No error encountered */
}
int GMTAPI_Destroy_Grid (struct GMTAPI_CTRL *API, struct GMT_GRID **G_obj)
{
	/* Delete the given grid resource. */

	if (!(*G_obj)) {	/* Probably not a good sign */
		//GMT_Report (API, GMT_MSG_DEBUG, "GMTAPI_Destroy_Grid: Passed NULL pointer - skipped\n");
		return (GMT_PTR_IS_NULL);
	}
	if ((*G_obj)->alloc_level != API->GMT->hidden.func_level) return (GMT_FREE_WRONG_LEVEL);	/* Not the right level */

	GMT_free_grid (API->GMT, G_obj, true);
	return GMT_OK;
}

/* return_address is a convenience function that, given type, calls the correct converter */
void *return_address (void *data, unsigned int type) {
	void *p = NULL;
	switch (type) {
		case GMT_IS_GRID:	p = gmt_get_grid_ptr (data);	break;
		case GMT_IS_DATASET:	p = gmt_get_dataset_ptr (data);	break;
		case GMT_IS_TEXTSET:	p = gmt_get_textset_ptr (data);	break;
		//case GMT_IS_CPT:	p = gmt_get_cpt_ptr (data);	break;
		//case GMT_IS_MATRIX:	p = gmt_get_matrix_ptr (data);	break;
		//case GMT_IS_VECTOR:	p = gmt_get_vector_ptr (data);	break;
		case GMT_IS_COORD:	p = gmt_get_coord_ptr ((double **)data);	break;
#ifdef HAVE_GDAL
		case GMT_IS_IMAGE:	p = gmt_get_image_ptr (data);	break;
#endif
	}
	return (p);
}

int GMTAPI_get_objectID_from_data_ptr (struct GMTAPI_CTRL *API, void *ptr)
{	/* Returns the ID of the first object whose data pointer matches *ptr.
 	 * This is necessary since many objects may have the same pointer
	 * but we only want to destroy the memory once.  This function is
	 * only used in GMT_Destroy_Data.
	 */
	unsigned int i;
	int object_ID = GMT_NOTSET;	/* Not found yet */
	void *data = NULL;

	for (i = 0; object_ID == GMT_NOTSET && i < API->n_objects; i++) {	/* Loop over all objects */
		if (!API->object[i]) continue;	/* Skip freed objects */
		data = return_address (ptr, API->object[i]->family);	/* Get void* pointer to resource */
		if (API->object[i]->data == data && object_ID == GMT_NOTSET) object_ID = API->object[i]->ID;	/* Found a matching data pointer */
	}
	return (object_ID);	/* Return ID or -1 if not found */
}
void GMT_free_dataset (struct GMT_CTRL *GMT, struct GMT_DATASET **data)
{	/* This takes pointer to data array and thus can return it as NULL */
	GMT_free_dataset_ptr (GMT, *data);
	GMT_free (GMT, *data);
}
int GMTAPI_Destroy_Dataset (struct GMTAPI_CTRL *API, struct GMT_DATASET **D_obj)
{
	/* Delete the given dataset resource. */

	if (!(*D_obj)) {	/* Probably not a good sign */
		//GMT_Report (API, GMT_MSG_DEBUG, "GMTAPI_Destroy_Dataset: Passed NULL pointer - skipped\n");
		return (GMT_PTR_IS_NULL);
	}
	if ((*D_obj)->alloc_level != API->GMT->hidden.func_level) return (GMT_FREE_WRONG_LEVEL);	/* Not the right level */

	GMT_free_dataset (API->GMT, D_obj);
	return GMT_OK;
}

void GMT_free_textset (struct GMT_CTRL *GMT, struct GMT_TEXTSET **data)
{	/* This takes pointer to data array and thus can return it as NULL */

	GMT_free_textset_ptr (GMT, *data);
	GMT_free (GMT, *data);
}
int GMTAPI_Destroy_Textset (struct GMTAPI_CTRL *API, struct GMT_TEXTSET **T_obj)
{
	/* Delete the given textset resource. */

	if (!(*T_obj)) {	/* Probably not a good sign */
		//GMT_Report (API, GMT_MSG_DEBUG, "GMTAPI_Destroy_Textset: Passed NULL pointer - skipped\n");
		return (GMT_PTR_IS_NULL);
	}
	if ((*T_obj)->alloc_level != API->GMT->hidden.func_level) return (GMT_FREE_WRONG_LEVEL);	/* Not the right level */

	GMT_free_textset (API->GMT, T_obj);
	return GMT_OK;
}

int GMTAPI_Destroy_Coord (struct GMTAPI_CTRL *API, double **ptr)
{
	GMT_free (API->GMT, *ptr);
	return GMT_OK;
}

int GMTAPI_Unregister_IO (struct GMTAPI_CTRL *API, int object_ID, unsigned int direction)
{	/* Remove specified object ID from active list of objects */
	int s_item;
	unsigned item;

	if (API == NULL) return (GMT_NOT_A_SESSION);		/* GMT_Create_Session has not been called */
	if (API->n_objects == 0) return (GMTAPI_report_error (API, GMT_NO_RESOURCES));	/* There are no known resources yet */

	/* Check if this is a valid ID and matches the direction */
	if ((s_item = GMTAPI_Validate_ID (API, GMT_NOTSET, object_ID, direction)) == GMT_NOTSET) return (GMTAPI_report_error (API, API->error));

	/* OK, now it is safe to remove the object; item >= 0 */

	item = s_item;
	//GMT_Report (API, GMT_MSG_DEBUG, "GMTAPI_Unregister_IO: Unregistering object no %d [n_objects = %d]\n", API->object[item]->ID, API->n_objects-1);
 	if (API->object[item]->data) //GMT_Report (API, GMT_MSG_DEBUG, "GMTAPI_Unregister_IO: Object no %d has non-NULL data pointer\n", API->object[item]->ID);
 	if (API->object[item]->resource) //GMT_Report (API, GMT_MSG_DEBUG, "GMTAPI_Unregister_IO: Object no %d has non-NULL resource pointer\n", API->object[item]->ID);

	if (API->object[item]->method == GMT_IS_FILE && API->object[item]->filename) free (API->object[item]->filename);	/* Free any strdup-allocated filenames */
	GMT_free (API->GMT, API->object[item]);		/* Free the current data object */
	API->n_objects--;				/* Tally of how many data sets are left */
	while (item < API->n_objects) {
		API->object[item] = API->object[item+1];	/* Shuffle pointers down one entry */
		item++;
	}

	/* All active resources are found consecutively from 0 to (API->n_objects-1); those with status == 0 (GMT_IS_UNUSED) are available for use. */
	return GMT_OK;
}

int GMT_Destroy_Data (void *V_API, void *object)
{
	/* Destroy a resource that is no longer needed.
	 * Returns the error code.
	 */
	int error, item, object_ID;
	enum GMT_enum_family family;
	struct GMTAPI_CTRL *API = NULL;

	if (V_API == NULL) return_error (V_API, GMT_NOT_A_SESSION);
	if (object == NULL) return (false);	/* Null address, quietly skip */
	if (!ptrvoid((char **)object)) return (false);	/* Null pointer, quietly skip */
	API = gmt_get_api_ptr ((struct GMTAPI_CTRL *)V_API);
	if ((object_ID = GMTAPI_get_objectID_from_data_ptr (API, object)) == GMT_NOTSET) return_error (API, GMT_OBJECT_NOT_FOUND);	/* Could not find it */
	if ((item = GMTAPI_Validate_ID (API, GMT_NOTSET, object_ID, GMT_NOTSET)) == GMT_NOTSET) return_error (API, API->error);

	family = (API->object[item]->actual_family) ? API->object[item]->actual_family : API->object[item]->family;	/* So if dataset via matrix we want family = matrix */
	switch (family) {	/* Standard 5 families, plus matrix/vector and coordinates */
		case GMT_IS_GRID:	/* GMT grid */
			error = GMTAPI_Destroy_Grid (API, (struct GMT_GRID **) object);
			break;
		case GMT_IS_DATASET:
			error = GMTAPI_Destroy_Dataset (API, (struct GMT_DATASET **)object);
			break;
		case GMT_IS_TEXTSET:
			error = GMTAPI_Destroy_Textset (API, (struct GMT_TEXTSET **)object);
			break;
		case GMT_IS_CPT:
			//error = GMTAPI_Destroy_CPT (API, object);
			break;
#ifdef HAVE_GDAL
		case GMT_IS_IMAGE:
			error = GMTAPI_Destroy_Image (API, object);
			break;
#endif

		/* Also allow destoying of intermediate vector and matrix containers */
		case GMT_IS_MATRIX:
			//error = GMTAPI_Destroy_Matrix (API, object);
			break;
		case GMT_IS_VECTOR:
			//error = GMTAPI_Destroy_Vector (API, object);
			break;
		case GMT_IS_COORD:
			error = GMTAPI_Destroy_Coord (API, (double **)object);
			break;
		default:
			return_error (API, GMT_NOT_A_VALID_FAMILY);
			break;
	}
	if (!error) {	/* We successfully freed the items, now remove from IO list */
		unsigned int j;
		void *address = API->object[item]->data;
		//GMT_Report (API, GMT_MSG_DEBUG, "GMT_Destroy_Data: freed memory for a %s for object %d\n", GMT_family[family], object_ID);
		if ((error = GMTAPI_Unregister_IO (API, object_ID, GMT_NOTSET))) return_error (API, error);	/* Did not find object */
		for (j = 0; j < API->n_objects; j++) {
			if (API->object[j]->data == address) API->object[j]->data = NULL;		/* Set repeated data references to NULL so we don't try to free twice */
			if (API->object[j]->resource == address) API->object[j]->resource = NULL;	/* Set matching resources to NULL so we don't try to read from there again */
		}
#ifdef DEBUG
		GMT_list_API (API, "GMT_Destroy_Data");
#endif
	}
	else {
		/* Quietly ignore these errors: GMT_PTR_IS_NULL, GMT_FREE_EXTERNAL_NOT_ALLOWED, GMT_FREE_WRONG_LEVEL as they are not considered errors here. */
		//GMT_Report (API, GMT_MSG_DEBUG, "GMT_Destroy_Data: Ignored warning %d for object %d\n", error, object_ID);
	}
	return (GMT_OK);	/* Returns number of items freed or an error */
}

int GMT_Register_IO (void *V_API, unsigned int family, unsigned int method, unsigned int geometry, unsigned int direction, double wesn[], void *resource)
{
	/* Adds a new data object to the list of registered objects and returns a unique object ID.
	 * Arguments are as listed for GMTAPI_Register_Im|Export (); see those for details.
	 * During the registration we make sure files exist and are readable.
	 *
	 * if direction == GMT_IN:
	 * A program uses this routine to pass information about input data to GMT.
	 * family:	Specifies the data type we are trying to import; select one of 5 families:
	 *   GMT_IS_CPT:	A GMT_PALETTE structure:
	 *   GMT_IS_DATASET:	A GMT_DATASET structure:
	 *   GMT_IS_TEXTSET:	A GMT_TEXTSET structure:
	 *   GMT_IS_GRID:	A GMT_GRID structure:
	 *   GMT_IS_IMAGE:	A GMT_IMAGE structure:
	 * method:	Specifies by what method we will import this data set:
	 *   GMT_IS_FILE:	A file name is given via input.  The program will read data from this file
	 *   GMT_IS_STREAM:	A file pointer to an open file is passed via input. --"--
	 *   GMT_IS_FDESC:	A file descriptor to an open file is passed via input. --"--
	 *   GMT_IS_DUPLICATE:	A pointer to a data set to be copied
	 *   GMT_IS_REFERENCE:	A pointer to a data set to be passed as is [we may reallocate sizes only if GMT-allocated]
	 * The following approaches can be added to the method for all but CPT:
	 *   GMT_VIA_MATRIX:	A 2-D user matrix is passed via input as a source for copying.
	 *			The GMT_MATRIX structure must have parameters filled out.
	 *   GMT_VIA_VECTOR:	An array of user column vectors is passed via input as a source for copying.
	 *			The GMT_VECTOR structure must have parameters filled out.
	 * geometry:	One of GMT_IS_{TEXT|POINT|LINE|POLY|SURF} (the last for GMT grids)
	 * input:	Pointer to the source filename, stream, handle, array position, etc.
	 * wesn:	Grid subset defined by 4 doubles; otherwise use NULL
	 * RETURNED:	Unique ID assigned to this input resouce, or GMT_NOTSET (-1) if error.
	 *
	 * An error status is returned if problems are encountered via API->error [GMT_OK].
	 *
	 * GMT_IS_GRID & GMT_VIA_MATRIX: Since GMT internally uses floats in C arrangement, anything else will be converted to float.
	 * GMT_IS_DATASET & GMT_VIA_MATRIX: Since GMT internally uses doubles in C arrangement, anything else will be converted to double.
	 *
	 * GMTAPI_Register_Import will allocate and populate a GMTAPI_DATA_OBJECT structure which
	 * is appended to the data list maintained by the GMTAPI_CTRL API structure.
	 *
	 * if direction == GMT_OUT:
	 * The main program uses this routine to pass information about output data from GMT.
	 * family:	Specifies the data type we are trying to export; select one of:
	 *   GMT_IS_CPT:	A GMT_PALETTE structure:
	 *   GMT_IS_DATASET:	A GMT_DATASET structure:
	 *   GMT_IS_TEXTSET:	A GMT_TEXTSET structure:
	 *   GMT_IS_IMAGE:	A GMT_IMAGE structure:
	 *   GMT_IS_GRID:	A GMT_GRID structure:
	 * method:	Specifies by what method we will export this data set:
	 *   GMT_IS_FILE:	A file name is given via output.  The program will write data to this file
	 *   GMT_IS_STREAM:	A file pointer to an open file is passed via output. --"--
	 *   GMT_IS_FDESC:	A file descriptor to an open file is passed via output. --"--
	 *   GMT_IS_DUPLICATE:	A pointer to a data set to be copied
	 *   GMT_IS_REFERENCE:	A pointer to a data set to be passed as is [we may reallocate sizes only if GMT-allocated]
	 * The following approaches can be added to the method for all but CPT:
	 *   GMT_VIA_MATRIX:	A 2-D user matrix is passed via input as a source for copying.
	 *			The GMT_MATRIX structure must have parameters filled out.
	 *   GMT_VIA_VECTOR:	An array of user column vectors is passed via input as a source for copying.
	 *			The GMT_VECTOR structure must have parameters filled out.
	 * geometry:	One of GMT_IS_{TEXT|POINT|LINE|POLY|SURF} (the last for GMT grids)
	 * output:	Pointer to the destination filename, stream, handle, array position, etc.
	 * wesn:	Grid subset defined by 4 doubles; otherwise use NULL
	 * RETURNED:	Unique ID assigned to this output resouce, or GMT_NOTSET (-1) if error.
	 *
	 * An error status is returned if problems are encountered via API->error [GMT_OK].
	 *
	 * GMTAPI_Register_Export will allocate and populate a GMTAPI_DATA_OBJECT structure which
	 * is appended to the data list maintained by the GMTAPI_CTRL API structure.
	 */
	int item, via = 0, m, object_ID;
	unsigned int mode = method & GMT_IO_RESET;	/* In case we wish to reuse this resource */
	char message[GMT_BUFSIZ];
	struct GMTAPI_DATA_OBJECT *S_obj = NULL;
	struct GMT_MATRIX *M_obj = NULL;
	struct GMT_VECTOR *V_obj = NULL;
	struct GMTAPI_CTRL *API = NULL;

	if (V_API == NULL) return_value (V_API, GMT_NOT_A_SESSION, GMT_NOTSET);
	API = gmt_get_api_ptr ((struct GMTAPI_CTRL *)V_API);
	if (GMTAPI_Validate_Geometry (API, family, geometry)) return_value (API, GMT_BAD_GEOMETRY, GMT_NOTSET);
	
	/* Check if this filename is an embedded API Object ID passed via the filename and of the right kind.  */
	if ((object_ID = GMTAPI_Memory_Registered (API, (enum GMT_enum_family) family, direction, resource)) != GMT_NOTSET)
	{
		printf("@@@@@@@@@@@GMT_Register_IO: %s %s %d   object_ID :%d \n",__FILE__ ,__func__,__LINE__,object_ID);
		return (object_ID);	/* OK, return the object ID */
	}
	if ((object_ID = GMTAPI_is_registered (API,(enum GMT_enum_family) family, geometry, direction, mode, NULL, resource)) != GMT_NOTSET) {	/* Registered before */
		printf("@@@@@@@@@@@GMT_Register_IO: %s %s %d   object_ID :%d \n",__FILE__ ,__func__,__LINE__,object_ID);
		if ((item = GMTAPI_Validate_ID (API, GMT_NOTSET, object_ID, direction)) == GMT_NOTSET) return_value (API, API->error, GMT_NOTSET);
		if ((family == GMT_IS_GRID || family == GMT_IS_IMAGE) && !full_region (wesn)) {	/* Update the subset region if given (for grids/images only) */
			S_obj = API->object[item];	/* Use S as shorthand */
			GMT_memcpy (S_obj->wesn, wesn, 4, double);
			S_obj->region = true;
		}
		return (object_ID);	/* Already registered so we are done */
	}
	method -= mode;	/* Remove GMT_IO_RESET if it was passed */
	if (method >= GMT_VIA_VECTOR) {
		via = (method / GMT_VIA_VECTOR) - 1;
		m = method - (via + 1) * GMT_VIA_VECTOR;	/* Array index that have any GMT_VIA_* removed */
	}
	else
		m = method;
	
	switch (method) {	/* Consider CPT, data, text, and grids, accessed via a variety of methods */
		case GMT_IS_FILE:	/* Registration via a single file name */
			/* No, so presumably it is a regular file name */
			if (direction == GMT_IN) {	/* For input we can check if the file exists and can be read. */
				char *p, *file = strdup ((const char *)resource);
				if ((family == GMT_IS_GRID || family == GMT_IS_IMAGE) && (p = strchr (file, '=')))
					{
					*p = '\0';	/* Chop off any =<stuff> for grids and images so access can work */
					}
				else if (family == GMT_IS_IMAGE && (p = strchr (file, '+'))) {
					char *c = strchr (file, '.');	/* The period before an extension */
					 /* PW 1/30/2014: Protect images with band requiest, e.g., my_image.jpg+b2 */
					if (c && c < p && p[1] == 'b' && isdigit (p[2])) {
						printf ("Truncating +b modifier for image filename %s\n", file);
						*p = '\0';	/* Chop off any +b<band> for images at end of extension so access can work */
					}
				}
				if (GMT_access (API->GMT, file, F_OK) && !GMT_check_url_name (file)) {	/* For input we can check if the file exists (except if via Web) */
					printf ("File %s not found\n", file);
					return_value (API, GMT_FILE_NOT_FOUND, GMT_NOTSET);
				}
				if (GMT_access (API->GMT, file, R_OK) && !GMT_check_url_name (file)) {	/* Found it but we cannot read. */
					printf ( "Not permitted to read file %s\n", file);
					return_value (API, GMT_BAD_PERMISSION, GMT_NOTSET);
				}
				free (file);
			}
			else if (resource == NULL) {	/* No file given [should this mean stdin/stdout?] */
				
				return_value (API, GMT_OUTPUT_NOT_SET, GMT_NOTSET);
			}
			/* Create a new data object and initialize variables */
			if ((S_obj = GMTAPI_Make_DataObject (API, (enum GMT_enum_family)family, method, geometry, NULL, direction)) == NULL) {
				
				return_value (API, GMT_MEMORY_ERROR, GMT_NOTSET);	/* No more memory */
			}
			if (strlen ((const char *)resource)) S_obj->filename = strdup ((const char *)resource);
			sprintf (message, "Object ID %%d : Registered %s %s %s as an %s resource with geometry %s [n_objects = %%d]\n", GMT_family[family], GMT_method[m], S_obj->filename, GMT_direction[direction], GMT_geometry[gmtry(geometry)]);
			
			break;

		case GMT_IS_STREAM:	/* Methods that indirectly involve a file */
		case GMT_IS_FDESC:
			if (resource == NULL) {	/* No file given [should this mean stdin/stdout?] */
				
				return_value (API, GMT_OUTPUT_NOT_SET, GMT_NOTSET);
			}
			if ((S_obj = GMTAPI_Make_DataObject (API, (enum GMT_enum_family) family, method, geometry, NULL, direction)) == NULL) {
				
				return_value (API, GMT_MEMORY_ERROR, GMT_NOTSET);	/* No more memory */
			}
			S_obj->fp = resource;	/* Pass the stream of fdesc onward */
			sprintf (message, "Object ID %%d : Registered %s %s "/*%" PRIxS*/ " as an %s resource with geometry %s [n_objects = %%d]\n", GMT_family[family], GMT_method[m], (size_t)resource, GMT_direction[direction], GMT_geometry[gmtry(geometry)]);
			break;

		case GMT_IS_DUPLICATE:
		case GMT_IS_REFERENCE:
#if 0
			if (direction == GMT_OUT && resource != NULL) {
				return_value (API, GMT_PTR_NOT_NULL, GMT_NOTSET);	/* Output registration of memory takes no resource */
			} else
#endif
			if (direction == GMT_IN && resource == NULL) {
				
				return_value (API, GMT_PTR_IS_NULL, GMT_NOTSET);	/* Input registration of memory takes a resource */
			}
			if ((S_obj = GMTAPI_Make_DataObject (API, (enum GMT_enum_family) family, method, geometry, resource, direction)) == NULL) {
			
				return_value (API, GMT_MEMORY_ERROR, GMT_NOTSET);	/* No more memory */
			}
			sprintf(message, "Object ID %%d : Registered %s %s "/*%" PRIxS*/" as an %s resource with geometry %s [n_objects = %%d]\n", GMT_family[family], GMT_method[m], (size_t)resource, GMT_direction[direction], GMT_geometry[gmtry(geometry)]);
			break;

		 case GMT_IS_DUPLICATE + GMT_VIA_MATRIX:	/* Here, a data grid is passed via a GMT_MATRIX structure */
		 case GMT_IS_REFERENCE + GMT_VIA_MATRIX:
			if ((M_obj = resource) == NULL) {
				
				return_value (API, GMT_PTR_IS_NULL, GMT_NOTSET);	/* Matrix container must be given for both input and output */
			}
			if (direction == GMT_IN) {	/* For input we can check if the GMT_MATRIX structure has proper parameters. */
				if (M_obj->n_rows <= 0 || M_obj->n_columns <= 0) {
					//GMT_Report (API, GMT_MSG_NORMAL, "Error in GMT_Register_IO (%s): Matrix dimensions not set.\n", GMT_direction[direction]);
					printf("@@@@@@@@@@@GMT_Register_IO: %s %s %d   object_ID :%d \n",__FILE__ ,__func__,__LINE__,object_ID);
					return_value (API, GMT_NO_PARAMETERS, GMT_NOTSET);
				}
			}
			if ((S_obj = GMTAPI_Make_DataObject (API, (enum GMT_enum_family) family, method, geometry, resource, direction)) == NULL) {
				printf("@@@@@@@@@@@GMT_Register_IO: %s %s %d   object_ID :%d \n",__FILE__ ,__func__,__LINE__,object_ID);
				return_value (API, GMT_MEMORY_ERROR, GMT_NOTSET);	/* No more memory */
			}
			API->GMT->common.b.active[direction] = true;
			sprintf (message, "Object ID %%d : Registered %s %s "/*%" PRIxS*/" via %s as an %s resource with geometry %s [n_objects = %%d]\n", GMT_family[family], GMT_method[m], (size_t)resource, GMT_via[via], GMT_direction[direction], GMT_geometry[gmtry(geometry)]);
			break;
		 case GMT_IS_DUPLICATE + GMT_VIA_VECTOR:	/* Here, some data vectors are passed via a GMT_VECTOR structure */
		 case GMT_IS_REFERENCE + GMT_VIA_VECTOR:
			if ((V_obj = resource) == NULL) {
				printf("@@@@@@@@@@@GMT_Register_IO: %s %s %d   object_ID :%d \n",__FILE__ ,__func__,__LINE__,object_ID);
				return_value (API, GMT_PTR_IS_NULL, GMT_NOTSET);	/* Vector container must be given for both input and output */
			}
			if (direction == GMT_IN) {	/* For input we can check if the GMT_VECTOR structure has proper parameters. */
				if (V_obj->n_rows <= 0 || V_obj->n_columns <= 0) {
					//GMT_Report (API, GMT_MSG_NORMAL, "Error in GMT_Register_IO (%s): Vector parameters not set.\n", GMT_direction[direction]);
					printf("@@@@@@@@@@@GMT_Register_IO: %s %s %d   object_ID :%d \n",__FILE__ ,__func__,__LINE__,object_ID);
					return_value (API, GMT_NO_PARAMETERS, GMT_NOTSET);
				}
			}
			if ((S_obj = GMTAPI_Make_DataObject (API, (enum GMT_enum_family) family, method, geometry, resource, direction)) == NULL) {
				printf("@@@@@@@@@@@GMT_Register_IO: %s %s %d   object_ID :%d \n",__FILE__ ,__func__,__LINE__,object_ID);
				return_value (API, GMT_MEMORY_ERROR, GMT_NOTSET);	/* No more memory */
			}
			API->GMT->common.b.active[direction] = true;
			sprintf (message, "Object ID %%d : Registered %s %s "/*%" PRIxS*/" via %s as an %s resource with geometry %s [n_objects = %%d]\n", GMT_family[family], GMT_method[m], (size_t)resource, GMT_via[via], GMT_direction[direction], GMT_geometry[gmtry(geometry)]);
			break;
		case GMT_IS_COORD:	/* Internal registration of coordinate arrays so that GMT_Destroy_Data can free them */
			if ((S_obj = GMTAPI_Make_DataObject (API, (enum GMT_enum_family) family, method, geometry, resource, direction)) == NULL) {
				printf("@@@@@@@@@@@GMT_Register_IO: %s %s %d   object_ID :%d \n",__FILE__ ,__func__,__LINE__,object_ID);
				return_value (API, GMT_MEMORY_ERROR, GMT_NOTSET);	/* No more memory */
			}
			sprintf (message, "Object ID %%d : Registered double array "/*%" PRIxS*/" as an %s resource [n_objects = %%d]\n", (size_t)resource, GMT_direction[direction]);
			break;
		default:
			//GMT_Report (API, GMT_MSG_NORMAL, "Error in GMT_Register_IO (%s): Unrecognized method %d\n", GMT_direction[direction], method);
			printf("@@@@@@@@@@@GMT_Register_IO: %s %s %d   object_ID :%d \n",__FILE__ ,__func__,__LINE__,object_ID);
			return_value (API, GMT_NOT_A_VALID_METHOD, GMT_NOTSET);
			break;
	}

	if (!full_region (wesn)) {	/* Copy the subset region if it was given (for grids) */
		GMT_memcpy (S_obj->wesn, wesn, 4, double);
		S_obj->region = true;
	}

	S_obj->alloc_level = API->GMT->hidden.func_level;	/* Object was allocated at this module nesting level */

	/* Here S is not NULL and no errors have occurred (yet) */

	if (method != GMT_IS_COORD) API->registered[direction] = true;	/* We have at least registered one item */
	object_ID = GMTAPI_Add_Data_Object (API, S_obj);
	//GMT_Report (API, GMT_MSG_DEBUG, message, object_ID, API->n_objects);
#ifdef DEBUG
	GMT_list_API (API, "GMT_Register_IO");
#endif
	
	return_value (API, API->error, object_ID);
}

struct GMT_DATASET * GMT_create_dataset (struct GMT_CTRL *GMT, uint64_t n_tables, uint64_t n_segments, uint64_t n_rows, uint64_t n_columns, unsigned int geometry, bool alloc_only)
{	/* Create an empty data set structure with the required number of empty tables, all set to hold n_segments with n_columns */
	uint64_t tbl;
	struct GMT_DATASET *D = NULL;

	D = GMT_memory (GMT, NULL, 1, struct GMT_DATASET);
	if (n_columns) {
		D->min = GMT_memory (GMT, NULL, n_columns, double);
		D->max = GMT_memory (GMT, NULL, n_columns, double);
	}
	D->n_columns = n_columns;
	D->geometry = geometry;
	if (n_tables) D->table = GMT_memory (GMT, NULL, n_tables, struct GMT_DATATABLE *);
	D->n_alloc = D->n_tables = n_tables;
	if (!alloc_only) D->n_segments = D->n_tables * n_segments;
	if (!alloc_only) D->n_records = D->n_segments * n_rows;
	for (tbl = 0; tbl < n_tables; tbl++) if ((D->table[tbl] = GMT_create_table (GMT, n_segments, n_rows, n_columns, alloc_only)) == NULL) return (NULL);
	D->alloc_level = GMT->hidden.func_level;	/* Must be freed at this level. */
	D->alloc_mode = GMT_ALLOCATED_BY_GMT;		/* So GMT_* modules can free this memory. */
	D->id = GMT->parent->unique_var_ID++;		/* Give unique identifier */

	return (D);
}

void * GMT_Create_Data (void *V_API, unsigned int family, unsigned int geometry, unsigned int mode, uint64_t dim[], double *range, double *inc, unsigned int registration, int pad, void *data)
{
	/* Create an empty container of the requested kind and allocate space for content.
	 * The known families are GMT_IS_{DATASET,TEXTSET,GRID,CPT,IMAGE}, but we
	 * also allow for creation of the containers for GMT_IS_{VECTOR,MATRIX}. Note
	 * that for VECTOR|MATRIX we dont allocate space to hold data as it is the users
	 * responsibility to hook their data pointers in.  The VECTOR allocates the array
	 * of column vector type and data pointers.
	 * geometry should reflect the resource, e.g. GMT_IS_SURFACE for grid, etc.
	 * There are two ways to define the dimensions needed to actually allocate memory:
	 * (A) Via uint64_t dim[]:
	 *   The dim array contains up to 4 dimensions for:
	 *	0: par[GMT_TBL] = number of tables,
	 *	1: par[GMT_SEG] = number of segments per table
	 *	2: par[GMT_ROW] = number of rows per segment.
	 *	3: par[GMT_COL] = number of columns per row [ignored for GMT_TEXTSET].
	 * The dim array is ignored for CPT and GMT grids.
	 *   For GMT_IS_VECTOR, par[0] holds the number of columns.
	 *   For GMT_IS_MATRIX, par[GMT_Z] = GMT[2] holds the number of layers (dim == NULL means just 1 layer).
	 * (B) Via range, inc, registration:
	 *   Convert user domain range, increments, and registration into dimensions
	 *   for the container.  For grids and images we fill out the GMT_GRID_HEADER;
	 *   for vectors and matrices we fill out their internal parameters.
	 *   For complex grids pass registration + GMT_GRID_IS_COMPLEX_{REAL|IMAG}
	 *   For GMT_IS_MATRIX, par[GMT_Z] = holds the number of layers (dim == NULL means just 1 layer).
	 * pad sets the padding for grids and images, ignored for other resources.
	 * Some default actions for grids:
	 * range = NULL: Select current -R setting if present.
	 * registration = -1: Gridline unless -r is in effect.
	 * Give -1 (GMT_NOTSET) to accept GMT default padding [2].
	 *
	 * For creating grids and images you can do it in one or two steps:
 	 * (A) Pass mode = GMT_GRID_ALL; this creates both header and allocates grid|image;
	 * (B) Call GMT_Create_Data twice:
	 * 	1. First with mode = GMT_GRID_HEADER_ONLY which creates header only
	 *	   and computes the dimensions based on the other arguments.
	 *	2. 2nd with mode = GMT_GRID_DATA_ONLY, which allocates the grid|image array
	 *	   based on the dimensions already set.  This time you pass NULL/0
	 *	   for dim, wesn,inc,registration,pad but let data be your grid|image returned
	 *	   to you after step 1.
	 * By default, the created resource is consider an input resource (direction == GMT_IN).
	 * However, you can change that by adding GMT_VIA_OUTPUT to the mode. This means you
	 * are passing in a blank container that a GMT module will fill in with results so
	 * that you can access this data outside of the module.
	 *
	 * Return: Pointer to resource, or NULL if an error (set via API->error).
	 */

	int error = GMT_OK, object_ID = GMT_NOTSET;
	uint64_t n_layers = 0;
	bool already_registered = false, has_ID = false;
	void *new_obj = NULL, *p_data = NULL;
	struct GMTAPI_CTRL *API = NULL;

	if (V_API == NULL) return_null (V_API, GMT_NOT_A_SESSION);
	API = gmt_get_api_ptr (V_API);

		p_data = data;	/* data can only be non-NULL for Grids/Images passing back G to get the data array */

	switch (family) {	/* dataset, cpt, text, grid , image, vector, matrix */
		case GMT_IS_GRID:	/* GMT grid, allocate header but not data array */
			if ((mode & GMT_GRID_DATA_ONLY) == 0) {	/* Create new grid unless we only ask for data only */
	 			if ((new_obj = GMT_create_grid (API->GMT)) == NULL)
					return_null (API, GMT_MEMORY_ERROR);	/* Allocation error */

				if (pad >= 0) GMT_set_pad (API->GMT, pad);	/* Change the default pad; give -1 to leave as is */
				
				error = GMTAPI_init_grid (API->GMT, NULL, range, inc, registration, mode, new_obj);
				if (pad >= 0) GMT_set_pad (API->GMT, API->pad);	/* Reset to the default pad */
			}
			else {	/* Already registered so has_ID must be false */
				if (has_ID || (new_obj = p_data) == NULL) return_null (API, GMT_PTR_IS_NULL);	/* Error if data is NULL */
				already_registered = true;
			}
			if ((mode & GMT_GRID_HEADER_ONLY) == 0) {	/* Allocate the grid array unless we asked for header only */
				if ((error = gmt_alloc_grid (API->GMT, new_obj)) != GMT_NOERROR) return_null (API, error);	/* Allocation error */
			}
			break;

		case GMT_IS_DATASET:	/* GMT dataset, allocate the requested tables, segments, rows, and columns */
			if (dim == NULL) return_null (API, GMT_PTR_IS_NULL);
			if (dim[GMT_TBL] > UINT_MAX || dim[GMT_ROW] > UINT_MAX) return_null (API, GMT_DIM_TOO_LARGE);
			if ((new_obj = GMT_create_dataset (API->GMT, dim[GMT_TBL], dim[GMT_SEG], dim[GMT_ROW], dim[GMT_COL], geometry, false)) == NULL) return_null (API, GMT_MEMORY_ERROR);	/* Allocation error */
			break;
		case GMT_IS_TEXTSET:	/* GMT text dataset, allocate the requested tables, segments, and rows */
			if (dim == NULL) return_null (API, GMT_PTR_IS_NULL);
			if (dim[GMT_TBL] > UINT_MAX) return_null (API, GMT_DIM_TOO_LARGE);
			if ((new_obj = GMT_create_textset (API->GMT, dim[GMT_TBL], dim[GMT_SEG], dim[GMT_ROW], false)) == NULL) return_null (API, GMT_MEMORY_ERROR);	/* Allocation error */
			break;
		//case GMT_IS_CPT:	/* GMT CPT table, allocate one with space for dim[0] color entries */
			//if (dim == NULL) return_null (API, GMT_PTR_IS_NULL);
		 //	if ((new_obj = GMT_create_palette (API->GMT, dim[0])) == NULL) return_null (API, GMT_MEMORY_ERROR);	/* Allocation error */
		//	break;
		case GMT_IS_MATRIX:	/* GMT matrix container, allocate one with the requested number of layers, rows & columns */
			//n_layers = (dim == NULL || dim[0] == 0) ? 1U : dim[GMT_Z];
		 	//new_obj = GMT_create_matrix (API->GMT, n_layers);
			//if (pad) //GMT_Report (API, GMT_MSG_VERBOSE, "Pad argument (%d) ignored in initialization of %s\n", pad, GMT_family[family]);
			//if ((API->error = GMTAPI_init_matrix (API, dim, range, inc, registration, mode, new_obj))) {	/* Failure, must free the object */
				//struct GMT_MATRIX *M = return_address (new_obj, GMT_IS_MATRIX);	/* Get pointer to resource */
				//GMT_free_matrix (API->GMT, &M, true);
			//}
			//break;
		case GMT_IS_VECTOR:	/* GMT vector container, allocate one with the requested number of columns & rows */
			//if (dim == NULL) return_null (API, GMT_PTR_IS_NULL);
			//if (dim[0] == 0) return_null (API, GMT_VALUE_NOT_SET);
	 		//new_obj = GMT_create_vector (API->GMT, dim[0]);
			//if (pad) //GMT_Report (API, GMT_MSG_VERBOSE, "Pad argument (%d) ignored in initialization of %s\n", pad, GMT_family[family]);
			//if ((API->error = GMTAPI_init_vector (API, dim, range, inc, registration, new_obj))) {	/* Failure, must free the object */
				//struct GMT_VECTOR *V = return_address (new_obj, GMT_IS_VECTOR);	/* Get pointer to resource */
				//GMT_free_vector (API->GMT, &V, true);
			//}
			//break;
		default:
			//API->error = GMT_NOT_A_VALID_FAMILY;
			break;
	}
	if (API->error) return_null (API, API->error);

	if (!already_registered) {	/* Now register this dataset so it can be deleted by GMT_Destroy_Data or GMT_Garbage_Collection */
		enum GMT_enum_via via = GMT_VIA_NONE;
		enum GMT_enum_family actual_family = family;
		int direction, def_direction = GMT_IN;	/* Default direction is GMT_IN unless mode & GMT_VIA_OUTPUT was passed */
		int item = GMT_NOTSET;
		if (mode & GMT_VIA_OUTPUT) def_direction = GMT_OUT;	/* Create item for output instead */
		direction = (object_ID == GMT_NOTSET) ? def_direction : GMT_NOTSET;	/* Do not consider direction if pre-registered */
		if (object_ID == GMT_NOTSET) {	/* Must register this object */
			if (family == GMT_IS_MATRIX) {
				if (geometry & GMT_IS_PLP) actual_family = GMT_IS_DATASET, via = GMT_VIA_MATRIX;
				if (geometry & GMT_IS_SURFACE) actual_family = GMT_IS_GRID, via = GMT_VIA_MATRIX;
			}
			else if (family == GMT_IS_VECTOR) {
				actual_family = GMT_IS_DATASET, via = GMT_VIA_VECTOR;
			}
			if ((object_ID = GMT_Register_IO (API, actual_family, GMT_IS_REFERENCE + via, geometry, def_direction, range, new_obj)) == GMT_NOTSET) return_null (API, API->error);	/* Failure to register */
		}
		if ((item = GMTAPI_Validate_ID (API, actual_family, object_ID, direction)) == GMT_NOTSET) return_null (API, API->error);
		API->object[item]->data = new_obj;	/* Retain pointer to the allocated data so we use garbage collection later */
		API->object[item]->actual_family = family;	/* So that if we got a matrix posing as dataset we can destroy the matrix later */
		if (mode & GMT_VIA_OUTPUT) API->object[item]->messenger = true;	/* We are passing a dummy container that should be destroyed before returning data */
		printf ( "Successfully created a new %s container\n", GMT_family[family]);
	}
	else
		printf ( "Successfully added data array to previously registered %s container\n", GMT_family[family]);

	return (new_obj);
}
/* Sorted array with information for all GMT core modules */

/* name, library, and purpose for each module */
struct Gmt_moduleinfo {
	const char *name;             /* Program name */
	const char *component;        /* Component (core, supplement, custom) */
	const char *purpose;          /* Program purpose */
#ifndef BUILD_SHARED_LIBS
	/* gmt module function pointer: */
	int (*p_func)(void*, int, void*);
#endif
};

struct Gmt_moduleinfo g_core_module[] = {
#ifdef BUILD_SHARED_LIBS
	{"blockmean", "core", "Block average (x,y,z) data tables by L2 norm"},
	{"blockmedian", "core", "Block average (x,y,z) data tables by L1 norm (spatial median)"},
	{"blockmode", "core", "Block average (x,y,z) data tables by mode estimation"},
	{"filter1d", "core", "Do time domain filtering of 1-D data tables"},
	{"fitcircle", "core", "Find mean position and best-fitting great- or small-circle to points on sphere"},
	{"gmt2kml", "core", "Convert GMT data tables to KML files for Google Earth"},
	{"gmtconnect", "core", "Connect individual lines whose end points match within tolerance"},
	{"gmtconvert", "core", "Convert, paste, or extract columns from data tables"},
	{"gmtdefaults", "core", "List current GMT default parameters"},
	{"gmtget", "core", "Get individual GMT default parameters"},
	{"gmtinfo", "core", "Get information about data tables"},
	{"gmtmath", "core", "Reverse Polish Notation (RPN) calculator for data tables"},
	{"gmtselect", "core", "Select data table subsets based on multiple spatial criteria"},
	{"gmtset", "core", "Change individual GMT default parameters"},
	{"gmtsimplify", "core", "Line reduction using the Douglas-Peucker algorithm"},
	{"gmtspatial", "core", "Do geospatial operations on lines and polygons"},
	{"gmtvector", "core", "Basic manipulation of Cartesian vectors"},
	{"gmtwhich", "core", "Find full path to specified files"},
	{"grd2cpt", "core", "Make linear or histogram-equalized color palette table from grid"},
	{"grd2rgb", "core", "Write r/g/b grid files from a grid file, a raw RGB file, or SUN rasterfile"},
	{"grd2xyz", "core", "Convert grid file to data table"},
	{"grdblend", "core", "Blend several partially over-lapping grids into one larger grid"},
	{"grdclip", "core", "Clip the range of grids"},
	{"grdcontour", "core", "Make contour map using a grid"},
	{"grdcut", "core", "Extract subregion from a grid"},
	{"grdedit", "core", "Modify header or content of a grid"},
	{"grdfft", "core", "Do mathematical operations on grids in the wavenumber (or frequency) domain"},
	{"grdfilter", "core", "Filter a grid in the space (or time) domain"},
	{"grdgradient", "core", "Compute directional gradients from a grid"},
	{"grdhisteq", "core", "Perform histogram equalization for a grid"},
	{"grdimage", "core", "Project grids or images and plot them on maps"},
	{"grdinfo", "core", "Extract information from grids"},
	{"grdlandmask", "core", "Create a \"wet-dry\" mask grid from shoreline data base"},
	{"grdmask", "core", "Create mask grid from polygons or point coverage"},
	{"grdmath", "core", "Reverse Polish Notation (RPN) calculator for grids (element by element)"},
	{"grdpaste", "core", "Join two grids along their common edge"},
	{"grdproject", "core", "Forward and inverse map transformation of grids"},
	{"grdraster", "core", "Extract subregion from a binary raster and save as a GMT grid"},
	{"grdreformat", "core", "Convert between different grid formats"},
	{"grdsample", "core", "Resample a grid onto a new lattice"},
	{"grdtrack", "core", "Sample grids at specified (x,y) locations"},
	{"grdtrend", "core", "Fit trend surface to grids and compute residuals"},
	{"grdvector", "core", "Plot vector field from two component grids"},
	{"grdview", "core", "Create 3-D perspective image or surface mesh from a grid"},
	{"grdvolume", "core", "Calculate grid volume and area constrained by a contour"},
	{"greenspline", "core", "Interpolate using Green's functions for splines in 1-3 dimensions"},
	{"kml2gmt", "core", "Extract GMT table data from Google Earth KML files"},
	{"makecpt", "core", "Make GMT color palette tables"},
	{"mapproject", "core", "Do forward and inverse map transformations, datum conversions and geodesy"},
	{"nearneighbor", "core", "Grid table data using a \"Nearest neighbor\" algorithm"},
	{"project", "core", "Project table data onto lines or great circles, generate tracks, or translate coordinates"},
	{"ps2raster", "core", "Convert [E]PS file(s) to other formats using GhostScript"},
	{"psbasemap", "core", "Plot PostScript base maps"},
	{"psclip", "core", "Initialize or terminate polygonal clip paths"},
	{"pscoast", "core", "Plot continents, countries, shorelines, rivers, and borders on maps"},
	{"pscontour", "core", "Contour table data by direct triangulation"},
	{"pshistogram", "core", "Calculate and plot histograms"},
	{"psimage", "core", "Place images or EPS files on maps"},
	{"pslegend", "core", "Plot legends on maps"},
	{"psmask", "core", "Use data tables to clip or mask map areas with no coverage"},
	{"psrose", "core", "Plot a polar histogram (rose, sector, windrose diagrams)"},
	{"psscale", "core", "Plot a gray-scale or color-scale on maps"},
	{"pstext", "core", "Plot or typeset text on maps"},
	{"pswiggle", "core", "Plot z = f(x,y) anomalies along tracks"},
	{"psxyz", "core", "Plot lines, polygons, and symbols in 3-D"},
	{"psxy", "core", "Plot lines, polygons, and symbols on maps"},
	{"read", "core", "Read GMT objects into external API"},
	{"sample1d", "core", "Resample 1-D table data using splines"},
	{"spectrum1d", "core", "Compute auto- [and cross-] spectra from one [or two] timeseries"},
	{"sph2grd", "core", "Compute grid from spherical harmonic coefficients"},
	{"sphdistance", "core", "Make grid of distances to nearest points on a sphere"},
	{"sphinterpolate", "core", "Spherical gridding in tension of data on a sphere"},
	{"sphtriangulate", "core", "Delaunay or Voronoi construction of spherical lon,lat data"},
	{"splitxyz", "core", "Split xyz[dh] data tables into individual segments"},
	{"surface", "core", "Grid table data using adjustable tension continuous curvature splines"},
	{"trend1d", "core", "Fit a [weighted] [robust] polynomial [or Fourier] model for y = f(x) to xy[w] data"},
	{"trend2d", "core", "Fit a [weighted] [robust] polynomial for z = f(x,y) to xyz[w] data"},
	{"triangulate", "core", "Do optimal (Delaunay) triangulation and gridding of Cartesian table data"},
	{"write", "core", "Write GMT objects from external API"},
	{"xyz2grd", "core", "Convert data table to a grid file"},
	{NULL, NULL, NULL} /* last element == NULL detects end of array */
#else

	//{"grd2cpt", "core", "Make linear or histogram-equalized color palette table from grid", &GMT_grd2cpt},
	//{"grd2rgb", "core", "Write r/g/b grid files from a grid file, a raw RGB file, or SUN rasterfile", &GMT_grd2rgb},
	//{"grd2xyz", "core", "Convert grid file to data table", &GMT_grd2xyz},
	//{"splitxyz", "core", "Split xyz[dh] data tables into individual segments", &GMT_splitxyz},
	{"surface", "core", "Grid table data using adjustable tension continuous curvature splines", &GMT_surface},
	//{"write", "core", "Write GMT objects from external API", &GMT_write},
	//{"xyz2grd", "core", "Convert data table to a grid file", &GMT_xyz2grd},
	{NULL, NULL, NULL, NULL}  /*last element == NULL detects end of array */
#endif
};
#if 0
int GMT_Call_Module (void *V_API, const char *module, int mode, void *args)
{	/* Call the specified shared module and pass it the mode and args.
 	 * mode can be one of the following:
	 * GMT_MODULE_EXIST [-3]:	Return GMT_NOERROR (0) if module exists, GMT_NOT_A_VALID_MODULE otherwise.
	 * GMT_MODULE_PURPOSE [-2]:	As GMT_MODULE_EXIST, but also print the module purpose.
	 * GMT_MODULE_OPT [-1]:		Args is a linked list of option structures.
	 * GMT_MODULE_CMD [0]:		Args is a single textstring with multiple options
	 * mode > 0:			Args is an array of text strings (argv[]).
	 */
	int status = GMT_NOERROR;
	unsigned int lib;
	struct GMTAPI_CTRL *API = NULL;
	char gmt_module[GMT_LEN32] = "GMT_";
	int (*p_func)(void*, int, void*) = NULL;       /* function pointer */

	if (V_API == NULL) return_error (V_API, GMT_NOT_A_SESSION);
	API = gmt_get_api_ptr ((struct GMTAPI_CTRL *)V_API);
	
	/* Here we call a named module */

	strncat (gmt_module, module, GMT_LEN32-5);		/* Concatenate GMT_-prefix and module name to get function name */
	for (lib = 0; lib < API->n_shared_libs; lib++) {	/* Look for gmt_module in any of the shared libs */
		*(void **) (&p_func) = gmt_get_module_func (API, gmt_module, lib);
		if (p_func) break;	/* Found it in this shared library */
	}
	if (p_func == NULL) {	/* Not in any of the shared libraries */
		//GMT_Report (API, GMT_MSG_VERBOSE, "Shared GMT module not found: %s \n", module);
		status = GMT_NOT_A_VALID_MODULE;
	}
	else if (mode == GMT_MODULE_EXIST)	/* Just wanted to know it is there */
		return (GMT_NOERROR);
	else	/* Call the function and return its return value */
		status = (*p_func) (V_API, mode, args);
	return (status);
}
#endif
void * GMT_Retrieve_Data (void *V_API, int object_ID)
{
	/* Function to return pointer to the container for a registered data set.
	 * Typically used when we wish a module to "write" its results to a memory
	 * location that we wish to access from the calling program.  The procedure
	 * is to use GMT_Register_IO with GMT_REF|COPY|READONLY and GMT_OUT but use
	 * NULL as the source/destination.  Data are "written" by GMT allocating a
	 * output container and updating the objects->resource pointer to this container.
	 * GMT_Retrieve_Data simply returns that pointer given the registered ID.
	 */

	int item;
	struct GMTAPI_CTRL *API = NULL;

	if (V_API == NULL) return_null (V_API, GMT_NOT_A_SESSION);

	/* Determine the item in the object list that matches this object_ID */
	API = gmt_get_api_ptr (V_API);
	if ((item = GMTAPI_Validate_ID (API, GMT_NOTSET, object_ID, GMT_NOTSET)) == GMT_NOTSET) {
		return_null (API, API->error);
	}
	/* Make sure the resource is present */
	if (API->object[item]->resource == NULL) {
		return_null (API, GMT_PTR_IS_NULL);
	}
	/* Make sure the data pointer has not been set */
	if (API->object[item]->data) {
		return_null (API, GMT_PTR_NOT_NULL);
	}
	/* Assign data from resource and wipe resource pointer */
	API->object[item]->data = API->object[item]->resource;
	API->object[item]->resource = NULL;

#ifdef DEBUG
	GMT_list_API (API, "GMT_Retrieve_Data");
#endif
	return (API->object[item]->data);	/* Return pointer to the data container */
}
/* Also used in gmt_io.c and prototyped in gmt_internals.h: */
char * GMT_create_header_item (struct GMTAPI_CTRL *API, unsigned int mode, void *arg)
{
	char *txt = (mode & GMT_COMMENT_IS_OPTION) ? GMT_Create_Cmd (API, arg) : (char *)arg;
	static char buffer[GMT_BUFSIZ];
	GMT_memset (buffer, GMT_BUFSIZ, char);
	if (mode & GMT_COMMENT_IS_TITLE) strcat (buffer, "  Title :");
	if (mode & GMT_COMMENT_IS_COMMAND) {
		strcat (buffer, " Command : ");
		strcat (buffer, API->GMT->init.module_name);
		strcat (buffer, " ");
	}
	if (mode & GMT_COMMENT_IS_REMARK) strcat (buffer, " Remark : ");
	strcat (buffer, txt);
	if (mode & GMT_COMMENT_IS_OPTION) GMT_free (API->GMT, txt);
	return (buffer);
}

