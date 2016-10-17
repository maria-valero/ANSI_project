#include <netcdf.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <unistd.h>
#include <math.h>
#include <float.h>
#include <dirent.h>
#include "gmt_macros.h"
#include "memory.h"
#include "gmt_unit.h"
#include "gmt_type.h"
#include "gmt_define.h"
#include "gmt_defaults.h"
#include "grd_io.h"
#include "gmt_grdio.h"
#include "common_byteswap.h"
/* Macro to test if filename is a special name indicating memory location */
#define GMT_LEN_UNITS2	"efkMnu"	/* Distances in meter, foot, survey foot, km, Mile, nautical mile */



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

//#define access _access
enum Grid_packing_mode {
	k_grd_pack = 0, /* scale and offset before writing to disk */
	k_grd_unpack    /* remove scale and offset after reading packed data */
};

	struct GRD_PAD {
		double wesn[4];
		unsigned int pad[4];
	};

	struct GMT_GRID *GMT_duplicate_grid (struct GMT_CTRL *GMT, struct GMT_GRID *G, unsigned int mode)
	{	/* Duplicates an entire grid, including data if requested. */
		struct GMT_GRID *Gnew = NULL;

		Gnew = GMT_create_grid (GMT);
		GMT_memcpy (Gnew->header, G->header, 1, struct GMT_GRID_HEADER);
		if ((mode & GMT_DUPLICATE_DATA) || (mode & GMT_DUPLICATE_ALLOC)) {	/* Also allocate and possiblhy duplicate data array */
			Gnew->data = GMT_memory_aligned (GMT, NULL, G->header->size, float);
			if (mode & GMT_DUPLICATE_DATA) GMT_memcpy (Gnew->data, G->data, G->header->size, float);
		}
		return (Gnew);
	}

int GMT_getincn (struct GMT_CTRL *GMT, char *line, double inc[], unsigned int n)
{
	unsigned int last, i, pos;
	char p[GMT_BUFSIZ];
	double scale = 1.0;

	/* Deciphers dx/dy/dz/dw/du/dv/... increment strings with n items */

	
	
	if (!line) { printf ("No argument given to GMT_getincn\n");
		GMT_exit (GMT, EXIT_FAILURE); // return EXIT_FAILURE; 
	}

	GMT_memset (inc, n, double);
	
	i = pos = GMT->current.io.inc_code[GMT_X] = GMT->current.io.inc_code[GMT_Y] = 0;

	
	while (i < n && (GMT_strtok (line, "/", &pos, p))) {
		last = (unsigned int)strlen (p) - 1;
		
		if (p[last] == '=') {	/* Let -I override -R */
			p[last] = 0;
			if (i < 2) GMT->current.io.inc_code[i] |= GMT_INC_IS_EXACT;
			last--;
		}
		else if (p[last] == '+' || p[last] == '!') {	/* Number of nodes given, determine inc from domain (! added since documentation mentioned this once... */
			p[last] = 0;
			if (i < 2) GMT->current.io.inc_code[i] |= GMT_INC_IS_NNODES;
			last--;
		}
		switch (p[last]) {
			case 'd':	/* Gave arc degree */
				p[last] = 0;
				break;
			//case 'm':	/* Gave arc minutes */
				//p[last] = 0;
				//scale = GMT_MIN2DEG;
				//break;
			case 'c':
				//if (GMT_compat_check (GMT, 4)) {	/* Warn and fall through */
					//GMT_Report (GMT->parent, GMT_MSG_COMPAT, "Warning: Second interval unit c is deprecated; use s instead\n");
				//}
				//else {
				//	scale = 1.0;
				//	break;
				//}
			case 's':	/* Gave arc seconds */
				//p[last] = 0;
				//scale = GMT_SEC2DEG;
				//break;
			case 'e':	/* Gave meters along mid latitude */
				p[last] = 0;
				if (i < 2) GMT->current.io.inc_code[i] |= GMT_INC_IS_M;
				break;
			case 'f':	/* Gave feet along mid latitude */
				p[last] = 0;
				if (i < 2) GMT->current.io.inc_code[i] |= GMT_INC_IS_FEET;
				break;
			case 'k':	/* Gave km along mid latitude */
				p[last] = 0;
				if (i < 2) GMT->current.io.inc_code[i] |= GMT_INC_IS_KM;
				break;
			case 'M':	/* Gave miles along mid latitude */
				p[last] = 0;
				if (i < 2) GMT->current.io.inc_code[i] |= GMT_INC_IS_MILES;
				break;
			case 'n':	/* Gave nautical miles along mid latitude */
				p[last] = 0;
				if (i < 2) GMT->current.io.inc_code[i] |= GMT_INC_IS_NMILES;
				break;
			case 'u':	/* Gave survey feet along mid latitude */
				//p[last] = 0;
				//if (i < 2) GMT->current.io.inc_code[i] |= GMT_INC_IS_SURVEY_FEET;
				//break;
			default:	/* No special flags or units */
				scale = 1.0;
				break;
		}
		if ((sscanf(p, "%lf", &inc[i])) != 1) {
			printf ( "Error: Unable to decode %s as a floating point number\n", p);
			GMT_exit (GMT, EXIT_FAILURE); return EXIT_FAILURE;
		}
		inc[i] *= scale;
		i++;	/* Goto next increment */
	}


	return (i);	/* Returns the number of increments found */
}

void GMT_check_lattice (struct GMT_CTRL *GMT, double *inc, unsigned int *registration, bool *active)
{	/* Uses provided settings to initialize the lattice settings from
	 * the -R<grdfile> if it was given; else it does nothing.
	 */
	if (!GMT->current.io.grd_info.active) return;	/* -R<grdfile> was not used; use existing settings */

	/* Here, -R<grdfile> was used and we will use the settings supplied by the grid file (unless overridden) */
	if (!active || *active == false) {	/* -I not set separately */
		GMT_memcpy (inc, GMT->current.io.grd_info.grd.inc, 2, double);
		inc[GMT_Y] = GMT->current.io.grd_info.grd.inc[GMT_Y];
	}
	if (registration) {	/* An pointer not NULL was passed that indicates grid registration */
		/* If a -r like option was set then toggle grid setting, else use grid setting */
		*registration = (*registration) ? !GMT->current.io.grd_info.grd.registration : GMT->current.io.grd_info.grd.registration;
	}
	if (active) *active = true;	/* When 4th arg is not NULL it is set to true (for Ctrl->active args) */
}

int GMT_parse_dash_option (struct GMT_CTRL *GMT, char *text)
{	/* parse any --PARAM[=value] arguments */
	int n;
	char *this_c = NULL, message[GMT_LEN128] = {""};
	if (!text)
		return (GMT_NOERROR);

	/* print version and exit */
	//if (strcmp (text, "version") == 0) {
		//sprintf (message, "%s\n", GMT_PACKAGE_VERSION_WITH_SVN_REVISION);
		//GMT->parent->print_func (stdout, message);
		/* cannot call GMT_Free_Options() from here, so we are leaking on exit.
		 * struct GMTAPI_CTRL *G = GMT->parent;
		 * if (GMT_Destroy_Session (G))
		 *   exit (EXIT_FAILURE); */
	//	exit (EXIT_SUCCESS);
	//}

	/* print GMT folders and exit */
	//if (strcmp (text, "show-datadir") == 0) {
	//	sprintf (message, "%s\n", GMT->session.SHAREDIR);
		//GMT->parent->print_func (stdout, message);
		/* leaking on exit same as above. */
		//exit (EXIT_SUCCESS);
	//}

	//if ((this_c = strchr (text, '='))) {
		/* Got --PAR=VALUE */
	//	this_c[0] = '\0';	/* Temporarily remove the '=' character */
	//	n = gmt_setparameter (GMT, text, &this_c[1]);
	//	this_c[0] = '=';	/* Put it back were it was */
	//}
	//else
		/* Got --PAR */
	//	n = gmt_setparameter (GMT, text, "true");
	//return (n);
}


bool GMT_getinc (struct GMT_CTRL *GMT, char *line, double inc[])
{	/* Special case of getincn use where n is two. */

	int n;

	/* Syntax: -I<xinc>[m|s|e|f|k|M|n|u|+|=][/<yinc>][m|s|e|f|k|M|n|u|+|=]
	 * Units: d = arc degrees
	 * 	  m = arc minutes
	 *	  s = arc seconds [was c]
	 *	  e = meter [Convert to degrees]
	 *	  f = feet [Convert to degrees]
	 *	  M = Miles [Convert to degrees]
	 *	  k = km [Convert to degrees]
	 *	  n = nautical miles [Convert to degrees]
	 *	  u = survey feet [Convert to degrees]
	 * Flags: = = Adjust -R to fit exact -I [Default modifies -I to fit -R]
	 *	  + = incs are actually nx/ny - convert to get xinc/yinc
	 */
	
	if (!line) { printf ("No argument given to GMT_getinc\n"); return (true); }

	n = GMT_getincn (GMT, line, inc, 2);
	
	if (n == 1) {	/* Must copy y info from x */
		inc[GMT_Y] = inc[GMT_X];
		GMT->current.io.inc_code[GMT_Y] = GMT->current.io.inc_code[GMT_X];	/* Use exact inc codes for both x and y */
	}
	
	if (GMT->current.io.inc_code[GMT_X] & GMT_INC_IS_NNODES && GMT->current.io.inc_code[GMT_X] & GMT_INC_UNITS) {
		printf ("Error: number of x nodes cannot have units\n");
		return (true);
	}
	
	if (GMT->current.io.inc_code[GMT_Y] & GMT_INC_IS_NNODES && GMT->current.io.inc_code[GMT_Y] & GMT_INC_UNITS) {
		printf ("Error: number of y nodes cannot have units\n");
		return (true);
	}
	
	return (false);
}

struct GMT_GRID * GMT_create_grid (struct GMT_CTRL *GMT)
{	/* Allocates space for a new grid container.  No space allocated for the float grid itself */
	struct GMT_GRID *G = NULL;

	G = GMT_memory (GMT, NULL, 1, struct GMT_GRID);
	G->header = GMT_memory (GMT, NULL, 1, struct GMT_GRID_HEADER);
	GMT_grd_init (GMT, G->header, NULL, false); /* Set default values */
	G->alloc_mode = GMT_ALLOCATED_BY_GMT;		/* Memory can be freed by GMT. */
	G->alloc_level = GMT->hidden.func_level;	/* Must be freed at this level. */
	G->id = GMT->parent->unique_var_ID++;		/* Give unique identifier */
	return (G);
}

/* There are many tools which requires grid x/y or cpt z to be in meters but the user may have these
 * data in km or miles.  Appending +u<unit> to the file addresses this conversion. */

char *GMT_file_unitscale (char *name)
{	/* Determine if this file ends in +u|U<unit>, with <unit> one of the valid Cartesian distance units */
	char *c = NULL;
	size_t len = strlen (name);					/* Get length of the file name */
	if (len < 4) return NULL;					/* Not enough space for name and modifier */
	c = &name[len-3];						/* c may be +u<unit>, +U<unit> or anything else */
	if (c[0] != '+') return NULL;					/* Does not start with + */
	if (! (c[1] == 'u' || c[1] == 'U')) return NULL;		/* Did not have the proper modifier u or U */
	if (strchr (GMT_LEN_UNITS2, c[2]) == NULL) return NULL;		/* Does no have a valid unit at the end */
	return c;							/* We passed, return c */
}

bool gmt_file_is_readable (struct GMT_CTRL *GMT, char *path)
{	/* Returns true if readable, otherwise give error and return false */
	if (!access (path, R_OK)) return (true);	/* Readable */
	/* Get here when found, but not readable */
	//GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "Unable to read %s (permissions?)\n", path);
	return (false);	/* Cannot read, give up */
}

bool gmt_traverse_dir (const char *file, char *path) {
	/* Look for file in the directory pointed to by path, recursively */
	DIR *D = NULL;
	struct dirent *F = NULL;
	int len, d_namlen;
	bool ok = false;
	char savedpath[GMT_BUFSIZ];

 	if ((D = opendir (path)) == NULL) return (false);	/* Unable to open directory listing */
	len = (int)strlen (file);
	strncpy (savedpath, path, GMT_BUFSIZ);	/* Make copy of current directory path */

	while (!ok && (F = readdir (D)) != NULL) {	/* For each directory entry until end or ok becomes true */
		d_namlen = (int)strlen (F->d_name);
		if (d_namlen == 1 && F->d_name[0] == '.') continue;				/* Skip current dir */
		if (d_namlen == 2 && F->d_name[0] == '.' && F->d_name[1] == '.') continue;	/* Skip parent dir */
#ifdef HAVE_SYS_DIR_H_
		if (F->d_type == DT_DIR) {	/* Entry is a directory; must search this directory recursively */
			sprintf (path, "%s/%s", savedpath, F->d_name);
			ok = gmt_traverse_dir (file, path);
		}
		else if (d_namlen == len && !strcmp (F->d_name, file)) {	/* Found the file in this dir (i.e., F_OK) */
			sprintf (path, "%s/%s", savedpath, file);
			ok = true;
		}
#endif /* HAVE_SYS_DIR_H_ */
	}
	(void)closedir (D);
	return (ok);	/* did or did not find file */
}

char *GMT_getdatapath (struct GMT_CTRL *GMT, const char *stem, char *path, int mode)
{
	/* stem is the name of the file, e.g., grid.img
	 * path is the full path to the file in question
	 * Returns full pathname if a workable path was found
	 * Looks for file stem in current directory and $GMT_{USER,DATA}DIR
	 * If the dir ends in / we traverse recursively [not under Windows].
	 */
	 
	unsigned int d, pos;
	size_t L;
	bool found;
	//GMT->session.USERDIR = strdup("/home/nishita/buildGMTUsingMake/");
	//GMT->session.DATADIR = strdup("data");
	
	char path1[GMT_BUFSIZ] ;//= "/home/nishita/buildGMTUsingMake/data";
		if (getcwd(path1, sizeof(path1)) == NULL)
		{
			printf("Error \n\n\n");
			return -1;		
		}

		/*********Added by Maria**********/
		FILE *f;
		char value[100],value2[100];
		strcat(path1,"/datapath.dat");
		f=fopen(path1,"r");
		fscanf(f,"%s",value);
		fclose(f);
		strcpy(value2,value);
		value2[strlen(value2)-6]='\0';	
		/*********************************/
		strcpy(path1,value2);			
		//strcpy(path1,"/home/labm/.core/TomoEK");	
	GMT->session.USERDIR = strdup(path1);
		strcpy(path1,value);		
		//strcpy(path1,"/home/labm/.core/TomoEK/data/");
	GMT->session.DATADIR = strdup(path1);
	//printf("file : %s line : %d func: %s  USERDIR :%s , DATADIR %s\n",__FILE__,__LINE__,__func__,GMT->session.USERDIR, GMT->session.DATADIR);
	char *udir[2] = {GMT->session.USERDIR, GMT->session.DATADIR}, dir[GMT_BUFSIZ];
	char path_separator[2] = {PATH_SEPARATOR, '\0'};

#ifdef HAVE_DIRENT_H_
	size_t N;
	bool gmt_traverse_dir (const char *file, char *path);
#endif /* HAVE_DIRENT_H_ */
	bool gmt_file_is_readable (struct GMT_CTRL *GMT, char *path);

	/* First look in the current working directory */
	if (!access (stem, F_OK)) {	/* Yes, found it */
		
		if (mode == F_OK || gmt_file_is_readable (GMT, (char *)stem)) {	/* Yes, found it or can read it */
			
			strcpy (path, stem);
			char path2[GMT_BUFSIZ] ;
			strcpy(path2, path);
			//printf("file : %s line : %d func: %s stem %s \n",__FILE__,__LINE__,__func__,stem);

			
				fflush(stdout);
			return (path2);
		}
		
		return (NULL);	/* Cannot read, give up */
	}
	
	/* If we got here and a full path is given, we give up ... unless it is one of those /vsi.../ files */
	if (stem[0] == '/') {
//#ifdef HAVE_GDAL
	//	if (GMT_check_url_name ((char *)stem))
		//	return ((char *)stem);			/* With GDAL all the /vsi-stuff is given existence credit */
//		else
//			return (NULL);
//#else
		return (NULL);
//#endif
	}

	/* Not found, see if there is a file in the GMT_{USER,DATA}DIR directories [if set] */

	for (d = 0; d < 2; d++) {	/* Loop over USER and DATA dirs */
		if (!udir[d]) continue;	/* This directory was not set */
		found = false;
		pos = 0;
		while (!found && (GMT_strtok (udir[d], path_separator, &pos, dir))) {
			L = strlen (dir);
#ifdef HAVE_DIRENT_H_
			if (dir[L-1] == '/' || (GMT_compat_check (GMT, 4) && dir[L-1] == '*')) {	/* Must search recursively from this dir */
				N = (dir[L-1] == '/') ? L - 1 : L - 2;
				strncpy (path, dir, N);	path[N] = 0;
				found = gmt_traverse_dir (stem, path);
			}
			else {
#endif /* HAVE_DIRENT_H_ */
				sprintf (path, "%s/%s", dir, stem);
				found = (!access (path, F_OK));
#ifdef HAVE_DIRENT_H_
			}
#endif /* HAVE_DIRENT_H_ */
		}
		if (found && gmt_file_is_readable (GMT, path)) return (path);	/* Yes, can read it */
	}

	return (NULL);	/* No file found, give up */
}

#include <unistd.h>
#include <fcntl.h>
int GMT_access (struct GMT_CTRL *GMT, const char* filename, int mode)
{	/* Like access but also checks the GMT_*DIR places */
	char file[GMT_BUFSIZ] = {""}, *c = NULL;
	
	if (GMT_File_Is_Memory (filename)) return (0);	/* Memory location always exists */
	//printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
	file[0] = '\0';		/* 'Initialize' it so we can test if it's still 'empty' after the sscanf below */
	if (!filename || !filename[0])
		return (-1);	/* No file given */
	sscanf (filename, "%[^=?]", file);	/* Exclude netcdf 3-D grid extensions to make sure we get a valid file name */
	if (file[0] == '\0')
		return (-1);		/* It happens for example when parsing grdmath args and it finds an isolated  "=" */
	//printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
	if ((c = GMT_file_unitscale (file))) c[0] = '\0';	/* Chop off any x/u unit specification */
	//printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
	if (mode == W_OK)
		return (access (file, mode));	/* When writing, only look in current directory */
	if (mode == R_OK || mode == F_OK) {	/* Look in special directories when reading or just checking for existance */
		
		char path[GMT_BUFSIZ];
		int fd;
		if(access (filename, mode) == -1)
			{
		//char path[GMT_BUFSIZ] ;//= "/home/nishita/buildGMTUsingMake/data";
		if (getcwd(path, sizeof(path)) == NULL)
		{
			printf("Error \n\n\n");
			return -1;		
		}
		/*********Added by Maria**********/
		FILE *f;
		char value[100];
		strcat(path,"/datapath.dat");
		f=fopen(path,"r");
		fscanf(f,"%s",value);
		fclose(f);
			
		/*********************************/
		strcpy(path,value);
		strcat(path,filename);
		fd = open(path, O_RDWR);
			}
		else
			{
		fd = open(filename, O_RDWR); //added by nishita
			}
    		if (fd == -1)
		{
			printf("file : %s not found !!  mode %d, file path %s \n",file,mode,path);
			return -1;
		}
		else
		{
			close(fd);
		}
		//return (GMT_getdatapath (GMT, file, path, mode) ? 0 : -1);
		return 0;
	}
	/* If we get here then mode is bad (X_OK)? */
	printf ( "GMT: Bad mode (%d) passed to GMT_access\n");
	return (-1);
}

int gmt_alloc_grid (struct GMT_CTRL *GMT, struct GMT_GRID *Grid)
{	/* Use information in Grid header to allocate the grid data.
	 * We assume gmt_init_grdheader has been called. */

	if (Grid->data) return (GMT_PTR_NOT_NULL);
	if (Grid->header->size == 0U) return (GMT_SIZE_IS_ZERO);
	if ((Grid->data = GMT_memory_aligned (GMT, NULL, Grid->header->size, float)) == NULL) return (GMT_MEMORY_ERROR);
	return (GMT_NOERROR);
}

bool GMT_check_url_name (char *fname) {
	/* File names starting as below should not be tested for existance or reading permissions as they
	   are either meant to be accessed on the fly (http & ftp) or they are compressed. So, if any of
	   the conditions holds true, returns true. All cases are read via GDAL support. */
	if ( !strncmp(fname,"http:",5)        || 
		!strncmp(fname,"https:",6)    || 
		!strncmp(fname,"ftp:",4)      || 
		!strncmp(fname,"/vsizip/",8)  || 
		!strncmp(fname,"/vsigzip/",9) || 
		!strncmp(fname,"/vsicurl/",9) ||
		!strncmp(fname,"/vsimem/",8)  || 
		!strncmp(fname,"/vsitar/",8) )

		return (true);
	else
		return (false);
}

int gmt_scanf_geo (char *s, double *val)
{
	/* Try to read a character string token stored in s, knowing that it should be a geographical variable.
	If successful, stores value in val and returns one of GMT_IS_FLOAT, GMT_IS_GEO, GMT_IS_LAT, GMT_IS_LON,
	whichever can be determined from the format of s.
	If unsuccessful, does not store anything in val and returns GMT_IS_NAN.
	This should have essentially the same functionality as the GMT3.4 GMT_scanf, except that the expectation
	is now used and returned, and this also permits a double precision format in the minutes or seconds,
	and does more error checking.  However, this is not optimized for speed (yet).  WHFS, 16 Aug 2001

	Note: Mismatch handling (e.g. this routine finds a lon but calling routine expected a lat) is not
	done here.
	*/

	int retval = GMT_IS_FLOAT, id, im;
	bool negate = false;
	unsigned int ncolons;
	size_t k;
	char scopy[GMT_LEN64] = {""}, suffix, *p = NULL, *p2 = NULL;
	double dd, dm, ds;

	k = strlen (s);
	if (k == 0) return (GMT_IS_NAN);
	if (!(isdigit ((int)s[k-1]))) {
		suffix = s[k-1];
		switch (suffix) {
			case 'W': case 'w':
				negate = true;
				retval = GMT_IS_LON;
				break;
			case 'E': case 'e':
				retval = GMT_IS_LON;
				break;
			case 'S': case 's':
				negate = true;
				retval = GMT_IS_LAT;
				break;
			case 'N': case 'n':
				retval = GMT_IS_LAT;
				break;
			case 'G': case 'g': case 'D': case 'd':
				retval = GMT_IS_GEO;
				break;
			case '.':	/* Decimal point without decimals, e.g., 123. */
				break;
			default:
				return (GMT_IS_NAN);
				break;
		}
		k--;
	}
	if (k >= GMT_LEN64) return (GMT_IS_NAN);
	strncpy (scopy, s, k);				/* Copy all but the suffix  */
	scopy[k] = 0;
	ncolons = 0;
	if ((p = strpbrk (scopy, "dD"))) {
		/* We found a D or d.  */
		if (strlen (p) == 1 || (strpbrk (&p[1], "dD:") ) ){
			/* It is at the end, or followed by a colon or another d or D.  */
			return (GMT_IS_NAN);
		}
		/* Map it to an e, permitting FORTRAN Double Precision formats.  */
		p[0] = 'e';
	}
	p = scopy;
	while ((p2 = strpbrk (p, ":"))) {
		if (strlen (p2) == 1) return (GMT_IS_NAN);	/* Shouldn't end with a colon  */
		ncolons++;
		if (ncolons > 2) return (GMT_IS_NAN);
		p = &p2[1];
	}

	if (ncolons && retval == GMT_IS_FLOAT) retval = GMT_IS_GEO;

	dd = 0.0;
	switch (ncolons) {
		case 0:
			if ((sscanf (scopy, "%lf", &dd)) != 1) return (GMT_IS_NAN);
			break;
		case 1:
			if ((sscanf (scopy, "%d:%lf", &id, &dm)) != 2) return (GMT_IS_NAN);
			dd = dm * GMT_MIN2DEG;
			if (id < 0)	/* Negative degrees present, subtract the fractional part */
				dd = id - dd;
			else if (id > 0)	/* Positive degrees present, add the fractional part */
				dd = id + dd;
			else {			/* degree part is 0; check if a leading sign is present */
				if (scopy[0] == '-') dd = -dd;	/* Make fraction negative */
			}
			break;
		case 2:
			if ((sscanf (scopy, "%d:%d:%lf", &id, &im, &ds)) != 3) return (GMT_IS_NAN);
			dd = im * GMT_MIN2DEG + ds * GMT_SEC2DEG;
			if (id < 0)	/* Negative degrees present, subtract the fractional part */
				dd = id - dd;
			else if (id > 0)	/* Positive degrees present, add the fractional part */
				dd = id + dd;
			else {			/* degree part is 0; check if a leading sign is present */
				if (scopy[0] == '-') dd = -dd;	/* Make fraction negative */
			}
			break;
	}
	*val = (negate) ? -dd : dd;
	return (retval);
}

int gmt_scanf_float (char *s, double *val)
{
	/* Try to decode a value from s and store
	in val.  s should not have any special format
	(neither geographical, with suffixes or
	separating colons, nor calendar nor clock).
	However, D and d are permitted to map to e
	if this would result in a success.  This
	allows Fortran Double Precision to be readable.

	On success, return GMT_IS_FLOAT and store val.
	On failure, return GMT_IS_NAN and do not touch val.
	*/

	char scopy[GMT_LEN64] = {""}, *p = NULL;
	double x;
	size_t j, k;

	x = strtod (s, &p);
	if (p[0] == 0) {	/* Success (non-Fortran).  */
		*val = x;
		return (GMT_IS_FLOAT);
	}
	if (p[0] != 'D' && p[0] != 'd') return (GMT_IS_NAN);
	k = strlen (p);
	if (k == 1) return (GMT_IS_NAN);	/* A string ending in e would be invalid  */
	/* Make a copy of s in scopy, mapping the d or D to an e */
	j = strlen (s);
	if (j > GMT_LEN64) return (GMT_IS_NAN);
	j -= k;
	strncpy (scopy, s, j);
	scopy[j] = 'e';
	strcpy (&scopy[j+1], &p[1]);
	x = strtod (scopy, &p);
	if (p[0] != 0) return (GMT_IS_NAN);
	*val = x;
	return (GMT_IS_FLOAT);
}
unsigned int GMT_unit_lookup (struct GMT_CTRL *GMT, int c, unsigned int unit)
{
	if (!isalpha ((int)c))	/* Not a unit modifier - just return the current default unit */
		return (unit);

	/* Now we check for the c-i-p units and barf otherwise */

	switch (c) {
		case 'c': case 'C':	/* Centimeters */
			unit = GMT_CM;
			break;
		case 'i': case 'I':	/* Inches */
			unit = GMT_INCH;
			break;
		case 'p': case 'P':	/* Points */
			unit = GMT_PT;
			break;
		default:
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: Length unit %c not supported - revert to default unit [%s]\n", (int)c, GMT->session.unit_name[unit]);
			break;
	}

	return (unit);
}

bool GMT_is_valid_number (char *t)
{
	int i, n;

	/* Checks if t fits the format [+|-][xxxx][.][yyyy][e|E[+|-]nn]. */

	if (!t) return (true);				/* Cannot be NULL */
	i = n = 0;
	if (t[i] == '+' || t[i] == '-') i++;		/* OK to have leading sign */
	while (isdigit ((int)t[i])) i++, n++;		/* OK to have numbers */
	if (t[i] == '.') {				/* Found a decimal */
		i++;	/* Go to next character */
		while (isdigit ((int)t[i])) i++, n++;	/* OK to have numbers following the decimal */
	}
	/* Here n must be > 0.  Also, we might find exponential notation */
	if (t[i] == 'e' || t[i] == 'E') {
		i++;
		if (t[i] == '+' || t[i] == '-') i++;	/* OK to have leading sign for exponent */
		while (isdigit ((int)t[i])) i++;	/* OK to have numbers for the exponent */
	}
	/* If all is well we should now have run out of characters in t and n > 0 - otherwise it is an error */
	return ((t[i] || n == 0) ? false : true);
}

double GMT_convert_units (struct GMT_CTRL *GMT, char *string, unsigned int default_unit, unsigned int target_unit)
{
	/* Converts the input string "value" to a float in units indicated by target_unit
	 * If value does not contain a unit (''c', 'i', or p') then the units indicated
	 * by default_unit will be used.
	 * Both target_unit and default_unit are either GMT_PT, GMT_CM, GMT_INCH or GMT_M.
	 */

	int c = 0, len, given_unit;
	bool have_unit = false;
	double value;

	if ((len = (int)strlen(string))) {
		c = string[len-1];
		if ((have_unit = isalpha ((int)c))) string[len-1] = '\0';	/* Temporarily remove unit */
	}

	/* So c is either 0 (meaning default unit) or any letter (even junk like z) */

	given_unit = GMT_unit_lookup (GMT, c, default_unit);	/* Will warn if c is not 0, 'c', 'i', 'p' */

	//if (!GMT_is_valid_number (string))
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: %s not a valid number and may not be decoded properly.\n", string);

	value = atof (string) * GMT->session.u2u[given_unit][target_unit];
	if (have_unit) string[len-1] = (char)GMT->session.unit_name[given_unit][0];	/* Put back the (implied) given unit */

	return (value);
}

int gmt_scanf_dim (struct GMT_CTRL *GMT, char *s, double *val)
{
	/* Try to decode a value from s and store
	in val.  s is a regular float with optional
	unit info, e.g., 8.5i or 7.5c.  If a valid unit
	is found we convert the number to inch.
	We also skip any trailing modifiers like +<mods>, e.g.
	vector specifications like 0.5i+jc+b+s

	We return GMT_IS_FLOAT and pass val.
	*/

	if (isalpha ((int)s[0]) || (s[1] == 0 && (s[0] == '-' || s[0] == '+')))	/* Probably a symbol character; quietly return 0 */
		*val = 0.0;
	else {	/* Probably a dimension with optional unit.  First check if there are modifiers to ignore here */
		char *p = NULL;
		if ((p = strchr (s, '+'))) { /* Found trailing +mod args */
			*p = 0;	/* Chop off modifier */
			*val = GMT_to_inch (GMT, s);	/* Get dimension */
			*p = '+';	/* Restore modifier */
		}
		else
			*val = GMT_to_inch (GMT, s);
	}
	return (GMT_IS_FLOAT);
}

int gmt_scanf_clock (struct GMT_CTRL *GMT, char *s, double *val)
{
	/* On failure, return -1.  On success, set val and return 0.

	Looks for apAP, but doesn't discover a failure if called with "11:13:15 Hello, Walter",
	because it will find an a.

	Doesn't check whether use of a or p matches stated intent to use twelve_hour_clock.

	ISO standard allows 24:00:00, so 86400 is not too big.
	If the day of this clock might be a day with a leap second, (this routine doesn't know that)
	then we should also allow 86401.  A value exceeding 86401 is an error.
	*/

	int k, hh, mm, add_noon = 0, hh_limit = 24;	/* ISO std allows 24:00:00  */
	double ss, x;
	char *p = NULL;

	if ( (p = strpbrk (s, "apAP") ) ) {
		switch (p[0]) {
			case 'a':
			case 'A':
				add_noon = 0;
				hh_limit = 12;
				break;
			case 'p':
			case 'P':
				add_noon = 43200;
				hh_limit = 12;
				break;
			default:
				return (-1);
				break;
		}
	}

	k = sscanf (s, GMT->current.io.clock_input.format, &hh, &mm, &ss);
	if (k == 0) return (-1);
	if (hh < 0 || hh > hh_limit) return (-1);

	x = (double)(add_noon + 3600*hh);
	if (k > 1) {
		if (mm < 0 || mm > 59) return (-1);
		x += 60*mm;
	}
	if (k > 2) {
		x += ss;
		if (x > 86401.0) return (-1);
	}
	*val = x;
	return (0);
}


int gmt_scanf_g_calendar (struct GMT_CTRL *GMT, char *s, int64_t *rd)
{
	/* Return -1 on failure.  Set rd and return 0 on success.

	For gregorian calendars.  */

	int i, k, ival[4];
	char month[16];

	if (GMT->current.io.date_input.day_of_year) {
		/* Calendar uses year and day of year format.  */
		if ( (k = sscanf (s, GMT->current.io.date_input.format,
			&ival[GMT->current.io.date_input.item_order[0]],
			&ival[GMT->current.io.date_input.item_order[1]]) ) == 0) return (-1);
		if (k < 2) {
			if (!GMT->current.io.date_input.truncated_cal_is_ok) return (-1);
			ival[1] = 1;	/* Set first day of year  */
		}
		if (GMT->current.io.date_input.Y2K_year) {
			if (ival[0] < 0 || ival[0] > 99) return (-1);
			ival[0] = GMT_y2_to_y4_yearfix (GMT, ival[0]);
		}
		k = (GMT_is_gleap (ival[0])) ? 366 : 365;
		if (ival[3] < 1 || ival[3] > k) return (-1);
		*rd = GMT_rd_from_gymd (GMT, ival[0], 1, 1) + ival[3] - 1;
		return (0);
	}

	/* Get here when calendar type has months and days of months.  */

	if (GMT->current.io.date_input.mw_text) {	/* Have month name abbreviation in data format */
		switch (GMT->current.io.date_input.item_pos[1]) {	/* Order of month in data string */
			case 0:	/* e.g., JAN-24-1987 or JAN-1987-24 */
				k = sscanf (s, GMT->current.io.date_input.format, month, &ival[GMT->current.io.date_input.item_order[1]], &ival[GMT->current.io.date_input.item_order[2]]);
				break;
			case 1:	/* e.g., 24-JAN-1987 or 1987-JAN-24 */
				k = sscanf (s, GMT->current.io.date_input.format, &ival[GMT->current.io.date_input.item_order[0]], month, &ival[GMT->current.io.date_input.item_order[2]]);
				break;
			case 2:	/* e.g., JAN-24-1987 ? */
				k = sscanf (s, GMT->current.io.date_input.format, month, &ival[GMT->current.io.date_input.item_order[1]], &ival[GMT->current.io.date_input.item_order[2]]);
				break;
			default:
				k = 0;
				return (-1);
				break;
		}
		GMT_str_toupper (month);
		for (i = ival[1] = 0; i < 12 && ival[1] == 0; i++) {
			if (!strcmp (month, GMT->current.time.language.month_name[3][i])) ival[1] = i + 1;
		}
		if (ival[1] == 0) return (-1);	/* No match for month name */
	}
	else if ((k = sscanf (s, GMT->current.io.date_input.format, &ival[GMT->current.io.date_input.item_order[0]], &ival[GMT->current.io.date_input.item_order[1]], &ival[GMT->current.io.date_input.item_order[2]])) == 0)
		return (-1);
	if (k < 3) {
		if (GMT->current.io.date_input.truncated_cal_is_ok) {
			ival[2] = 1;	/* Set first day of month  */
			if (k == 1) ival[1] = 1;	/* Set first month of year */
		}
		else
			return (-1);
	}
	if (GMT->current.io.date_input.Y2K_year) {
		if (ival[0] < 0 || ival[0] > 99) return (-1);
		ival[0] = GMT_y2_to_y4_yearfix (GMT, ival[0]);
	}

	if (GMT_g_ymd_is_bad (ival[0], ival[1], ival[2]) ) return (-1);

	*rd = GMT_rd_from_gymd (GMT, ival[0], ival[1], ival[2]);
	return (0);
}

int gmt_scanf_ISO_calendar (struct GMT_CTRL *GMT, char *s, int64_t *rd) {

	/* On failure, return -1.  On success, set rd and return 0.
	Assumes that year, week of year, day of week appear in that
	order only, and that the format string can handle the W.
	Assumes also that it is always OK to fill in missing bits.  */

	int k, n, ival[3];

	if ((n = sscanf (s, GMT->current.io.date_input.format, &ival[0], &ival[1], &ival[2])) == 0) return (-1);

	/* Handle possible missing bits */
	for (k = n; k < 3; k++) ival[k] = 1;

	if (ival[1] < 1 || ival[1] > 53) return (-1);
	if (ival[2] < 1 || ival[2] > 7) return (-1);
	if (GMT->current.io.date_input.Y2K_year) {
		if (ival[0] < 0 || ival[0] > 99) return (-1);
		ival[0] = GMT_y2_to_y4_yearfix (GMT, ival[0]);
	}
	*rd = GMT_rd_from_iywd (GMT, ival[0], ival[1], ival[2]);
	return (0);
}

int gmt_scanf_calendar (struct GMT_CTRL *GMT, char *s, int64_t *rd)
{
	/* On failure, return -1.  On success, set rd and return 0 */
	if (GMT->current.io.date_input.iso_calendar) return (gmt_scanf_ISO_calendar (GMT, s, rd));
	return (gmt_scanf_g_calendar (GMT, s, rd));
}


int GMT_scanf (struct GMT_CTRL *GMT, char *s, unsigned int expectation, double *val)
{
	/* Called with s pointing to a char string, expectation
	indicating what is known/required/expected about the
	format of the string.  Attempts to decode the string to
	find a double value.  Upon success, loads val and
	returns type found.  Upon failure, does not touch val,
	and returns GMT_IS_NAN.  Expectations permitted on call
	are
		GMT_IS_FLOAT	we expect an uncomplicated float.
	*/

	char calstring[GMT_LEN64] = {""}, clockstring[GMT_LEN64] = {""}, *p = NULL;
	double x;
	int64_t rd;
	size_t callen, clocklen;
	
	if (s[0] == 'T') {	/* Numbers cannot start with letters except for clocks, e.g., T07:0 */
		if (!isdigit ((int)s[1])) return (GMT_IS_NAN);	/* Clocks must have T followed by digit, e.g., T07:0 otherwise junk*/
	}
	else if (isalpha ((int)s[0])) return (GMT_IS_NAN);	/* Numbers cannot start with letters */

	if (expectation & GMT_IS_GEO) {
		/* True if either a lat or a lon is expected  */
		return (gmt_scanf_geo (s, val));
	}

	else if (expectation == GMT_IS_FLOAT) {
		/* True if no special format is expected or allowed  */
		return (gmt_scanf_float (s, val));
	}

	else if (expectation == GMT_IS_DIMENSION) {
		/* True if units might be appended, e.g. 8.4i  */
		return (gmt_scanf_dim (GMT, s, val));
	}

	else if (expectation == GMT_IS_RELTIME) {
		/* True if we expect to read a float with no special
		formatting (except for an optional trailing 't'), and then
		assume it is relative time in user's units since epoch.  */
		callen = strlen (s) - 1;
		if (s[callen] == 't') s[callen] = '\0';
		if ((gmt_scanf_float (s, val)) == GMT_IS_NAN) return (GMT_IS_NAN);
		return (GMT_IS_ABSTIME);
	}

	else if (expectation == GMT_IS_ABSTIME) {
		/* True when we expect to read calendar and/or
		clock strings in user-specified formats.  If both
		are present, they must be in the form
		<calendar_string>T<clock_string>.
		If only a calendar string is present, then either
		<calendar_string> or <calendar_string>T are valid.
		If only a clock string is present, then it must
		be preceded by a T:  T<clock_string>, and the time
		will be treated as if on day one of our calendar.  */

		callen = strlen (s);
		if (callen < 2) return (GMT_IS_NAN);	/* Maybe should be more than 2  */

		//if ((p = strchr ( s, (int)('T'))) == NULL) {	/* This was too naive, being tricked by data like 12-OCT-20 (no trailing T, so OCT was it) */
		if (s[0] == 'T') {	/* Got T<clock> presumably */
			strncpy (clockstring, &s[1], GMT_LEN64);
			clocklen = callen - 1;
			callen = 0;
		}
		else if (callen <= GMT->current.io.date_input.T_pos) {	/* There is no trailing T.  Put all of s in calstring.  */
			clocklen = 0;
			strncpy (calstring, s, GMT_LEN64);
		}
		else if (callen == (GMT->current.io.date_input.T_pos+1)) {	/* Just a trailing T but no clock */
			clocklen = 0;
			strncpy (calstring, s, GMT->current.io.date_input.T_pos);
			calstring[GMT->current.io.date_input.T_pos] = 0;
		}
		else {	/* Have something following the T */
			p = &s[GMT->current.io.date_input.T_pos+1];
			if (p[0] == 'T') ++p;	/* Probably negative year pushed everything off by one */
			clocklen = strlen (p);
			callen -= (clocklen + 1);	/* The one is for the 'T' */
			strncpy (calstring, s, callen);
			calstring[callen] = 0;
			strncpy (clockstring, p, GMT_LEN64);
			if (clocklen) clocklen--;
		}
		x = 0.0;	/* Default to 00:00:00 if no clock is given */
		if (clocklen && gmt_scanf_clock (GMT, clockstring, &x)) return (GMT_IS_NAN);
		rd = GMT->current.time.today_rata_die;	/* Default to today if no date is given */
		if (callen && gmt_scanf_calendar (GMT, calstring, &rd)) return (GMT_IS_NAN);
		*val = GMT_rdc2dt (GMT, rd, x);
		if (GMT->current.setting.time_is_interval) {	/* Must truncate and center on time interval */
			GMT_moment_interval (GMT, &GMT->current.time.truncate.T, *val, true);	/* Get the current interval */
			if (GMT->current.time.truncate.direction) {	/* Actually need midpoint of previous interval... */
				x = GMT->current.time.truncate.T.dt[0] - 0.5 * (GMT->current.time.truncate.T.dt[1] - GMT->current.time.truncate.T.dt[0]);
				GMT_moment_interval (GMT, &GMT->current.time.truncate.T, x, true);	/* Get the current interval */
			}
			/* Now get half-point of interval */
			*val = 0.5 * (GMT->current.time.truncate.T.dt[1] + GMT->current.time.truncate.T.dt[0]);
		}
		return (GMT_IS_ABSTIME);
	}
	
	else if (expectation == GMT_IS_ARGTIME) {
		printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
		//fflush(stdout);
		//return (gmt_scanf_argtime (GMT, s, val));
	}

	else if (expectation & GMT_IS_UNKNOWN) {
		/* True if we dont know but must try both geographic or float formats  */
		int type = gmt_scanf_geo (s, val);
		if ((type == GMT_IS_LON) && GMT->current.io.warn_geo_as_cartesion) {
			printf ("GMT: Longitude input data detected and successfully converted but will be considered Cartesian coordinates.\n");
			printf ( "GMT: If you need longitudes to be processed as periodic in 360 degrees then you must use -fg.\n");
			GMT->current.io.warn_geo_as_cartesion = false;	/* OK, done with the warning */
		}
		return (type);
	}

	else {
		printf("file : %s line : %d func: %s\n",__FILE__,__LINE__,__func__);
		printf ( "GMT_LOGIC_BUG: GMT_scanf() called with invalid expectation.\n");
		return (GMT_IS_NAN);
	}
}


/* Various functions to support {grd2xyz,xyz2grd}_func.c */

/* NOTE: In the following we check GMT->current.io.col_type[GMT_IN][2] and GMT->current.io.col_type[GMT_OUT][2] for formatting help for the first column.
 * We use column 3 ([2] or GMT_Z) instead of the first ([0]) since we really are dealing with the z in z (x,y) here
 * and the x,y are implicit from the -R -I arguments.
 */

int gmt_A_read (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{ /* Can read one or more items from input records. Limitation is
	 * that they must be floating point values (no dates or ddd:mm:ss) */
	uint64_t i;
	for (i = 0; i < n; ++i) {
		if (fscanf (fp, "%lg", &d[i]) <= 0)
			/* Read was unsuccessful */
			return (GMT_DATA_READ_ERROR);
	}
	return (GMT_OK);
}

int gmt_a_read (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{ /* Only reads one item regardless of *n */
	char line[GMT_LEN64] = {""}, *p;
	if (!fgets (line, GMT_LEN64, fp)) {
		/* Read was unsuccessful */
		GMT->current.io.status = GMT_IO_EOF;
		return (GMT_DATA_READ_ERROR);
	}
	/* Find end of string */
	p = line;
	while (*p)
		++p;
	/* Remove trailing whitespace */
	while ((--p != line) && strchr (" \t,\r\n", (int)*p));
	*(p + 1) = '\0';
	/* Convert whatever it is to double */
	GMT_scanf (GMT, line, GMT->current.io.col_type[GMT_IN][GMT_Z], d);
	return (GMT_OK);
}

int gmt_c_read (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* read int8_t aka char */
	uint64_t i;
	int8_t s;
	for (i = 0; i < n; ++i) {
		if (!GMT_fread (&s, sizeof (int8_t), 1U, fp)) {
			/* Read was unsuccessful */
			GMT->current.io.status = GMT_IO_EOF;
			return (GMT_DATA_READ_ERROR);
		}
		d[i] = (double) s;
	}
	return (GMT_OK);
}

int gmt_u_read (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* read uint8_t aka unsigned char */
	uint64_t i;
	uint8_t u;
	for (i = 0; i < n; ++i) {
		if (!GMT_fread (&u, sizeof (uint8_t), 1U, fp)) {
			GMT->current.io.status = GMT_IO_EOF;
			return (GMT_DATA_READ_ERROR);
		}
		d[i] = (double) u;
	}
	return (GMT_OK);
}

int gmt_h_read (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* read int16_t */
	uint64_t i;
	int16_t s;
	for (i = 0; i < n; ++i) {
		if (!GMT_fread (&s, sizeof (int16_t), 1U, fp)) {
			GMT->current.io.status = GMT_IO_EOF;
			return (GMT_DATA_READ_ERROR);
		}
		d[i] = (double) s;
	}
	return (GMT_OK);
}

int gmt_h_read_swab (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* read byteswapped int16_t */
	uint64_t i;
	uint16_t u;
	int16_t *s = (int16_t *)&u;
	for (i = 0; i < n; ++i) {
		if (!GMT_fread (&u, sizeof (uint16_t), 1U, fp)) {
			GMT->current.io.status = GMT_IO_EOF;
			return (GMT_DATA_READ_ERROR);
		}
		u = bswap16 (u);
		d[i] = (double) *s;
	}
	return (GMT_OK);
}

int gmt_H_read (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* read uint16_t */
	uint64_t i;
	uint16_t u;
	for (i = 0; i < n; ++i) {
		if (!GMT_fread (&u, sizeof (uint16_t), 1U, fp)) {
			GMT->current.io.status = GMT_IO_EOF;
			return (GMT_DATA_READ_ERROR);
		}
		d[i] = (double) u;
	}
	return (GMT_OK);
}

int gmt_H_read_swab (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* read byteswapped uint16_t */
	uint64_t i;
	uint16_t u;
	for (i = 0; i < n; ++i) {
		if (!GMT_fread (&u, sizeof (uint16_t), 1U, fp)) {
			GMT->current.io.status = GMT_IO_EOF;
			return (GMT_DATA_READ_ERROR);
		}
		d[i] = (double) bswap16 (u);
	}
	return (GMT_OK);
}

int gmt_i_read (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* read int32_t */
	uint64_t i;
	int32_t s;
	for (i = 0; i < n; ++i) {
		if (!GMT_fread (&s, sizeof (int32_t), 1U, fp)) {
			GMT->current.io.status = GMT_IO_EOF;
			return (GMT_DATA_READ_ERROR);
		}
		d[i] = (double) s;
	}
	return (GMT_OK);
}

int gmt_i_read_swab (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* read byteswapped int32_t */
	uint64_t i;
	uint32_t u;
	int32_t *s = (int32_t *)&u;
	for (i = 0; i < n; ++i) {
		if (!GMT_fread (&u, sizeof (uint32_t), 1U, fp)) {
			GMT->current.io.status = GMT_IO_EOF;
			return (GMT_DATA_READ_ERROR);
		}
		u = bswap32 (u);
		d[i] = (double) *s;
	}
	return (GMT_OK);
}

int gmt_I_read (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* read uint32_t */
	uint64_t i;
	uint32_t u;
	for (i = 0; i < n; ++i) {
		if (!GMT_fread (&u, sizeof (uint32_t), 1U, fp)) {
			GMT->current.io.status = GMT_IO_EOF;
			return (GMT_DATA_READ_ERROR);
		}
		d[i] = (double) u;
	}
	return (GMT_OK);
}

int gmt_I_read_swab (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* read byteswapped uint32_t */
	uint64_t i;
	uint32_t u;
	for (i = 0; i < n; ++i) {
		if (!GMT_fread (&u, sizeof (uint32_t), 1U, fp)) {
			GMT->current.io.status = GMT_IO_EOF;
			return (GMT_DATA_READ_ERROR);
		}
		d[i] = (double) bswap32 (u);
	}
	return (GMT_OK);
}

int gmt_l_read (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* read int64_t */
	uint64_t i;
	int64_t s;

	for (i = 0; i < n; ++i) {
		if (!GMT_fread (&s, sizeof (int64_t), 1U, fp)) {
			GMT->current.io.status = GMT_IO_EOF;
			return (GMT_DATA_READ_ERROR);
		}
		d[i] = (double) s;
	}
	return (GMT_OK);
}

int gmt_l_read_swab (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* read byteswapped int64_t */
	uint64_t i;
	uint64_t u;
	int64_t *s = (int64_t *)&u;
	for (i = 0; i < n; ++i) {
		if (!GMT_fread (&u, sizeof (uint64_t), 1U, fp)) {
			GMT->current.io.status = GMT_IO_EOF;
			return (GMT_DATA_READ_ERROR);
		}
		u = bswap64(u);
		d[i] = (double) *s;
	}
	return (GMT_OK);
}

int gmt_L_read (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* read uint64_t */
	uint64_t i;
	uint64_t u;

	for (i = 0; i < n; ++i) {
		if (!GMT_fread (&u, sizeof (uint64_t), 1U, fp)) {
			GMT->current.io.status = GMT_IO_EOF;
			return (GMT_DATA_READ_ERROR);
		}
		d[i] = (double) u;
	}
	return (GMT_OK);
}

int gmt_L_read_swab (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* read byteswapped uint64_t */
	uint64_t i;
	uint64_t u;

	for (i = 0; i < n; ++i) {
		if (!GMT_fread (&u, sizeof (uint64_t), 1U, fp)) {
			GMT->current.io.status = GMT_IO_EOF;
			return (GMT_DATA_READ_ERROR);
		}
		d[i] = (double) bswap64(u);
	}
	return (GMT_OK);
}

int gmt_f_read (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* read float */
	uint64_t i;
	float f;
	for (i = 0; i < n; ++i) {
		if (!GMT_fread (&f, sizeof (float), 1U, fp)) {
			GMT->current.io.status = GMT_IO_EOF;
			return (GMT_DATA_READ_ERROR);
		}
		d[i] = (double) f;
	}
	return (GMT_OK);
}

int gmt_f_read_swab (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* read byteswapped float */
	uint64_t i;
	union {
		float f;
		uint32_t bits;
	} u;
	for (i = 0; i < n; ++i) {
		if (!GMT_fread (&u.bits, sizeof (uint32_t), 1U, fp)) {
			GMT->current.io.status = GMT_IO_EOF;
			return (GMT_DATA_READ_ERROR);
		}
		u.bits = bswap32 (u.bits);
		d[i] = (double) u.f;
	}
	return (GMT_OK);
}

int gmt_d_read (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* read double */
	uint64_t i;
	for (i = 0; i < n; ++i) {
		if (!GMT_fread (&d[i], sizeof (double), 1U, fp)) {
			GMT->current.io.status = GMT_IO_EOF;
			return (GMT_DATA_READ_ERROR);
		}
	}
	return (GMT_OK);
}

int gmt_d_read_swab (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* read byteswapped double */
	uint64_t i;
	union {
		double d;
		uint64_t bits;
	} u;
	for (i = 0; i < n; ++i) {
		if (!GMT_fread (&u.bits, sizeof (uint64_t), 1U, fp)) {
			GMT->current.io.status = GMT_IO_EOF;
			return (GMT_DATA_READ_ERROR);
		}
		u.bits = bswap64 (u.bits);
		d[i] = u.d;
	}
	return (GMT_OK);
}

bool GMT_geo_to_dms (double val, int n_items, double fact, int *d, int *m,  int *s,  int *ix)
{
	/* Convert floating point degrees to dd:mm[:ss][.xxx].  Returns true if d = 0 and val is negative */
	bool minus;
	int isec, imin;
	double sec, fsec, min, fmin, step;

	minus = (val < 0.0);
	step = (fact == 0.0) ? GMT_CONV_LIMIT : 0.5 / fact;  	/* Precision desired in seconds (or minutes); else just deal with roundoff */

	if (n_items == 3) {		/* Want dd:mm:ss[.xxx] format */
		sec = GMT_DEG2SEC_F * fabs (val) + step;	/* Convert to seconds */
		isec = irint (floor (sec));			/* Integer seconds */
		fsec = sec - (double)isec;  			/* Leftover fractional second */
		*d = isec / GMT_DEG2SEC_I;			/* Integer degrees */
		isec -= ((*d) * GMT_DEG2SEC_I);			/* Left-over seconds in the last degree */
		*m = isec / GMT_MIN2SEC_I;			/* Integer minutes */
		isec -= ((*m) * GMT_MIN2SEC_I);			/* Leftover seconds in the last minute */
		*s = isec;					/* Integer seconds */
		*ix = irint (floor (fsec * fact));		/* Fractional seconds scaled to integer */
	}
	else if (n_items == 2) {		/* Want dd:mm[.xxx] format */
		min = GMT_DEG2MIN_F * fabs (val) + step;	/* Convert to minutes */
		imin = irint (floor (min));			/* Integer minutes */
		fmin = min - (double)imin;  			/* Leftover fractional minute */
		*d = imin / GMT_DEG2MIN_I;			/* Integer degrees */
		imin -= ((*d) * GMT_DEG2MIN_I);			/* Left-over seconds in the last degree */
		*m = imin;					/* Integer minutes */
		*s = 0;						/* No seconds */
		*ix = irint (floor (fmin * fact));		/* Fractional minutes scaled to integer */
	}
	else {		/* Want dd[.xxx] format */
		min = fabs (val) + step;			/* Convert to degrees */
		imin = irint (floor (min));			/* Integer degrees */
		fmin = min - (double)imin;  			/* Leftover fractional degree */
		*d = imin;					/* Integer degrees */
		*m = 0;						/* Integer minutes */
		*s = 0;						/* No seconds */
		*ix = irint (floor (fmin * fact));		/* Fractional degrees scaled to integer */
	}
	if (minus) {	/* OK, change sign, but watch for *d = 0 */
		if (*d)	/* Non-zero degree term is easy */
			*d = -(*d);
		else	/* Cannot change 0 to -0, so pass flag back to calling function */
			return (true);
	}
	return (false);
}

void gmt_format_geo_output (struct GMT_CTRL *GMT, bool is_lat, double geo, char *text)
{
	int k, n_items, d, m, s, m_sec, h_pos = 0;
	bool minus;
	char hemi[3] = {""}, *f = NULL;

	if (is_lat) {	/* Column is supposedly latitudes */
		if (fabs (geo) > 90.0) {
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Column selected for latitude-formatting has values that exceed +/- 90; set to NaN\n");
			sprintf (text, "NaN");
			return;
		}
	}
	else GMT_lon_range_adjust (GMT->current.io.geo.range, &geo);	/* Adjust longitudes */
	if (GMT->current.io.geo.decimal) {	/* Easy */
		f = (GMT->current.io.o_format[is_lat]) ? GMT->current.io.o_format[is_lat] : GMT->current.setting.format_float_out;
		sprintf (text, f, geo);
		return;
	}

	if (GMT->current.io.geo.wesn) {	/* Trailing WESN */
		if (GMT->current.io.geo.wesn == 2) hemi[h_pos++] = ' ';	/* Want space between numbers and hemisphere letter */
		if (is_lat)
			hemi[h_pos] = (GMT_IS_ZERO (geo)) ? 0 : ((geo < 0.0) ? 'S' : 'N');
		else
			hemi[h_pos] = (GMT_IS_ZERO (geo) || doubleAlmostEqual (geo, 180.0)) ? 0 : ((geo < 0.0) ? 'W' : 'E');
		geo = fabs (geo);
		if (hemi[h_pos] == 0) hemi[0] = 0;
	}

	for (k = n_items = 0; k < 3; k++) if (GMT->current.io.geo.order[k] >= 0) n_items++;	/* How many of d, m, and s are requested as integers */
	minus = GMT_geo_to_dms (geo, n_items, GMT->current.io.geo.f_sec_to_int, &d, &m, &s, &m_sec);	/* Break up into d, m, s, and remainder */
	if (minus) text[0] = '-';	/* Must manually insert leading minus sign when degree == 0 */
	if (GMT->current.io.geo.n_sec_decimals) {		/* Wanted fraction printed */
		if (n_items == 3)
			sprintf (&text[minus], GMT->current.io.geo.y_format, d, m, s, m_sec, hemi);
		else if (n_items == 2)
			sprintf (&text[minus], GMT->current.io.geo.y_format, d, m, m_sec, hemi);
		else
			sprintf (&text[minus], GMT->current.io.geo.y_format, d, m_sec, hemi);
	}
	else if (n_items == 3)
		sprintf (&text[minus], GMT->current.io.geo.y_format, d, m, s, hemi);
	else if (n_items == 2)
		sprintf (&text[minus], GMT->current.io.geo.y_format, d, m, hemi);
	else
		sprintf (&text[minus], GMT->current.io.geo.y_format, d, hemi);
}

void gmt_format_abstime_output (struct GMT_CTRL *GMT, double dt, char *text)
{
	char date[GMT_LEN16] = {""}, tclock[GMT_LEN16] = {""};

	GMT_format_calendar (GMT, date, tclock, &GMT->current.io.date_output, &GMT->current.io.clock_output, false, 1, dt);
	if (date[0] == '\0')	/* No date wanted hence dont use T */
		sprintf (text, "%s", tclock);
	else if (tclock[0] == '\0')	/* No clock wanted hence dont use T */
		sprintf (text, "%s", date);
	else	/* ISO format */
		sprintf (text, "%sT%s", date, tclock);
}

void GMT_ascii_format_col (struct GMT_CTRL *GMT, char *text, double x, unsigned int direction, uint64_t col)
{	/* Format based on column position in in or out direction */

	//printf("hi .................\n");
	if (GMT_is_dnan (x)) {	/* NaN, just write it as a string */
		sprintf (text, "NaN");
		return;
	}
	switch (GMT->current.io.col_type[direction][col]) {
		case GMT_IS_LON:
			gmt_format_geo_output (GMT, false, x, text);
			break;
		case GMT_IS_LAT:
			gmt_format_geo_output (GMT, true, x, text);
			break;
		case GMT_IS_ABSTIME:
			gmt_format_abstime_output (GMT, x, text);
			break;
		default:	/* Floating point */
			if (GMT->current.io.o_format[col])	/* Specific to this column */
				{
				sprintf (text, GMT->current.io.o_format[col], x);
				//printf("hi ..........%f",x);
				}
			else	/* Use the general float format */
				{
				//sprintf (text, GMT->current.setting.format_float_out, x);
				//printf("text .........X: .%lf",x);
				sprintf (text,"%.9f ", x);//nishita
				//sprintf (text, "%f ", x);
				
				}
			break;
	}
}

int GMT_ascii_output_col (struct GMT_CTRL *GMT, FILE *fp, double x, uint64_t col)
{	/* Formats x according to to output column number */
	char text[GMT_LEN256] = {""};

	//nishita start

	//end

	GMT_ascii_format_col (GMT, text, x, GMT_OUT, col);
	
	return (fprintf (fp, "%s", text));
}

int gmt_a_write (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* write ascii */
	uint64_t i;
	for (i = 0; i < (n - 1); ++i) {
		GMT_ascii_output_col (GMT, fp, d[i], GMT_Z);
		fprintf (fp, "\t");
	}
	/* last col */
	GMT_ascii_output_col (GMT, fp, d[i], GMT_Z);
	fprintf (fp, "\n");
	return (GMT_OK);
}

int gmt_c_write (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* write int8_t aka char */
	uint64_t i;
	int8_t s;
	for (i = 0; i < n; ++i) {
		s = (int8_t) d[i];
		if (GMT_fwrite (&s, sizeof (int8_t), 1U, fp) != 1U)
			return (GMT_DATA_WRITE_ERROR);
	}
	return (GMT_OK);
}

int gmt_u_write (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* write uint8_t aka unsigned char */
	uint64_t i;
	uint8_t u;
	for (i = 0; i < n; ++i) {
		u = (uint8_t) d[i];
		if (GMT_fwrite (&u, sizeof (uint8_t), 1U, fp) != 1U)
			return (GMT_DATA_WRITE_ERROR);
	}
	return (GMT_OK);
}

int gmt_h_write (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* write int16_t */
	uint64_t i;
	int16_t s;
	for (i = 0; i < n; ++i) {
		s = (int16_t) d[i];
		if (GMT_fwrite (&s, sizeof (int16_t), 1U, fp) != 1U)
			return (GMT_DATA_WRITE_ERROR);
	}
	return (GMT_OK);
}

int gmt_h_write_swab (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* write byteswapped int16_t */
	uint64_t i;
	uint16_t u;
	int16_t *s = (int16_t *)&u;
	for (i = 0; i < n; ++i) {
		*s = (int16_t) d[i];
		u = bswap16 (u);
		if (GMT_fwrite (&u, sizeof (uint16_t), 1U, fp) != 1U)
			return (GMT_DATA_WRITE_ERROR);
	}
	return (GMT_OK);
}

int gmt_H_write (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* write uint16_t */
	uint64_t i;
	uint16_t u;
	for (i = 0; i < n; ++i) {
		u = (uint16_t) d[i];
		if (GMT_fwrite (&u, sizeof (uint16_t), 1U, fp) != 1U)
			return (GMT_DATA_WRITE_ERROR);
	}
	return (GMT_OK);
}

int gmt_H_write_swab (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* write byteswapped uint16_t */
	uint64_t i;
	uint16_t u;
	for (i = 0; i < n; ++i) {
		u = bswap16 ((uint16_t) d[i]);
		if (GMT_fwrite (&u, sizeof (uint16_t), 1U, fp) != 1U)
			return (GMT_DATA_WRITE_ERROR);
	}
	return (GMT_OK);
}

int gmt_i_write (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* write int32_t */
	uint64_t i;
	int32_t s;
	for (i = 0; i < n; ++i) {
		s = (int32_t) d[i];
		if (GMT_fwrite (&s, sizeof (int32_t), 1U, fp) != 1U)
			return (GMT_DATA_WRITE_ERROR);
	}
	return (GMT_OK);
}

int gmt_i_write_swab (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* write byteswapped int32_t */
	uint64_t i;
	uint32_t u;
	int32_t *s = (int32_t *)&u;
	for (i = 0; i < n; ++i) {
		*s = (int32_t) d[i];
		u = bswap32 (u);
		if (GMT_fwrite (&u, sizeof (uint32_t), 1U, fp) != 1U)
			return (GMT_DATA_WRITE_ERROR);
	}
	return (GMT_OK);
}

int gmt_I_write (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* write uint32_t */
	uint64_t i;
	uint32_t u;
	for (i = 0; i < n; ++i) {
		u = (uint32_t) d[i];
		if (GMT_fwrite (&u, sizeof (uint32_t), 1U, fp) != 1U)
			return (GMT_DATA_WRITE_ERROR);
	}
	return (GMT_OK);
}

int gmt_I_write_swab (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* write byteswapped uint32_t */
	uint64_t i;
	uint32_t u;
	for (i = 0; i < n; ++i) {
		u = bswap32 ((uint32_t) d[i]);
		if (GMT_fwrite (&u, sizeof (uint32_t), 1U, fp) != 1U)
			return (GMT_DATA_WRITE_ERROR);
	}
	return (GMT_OK);
}

int gmt_l_write (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* write int64_t */
	uint64_t i;
	int64_t s;
	for (i = 0; i < n; ++i) {
		s = (int64_t) d[i];
		if (GMT_fwrite (&s, sizeof (int64_t), 1U, fp) != 1U)
			return (GMT_DATA_WRITE_ERROR);
	}
	return (GMT_OK);
}

int gmt_l_write_swab (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* write byteswapped int64_t */
	uint64_t i;
	uint64_t u;
	int64_t *s = (int64_t *)&u;
	for (i = 0; i < n; ++i) {
		*s = (int64_t) d[i];
		u = bswap64(u);
		if (GMT_fwrite (&u, sizeof (uint64_t), 1U, fp) != 1U)
			return (GMT_DATA_WRITE_ERROR);
	}
	return (GMT_OK);
}

int gmt_L_write (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* write uint64_t */
	uint64_t i;
	uint64_t u;
	for (i = 0; i < n; ++i) {
		u = (uint64_t) d[i];
		if (GMT_fwrite (&u, sizeof (int64_t), 1U, fp) != 1U)
			return (GMT_DATA_WRITE_ERROR);
	}
	return (GMT_OK);
}

int gmt_L_write_swab (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* write byteswapped uint64_t */
	uint64_t i;
	uint64_t u;
	for (i = 0; i < n; ++i) {
		u = bswap64((uint64_t) d[i]);
		if (GMT_fwrite (&u, sizeof (uint64_t), 1U, fp) != 1U)
			return (GMT_DATA_WRITE_ERROR);
	}
	return (GMT_OK);
}

int gmt_f_write (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* write float */
	uint64_t i;
	for (i = 0; i < n; ++i) {
		float f = (float) d[i];
		if (GMT_fwrite (&f, sizeof (float), 1U, fp) != 1U)
			return (GMT_DATA_WRITE_ERROR);
	}
	return (GMT_OK);
}

int gmt_f_write_swab (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* write byteswapped float */
	uint64_t i;
	union {
		float f;
		uint32_t bits;
	} u;
	for (i = 0; i < n; ++i) {
		u.f = (float) d[i];
		u.bits = bswap32(u.bits);
		if (GMT_fwrite (&u.bits, sizeof (uint32_t), 1U, fp) != 1U)
			return (GMT_DATA_WRITE_ERROR);
	}
	return (GMT_OK);
}

int gmt_d_write (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* write double */
	if (GMT_fwrite (d, sizeof (double), n, fp) != n)
		return (GMT_DATA_WRITE_ERROR);
	return (GMT_OK);
}

int gmt_d_write_swab (struct GMT_CTRL *GMT, FILE *fp, uint64_t n, double *d)
{
	/* write byteswapped double */
	uint64_t i;
	union {
		double d;
		uint64_t bits;
	} u;
	for (i = 0; i < n; ++i) {
		u.d = d[i];
		u.bits = bswap64 (u.bits);
		if (GMT_fwrite (&u.bits, sizeof (uint64_t), 1U, fp) != 1U)
			return (GMT_DATA_WRITE_ERROR);
	}
	return (GMT_OK);
}

p_to_io_func GMT_get_io_ptr (struct GMT_CTRL *GMT, int direction, enum GMT_swap_direction swap, char type)
{	/* Return pointer to read or write function for this data type */
	/* swap is 0 for no swap, 1 for swap input, 2 for swap output, 3 for swap both */
	p_to_io_func p = NULL;

	switch (type) {	/* Set read pointer depending on data format */
		case 'A':	/* ASCII with more than one per record */
			p = (direction == GMT_IN) ? &gmt_A_read : &gmt_a_write;
			break;
		case 'a':	/* ASCII */
			p = (direction == GMT_IN) ? &gmt_a_read : &gmt_a_write;
			break;
		case 'c':	/* Binary int8_t */
			p = (direction == GMT_IN) ? &gmt_c_read : &gmt_c_write;
			break;
		case 'u':	/* Binary uint8_t */
			p = (direction == GMT_IN) ? &gmt_u_read : &gmt_u_write;
			break;
		case 'h':	/* Binary int16_t */
			if (direction == GMT_IN)
				p = (swap & k_swap_in) ? &gmt_h_read_swab : &gmt_h_read;
			else
				p = (swap & k_swap_out) ? &gmt_h_write_swab : &gmt_h_write;
			break;
		case 'H':	/* Binary uint16_t */
			if (direction == GMT_IN)
				p = (swap & k_swap_in) ? &gmt_H_read_swab : &gmt_H_read;
			else
				p = (swap & k_swap_out) ? &gmt_H_write_swab : &gmt_H_write;
			break;
		case 'i':	/* Binary int32_t */
			if (direction == GMT_IN)
				p = (swap & k_swap_in) ? &gmt_i_read_swab : &gmt_i_read;
			else
				p = (swap & k_swap_out) ? &gmt_i_write_swab : &gmt_i_write;
			break;
		case 'I':	/* Binary uint32_t */
			if (direction == GMT_IN)
				p = (swap & k_swap_in) ? &gmt_I_read_swab : &gmt_I_read;
			else
				p = (swap & k_swap_out) ? &gmt_I_write_swab : &gmt_I_write;
			break;
		case 'l':	/* Binary int64_t */
			if (direction == GMT_IN)
				p = (swap & k_swap_in) ? &gmt_l_read_swab : &gmt_l_read;
			else
				p = (swap & k_swap_out) ? &gmt_l_write_swab : &gmt_l_write;
			break;
		case 'L':	/* Binary uint64_t */
			if (direction == GMT_IN)
				p = (swap & k_swap_in) ? &gmt_L_read_swab : &gmt_L_read;
			else
				p = (swap & k_swap_out) ? &gmt_L_write_swab : &gmt_L_write;
			break;
		case 'f':	/* Binary 4-byte float */
			if (direction == GMT_IN)
				p = (swap & k_swap_in) ? &gmt_f_read_swab : &gmt_f_read;
			else
				p = (swap & k_swap_out) ? &gmt_f_write_swab : &gmt_f_write;
			break;
		case 'd':	/* Binary 8-byte double */
			if (direction == GMT_IN)
				p = (swap & k_swap_in) ? &gmt_d_read_swab : &gmt_d_read;
			else
				p = (swap & k_swap_out) ? &gmt_d_write_swab : &gmt_d_write;
			break;
		case 'x':
			break;	/* Binary skip */

		default:
			printf ("%c not a valid data type!\n", type);
			GMT_exit (GMT, EXIT_FAILURE); return NULL;
			break;
	}

	return (p);
}

int GMT_get_io_type (struct GMT_CTRL *GMT, char type)
{
	int t = -1;
	switch (type) {
		/* Set read pointer depending on data format */
		case 'a': case 'A':          break; /* ASCII */
		case 'c': t = GMT_CHAR;   break; /* Binary int8_t */
		case 'u': t = GMT_UCHAR;  break; /* Binary uint8_t */
		case 'h': t = GMT_SHORT;  break; /* Binary int16_t */
		case 'H': t = GMT_USHORT; break; /* Binary uint16_t */
		case 'i': t = GMT_INT;    break; /* Binary int32_t */
		case 'I': t = GMT_UINT;   break; /* Binary uint32_t */
		case 'l': t = GMT_LONG;   break; /* Binary int64_t */
		case 'L': t = GMT_ULONG;  break; /* Binary uint64_t */
		case 'f': t = GMT_FLOAT;  break; /* Binary 4-byte float */
		case 'd': t = GMT_DOUBLE; break; /* Binary 8-byte double */
		default:
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "%c not a valid data type!\n", type);
			GMT_exit (GMT, EXIT_FAILURE); return EXIT_FAILURE;
			break;
	}
	return (t+1);	/* Since 0 means not set */
}

int GMT_set_cols (struct GMT_CTRL *GMT, unsigned int direction, uint64_t expected)
{	/* Initializes the internal GMT->common.b.ncol[] settings.
	 * direction is either GMT_IN or GMT_OUT.
	 * expected is the expected or known number of columns.  Use 0 if not known.
	 * For binary input or output the number of columns must be specified.
	 * For ascii output the number of columns must also be specified.
	 * For ascii input the i/o machinery will set this automatically so expected is ignored.
	 * Programs that need to read an input record in order to determine how
	 * many columns on output should call this function after returning the
	 * first data record; otherwise, call it before registering the resource.
	 */
	static char *mode[2] = {"input", "output"};

	if (! (direction == GMT_IN || direction == GMT_OUT)) return (GMT_NOT_A_VALID_DIRECTION);

	if (direction == GMT_IN && GMT->common.b.ncol[direction]) return (GMT_OK);	/* Already set once by -bi */

	if (expected == 0 && (direction == GMT_OUT || GMT->common.b.active[direction])) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Number of %s columns has not been set\n", mode[direction]);

		return (GMT_N_COLS_NOT_SET);
	}
	/* Here we may set the number of data columns */
	if (GMT->common.b.active[direction]) {	/* Must set uninitialized input/output pointers */
		uint64_t col;
		char type = (GMT->common.b.type[direction]) ? GMT->common.b.type[direction] : 'd';
		for (col = GMT->common.b.ncol[direction]; col < expected; col++) {
			GMT->current.io.fmt[direction][col].io = GMT_get_io_ptr (GMT, direction, GMT->common.b.swab[direction], type);
			GMT->current.io.fmt[direction][col].type = GMT_get_io_type (GMT, type);
		}
		GMT->common.b.ncol[direction] = expected;
	}
	else
		GMT->common.b.ncol[direction] = (direction == GMT_IN && expected == 0) ? GMT_MAX_COLUMNS : expected;
	if (direction == GMT_OUT && GMT->common.b.o_delay) {	/* Issue delayed message (see GMT_io_banner) */
		//GMT_io_banner (GMT, direction);
		GMT->common.b.o_delay = false;
	}
	if (direction == GMT_IN && GMT->common.i.active && GMT->common.i.n_cols > expected)
		printf ("Warning: Number of %s columns required [%d] is less that implied by -i [%d]\n", mode[GMT_IN], expected, GMT->common.i.n_cols);
	return (GMT_OK);
}

FILE *GMT_fopen (struct GMT_CTRL *GMT, const char *filename, const char *mode)
{
	char path[GMT_BUFSIZ];
	FILE *fd = NULL;
	
	if (mode[0] != 'r')	/* Open file for writing (so cannot be netCDF) */
		{
		
		//added by nishita

			char path1[GMT_BUFSIZ] ;//= "/home/nishita/buildGMTUsingMake/data";
		if (getcwd(path1, sizeof(path1)) == NULL)
		{
			printf("Error \n\n\n");
			return -1;		
		}
		
		/*********Added by Maria**********/
		FILE *f;
		char value[100];
		strcat(path1,"/datapath.dat");
		f=fopen(path1,"r");
		fscanf(f,"%s",value);
		fclose(f);
		/*********************************/
		strcpy(path1,value);			

	//	strcpy(path1,"/home/labm/.core/TomoEK/data/");
		strcat(path1,filename);
		//end
		//return (fopen (filename, mode)); //changed by nishita
		return (fopen (path1, mode)); 
		}
	else if (GMT->common.b.active[GMT_IN])	/* Definitely not netCDF */
		{
		
		return (fopen (GMT_getdatapath(GMT, filename, path, R_OK), mode));
		}
	//else if (GMT_compat_check (GMT, 4) && GMT->common.b.varnames[0])	/* Definitely netCDF */
		//return (gmt_nc_fopen (GMT, filename, mode));
	//else if (strchr (filename, '?'))	/* Definitely netCDF */
		//return (gmt_nc_fopen (GMT, filename, mode));
//#ifdef WIN32
//	else if (!strcmp (filename, "NUL"))	/* Special case of /dev/null under Windows */
//#else
	else if (!strcmp (filename, "/dev/null"))	/* The Unix null device; catch here to avoid gmt_nc_fopen */
//#endif
	{
		
		return (fopen (GMT_getdatapath(GMT, filename, path, R_OK), mode));
	}
	else {	/* Maybe netCDF */
		
		fd = gmt_nc_fopen (GMT, filename, mode);
		if (!fd) {
			char *c;
			if ((c = GMT_getdatapath(GMT, filename, path, R_OK)) != NULL) fd = fopen(c, mode);
		}
		return (fd);
	}
}

void GMT_io_binary_header (struct GMT_CTRL *GMT, FILE *fp, unsigned int dir)
{
	uint64_t k;
	char c = ' ';
	if (dir == GMT_IN) {	/* Use fread since we dont know if input is a stream or a file */
		for (k = 0; k < GMT->current.setting.io_n_header_items; k++) (void)GMT_fread (&c, sizeof (char), 1U, fp);
	}
	else {
		for (k = 0; k < GMT->current.setting.io_n_header_items; k++) GMT_fwrite (&c, sizeof (char), 1U, fp);
	}
}

int GMT_alloc_univector (struct GMT_CTRL *GMT, union GMT_UNIVECTOR *u, unsigned int type, uint64_t n_rows)
{
	/* Allocate space for one univector according to data type */
	int error = GMT_OK;
	switch (type) {
		case GMT_UCHAR:  u->uc1 = GMT_memory (GMT, u->uc1, n_rows, uint8_t);   if (u->uc1 == NULL) error = GMT_MEMORY_ERROR; break;
		case GMT_CHAR:   u->sc1 = GMT_memory (GMT, u->sc1, n_rows, int8_t);    if (u->sc1 == NULL) error = GMT_MEMORY_ERROR; break;
		case GMT_USHORT: u->ui2 = GMT_memory (GMT, u->ui2, n_rows, uint16_t);  if (u->ui2 == NULL) error = GMT_MEMORY_ERROR; break;
		case GMT_SHORT:  u->si2 = GMT_memory (GMT, u->si2, n_rows, int16_t);   if (u->si2 == NULL) error = GMT_MEMORY_ERROR; break;
		case GMT_UINT:   u->ui4 = GMT_memory (GMT, u->ui4, n_rows, uint32_t);  if (u->ui4 == NULL) error = GMT_MEMORY_ERROR; break;
		case GMT_INT:    u->si4 = GMT_memory (GMT, u->si4, n_rows, int32_t);   if (u->si4 == NULL) error = GMT_MEMORY_ERROR; break;
		case GMT_ULONG:  u->ui8 = GMT_memory (GMT, u->ui8, n_rows, uint64_t);  if (u->ui8 == NULL) error = GMT_MEMORY_ERROR; break;
		case GMT_LONG:   u->si8 = GMT_memory (GMT, u->si8, n_rows, int64_t);   if (u->si8 == NULL) error = GMT_MEMORY_ERROR; break;
		case GMT_FLOAT:  u->f4  = GMT_memory (GMT, u->f4,  n_rows, float);     if (u->f4  == NULL) error = GMT_MEMORY_ERROR; break;
		case GMT_DOUBLE: u->f8  = GMT_memory (GMT, u->f8,  n_rows, double);    if (u->f8  == NULL) error = GMT_MEMORY_ERROR; break;
	}
	return (error);
}

int gmt_alloc_vectors (struct GMT_CTRL *GMT, struct GMT_VECTOR *V)
{	/* Allocate space for each column according to data type */
	uint64_t col;
	int error;

	if (!V) return (GMT_PTR_IS_NULL);			/* Nothing to allocate to */
	if (V->n_columns == 0) return (GMT_PTR_IS_NULL);	/* No columns specified */
	if (V->n_rows == 0) return (GMT_N_COLS_NOT_SET);	/* No rows specified */
	if (V->data) return (GMT_PTR_IS_NULL);			/* Array of columns have not been allocated */
	for (col = 0; col < V->n_columns; col++) {
		if ((error = GMT_alloc_univector (GMT, &V->data[col], V->type[col],  V->n_rows)) != GMT_OK) return (error);
	}
	return (GMT_OK);
}

int GMT_fclose (struct GMT_CTRL *GMT, FILE *stream)
{
	if (!stream || stream == NULL)  return (0);
	/* First skip any stream related to the three Unix i/o descriptors */
	if (stream == GMT->session.std[GMT_IN])  return (0);
	if (stream == GMT->session.std[GMT_OUT]) return (0);
	if (stream == GMT->session.std[GMT_ERR]) return (0);
	if ((size_t)stream == (size_t)-GMT->current.io.ncid) {
		/* Special treatment for netCDF files */
		nc_close (GMT->current.io.ncid);
		GMT_free (GMT, GMT->current.io.varid);
		GMT_free (GMT, GMT->current.io.add_offset);
		GMT_free (GMT, GMT->current.io.scale_factor);
		GMT_free (GMT, GMT->current.io.missing_value);
		GMT->current.io.ncols = 0;
		GMT->current.io.ncid = GMT->current.io.nvars = 0;
		GMT->current.io.ndim = GMT->current.io.nrec = 0;
		GMT->current.io.input = GMT->session.input_ascii;
		return (0);
	}
	/* Regular file */
	return (fclose (stream));
}

		/* isspace, isalpha, ...: avoid assert (only happens with debug CRT) 
		   if passed a parameter that isn't EOF or in the range of 0 through 0xFF. */
#		ifdef _DEBUG
#			define isspace(c) (c > 0 && c < 0xFF && isspace(c))
#			define isalpha(c) (c > 0 && c < 0xFF && isalpha(c))
#		endif /* _DEBUG */

void GMT_strstrip(char *string, bool strip_leading) {
	/* Strip leading and trailing whitespace from string */
	char *start = string;
	char *end;

	assert (string != NULL); /* NULL pointer */

	if (strip_leading) {
		/* Skip over leading whitespace */
		while ((*start) && isspace(*start))
			++start;
		/* Is string just whitespace? */
		if (!(*start)) {
			*string = '\0'; /* Truncate entire string */
			return;
		}
	}

	/* Find end of string */
	end = start;
	while (*end)
		++end;

	/* Step backward until first non-whitespace */
	while ((--end != start) && isspace(*end));

	/* Chop off trailing whitespace */
	*(end + 1) = '\0';

	/* If leading whitespace, then move entire string back */
	if (string != start)
		memmove(string, start, end-start+2);
}

char *GMT_trim_segheader (struct GMT_CTRL *GMT, char *line) {
	/* Trim trailing junk and return pointer to first non-space/tab/> part of segment header
	 * Do not try to free the returned pointer!
	 */
	GMT_strstrip (line, false); /* Strip trailing whitespace */
	/* Skip over leading whitespace and segment marker */
	while (*line && (isspace(*line) || *line == GMT->current.setting.io_seg_marker[GMT_IN]))
		++line;
	/* Return header string */
	return (line);
}

void GMT_free_ogr (struct GMT_CTRL *GMT, struct GMT_OGR **G, unsigned int mode)
{	/* Free up GMT/OGR structure, if used */
	unsigned int k;
	if (!(*G)) return;	/* Nothing to do */
	/* mode = 0 only frees the aspatial data value array, while mode = 1 frees the entire struct and contents */
	for (k = 0; k < (*G)->n_aspatial; k++) {
		if (mode == 1 && (*G)->name && (*G)->name[k]) free ((*G)->name[k]);
		if ((*G)->tvalue && (*G)->tvalue[k]) free ((*G)->tvalue[k]);
	}
	if ((*G)->tvalue) GMT_free (GMT, (*G)->tvalue);
	if ((*G)->dvalue) GMT_free (GMT, (*G)->dvalue);
	if (mode == 0) return;	/* That's all we do for now */
	/* Here we free up everything */
	GMT_free (GMT, (*G)->name);
	GMT_free (GMT, (*G)->type);
	if ((*G)->region) free ((*G)->region);
	for (k = 0; k < 4; k++) free ((*G)->proj[k]);
	GMT_free (GMT, (*G));
}

void GMT_assign_segment (struct GMT_CTRL *GMT, struct GMT_DATASEGMENT *S, uint64_t n_rows, uint64_t n_columns)
{	/* Allocates and memcpy over vectors from GMT->hidden.mem_coord.
  	 * If n_rows > GMT_INITIAL_MEM_ROW_ALLOC then we pass the arrays and reset the tmp arrays to NULL
	 */
	uint64_t col;
	if (n_rows == 0) return;	/* Nothing to do */
	/* First allocate struct member arrays */
	S->coord = GMT_memory (GMT, NULL, n_columns, double *);
	S->min   = GMT_memory (GMT, NULL, n_columns, double);
	S->max   = GMT_memory (GMT, NULL, n_columns, double);

	if (n_rows > GMT_INITIAL_MEM_ROW_ALLOC) {	/* Large segment, just pass allocated pointers and start over with new tmp vectors later */
		//GMT_Report (GMT->parent, GMT_MSG_DEBUG, "GMT_assign_segment: Pass %" PRIu64 " large arrays with length = %" PRIu64 " off and get new tmp arrays\n", n_columns, n_rows);
		for (col = 0; col < n_columns; col++) {	/* Initialize the min/max array */
			if (n_rows < GMT->hidden.mem_rows) GMT->hidden.mem_coord[col] = GMT_memory (GMT, GMT->hidden.mem_coord[col], n_rows, double);	/* Trim back */
			S->coord[col] = GMT->hidden.mem_coord[col];	/* Pass the pointer */
			GMT->hidden.mem_coord[col] = NULL;		/* Null this out to start over for next segment */
		}
		GMT->hidden.mem_cols = 0;	/* Flag that we need to reallocate new temp arrays for next segment, if any */
	}
	else {	/* Small segments, allocate and memcpy, leave tmp array as is for further use */
		for (col = 0; col < n_columns; col++) {	/* Initialize the min/max array */
			//S->coord[col] = GMT_memory (GMT, NULL, n_rows, double);
			S->coord[col] = GMT_memory (GMT, S->coord[col], n_rows, double);
			GMT_memcpy (S->coord[col], GMT->hidden.mem_coord[col], n_rows, double);
		}
	}
	S->n_rows = n_rows;
	S->n_columns = n_columns;
}

void GMT_quad_reset (struct GMT_CTRL *GMT, struct GMT_QUAD *Q, uint64_t n_items)
{	/* Allocate and initialize the QUAD struct needed to find min/max of a set of longitudes */
	uint64_t i;

	GMT_memset (Q, n_items, struct GMT_QUAD);	/* Set all to NULL/0 */
	for (i = 0; i < n_items; i++) {
		Q[i].min[0] = Q[i].min[1] = +DBL_MAX;
		Q[i].max[0] = Q[i].max[1] = -DBL_MAX;
		Q[i].range[0] = GMT_IS_M180_TO_P180_RANGE;
		Q[i].range[1] = GMT_IS_0_TO_P360_RANGE;
	}
}

struct GMT_QUAD * GMT_quad_init (struct GMT_CTRL *GMT, uint64_t n_items)
{	/* Allocate an initialize the QUAD struct needed to find min/max of longitudes */
	struct GMT_QUAD *Q = GMT_memory (GMT, NULL, n_items, struct GMT_QUAD);

	GMT_quad_reset (GMT, Q, n_items);

	return (Q);
}

void GMT_lon_range_adjust (unsigned int range, double *lon)
{
	switch (range) {	/* Adjust to the desired range */
		case GMT_IS_0_TO_P360_RANGE:		/* Make 0 <= lon <= 360 */
			while ((*lon) < 0.0) (*lon) += 360.0;
			while ((*lon) > 360.0) (*lon) -= 360.0;
			break;
		case GMT_IS_0_TO_P360:		/* Make 0 <= lon < 360 */
			while ((*lon) < 0.0) (*lon) += 360.0;
			while ((*lon) >= 360.0) (*lon) -= 360.0;
			break;
		case GMT_IS_M360_TO_0_RANGE:		/* Make -360 <= lon <= 0 */
			while ((*lon) < -360.0) (*lon) += 360.0;
			while ((*lon) > 0) (*lon) -= 360.0;
			break;
		case GMT_IS_M360_TO_0:		/* Make -360 < lon <= 0 */
			while ((*lon) <= -360.0) (*lon) += 360.0;
			while ((*lon) > 0) (*lon) -= 360.0;
			break;
		case GMT_IS_M180_TO_P180_RANGE:	/* Make -180 <= lon <= +180 */
			while ((*lon) < -180.0) (*lon) += 360.0;
			while ((*lon) > 180.0) (*lon) -= 360.0;
			break;
		case GMT_IS_M180_TO_P180:	/* Make -180 <= lon < +180 [Special case where +180 is not desired] */
			while ((*lon) < -180.0) (*lon) += 360.0;
			while ((*lon) >= 180.0) (*lon) -= 360.0;
			break;
		case GMT_IS_M180_TO_P270_RANGE:	/* Make -180 <= lon < +270 [Special case for GSHHG only] */
			while ((*lon) < -180.0) (*lon) += 360.0;
			while ((*lon) >= 270.0) (*lon) -= 360.0;
			break;
		default:	/* Do nothing */
			break;
	}
}

void GMT_quad_add (struct GMT_CTRL *GMT, struct GMT_QUAD *Q, double x)
{	/* Update quad array for this longitude x */
	unsigned int way, quad_no;
	if (GMT_is_dnan (x)) return;	/* Cannot handle a NaN */
	for (way = 0; way < 2; way++) {
		GMT_lon_range_adjust (Q->range[way], &x);	/* Set -180/180, then 0-360 range */
		Q->min[way] = MIN (x, Q->min[way]);
		Q->max[way] = MAX (x, Q->max[way]);
	}
	quad_no = urint (floor (x / 90.0));	/* Now x is 0-360; this yields quadrants 0-3 */
	if (quad_no == 4) quad_no = 0;		/* When x == 360.0 */
	Q->quad[quad_no] = true;		/* Our x fell in this quadrant */
}

unsigned int GMT_quad_finalize (struct GMT_CTRL *GMT, struct GMT_QUAD *Q)
{
	/* Finalize longitude range settings */
	uint64_t n_quad;
	unsigned int way;

	n_quad = Q->quad[0] + Q->quad[1] + Q->quad[2] + Q->quad[3];		/* How many quadrants had data */
	if (Q->quad[0] && Q->quad[3])		/* Longitudes on either side of Greenwich only, must use -180/+180 notation */
		way = 0;
	else if (Q->quad[1] && Q->quad[2])	/* Longitudes on either side of the date line, must user 0/360 notation */
		way = 1;
	else if (n_quad == 2 && ((Q->quad[0] && Q->quad[2]) || (Q->quad[1] && Q->quad[3])))	/* Funny quadrant gap, pick shortest longitude extent */
		way = ((Q->max[0] - Q->min[0]) < (Q->max[1] - Q->min[1])) ? 0 : 1;
	else					/* Either will do, use default settings */
		way = (GMT->current.io.geo.range == GMT_IS_0_TO_P360_RANGE) ? 1 : 0;
	/* Final adjustments */
	if (Q->min[way] > Q->max[way]) Q->min[way] -= 360.0;
	if (Q->min[way] < 0.0 && Q->max[way] < 0.0) Q->min[way] += 360.0, Q->max[way] += 360.0;
	return (way);
}

void GMT_get_lon_minmax (struct GMT_CTRL *GMT, double *lon, uint64_t n_rows, double *min, double *max)
{	/* Return the min/max longitude in array lon using clever quadrant checking. */
	unsigned int way;
	uint64_t row;
	struct GMT_QUAD *Q = GMT_quad_init (GMT, 1);	/* Allocate and initialize one QUAD structure */

	/* We must keep separate min/max for both Dateline and Greenwich conventions */
	for (row = 0; row < n_rows; row++) GMT_quad_add (GMT, Q, lon[row]);

	/* Finalize longitude range settings */
	way = GMT_quad_finalize (GMT, Q);
	*min = Q->min[way];		*max = Q->max[way];
	GMT_free (GMT, Q);
}

void GMT_set_seg_minmax (struct GMT_CTRL *GMT, struct GMT_DATASEGMENT *S)
{	/* Determine the min/max values for each column in the segment */
	uint64_t row, col;

	for (col = 0; col < S->n_columns; col++) {
		if (GMT->current.io.col_type[GMT_IN][col] == GMT_IS_LON) /* Requires separate quandrant assessment */
			GMT_get_lon_minmax (GMT, S->coord[col], S->n_rows, &(S->min[col]), &(S->max[col]));
		else {	/* Simple Cartesian-like arrangement */
			S->min[col] = S->max[col] = S->coord[col][0];
			for (row = 1; row < S->n_rows; row++) {
				if (S->coord[col][row] < S->min[col]) S->min[col] = S->coord[col][row];
				if (S->coord[col][row] > S->max[col]) S->max[col] = S->coord[col][row];
			}
		}
	}
}

void GMT_set_tbl_minmax (struct GMT_CTRL *GMT, struct GMT_DATATABLE *T)
{	/* Update the min/max of all segments and the entire table */
	uint64_t seg, col;
	struct GMT_DATASEGMENT *S = NULL;

	if (!T) return;	/* No table given */
	if (!T->n_columns) return;	/* No columns given */
	if (!T->min) T->min = GMT_memory (GMT, NULL, T->n_columns, double);
	if (!T->max) T->max = GMT_memory (GMT, NULL, T->n_columns, double);
	for (col = 0; col < T->n_columns; col++) {	/* Initialize */
		T->min[col] = DBL_MAX;
		T->max[col] = -DBL_MAX;
	}
	for (seg = 0; seg < T->n_segments; seg++) {
		S = T->segment[seg];
		GMT_set_seg_minmax (GMT, S);
		for (col = 0; col < T->n_columns; col++) {
			if (S->min[col] < T->min[col]) T->min[col] = S->min[col];
			if (S->max[col] > T->max[col]) T->max[col] = S->max[col];
		}
	}
}

void GMT_set_dataset_minmax (struct GMT_CTRL *GMT, struct GMT_DATASET *D)
{
	uint64_t tbl, col;
	struct GMT_DATATABLE *T = NULL;
	if (!D) return;	/* No dataset given */
	if (!D->n_columns) return;	/* No columns given */
	if (!D->min) D->min = GMT_memory (GMT, NULL, D->n_columns, double);
	if (!D->max) D->max = GMT_memory (GMT, NULL, D->n_columns, double);
	for (col = 0; col < D->n_columns; col++) {	/* Initialize */
		D->min[col] = DBL_MAX;
		D->max[col] = -DBL_MAX;
	}
	for (tbl = 0; tbl < D->n_tables; tbl++) {
		T = D->table[tbl];
		for (col = 0; col < D->n_columns; col++) {
			if (T->min[col] < D->min[col]) D->min[col] = T->min[col];
			if (T->max[col] > D->max[col]) D->max[col] = T->max[col];
		}
	}
	
}

int GMT_alloc_segment (struct GMT_CTRL *GMT, struct GMT_DATASEGMENT *S, uint64_t n_rows, uint64_t n_columns, bool first)
{	/* (re)allocates memory for a segment of given dimensions.
 	 * If n_rows is 0 then we do not set S->n_rows.  */
	uint64_t col;
	if (first && n_columns) {	/* First time we allocate the number of columns needed */
		S->coord = GMT_memory (GMT, NULL, n_columns, double *);
		S->min = GMT_memory (GMT, NULL, n_columns, double);
		S->max = GMT_memory (GMT, NULL, n_columns, double);
		S->n_columns = n_columns;
		for (col = 0; col < n_columns; col++) {	/* Initialize the min/max array */
			S->min[col] = +DBL_MAX;
			S->max[col] = -DBL_MAX;
		}
	}
	if (n_rows) S->n_rows = n_rows;
	S->n_alloc = n_rows;
	if (n_rows) for (col = 0; col < n_columns; col++) S->coord[col] = GMT_memory (GMT, S->coord[col], n_rows, double);
	return (GMT_OK);
}

struct GMT_DATATABLE * gmt_alloc_table (struct GMT_CTRL *GMT, struct GMT_DATATABLE *Tin, uint64_t n_columns, uint64_t n_rows)
{
	/* Allocate the new Table structure with same # of segments and rows/segment as input table.
	 * However, n_columns is given separately and could differ.
	 * If n_rows is > 0 we well override the Tin rows counts by using n_rows instead.  */
	unsigned int hdr;
	uint64_t seg, nr;
	struct GMT_DATATABLE *T = GMT_memory (GMT, NULL, 1, struct GMT_DATATABLE);

	T->n_segments = T->n_alloc = Tin->n_segments;	/* Same number of segments as input table */
	T->n_headers  = Tin->n_headers;
	T->n_columns  = n_columns;		/* Separately specified n_columns */
	T->min = GMT_memory (GMT, NULL, n_columns, double);
	T->max = GMT_memory (GMT, NULL, n_columns, double);
	if (T->n_headers) {
		T->header = GMT_memory (GMT, NULL, Tin->n_headers, char *);
		for (hdr = 0; hdr < T->n_headers; hdr++) T->header[hdr] = strdup (Tin->header[hdr]);
	}
	T->segment = GMT_memory (GMT, NULL, Tin->n_segments, struct GMT_DATASEGMENT *);
	for (seg = 0; seg < T->n_segments; seg++) {
		T->segment[seg] = GMT_memory (GMT, NULL, 1, struct GMT_DATASEGMENT);
		nr = (n_rows) ? n_rows : Tin->segment[seg]->n_rows;
		GMT_alloc_segment (GMT, T->segment[seg], nr, n_columns, true);
		T->segment[seg]->n_rows = nr;
		T->segment[seg]->n_columns = n_columns;
		T->n_records += nr;
		if (Tin->segment[seg]->header) T->segment[seg]->header = strdup (Tin->segment[seg]->header);
		if (Tin->segment[seg]->label) T->segment[seg]->label = strdup (Tin->segment[seg]->label);
	}
	return (T);
}

struct GMT_DATASET * GMT_alloc_dataset (struct GMT_CTRL *GMT, struct GMT_DATASET *Din, uint64_t n_rows, uint64_t n_columns, unsigned int mode)
{
	/* Allocate new dataset structure with same # of tables, segments and rows/segment as input data set.
	 * However, n_columns is given separately and could differ.  Also, if n_rows > 0 we let that override the segment row counts.
	 * We copy over headers and segment headers.
	 * mode controls how the new dataset is to be allocated;
	 * mode = GMT_ALLOC_NORMAL means we replicate the number of tables and the layout of the Din dataset
	 * mode = GMT_ALLOC_VERTICAL means we concatenate all the tables in Din into a single table for Dout
	 * mode = GMT_ALLOC_HORIZONTAL means we base the Dout size only on the first Din table
	 *	(# of segments, # of rows/segment) because tables will be pasted horizontally and not vertically.
	 */
	unsigned int hdr;
	size_t len;
	uint64_t nr, tbl, seg, n_seg, seg_in_tbl;
	struct GMT_DATASET *D = GMT_memory (GMT, NULL, 1, struct GMT_DATASET);

	D->n_columns = (n_columns) ? n_columns : Din->n_columns;
	D->geometry = Din->geometry;
	D->min = GMT_memory (GMT, NULL, D->n_columns, double);
	D->max = GMT_memory (GMT, NULL, D->n_columns, double);
	if (mode) {	/* Pack everything into a single table */
		D->n_alloc = D->n_tables = 1;
		if (mode == GMT_ALLOC_VERTICAL)
			for (tbl = n_seg = 0; tbl < Din->n_tables; tbl++) n_seg += Din->table[tbl]->n_segments;
		else	/* mode == GMT_ALLOC_HORIZONTAL */
			n_seg = Din->table[0]->n_segments;
		D->table = GMT_memory (GMT, NULL, 1, struct GMT_DATATABLE *);
		D->table[0] = GMT_memory (GMT, NULL, 1, struct GMT_DATATABLE);

		/* As for file headers we concatenate the headers from all tables */
		D->table[0]->n_headers  = Din->table[0]->n_headers;
		if (D->table[0]->n_headers) D->table[0]->header = GMT_memory (GMT, NULL, D->table[0]->n_headers, char *);
		for (hdr = 0; hdr < D->table[0]->n_headers; hdr++) {	/* Concatenate headers */
			for (tbl = len = 0; tbl < Din->n_tables; tbl++) len += (strlen (Din->table[tbl]->header[hdr]) + 2);
			D->table[0]->header[hdr] = calloc (len, sizeof (char));
			strncpy (D->table[0]->header[hdr], Din->table[0]->header[hdr], len);
			if (Din->n_tables > 1) GMT_chop (D->table[0]->header[hdr]);	/* Remove newline */
			for (tbl = 1; tbl < Din->n_tables; tbl++) {	/* Now go across tables to paste */
				if (tbl < (Din->n_tables - 1)) GMT_chop (Din->table[tbl]->header[hdr]);
				strcat (D->table[0]->header[hdr], "\t");
				strcat (D->table[0]->header[hdr], Din->table[tbl]->header[hdr]);
			}
		}

		D->n_segments = D->table[0]->n_segments = D->table[0]->n_alloc = n_seg;
		D->table[0]->n_columns = D->n_columns;
		D->table[0]->segment = GMT_memory (GMT, NULL, n_seg, struct GMT_DATASEGMENT *);
		D->table[0]->min = GMT_memory (GMT, NULL, D->n_columns, double);
		D->table[0]->max = GMT_memory (GMT, NULL, D->n_columns, double);
		for (seg = tbl = seg_in_tbl = 0; seg < D->n_segments; seg++) {
			if (seg == Din->table[tbl]->n_segments) { tbl++; seg_in_tbl = 0; }	/* Go to next table */
			D->table[0]->segment[seg] = GMT_memory (GMT, NULL, 1, struct GMT_DATASEGMENT);
			nr = (n_rows) ? n_rows : Din->table[tbl]->segment[seg_in_tbl]->n_rows;
			D->table[0]->segment[seg]->n_rows = nr;
			GMT_alloc_segment (GMT, D->table[0]->segment[seg], nr, D->n_columns, true);
			D->table[0]->segment[seg]->n_columns = D->n_columns;
			if (mode != GMT_ALLOC_HORIZONTAL && Din->table[tbl]->segment[seg_in_tbl]->header) D->table[0]->segment[seg]->header = strdup (Din->table[tbl]->segment[seg_in_tbl]->header);
			D->n_records += nr;
			seg_in_tbl++;
		}
	}
	else {	/* Just copy over the same dataset layout except for columns */
		D->n_alloc  = D->n_tables = Din->n_tables;		/* Same number of tables as input dataset */
		D->n_segments  = Din->n_segments;	/* Same number of segments as input dataset */
		D->n_records  = Din->n_records;		/* Same number of records as input dataset */
		D->table = GMT_memory (GMT, NULL, D->n_tables, struct GMT_DATATABLE *);
		for (tbl = 0; tbl < D->n_tables; tbl++) {
			D->table[tbl] = gmt_alloc_table (GMT, Din->table[tbl], D->n_columns, n_rows);
		}
	}
	D->alloc_level = GMT->hidden.func_level;	/* Must be freed at this level. */
	D->alloc_mode = GMT_ALLOCATED_BY_GMT;		/* So GMT_* modules can free this memory. */
	D->id = GMT->parent->unique_var_ID++;		/* Give unique identifier */
	return (D);
}

void GMT_copy_segment (struct GMT_CTRL *GMT, struct GMT_DATASEGMENT *Sout, struct GMT_DATASEGMENT *Sin)
{	/* Duplicates the segment */
	uint64_t col;
	for (col = 0; col < Sin->n_columns; col++) GMT_memcpy (Sout->coord[col], Sin->coord[col], Sin->n_rows, double);
	GMT_memcpy (Sout->min, Sin->min, Sin->n_columns, double);
	GMT_memcpy (Sout->max, Sin->max, Sin->n_columns, double);
	Sout->n_rows = Sin->n_rows;
}

struct GMT_DATASET * GMT_duplicate_dataset (struct GMT_CTRL *GMT, struct GMT_DATASET *Din, unsigned int mode, unsigned int *geometry)
{	/* Make an exact replica, return geometry if not NULL */
	uint64_t tbl, seg;
	struct GMT_DATASET *D = NULL;
	D = GMT_alloc_dataset (GMT, Din, 0, Din->n_columns, mode);
	GMT_memcpy (D->min, Din->min, Din->n_columns, double);
	GMT_memcpy (D->max, Din->max, Din->n_columns, double);
	for (tbl = 0; tbl < Din->n_tables; tbl++) {
		for (seg = 0; seg < Din->table[tbl]->n_segments; seg++) {
			GMT_copy_segment (GMT, D->table[tbl]->segment[seg], Din->table[tbl]->segment[seg]);
		}
		GMT_memcpy (D->table[tbl]->min, Din->table[tbl]->min, Din->table[tbl]->n_columns, double);
		GMT_memcpy (D->table[tbl]->max, Din->table[tbl]->max, Din->table[tbl]->n_columns, double);
	}
	if (geometry) *geometry = D->geometry;
	return (D);
}

struct GMT_TEXTTABLE * gmt_alloc_texttable (struct GMT_CTRL *GMT, struct GMT_TEXTTABLE *Tin)
{
	/* Allocate the new Text Table structure with same # of segments and rows/segment as input table. */
	uint64_t seg;
	unsigned int hdr;
	struct GMT_TEXTTABLE *T = GMT_memory (GMT, NULL, 1, struct GMT_TEXTTABLE);

	T->n_segments = T->n_alloc = Tin->n_segments;	/* Same number of segments as input table */
	T->n_records  = Tin->n_records;		/* Same number of records as input table */
	T->n_headers  = Tin->n_headers;
	if (T->n_headers) {
		T->header = GMT_memory (GMT, NULL, Tin->n_headers, char *);
		for (hdr = 0; hdr < T->n_headers; hdr++) T->header[hdr] = strdup (Tin->header[hdr]);
	}
	T->segment = GMT_memory (GMT, NULL, Tin->n_segments, struct GMT_TEXTSEGMENT *);
	for (seg = 0; seg < T->n_segments; seg++) {
		T->segment[seg] = GMT_memory (GMT, NULL, 1, struct GMT_TEXTSEGMENT);
		T->segment[seg]->record = GMT_memory (GMT, NULL, Tin->segment[seg]->n_rows, char *);
		T->segment[seg]->n_rows = T->segment[seg]->n_alloc = Tin->segment[seg]->n_rows;
		if (Tin->segment[seg]->header) T->segment[seg]->header = strdup (Tin->segment[seg]->header);
	}
	return (T);
}

struct GMT_TEXTSET * GMT_alloc_textset (struct GMT_CTRL *GMT, struct GMT_TEXTSET *Din, unsigned int mode)
{
	/* Allocate new textset structure with same # of tables, segments and rows/segment as input data set.
	 * We copy over headers and segment headers.
	 * mode controls how the new dataset is to be allocated;
	 * mode = GMT_ALLOC_NORMAL means we replicate the number of tables and the layout of the Din dataset
	 * mode = GMT_ALLOC_VERTICAL means we concatenate all the tables in Din into a single table for Dout
	 * mode = GMT_ALLOC_HORIZONTAL means we base the Dout size only on the first Din table
	 *	(# of segments, # of rows/segment) because tables will be pasted horizontally and not vertically.
	 */
	unsigned int hdr;
	uint64_t tbl, seg, n_seg, seg_in_tbl;
	size_t len;
	struct GMT_TEXTSET *D = GMT_memory (GMT, NULL, 1, struct GMT_TEXTSET);

	if (mode) {	/* Pack everything into a single table */
		D->n_alloc = D->n_tables = 1;
		if (mode == GMT_ALLOC_VERTICAL)
			for (n_seg = tbl = 0; tbl < Din->n_tables; tbl++) n_seg += Din->table[tbl]->n_segments;
		else /* mode == GMT_ALLOC_HORIZONTAL */
			n_seg = Din->table[0]->n_segments;
		D->table = GMT_memory (GMT, NULL, 1, struct GMT_TEXTTABLE *);
		D->table[0] = GMT_memory (GMT, NULL, 1, struct GMT_TEXTTABLE);

		/* As for file headers we concatenate the headers from all tables */
		D->table[0]->n_headers  = Din->table[0]->n_headers;
		if (D->table[0]->n_headers) D->table[0]->header = GMT_memory (GMT, NULL, D->table[0]->n_headers, char *);
		for (hdr = 0; hdr < D->table[0]->n_headers; hdr++) {	/* Concatenate headers */
			for (len = tbl = 0; tbl < Din->n_tables; tbl++) len += (strlen (Din->table[tbl]->header[hdr]) + 2);
			D->table[0]->header[hdr] = calloc (len, sizeof (char));
			strncpy (D->table[0]->header[hdr], Din->table[0]->header[hdr], GMT_BUFSIZ);
			if (Din->n_tables > 1) GMT_chop (D->table[0]->header[hdr]);	/* Remove newline */
			for (tbl = 1; tbl < Din->n_tables; tbl++) {	/* Now go across tables to paste */
				if (tbl < (Din->n_tables - 1)) GMT_chop (Din->table[tbl]->header[hdr]);
				strcat (D->table[0]->header[hdr], "\t");
				strcat (D->table[0]->header[hdr], Din->table[tbl]->header[hdr]);
			}
		}

		D->n_segments = D->table[0]->n_segments = D->table[0]->n_alloc = n_seg;
		D->table[0]->segment = GMT_memory (GMT, NULL, n_seg, struct GMT_TEXTSEGMENT *);
		for (seg = tbl = seg_in_tbl = 0; seg < D->n_segments; seg++) {
			if (seg == Din->table[tbl]->n_segments) { tbl++; seg_in_tbl = 0; }	/* Go to next table */
			D->table[0]->segment[seg] = GMT_memory (GMT, NULL, 1, struct GMT_TEXTSEGMENT);
			D->table[0]->segment[seg]->n_rows = Din->table[tbl]->segment[seg_in_tbl]->n_rows;
			D->table[0]->segment[seg]->record = GMT_memory (GMT, NULL, D->table[0]->segment[seg]->n_rows, char *);
			if (mode == GMT_ALLOC_VERTICAL && Din->table[tbl]->segment[seg_in_tbl]->header) D->table[0]->segment[seg]->header = strdup (Din->table[tbl]->segment[seg_in_tbl]->header);
			seg_in_tbl++;
		}
	}
	else {	/* Just copy over the same dataset layout except for columns */
		D->n_alloc = D->n_tables = Din->n_tables;		/* Same number of tables as input dataset */
		D->n_segments  = Din->n_segments;	/* Same number of segments as input dataset */
		D->n_records  = Din->n_records;		/* Same number of records as input dataset */
		D->table = GMT_memory (GMT, NULL, D->n_tables, struct GMT_TEXTTABLE *);
		for (tbl = 0; tbl < D->n_tables; tbl++) D->table[tbl] = gmt_alloc_texttable (GMT, Din->table[tbl]);
	}
	D->geometry = Din->geometry;

	return (D);
}

struct GMT_TEXTSET * GMT_duplicate_textset (struct GMT_CTRL *GMT, struct GMT_TEXTSET *Din, unsigned int mode)
{
	uint64_t tbl, row, seg;
	struct GMT_TEXTSET *D = NULL;
	D = GMT_alloc_textset (GMT, Din, mode);
	for (tbl = 0; tbl < Din->n_tables; tbl++) for (seg = 0; seg < Din->table[tbl]->n_segments; seg++) {
		for (row = 0; row < Din->table[tbl]->segment[seg]->n_rows; row++) D->table[tbl]->segment[seg]->record[row] = strdup (Din->table[tbl]->segment[seg]->record[row]);
	}
	return (D);
}

void gmt_free_ogr_seg (struct GMT_CTRL *GMT, struct GMT_DATASEGMENT *S)
{	/* Frees the OGR structure for a given segment */
	unsigned int k, n;
	n = (GMT->current.io.OGR) ? GMT->current.io.OGR->n_aspatial : GMT->common.a.n_aspatial;
	if (n) {
		for (k = 0; S->ogr->tvalue && k < n; k++) if (S->ogr->tvalue[k]) free (S->ogr->tvalue[k]);
		GMT_free (GMT, S->ogr->tvalue);
		GMT_free (GMT, S->ogr->dvalue);
	}
	GMT_free (GMT, S->ogr);
}

void GMT_free_segment (struct GMT_CTRL *GMT, struct GMT_DATASEGMENT **S, enum GMT_enum_alloc alloc_mode)
{
	/* Free memory allocated by GMT_read_table */

	unsigned int k;
	uint64_t col;
	struct GMT_DATASEGMENT *segment = *S;
	if (!segment) return;	/* Do not try to free NULL pointer */
	if (alloc_mode == GMT_ALLOCATED_BY_GMT) {	/* Free data GMT allocated */
		for (col = 0; col < segment->n_columns; col++) if (segment->coord[col]) GMT_free (GMT, segment->coord[col]);
	}
	if (segment->coord) GMT_free (GMT, segment->coord);
	if (segment->min) GMT_free (GMT, segment->min);
	if (segment->max) GMT_free (GMT, segment->max);
	if (segment->label) free ( segment->label);
	if (segment->header) free ( segment->header);
	for (k = 0; k < 2; k++) if (segment->file[k]) free (segment->file[k]);
	if (segment->ogr) gmt_free_ogr_seg (GMT, segment);	/* OGR metadata */
	GMT_free (GMT, segment);
	*S = NULL;
}

void GMT_free_table (struct GMT_CTRL *GMT, struct GMT_DATATABLE *table, enum GMT_enum_alloc alloc_mode)
{
	unsigned int k;
	if (!table) return;		/* Do not try to free NULL pointer */
	for (k = 0; k < table->n_headers; k++) free (table->header[k]);
	if (table->n_headers) GMT_free (GMT, table->header);
	if (table->min) GMT_free (GMT, table->min);
	if (table->max) GMT_free (GMT, table->max);
	for (k = 0; k < 2; k++) if (table->file[k]) free (table->file[k]);
	GMT_free_ogr (GMT, &(table->ogr), 1);
	if (table->segment) {	/* Free segments */
		uint64_t seg;
		for (seg = 0; seg < table->n_segments; seg++) GMT_free_segment (GMT, &(table->segment[seg]), alloc_mode);
		GMT_free (GMT, table->segment);
	}
	GMT_free (GMT, table);
}
void gmt_grd_parse_xy_units (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *h, char *file, unsigned int direction)
{	/* Decode the optional +u|U<unit> and determine scales */
	enum GMT_enum_units u_number;
	unsigned int mode = 0;
	char *c = NULL, *name = NULL;
	
	if (GMT_is_geographic (GMT, direction)) return;	/* Does not apply to geographic data */
	name = (file) ? file : h->name;
	if ((c = GMT_file_unitscale (name)) == NULL) return;	/* Did not find any modifier */
	mode = (c[1] == 'u') ? 0 : 1;
	u_number = GMT_get_unit_number (GMT, c[2]);		/* Convert char unit to enumeration constant for this unit */
	if (u_number == GMT_IS_NOUNIT) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Grid file x/y unit specification %s was unrecognized (part of file name?) and is ignored.\n", c);
		return;
	}
	/* Got a valid unit */
	h->xy_unit_to_meter[direction] = GMT->current.proj.m_per_unit[u_number];	/* Converts unit to meters */
	if (mode) h->xy_unit_to_meter[direction] = 1.0 / h->xy_unit_to_meter[direction];	/* Wanted the inverse */
	h->xy_unit[direction] = u_number;	/* Unit ID */
	h->xy_adjust[direction] |= 1;		/* Says we have successfully parsed and readied the x/y scaling */
	h->xy_mode[direction] = mode;
	c[0] = '\0';	/* Chop off the unit specification from the file name */
}


void gmt_expand_filename (struct GMT_CTRL *GMT, char *file, char *fname)
{
	bool found;
	unsigned int i;
	size_t f_length, length;

	if (GMT->current.setting.io_gridfile_shorthand) {	/* Look for matches */
		f_length = strlen (file);
		for (i = 0, found = false; !found && i < GMT->session.n_shorthands; ++i) {
			length = strlen (GMT->session.shorthand[i].suffix);
			found = (length > f_length) ? false : !strncmp (&file[f_length - length], GMT->session.shorthand[i].suffix, length);
		}
		if (found) {	/* file ended in a recognized shorthand extension */
			--i;
			sprintf (fname, "%s=%s", file, GMT->session.shorthand[i].format);
		}
		else
			strcpy (fname, file);
	}
	else	/* Simply copy the full name */
		strcpy (fname, file);
}

int GMT_grd_format_decoder (struct GMT_CTRL *GMT, const char *code, unsigned int *type_id) {
	/* Returns the integer grid format ID that goes with the specified 2-character code */
	if (isdigit ((int)code[0])) {
		/* File format number given, look for old code */
		unsigned id_candidate = (unsigned) abs (atoi (code));
		if (id_candidate > 0 && id_candidate < GMT_N_GRD_FORMATS) {
			*type_id = id_candidate;
			return GMT_NOERROR;
		}
	}
	else {
		/* Character code given */
		unsigned i;
		for (i = 1; i < GMT_N_GRD_FORMATS; i++) {
			
			if (strncmp (GMT->session.grdformat[i], code, 2) == 0) {
				*type_id = i;
				return GMT_NOERROR;
			}
		}
	}
	return GMT_GRDIO_UNKNOWN_ID;
}

int parse_grd_format_scale (struct GMT_CTRL *Ctrl, struct GMT_GRID_HEADER *header, char *format) {
	/* parses format string after =-suffix: ff/scale/offset/invalid
	 * ff:      can be one of [abcegnrs][bsifd]
	 * scale:   can be any non-zero normalized number or 'a' for scale and
	 *          offset auto-adjust, defaults to 1.0 if omitted
	 * offset:  can be any finite number or 'a' for offset auto-adjust, defaults to 0 if omitted
	 * invalid: can be any finite number, defaults to NaN if omitted
	 * scale and offset may be left empty (e.g., ns//a will auto-adjust the offset only)
	 */

	char type_code[3];
	char *p;
	int err; /* GMT_err_trap */

	/* decode grid type */
	strncpy (type_code, format, 2);
	type_code[2] = '\0';
	err = GMT_grd_format_decoder (Ctrl, type_code, &header->type); /* update header type id */
	if (err != GMT_NOERROR)
		return err;

	/* parse scale/offset/invalid if any */
	p = strchr (format, '/');
	if (p != NULL && *p) {
		++p;
		/* parse scale */
		if (*p == 'a')
			header->z_scale_autoadust = header->z_offset_autoadust = true;
		else
			sscanf (p, "%lf", &header->z_scale_factor);
	}
	else
		return GMT_NOERROR;

	p = strchr (p, '/');
	if (p != NULL && *p) {
		++p;
		/* parse offset */
		if (*p != 'a') {
			header->z_offset_autoadust = false;
			sscanf (p, "%lf", &header->z_add_offset);
		}
		else
			header->z_offset_autoadust = true;
	}
	else
		return GMT_NOERROR;

	p = strchr (p, '/');
	if (p != NULL && *p) {
		++p;
		/* parse invalid value */
		sscanf (p, "%f", &header->nan_value);

		/* header->nan_value should be of same type as (float)*grid to avoid
		 * round-off errors. For example, =gd///-3.4028234e+38:gtiff, would fail
		 * otherwise because the GTiff'd NoData values are of type double but the
		 * grid is truncated to float.
		 * Don't allow infitiy: */
		if (!isfinite (header->nan_value))
			header->nan_value = (float)NAN;
	}

	return GMT_NOERROR;
}

int gmt_scanf_epoch (struct GMT_CTRL *GMT, char *s, int64_t *rata_die, double *t0) {

	/* Read a string which must be in one of these forms:
		[-]yyyy-mm-dd[T| [hh:mm:ss.sss]]
		[-]yyyy-Www-d[T| [hh:mm:ss.sss]]
	   Hence, data and clock can be separated by 'T' or ' ' (space), and the clock string is optional.
	   In fact, seconds can be decimal or integer, or missing. Minutes and hour are optional too.
	   Examples: 2000-01-01, 2000-01-01T, 2000-01-01 00:00, 2000-01-01T00, 2000-01-01T00:00:00.000
	*/

	double ss = 0.0;
	int i, yy, mo, dd, hh = 0, mm = 0;
	int64_t rd;
	char tt[8];

	i = 0;
	while (s[i] && s[i] == ' ') i++;
	if (!(s[i])) return (-1);
	if (strchr (&s[i], 'W') ) {	/* ISO calendar string, date with or without clock */
		if (sscanf (&s[i], "%5d-W%2d-%1d%[^0-9:-]%2d:%2d:%lf", &yy, &mo, &dd, tt, &hh, &mm, &ss) < 3) return (-1);
		if (GMT_iso_ywd_is_bad (yy, mo, dd) ) return (-1);
		rd = GMT_rd_from_iywd (GMT, yy, mo, dd);
	}
	else {				/* Gregorian calendar string, date with or without clock */
		if (sscanf (&s[i], "%5d-%2d-%2d%[^0-9:-]%2d:%2d:%lf", &yy, &mo, &dd, tt, &hh, &mm, &ss) < 3) return (-1);
		if (GMT_g_ymd_is_bad (yy, mo, dd) ) return (-1);
		rd = GMT_rd_from_gymd (GMT, yy, mo, dd);
	}
	if (GMT_hms_is_bad (hh, mm, ss)) return (-1);

	*rata_die = rd;								/* Rata day number of epoch */
	*t0 =  (GMT_HR2SEC_F * hh + GMT_MIN2SEC_F * mm + ss) * GMT_SEC2DAY;	/* Fractional day (0<= t0 < 1) since rata_die of epoch */
	return (GMT_NOERROR);
}

int GMT_init_time_system_structure (struct GMT_CTRL *GMT, struct GMT_TIME_SYSTEM *time_system) {
	/* Processes strings time_system.unit and time_system.epoch to produce a time system scale
	   (units in seconds), inverse scale, and rata die number and fraction of the epoch (days).
	   Return values: 0 = no error, 1 = unit error, 2 = epoch error, 3 = unit and epoch error.
	*/
	int error = GMT_NOERROR;

	/* Check the unit sanity */
	switch (time_system->unit) {
		case 'y':
		case 'Y':
			/* This is a kludge: we assume all years are the same length, thinking that a user
			with decimal years doesn't care about precise time.  To do this right would
			take an entirely different scheme, not a simple unit conversion. */
			time_system->scale = GMT_YR2SEC_F;
			break;
		case 'o':
		case 'O':
			/* This is also a kludge: we assume all months are the same length, thinking that a user
			with decimal years doesn't care about precise time.  To do this right would
			take an entirely different scheme, not a simple unit conversion. */
			time_system->scale = GMT_MON2SEC_F;
			break;
		case 'd':
		case 'D':
			time_system->scale = GMT_DAY2SEC_F;
			break;
		case 'h':
		case 'H':
			time_system->scale = GMT_HR2SEC_F;
			break;
		case 'm':
		case 'M':
			time_system->scale = GMT_MIN2SEC_F;
			break;
		case 's':
		case 'S':
			time_system->scale = 1.0;
			break;
		case 'c':
		case 'C':
			if (GMT_compat_check (GMT, 4)) {
				//GMT_Report (GMT->parent, GMT_MSG_COMPAT, "Warning: Unit c (seconds) is deprecated; use s instead.\n");
				time_system->scale = 1.0;
			}
			else
				error ++;
			break;
		default:
			error ++;
			break;
	}

	/* Set inverse scale and store it to avoid divisions later */
	time_system->i_scale = 1.0 / time_system->scale;

	/* Now convert epoch into rata die number and fraction */
	if (gmt_scanf_epoch (GMT, time_system->epoch, &time_system->rata_die, &time_system->epoch_t0)) error += 2;

	if (error & 1) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: TIME_UNIT is invalid.  Default assumed.\n");
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Choose one only from y o d h m s\n");
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Corresponding to year month day hour minute second\n");
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Note year and month are simply defined (365.2425 days and 1/12 of a year)\n");
	}
	if (error & 2) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: TIME_EPOCH format is invalid.  Default assumed.\n");
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "    A correct format has the form [-]yyyy-mm-ddThh:mm:ss[.xxx]\n");
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "    or (using ISO weekly calendar)   yyyy-Www-dThh:mm:ss[.xxx]\n");
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "    An example of a correct format is:  2000-01-01T12:00:00\n");
	}
	return (error);
}

bool GMT_get_time_system (struct GMT_CTRL *GMT, char *name, struct GMT_TIME_SYSTEM *time_system)
{
	/* Convert TIME_SYSTEM into TIME_EPOCH and TIME_UNIT.
	   TIME_SYSTEM can be one of the following: j2000, jd, mjd, s1985, unix, dr0001, rata
	   or any string in the form "TIME_UNIT since TIME_EPOCH", like "seconds since 1985-01-01".
	   This function only splits the strings, no validation or analysis is done.
	   See GMT_init_time_system_structure for that.
	   TIME_SYSTEM = other is completely ignored.
	*/
	char *epoch = NULL;

	if (!strcmp (name, "j2000")) {
		strcpy (time_system->epoch, "2000-01-01T12:00:00");
		time_system->unit = 'd';
	}
	else if (!strcmp (name, "jd")) {
		strcpy (time_system->epoch, "-4713-11-24T12:00:00");
		time_system->unit = 'd';
	}
	else if (!strcmp (name, "mjd")) {
		strcpy (time_system->epoch, "1858-11-17T00:00:00");
		time_system->unit = 'd';
	}
	else if (!strcmp (name, "s1985")) {
		strcpy (time_system->epoch, "1985-01-01T00:00:00");
		time_system->unit = 's';
	}
	else if (!strcmp (name, "unix")) {
		strcpy (time_system->epoch, "1970-01-01T00:00:00");
		time_system->unit = 's';
	}
	else if (!strcmp (name, "dr0001")) {
		strcpy (time_system->epoch, "0001-01-01T00:00:00");
		time_system->unit = 's';
	}
	else if (!strcmp (name, "rata")) {
		strcpy (time_system->epoch, "0000-12-31T00:00:00");
		time_system->unit = 'd';
	}
	else if (!strcmp (name, "other")) {
		/* Ignore completely */
	}
	else if ((epoch = strstr (name, "since"))) {
		epoch += 6;
		strncpy (time_system->epoch, epoch, GMT_LEN64);
		time_system->unit = name[0];
		if (!strncmp (name, "mon", 3U)) time_system->unit = 'o';
	}
	else
		return (true);
	return (false);
}

void gmt_grd_get_units (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *header)
{
	/* Set input data types for columns 0, 1 and 2 based on unit strings for
	   grid coordinates x, y and z.
	   When "Time": transform the data scale and offset to match the current time system.
	*/
	unsigned int i;
	char string[3][GMT_LEN256], *units = NULL;
	double scale = 1.0, offset = 0.0;
	struct GMT_TIME_SYSTEM time_system;

	/* Copy unit strings */
	strncpy (string[0], header->x_units, GMT_GRID_UNIT_LEN80);
	strncpy (string[1], header->y_units, GMT_GRID_UNIT_LEN80);
	strncpy (string[2], header->z_units, GMT_GRID_UNIT_LEN80);

	/* Parse the unit strings one by one */
	for (i = 0; i < 3; i++) {
		/* Skip parsing when input data type is already set */
		if (GMT->current.io.col_type[GMT_IN][i] & GMT_IS_GEO) continue;
		if (GMT->current.io.col_type[GMT_IN][i] & GMT_IS_RATIME) {
			GMT->current.proj.xyz_projection[i] = GMT_TIME;
			continue;
		}

		/* Change name of variable and unit to lower case for comparison */
		GMT_str_tolower (string[i]);

		if ((!strncmp (string[i], "longitude", 9U) || strstr (string[i], "degrees_e")) && (header->wesn[XLO] > -360.0 && header->wesn[XHI] <= 360.0)) {
			/* Input data type is longitude */
			GMT->current.io.col_type[GMT_IN][i] = GMT_IS_LON;
		}
		else if ((!strncmp (string[i], "latitude", 8U) || strstr (string[i], "degrees_n")) && (header->wesn[YLO] >= -90.0 && header->wesn[YHI] <= 90.0)) {
			/* Input data type is latitude */
			GMT->current.io.col_type[GMT_IN][i] = GMT_IS_LAT;
		}
		else if (!strcmp (string[i], "time") || !strncmp (string[i], "time [", 6U)) {
			/* Input data type is time */
			GMT->current.io.col_type[GMT_IN][i] = GMT_IS_RELTIME;
			GMT->current.proj.xyz_projection[i] = GMT_TIME;

			/* Determine coordinates epoch and units (default is internal system) */
			GMT_memcpy (&time_system, &GMT->current.setting.time_system, 1, struct GMT_TIME_SYSTEM);
			units = strchr (string[i], '[');
			//if (!units || GMT_get_time_system (GMT, ++units, &time_system) || GMT_init_time_system_structure (GMT, &time_system))
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: Time units [%s] in grid not recognised, defaulting to gmt.conf.\n", units);

			/* Determine scale between grid and internal time system, as well as the offset (in internal units) */
			scale = time_system.scale * GMT->current.setting.time_system.i_scale;
			offset = (time_system.rata_die - GMT->current.setting.time_system.rata_die) + (time_system.epoch_t0 - GMT->current.setting.time_system.epoch_t0);
			offset *= GMT_DAY2SEC_F * GMT->current.setting.time_system.i_scale;

			/* Scale data scale and extremes based on scale and offset */
			if (i == 0) {
				header->wesn[XLO] = header->wesn[XLO] * scale + offset;
				header->wesn[XHI] = header->wesn[XHI] * scale + offset;
				header->inc[GMT_X] *= scale;
			}
			else if (i == 1) {
				header->wesn[YLO] = header->wesn[YLO] * scale + offset;
				header->wesn[YHI] = header->wesn[YHI] * scale + offset;
				header->inc[GMT_Y] *= scale;
			}
			else {
				header->z_add_offset = header->z_add_offset * scale + offset;
				header->z_scale_factor *= scale;
			}
		}
	}
}

void gmt_grd_set_units (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *header)
{
	/* Set unit strings for grid coordinates x, y and z based on
	   output data types for columns 0, 1, and 2.
	*/
	unsigned int i;
	char *string[3] = {NULL, NULL, NULL}, unit[GMT_GRID_UNIT_LEN80] = {""};
	char date[GMT_LEN16] = {""}, clock[GMT_LEN16] = {""};

	/* Copy pointers to unit strings */
	string[0] = header->x_units;
	string[1] = header->y_units;
	string[2] = header->z_units;

	/* Use input data type as backup fr output data type */
	for (i = 0; i < 3; i++) 
		if (GMT->current.io.col_type[GMT_OUT][i] == GMT_IS_UNKNOWN) GMT->current.io.col_type[GMT_OUT][i] = GMT->current.io.col_type[GMT_IN][i];

	/* Catch some anomalies */
	if (GMT->current.io.col_type[GMT_OUT][GMT_X] == GMT_IS_LAT) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Output type for X-coordinate of grid %s is LAT. Replaced by LON.\n", header->name);
		GMT->current.io.col_type[GMT_OUT][GMT_X] = GMT_IS_LON;
	}
	if (GMT->current.io.col_type[GMT_OUT][GMT_Y] == GMT_IS_LON) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Output type for Y-coordinate of grid %s is LON. Replaced by LAT.\n", header->name);
		GMT->current.io.col_type[GMT_OUT][GMT_Y] = GMT_IS_LAT;
	}

	/* Set unit strings one by one based on output type */
	for (i = 0; i < 3; i++) {
		switch (GMT->current.io.col_type[GMT_OUT][i]) {
		case GMT_IS_LON:
			strcpy (string[i], "longitude [degrees_east]"); break;
		case GMT_IS_LAT:
			strcpy (string[i], "latitude [degrees_north]"); break;
		case GMT_IS_ABSTIME:
		case GMT_IS_RELTIME:
		case GMT_IS_RATIME:
			/* Determine time unit */
			switch (GMT->current.setting.time_system.unit) {
			case 'y':
				strcpy (unit, "years"); break;
			case 'o':
				strcpy (unit, "months"); break;
			case 'd':
				strcpy (unit, "days"); break;
			case 'h':
				strcpy (unit, "hours"); break;
			case 'm':
				strcpy (unit, "minutes"); break;
			default:
				strcpy (unit, "seconds"); break;
			}
			GMT_format_calendar (GMT, date, clock, &GMT->current.io.date_output, &GMT->current.io.clock_output, false, 1, 0.0);
			sprintf (string[i], "time [%s since %s %s]", unit, date, clock);
			/* Warning for non-double grids */
			if (i == 2 && GMT->session.grdformat[header->type][1] != 'd')
				//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: Use double precision output grid to avoid loss of significance of time coordinate.\n");
			break;
		}
	}
}

int GMT_update_grd_info (struct GMT_CTRL *GMT, char *file, struct GMT_GRID_HEADER *header)
{	/* file:	- IGNORED -
	 * header:	grid structure header
	 */

	/* pack z-range: */
	header->z_min = (header->z_min - header->z_add_offset) / header->z_scale_factor;
	header->z_max = (header->z_max - header->z_add_offset) / header->z_scale_factor;
	gmt_grd_set_units (GMT, header);
	return ((*GMT->session.updateinfo[header->type]) (GMT, header));
}
int GMT_grd_get_format (struct GMT_CTRL *GMT, char *file, struct GMT_GRID_HEADER *header, bool magic)
{
	/* This functions does a couple of things:
	 * 1. It tries to determine what kind of grid file this is. If a file is openeed for
	 *    reading we see if (a) a particular format has been specified with
	 *    the =<code> suffix, or (b) we are able to guess the format based on known
	 *    characteristics of various formats, or (c) assume the default grid format.
	 *    If a file is opened for writing, only option (a) and (c) apply.
	 *    If we cannot obtain the format we return an error.
	 * 2. We strip the suffix off. The relevant info is stored in the header struct.
	 * 3. In case of netCDF grids, the optional ?<varname> is stripped off as well.
	 *    The info is stored in header->varname.
	 * 4. If a file is open for reading, we set header->name to the full path of the file
	 *    by seaching in current dir and the various GMT_*DIR paths.
	 */

	size_t i = 0, j;
	int val;
	unsigned int direction = (magic) ? GMT_IN : GMT_OUT;
	char tmp[GMT_BUFSIZ];

	gmt_grd_parse_xy_units (GMT, header, file, direction);	/* Parse and strip xy scaling via +u<unit> modifier */

	gmt_expand_filename (GMT, file, header->name);	/* May append a suffix to header->name */

	/* Must reset scale and invalid value because sometimes headers from input grids are
	 * 'recycled' and used for output grids that may have a different type and z-range: */
	header->z_scale_factor = 1.0;
	header->z_add_offset   = 0.0;
	header->nan_value      = (float)NAN;

	i = strcspn (header->name, "="); /* get number of chars until first '=' or '\0' */

	if (header->name[i]) {	/* Reading or writing when =suffix is present: get format type, scale, offset and missing value */
		printf("did nt use \n");
	} /* if (header->name[i]) */
	else if (magic)
		{	/* Reading: determine file format automatically based on grid content */
		
		int choice = 0;
		sscanf (header->name, "%[^?]?%s", tmp, header->varname);    /* Strip off variable name */
#ifdef HAVE_GDAL
		/* Check if file is an URL */
		if (GMT_check_url_name(header->name)) {
			/* Then check for GDAL grid */
			if (GMT_is_gdal_grid (GMT, header) == GMT_NOERROR)
				return (GMT_NOERROR);
		}
#endif
		if (!GMT_getdatapath (GMT, tmp, header->name, R_OK))
			return (GMT_GRDIO_FILE_NOT_FOUND);	/* Possibly prepended a path from GMT_[GRID|DATA|IMG]DIR */
		/* First check if we have a netCDF grid. This MUST be first, because ?var needs to be stripped off. */
		if ((val = GMT_is_nc_grid (GMT, header)) == GMT_NOERROR)
			return (GMT_NOERROR);
		return (GMT_GRDIO_UNKNOWN_FORMAT);	/* No supported format found */
	}
	else {			/* Writing: get format type, scale, offset and missing value from GMT->current.setting.io_gridfile_format */
		if (sscanf (header->name, "%[^?]?%s", tmp, header->varname) > 1)
			strncpy (header->name, tmp, GMT_LEN256);    /* Strip off variable name */
		/* parse grid format string: */
		strcpy(GMT->current.setting.io_gridfile_format,"nf");
		//strcpy(GMT->current.setting.io_gridfile_format,"nc");
		//GMT->current.setting.io_gridfile_format = "nc";//hardcoded
		if ((val = parse_grd_format_scale (GMT, header, GMT->current.setting.io_gridfile_format)) != GMT_NOERROR)
		{
			
			return val;
		}
	}
	if (header->type == GMT_GRID_IS_AF)
		header->nan_value = 0.0f; /* NaN value for AGC format */
	
	return (GMT_NOERROR);
}
int GMT_write_grd_info (struct GMT_CTRL *GMT, char *file, struct GMT_GRID_HEADER *header)
{	/* file:	File name
	 * header:	grid structure header
	 */

	int err;	/* Implied by GMT_err_trap */

	//GMT_err_trap (GMT_grd_get_format (GMT, file, header, false));
	err = GMT_grd_get_format (GMT, file, header, false);
	if (err != GMT_NOERROR) return(err);
	gmt_grd_xy_scale (GMT, header, GMT_OUT);	/* Possibly scale wesn,inc */
	/* pack z-range: */
	header->z_min = (header->z_min - header->z_add_offset) / header->z_scale_factor;
	header->z_max = (header->z_max - header->z_add_offset) / header->z_scale_factor;
	gmt_grd_set_units (GMT, header);
	return ((*GMT->session.writeinfo[header->type]) (GMT, header));
}

size_t GMT_grd_data_size (struct GMT_CTRL *GMT, unsigned int format, float *nan_value)
{
	/* Determine size of data type and set NaN value, if not yet done so (integers only) */

	switch (GMT->session.grdformat[format][1]) {
		case 'b':
			if (isnan (*nan_value)) *nan_value = CHAR_MIN;
			return (sizeof (char));
			break;
		case 's':
			if (isnan (*nan_value)) *nan_value = SHRT_MIN;
			return (sizeof (int16_t));
			break;
		case 'i':
			if (isnan (*nan_value)) *nan_value = INT_MIN;
		case 'm':
			return (sizeof (int32_t));
			break;
		case 'f':
			return (sizeof (float));
			break;
		case 'd':
			return (sizeof (double));
			break;
		default:
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Unknown grid data type: %c\n", GMT->session.grdformat[format][1]);
			return (GMT_GRDIO_UNKNOWN_TYPE);
	}
}

bool GMT_grd_pad_status (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *header, unsigned int *pad)
{	/* Determines if this grid has padding at all (pad = NULL) OR
	 * if pad is given, determines if the pads are different.
	 * Return codes are:
	 * 1) If pad == NULL:
	 *    false: Grid has zero padding.
	 *    true:  Grid has non-zero padding.
	 * 2) If pad contains the desired pad:
	 *    true:  Grid padding matches pad exactly.
	 *    false: Grid padding failed to match pad exactly.
	 */
	unsigned int side;

	if (pad) {	/* Determine if the grid's pad differ from given pad (false) or not (true) */
		for (side = 0; side < 4; side++) if (header->pad[side] != pad[side]) return (false);	/* Pads differ */
		return (true);	/* Pads match */
	}
	else {	/* We just want to determine if the grid has padding already (true) or not (false) */
		for (side = 0; side < 4; side++) if (header->pad[side]) return (true);	/* Grid has a pad */
		return (false);	/* Grid has no pad */
	}
}

void grd_pad_off_sub (struct GMT_GRID *G, float *data)
{
	uint64_t ijp, ij0;
	unsigned int row;
	for (row = 0; row < G->header->ny; row++) {
		ijp = GMT_IJP (G->header, row, 0);	/* Index of start of this row's first column in padded grid  */
		ij0 = GMT_IJ0 (G->header, row, 0);	/* Index of start of this row's first column in unpadded grid */
		GMT_memcpy (&(data[ij0]), &(data[ijp]), G->header->nx, float);	/* Only copy the nx data values */
	}
}

void GMT_grd_pad_off (struct GMT_CTRL *GMT, struct GMT_GRID *G)
{	/* Shifts the grid contents so there is no pad.  The remainder of
	 * the array is not reset and should not be addressed, but
	 * we set it to zero just in case.
	 * If pad is zero then we do nothing.
	 */
	bool is_complex;
	uint64_t nm;
	if (G->header->arrangement == GMT_GRID_IS_INTERLEAVED) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Calling GMT_grd_pad_off on interleaved complex grid! Programming error?\n");
		return;
	}
	if (!GMT_grd_pad_status (GMT, G->header, NULL)) return;	/* No pad so nothing to do */

	/* Here, G has a pad which we need to eliminate */
	is_complex = (G->header->complex_mode & GMT_GRID_IS_COMPLEX_MASK);
	if (!is_complex || (G->header->complex_mode & GMT_GRID_IS_COMPLEX_REAL))
		grd_pad_off_sub (G, G->data);	/* Remove pad around real component only or entire normal grid */
	if (is_complex && (G->header->complex_mode & GMT_GRID_IS_COMPLEX_IMAG))
		grd_pad_off_sub (G, &G->data[G->header->size/2]);	/* Remove pad around imaginary component */
	nm = G->header->nm;	/* Number of nodes in one component */
	if (is_complex) nm *= 2;	/* But there might be two */
	if (G->header->size > nm) {	/* Just wipe the remaineder of the array to be sure */
		size_t n_to_cleen = G->header->size - nm;
		GMT_memset (&(G->data[nm]), n_to_cleen, float);	/* nm is 1st position after last row */
	}
	GMT_memset (G->header->pad, 4, int);	/* Pad is no longer active */
	GMT_set_grddim (GMT, G->header);		/* Update all dimensions to reflect the padding */
}

void gmt_grd_wipe_pad (struct GMT_CTRL *GMT, struct GMT_GRID *G)
{	/* Reset padded areas to 0. */
	unsigned int row;
	size_t ij0;
	
	if (G->header->arrangement == GMT_GRID_IS_INTERLEAVED) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Calling gmt_grd_wipe_pad on interleaved complex grid! Programming error?\n");
		return;
	}
	if (G->header->pad[YHI]) GMT_memset (G->data, G->header->mx * G->header->pad[YHI], float);	/* Wipe top pad */
	if (G->header->pad[YLO]) {	/* Wipe bottom pad */
		ij0 = GMT_IJ (G->header, G->header->my - G->header->pad[YLO], 0);	/* Index of start of bottom pad */
		GMT_memset (&(G->data[ij0]), G->header->mx * G->header->pad[YLO], float);
	}
	if (G->header->pad[XLO] == 0 && G->header->pad[XHI] == 0) return;	/* Nothing to do */
	for (row = G->header->pad[YHI]; row < G->header->my - G->header->pad[YLO]; row++) {	/* Wipe left and right pad which is trickier */
		ij0 = GMT_IJ (G->header, row, 0);	/* Index of this row's left column (1st entry in west pad) */
		if (G->header->pad[XLO]) GMT_memset (&(G->data[ij0]), G->header->pad[XLO], float);
		ij0 += (G->header->mx - G->header->pad[XHI]);	/* Start of this rows east pad's 1st column */
		if (G->header->pad[XHI]) GMT_memset (&(G->data[ij0]), G->header->pad[XHI], float);
	}
}

void grd_pad_on_sub (struct GMT_CTRL *GMT, struct GMT_GRID *G, struct GMT_GRID_HEADER *h_old, float *data)
{
	unsigned int row;
	uint64_t ij_new, ij_old;
	for (row = G->header->ny; row > 0; row--) {
		ij_new = GMT_IJP (G->header, row-1, 0);	/* Index of start of this row's first column in new padded grid  */
		ij_old = GMT_IJP (h_old, row-1, 0);	/* Index of start of this row's first column in old padded grid */
		GMT_memcpy (&(G->data[ij_new]), &(G->data[ij_old]), G->header->nx, float);
	}
	gmt_grd_wipe_pad (GMT, G);	/* Set pad areas to 0 */
}

void GMT_grd_zminmax (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *h, float *z)
{	/* Reset the xmin/zmax values in the header */
	unsigned int row, col;
	uint64_t node, n = 0;

	h->z_min = DBL_MAX;	h->z_max = -DBL_MAX;
	for (row = 0; row < h->ny; row++) {
		for (col = 0, node = GMT_IJP (h, row, 0); col < h->nx; col++, node++) {
			if (isnan (z[node]))
				continue;
			/* Update z_min, z_max */
			h->z_min = MIN (h->z_min, (double)z[node]);
			h->z_max = MAX (h->z_max, (double)z[node]);
			n++;
		}
	}
	if (n == 0) h->z_min = h->z_max = GMT->session.d_NaN;	/* No non-NaNs in the entire grid */
}

struct GMT_GRID_HEADER *GMT_duplicate_gridheader (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *h)
{	/* Duplicates a grid header. */
	struct GMT_GRID_HEADER *hnew = NULL;

	hnew = GMT_memory (GMT, NULL, 1, struct GMT_GRID_HEADER);
	GMT_memcpy (hnew, h, 1, struct GMT_GRID_HEADER);
	return (hnew);
}

void GMT_grd_pad_on (struct GMT_CTRL *GMT, struct GMT_GRID *G, unsigned int *pad)
{ /* Shift grid content from a non-padded (or differently padded) to a padded organization.
	 * We check that the grid size can handle this and allocate more space if needed.
	 * If pad matches the grid's pad then we do nothing.
	 */
	bool is_complex;
	size_t size;
	struct GMT_GRID_HEADER *h = NULL;

	if (G->header->arrangement == GMT_GRID_IS_INTERLEAVED) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Calling GMT_grd_pad_off on interleaved complex grid! Programming error?\n");
		return;
	}
	if (GMT_grd_pad_status (GMT, G->header, pad)) return;	/* Already padded as requested so nothing to do */
	if (pad[XLO] == 0 && pad[XHI] == 0 && pad[YLO] == 0 && pad[YHI] == 0) {	/* Just remove the existing pad entirely */
		GMT_grd_pad_off (GMT, G);
		return;
	}
	/* Here the pads differ (or G has no pad at all) */
	is_complex = (G->header->complex_mode & GMT_GRID_IS_COMPLEX_MASK);
	size = gmt_grd_get_nxpad (G->header, pad) * gmt_grd_get_nypad (G->header, pad);	/* New array size after pad is added */
	if (is_complex) size *= 2;	/* Twice the space for complex grids */
	if (size > G->header->size) {	/* Must allocate more space, but since no realloc for aligned memory we must do it the hard way */
		float *f = NULL;
		//GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "Extend grid via copy onto larger grid\n");
		f = GMT_memory_aligned (GMT, NULL, size, float);		/* New, larger grid size */
		GMT_memcpy (f, G->data, G->header->size, float);	/* Copy over previous grid values */
		GMT_free_aligned (GMT, G->data);				/* Free previous aligned grid memory */
		G->data = f;						/* Attach the new, larger aligned memory */
		G->header->size = size;					/* Update the size */
	}
	/* Because G may have a pad that is nonzero (but different from pad) we need a different header structure in the macros below */
	h = GMT_duplicate_gridheader (GMT, G->header);

	GMT_grd_setpad (GMT, G->header, pad);	/* Pad is now active and set to specified dimensions */
	GMT_set_grddim (GMT, G->header);		/* Update all dimensions to reflect the padding */
	if (is_complex && (G->header->complex_mode & GMT_GRID_IS_COMPLEX_IMAG))
		grd_pad_on_sub (GMT, G, h, &G->data[size/2]);	/* Add pad around imaginary component first */
	if (!is_complex || (G->header->complex_mode & GMT_GRID_IS_COMPLEX_REAL))
		grd_pad_on_sub (GMT, G, h, G->data);	/* Add pad around real component */
	GMT_free (GMT, h);	/* Done with this header */
}

int GMT_duplicate_univector (struct GMT_CTRL *GMT, union GMT_UNIVECTOR *u_out, union GMT_UNIVECTOR *u_in, unsigned int type, uint64_t n_rows)
{
	/* Allocate space for one univector according to data type */
	switch (type) {
		case GMT_UCHAR:  GMT_memcpy (u_out->uc1, u_in->uc1, n_rows, uint8_t);   break;
		case GMT_CHAR:   GMT_memcpy (u_out->sc1, u_in->sc1, n_rows, int8_t);    break;
		case GMT_USHORT: GMT_memcpy (u_out->ui2, u_in->ui2, n_rows, uint16_t);  break;
		case GMT_SHORT:  GMT_memcpy (u_out->si2, u_in->si2, n_rows, int16_t);   break;
		case GMT_UINT:   GMT_memcpy (u_out->ui4, u_in->ui4, n_rows, uint32_t);  break;
		case GMT_INT:    GMT_memcpy (u_out->si4, u_in->si4, n_rows, int32_t);   break;
		case GMT_ULONG:  GMT_memcpy (u_out->ui8, u_in->ui8, n_rows, uint64_t);  break;
		case GMT_LONG:   GMT_memcpy (u_out->si8, u_in->si8, n_rows, int64_t);   break;
		case GMT_FLOAT:  GMT_memcpy (u_out->f4,  u_in->f4,  n_rows, float);     break;
		case GMT_DOUBLE: GMT_memcpy (u_out->f8,  u_in->f8,  n_rows, double);    break;
	}
	return (GMT_OK);
}

struct GMT_MATRIX * GMT_duplicate_matrix (struct GMT_CTRL *GMT, struct GMT_MATRIX *M_in, bool duplicate_data)
{	/* Duplicates a matrix container - optionally duplicates the data array */
	struct GMT_MATRIX *M = NULL;
	M = GMT_memory (GMT, NULL, 1, struct GMT_MATRIX);
	GMT_memcpy (M, M_in, 1, struct GMT_MATRIX);
	GMT_memset (&M->data, 1, union GMT_UNIVECTOR);
	if (duplicate_data) {
		size_t size = M->n_rows * M->n_columns;
		if (GMT_alloc_univector (GMT, &(M->data), M->type, size)) return (NULL);
		GMT_duplicate_univector (GMT, &M->data, &M_in->data, M->type, size);
	}
	return (M);
}



void GMT_scale_and_offset_f (struct GMT_CTRL *GMT, float *data, size_t length, double scale, double offset) {
	/* Routine that does the data conversion and sanity checking before
	 * calling scale_and_offset_f() to scale and offset the data in a grid */
	float scale_f  = (float)scale;
	float offset_f = (float)offset;

	if (scale_f == 1.0 && offset_f == 0.0)
		return; /* No work needed */

	/* Sanity checks */
	if (!isnormal (scale)) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Scale must be a non-zero normalized number (%g).\n", scale);
		printf("Scale must be a non-zero normalized number (%g).\n", scale);
		scale_f = 1.0f;
	}
	if (!isfinite (offset)) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Offset must be a finite number (%g).\n", offset);
		printf( "Offset must be a finite number (%g).\n", offset);
		offset_f = 0.0f;
	}

	/* Call workhorse */
	scale_and_offset_f (data, length, scale_f, offset_f);
}
void GMT_pack_grid (struct GMT_CTRL *Ctrl, struct GMT_GRID_HEADER *header, float *grid, unsigned pack_mode) {
	size_t n_representations = 0; /* number of distinct values >= 0 that a signed integral type can represent */

	if (pack_mode == k_grd_pack && (header->z_scale_autoadust || header->z_offset_autoadust)) {
		switch (Ctrl->session.grdformat[header->type][1]) {
			case 'b':
				n_representations = 128;         /* exp2 (8 * sizeof (int8_t)) / 2 */
				break;
			case 's':
				n_representations = 32768;       /* exp2 (8 * sizeof (int16_t)) / 2 */
				break;
			case 'i':
				/* A single precision float's significand has a precision of 24 bits.
				 * In order to avoid round-off errors, we must not use all 2^32
				 * n_representations of an int32_t. */
				n_representations = 0x1000000;   /* exp2 (24) */
				break;
			/* default: do not auto-scale floating point */
		}
	}

	if (n_representations != 0) {
		/* Calculate auto-scale and offset */
		GMT_grd_zminmax (Ctrl, header, grid); /* Calculate z_min/z_max */
		if (header->z_offset_autoadust) {
			/* shift to center values around 0 but shift only by integral value */
			double z_range = header->z_max - header->z_min;
			if (isfinite (z_range))
				header->z_add_offset = rint(z_range / 2.0 + header->z_min);
		}
		if (header->z_scale_autoadust) {
			/* scale z-range to use all n_representations */
			double z_max = header->z_max - header->z_add_offset;
			double z_min = fabs(header->z_min - header->z_add_offset);
			double z_0_n_range = MAX (z_max, z_min); /* use [0,n] range because of signed int */
			--n_representations;                     /* subtract 1 for NaN value */
			if (isnormal (z_0_n_range))
				header->z_scale_factor = z_0_n_range / n_representations;
		}
	}

	/* Do actual packing/unpacking: */
	switch (pack_mode) {
		case k_grd_unpack:
			GMT_scale_and_offset_f (Ctrl, grid, header->size, header->z_scale_factor, header->z_add_offset);
			/* Adjust z-range in header: */
			header->z_min = header->z_min * header->z_scale_factor + header->z_add_offset;
			header->z_max = header->z_max * header->z_scale_factor + header->z_add_offset;
			break;
		case k_grd_pack:
			GMT_scale_and_offset_f (Ctrl, grid, header->size, 1.0/header->z_scale_factor, -header->z_add_offset/header->z_scale_factor);
			break;
		default:
			assert (false); /* GMT_pack_grid() called with illegal pack_mode */
	}
}

void GMT_grd_mux_demux (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *header, float *data, unsigned int desired_mode)
{	/* Multiplex and demultiplex complex grids.
 	 * Complex grids are read/written by dealing with just one component: real or imag.
	 * Thus, at read|write the developer must specify which component (GMT_CMPLX_REAL|IMAG).
	 * For fast disk i/o we read complex data in serial.  I.e., if we ask for GMT_CMPLX_REAL
	 * then the array will contain RRRRR....________, where ______ is unused space for the
	 * imaginary components.  Likewise, if we requested GMT_CMPLX_IMAG then the array will
	 * be returned as _______...IIIIIII....
	 * Operations like FFTs typically required the data to be interleaved, i.e., in the
	 * form RIRIRIRI.... Then, when the FFT work is done and we wish to write out the
	 * result we will need to demultiplex the array back to its serial RRRRR....IIIII
	 * format before writing takes place.
	 * GMT_grd_mux_demux performs either multiplex or demultiplex, depending on desired_mode.
	 * If grid is not complex then we just return doing nothing.
	 * Note: At this point the grid is mx * my and we visit all the nodes, including the pads.
	 * hence we use header->mx/my and GMT_IJ below.
	 */
	uint64_t row, col, col_1, col_2, left_node_1, left_node_2, offset, ij, ij2;
	float *array = NULL;
	
	if (! (desired_mode == GMT_GRID_IS_INTERLEAVED || desired_mode == GMT_GRID_IS_SERIAL)) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "GMT_grd_mux_demux called with inappropriate mode - skipped.\n");
		printf("GMT_grd_mux_demux called with inappropriate mode - skipped.\n");
		return;
	}
	if ((header->complex_mode & GMT_GRID_IS_COMPLEX_MASK) == 0) return;	/* Nuthin' to do */
	if (header->arrangement == desired_mode) return;				/* Already has the right layout */
	
	/* In most cases we will actually have RRRRR...______ or _____...IIIII..
	 * which means half the array is empty and it is easy to shuffle. However,
	 * in the case with actual RRRR...IIII there is no simple way to do this inplace;
	 * see http://stackoverflow.com/questions/1777901/array-interleaving-problem */
	
	if (desired_mode == GMT_GRID_IS_INTERLEAVED) {	/* Must multiplex the grid */
		if ((header->complex_mode & GMT_GRID_IS_COMPLEX_MASK) == GMT_GRID_IS_COMPLEX_MASK) {
			/* Transform from RRRRR...IIIII to RIRIRIRIRI... */
			/* Implement properly later; for now waste memory by duplicating [memory is cheap and plentiful] */
			/* One advantage is that the padding is all zero by virtue of the allocation */
			array = GMT_memory_aligned (GMT, NULL, header->size, float);
			offset = header->size / 2;	/* Position of 1st row in imag portion of RRRR...IIII... */
			for (row = 0; row < header->my; row++) {	/* Going from first to last row */
				for (col = 0; col < header->mx; col++) {
					ij = GMT_IJ (header, row, col);	/* Position of an 'R' in the RRRRR portion */
					ij2 = 2 * ij;
					array[ij2++] = data[ij];
					array[ij2] = data[ij+offset];
				}
			}
			GMT_memcpy (data, array, header->size, float);	/* Overwrite serial array with interleaved array */
			GMT_free (GMT, array);
		}
		else if (header->complex_mode & GMT_GRID_IS_COMPLEX_REAL) {
			/* Here we have RRRRRR..._________ and want R_R_R_R_... */
			for (row = header->my; row > 0; row--) {	/* Going from last to first row */
				left_node_1 = GMT_IJ (header, row-1, 0);	/* Start of row in RRRRR layout */
				left_node_2 = 2 * left_node_1;			/* Start of same row in R_R_R_ layout */
				for (col = header->mx, col_1 = col - 1, col_2 = 2*col - 1; col > 0; col--, col_1--) { /* Go from right to left */
					data[left_node_2+col_2] = 0.0f;	col_2--;	/* Set the Imag component to zero */
					data[left_node_2+col_2] = data[left_node_1+col_1];	col_2--;
				}
			}
		}
		else {
			/* Here we have _____...IIIII and want _I_I_I_I */
			offset = header->size / 2;	/* Position of 1st row in imag portion of ____...IIII... */
			for (row = 0; row < header->my; row++) {	/* Going from first to last row */
				left_node_1 = GMT_IJ (header, row, 0);		/* Start of row in _____IIII layout not counting ____*/
				left_node_2 = 2 * left_node_1;			/* Start of same row in _I_I_I... layout */
				left_node_1 += offset;				/* Move past length of all ____... */
				for (col_1 = 0, col_2 = 1; col_1 < header->mx; col_1++, col_2 += 2) {
					data[left_node_2+col_2] = data[left_node_1+col_1];
					data[left_node_1+col_1] = 0.0f;	/* Set the Real component to zero */
				}
			}
		}
	}
	else if (desired_mode == GMT_GRID_IS_SERIAL) {	/* Must demultiplex the grid */
		if ((header->complex_mode & GMT_GRID_IS_COMPLEX_MASK) == GMT_GRID_IS_COMPLEX_MASK) {
			/* Transform from RIRIRIRIRI... to RRRRR...IIIII  */
			/* Implement properly later; for now waste memory by duplicating [memory is cheap and plentiful] */
			/* One advantage is that the padding is all zero by virtue of the allocation */
			array = GMT_memory_aligned (GMT, NULL, header->size, float);
			offset = header->size / 2;	/* Position of 1st row in imag portion of RRRR...IIII... */
			for (row = 0; row < header->my; row++) {	/* Going from first to last row */
				for (col = 0; col < header->mx; col++) {
					ij = GMT_IJ (header, row, col);	/* Position of an 'R' in the RRRRR portion */
					ij2 = 2 * ij;
					array[ij] = data[ij2++];
					array[ij+offset] = data[ij2];
				}
			}
			GMT_memcpy (data, array, header->size, float);	/* Overwrite interleaved array with serial array */
			GMT_free (GMT, array);
		}
		else if (header->complex_mode & GMT_GRID_IS_COMPLEX_REAL) {
			/* Here we have R_R_R_R_... and want RRRRRR..._______  */
			for (row = 0; row < header->my; row++) {	/* Doing from first to last row */
				left_node_1 = GMT_IJ (header, row, 0);	/* Start of row in RRRRR... */
				left_node_2 = 2 * left_node_1;		/* Start of same row in R_R_R_R... layout */
				for (col_1 = col_2 = 0; col_1 < header->mx; col_1++, col_2 += 2) {
					data[left_node_1+col_1] = data[left_node_2+col_2];
				}
			}
			offset = header->size / 2;			/* Position of 1st _ in RRRR...____ */
			GMT_memset (&data[offset], offset, float);	/* Wipe _____ portion clean */
		}
		else {	/* Here we have _I_I_I_I and want _____...IIIII */
			offset = header->size / 2;	/* Position of 1st row in imag portion of ____...IIII... */
			for (row = header->my; row > 0; row--) {	/* Going from last to first row */
				left_node_1 = GMT_IJ (header, row, 0);	/* Start of row in _____IIII layout not counting ____*/
				left_node_2 = 2 * left_node_1;		/* Start of same row in _I_I_I... layout */
				left_node_1 += offset;			/* Move past length of all ____... */
				for (col = header->mx, col_1 = col - 1, col_2 = 2*col - 1; col > 0; col--, col_1--, col_2 -= 2) { /* Go from right to left */
					data[left_node_1+col_1] = data[left_node_2+col_2];
				}
			}
			GMT_memset (data, offset, float);	/* Wipe leading _____ portion clean */
		}
	}
	header->arrangement = desired_mode;
}

int gmt_grd_layout (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *header, float *grid, unsigned int complex_mode, unsigned int direction)
{	/* Checks or sets the array arrangement for a complex array */
	size_t needed_size;	/* Space required to hold both components of a complex grid */
	
	if ((complex_mode & GMT_GRID_IS_COMPLEX_MASK) == 0) return GMT_OK;	/* Regular, non-complex grid, nothing special to do */
	
	needed_size = 2ULL * header->mx * header->my;	/* For the complex array */
	if (header->size < needed_size) {
		//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Internal Error: Complex grid not large enough to hold both components!.\n");
		printf("Internal Error: Complex grid not large enough to hold both components!.\n");
		GMT_exit (GMT, EXIT_FAILURE); return EXIT_FAILURE;
	}
	if (direction == GMT_IN) {	/* About to read in a complex component; another one might have been read in earlier */
		if (header->arrangement == GMT_GRID_IS_INTERLEAVED) {
			//GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "Demultiplexing complex grid before reading can take place.\n");
			printf("Demultiplexing complex grid before reading can take place.\n");
			GMT_grd_mux_demux (GMT, header, grid, GMT_GRID_IS_SERIAL);
		}
		if ((header->complex_mode & GMT_GRID_IS_COMPLEX_MASK) == GMT_GRID_IS_COMPLEX_MASK) {	/* Already have both component; this will overwrite one of them */
			unsigned int type = (complex_mode & GMT_GRID_IS_COMPLEX_REAL) ? 0 : 1;
			char *kind[2] = {"read", "imaginary"};
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Overwriting previously stored %s component in complex grid.\n", kind[type]);
			printf("Overwriting previously stored %s component in complex grid.\n", kind[type]);
		}
		header->complex_mode |= complex_mode;	/* Update the grids complex mode */
	}
	else {	/* About to write out a complex component */
		if ((header->complex_mode & GMT_GRID_IS_COMPLEX_MASK) == 0) {	/* Not a complex grid */
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Internal Error: Asking to write out complex components from a non-complex grid.\n");
			printf("Internal Error: Asking to write out complex components from a non-complex grid.\n");
			GMT_exit (GMT, EXIT_FAILURE); return EXIT_FAILURE;
		}
		if ((complex_mode & GMT_GRID_IS_COMPLEX_REAL) && (header->complex_mode & GMT_GRID_IS_COMPLEX_REAL) == 0) {
			/* Programming error: Requesting to write real components when there are none */
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Internal Error: Complex grid has no real components that can be written to file.\n");
			printf("Internal Error: Complex grid has no real components that can be written to file.\n");
			GMT_exit (GMT, EXIT_FAILURE); return EXIT_FAILURE;
		}
		else if ((complex_mode & GMT_GRID_IS_COMPLEX_IMAG) && (header->complex_mode & GMT_GRID_IS_COMPLEX_IMAG) == 0) {
			/* Programming error: Requesting to write imag components when there are none */
			//GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Internal Error: Complex grid has no imaginary components that can be written to file.\n");
			printf("Internal Error: Complex grid has no imaginary components that can be written to file.\n");
			GMT_exit (GMT, EXIT_FAILURE); return EXIT_FAILURE;
		}
		if (header->arrangement == GMT_GRID_IS_INTERLEAVED) {	/* Must first demultiplex the grid */
			//GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "Demultiplexing complex grid before writing can take place.\n");
			printf("Demultiplexing complex grid before writing can take place.\n");
			GMT_grd_mux_demux (GMT, header, grid, GMT_GRID_IS_SERIAL);
		}
	}
	/* header->arrangment might now have been changed accordingly */
	return GMT_OK;
}
void gmt_grd_check_consistency (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *header, float *grid)
{	/* Enforce before writing a grid that periodic grids with repeating columns
	 * agree on the node values in those columns; if different replace with average.
	 * This only affects geographic grids of 360-degree extent with gridline registration.
	 * Also, if geographic grid with gridline registration, if the N or S pole row is present
	 * we ensure that they all have identical values, otherwise replace by mean value */
	unsigned int row = 0, col = 0;
	unsigned int we_conflicts = 0, p_conflicts = 0;
	uint64_t left = 0, right = 0, node = 0;

	if (header->registration == GMT_GRID_PIXEL_REG) return;	/* Not gridline registered */
	if (!GMT_is_geographic (GMT, GMT_OUT)) return;		/* Not geographic */
	if (header->wesn[YLO] == -90.0) {	/* Check consistency of S pole duplicates */
		double sum;
		node = GMT_IJP (header, 0, 0);	/* First node at S pole */
		sum = grid[node++];
		p_conflicts = 0;
		for (col = 1; col < header->nx; col++, node++) {
			if (grid[node] != grid[node-1]) p_conflicts++;
			sum += grid[node];
		}
		if (p_conflicts) {
			float f_value = (float)(sum / header->nx);
			//GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "Warning: detected %u inconsistent values at south pole. Values fixed by setting all to average row value.\n", p_conflicts);
			printf("Warning: detected %u inconsistent values at south pole. Values fixed by setting all to average row value.\n", p_conflicts);
			node = GMT_IJP (header, 0, 0);	/* First node at S pole */
			for (col = 0; col < header->nx; col++, node++) grid[node] = f_value;
		}
	}
	if (header->wesn[YHI] == +90.0) {	/* Check consistency of N pole duplicates */
		double sum;
		node = GMT_IJP (header, header->ny-1, 0);	/* First node at N pole */
		sum = grid[node++];
		p_conflicts = 0;
		for (col = 1; col < header->nx; col++, node++) {
			if (grid[node] != grid[node-1]) p_conflicts++;
			sum += grid[node];
		}
		if (p_conflicts) {
			float f_value = (float)(sum / header->nx);
			//GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "Warning: detected %u inconsistent values at north pole. Values fixed by setting all to average row value.\n", p_conflicts);
			printf("Warning: detected %u inconsistent values at north pole. Values fixed by setting all to average row value.\n", p_conflicts);
			node = GMT_IJP (header, header->ny-1, 0);	/* First node at N pole */
			for (col = 0; col < header->nx; col++, node++) grid[node] = f_value;
		}
	}
	if (!GMT_360_RANGE (header->wesn[XLO], header->wesn[XHI])) return;	/* Not 360-degree range */
	
	for (row = 0; row < header->ny; row++) {
		left = GMT_IJP (header, row, 0);	/* Left node */
		right = left + header->nx - 1;		/* Right node */
		if (grid[right] != grid[left]) {
			grid[right] = grid[left];
			we_conflicts++;
		}
	}
	if (we_conflicts)
	{
		//GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "Warning: detected %u inconsistent values along periodic east boundary of grid. Values fixed by duplicating west boundary.\n", we_conflicts);
		printf("Warning: detected %u inconsistent values along periodic east boundary of grid. Values fixed by duplicating west boundary.\n", we_conflicts);
	}
}

int GMT_write_grd (struct GMT_CTRL *GMT, char *file, struct GMT_GRID_HEADER *header, float *grid, double *wesn, unsigned int *pad, int complex_mode)
{	/* file:	File name
	 * header:	grid structure header
	 * grid:	array with final grid
	 * wesn:	Sub-region to write out  [Use entire file if NULL or contains 0,0,0,0]
	 * padding:	# of empty rows/columns to add on w, e, s, n of grid, respectively
	 * complex_mode:	&1 | &2 if complex array is to hold real (1) and imaginary (2) parts (otherwise read as real only)
	 *		Note: The file has only real values, we simply allow space in the array
	 *		for imaginary parts when processed by grdfft etc.
	 */

	int err;	/* Implied by GMT_err_trap */
	//if()
	//printf("iytqwer ");
	GMT_err_trap (GMT_grd_get_format (GMT, file, header, false));
	
	gmt_grd_set_units (GMT, header);
	
	GMT_pack_grid (GMT, header, grid, k_grd_pack); /* scale and offset */
	
	gmt_grd_xy_scale (GMT, header, GMT_OUT);	/* Possibly scale wesn,inc */
	
	gmt_grd_layout (GMT, header, grid, complex_mode, GMT_OUT);	/* Deal with complex layout */
	;
	gmt_grd_check_consistency (GMT, header, grid);			/* Fix east repeating columns and polar values */
	
	err = (*GMT->session.writegrd[header->type]) (GMT, header, grid, wesn, pad, complex_mode);
	if (GMT->parent->leave_grid_scaled == 0) GMT_pack_grid (GMT, header, grid, k_grd_unpack); /* revert scale and offset to leave grid as it was before writing unless session originated from gm*/
	
	return (err);
}

int GMT_grd_prep_io (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *header, double wesn[], unsigned int *width, unsigned int *height, int *first_col, int *last_col, int *first_row, int *last_row, unsigned int **index)
{
	/* Determines which rows and columns to extract to extract from a grid, based on w,e,s,n.
	 * This routine first rounds the w,e,s,n boundaries to the nearest gridlines or pixels,
	 * then determines the first and last columns and rows, and the width and height of the subset (in cells).
	 * The routine also returns and array of the x-indices in the source grid to be used in the target (subset) grid.
	 * All integers represented positive definite items hence unsigned variables.
	 */

	bool geo = false;
	unsigned int one_or_zero, col, *actual_col = NULL;	/* Column numbers */
	double small = 0.1, half_or_zero, x;
	//GMT_Report (GMT->parent, GMT_MSG_DEBUG, "region: %g %g, grid: %g %g\n", wesn[XLO], wesn[XHI], header->wesn[XLO], header->wesn[XHI]);

	half_or_zero = (header->registration == GMT_GRID_PIXEL_REG) ? 0.5 : 0.0;

	if (!GMT_is_subset (GMT, header, wesn)) {	/* Get entire file */
		*width  = header->nx;
		*height = header->ny;
		*first_col = *first_row = 0;
		*last_col  = header->nx - 1;
		*last_row  = header->ny - 1;
		GMT_memcpy (wesn, header->wesn, 4, double);
	}
	else {				/* Must deal with a subregion */
		if (GMT_x_is_lon (GMT, GMT_IN))
			geo = true;	/* Geographic data for sure */
		else if (wesn[XLO] < header->wesn[XLO] || wesn[XHI] > header->wesn[XHI])
			geo = true;	/* Probably dealing with periodic grid */

		if (wesn[YLO] < header->wesn[YLO] || wesn[YHI] > header->wesn[YHI]) return (GMT_GRDIO_DOMAIN_VIOLATION);	/* Calling program goofed... */

		one_or_zero = (header->registration == GMT_GRID_PIXEL_REG) ? 0 : 1;

		/* Make sure w,e,s,n are proper multiples of x_inc,y_inc away from x_min,y_min */

		//GMT_err_pass (GMT, GMT_adjust_loose_wesn (GMT, wesn, header), header->name);
		GMT_adjust_loose_wesn (GMT, wesn, header);
		/* Get dimension of subregion */

		*width  = urint ((wesn[XHI] - wesn[XLO]) * header->r_inc[GMT_X]) + one_or_zero;
		*height = urint ((wesn[YHI] - wesn[YLO]) * header->r_inc[GMT_Y]) + one_or_zero;

		/* Get first and last row and column numbers */

		*first_col = irint (floor ((wesn[XLO] - header->wesn[XLO]) * header->r_inc[GMT_X] + small));
		*last_col  = irint (ceil  ((wesn[XHI] - header->wesn[XLO]) * header->r_inc[GMT_X] - small)) - 1 + one_or_zero;
		*first_row = irint (floor ((header->wesn[YHI] - wesn[YHI]) * header->r_inc[GMT_Y] + small));
		*last_row  = irint (ceil  ((header->wesn[YHI] - wesn[YLO]) * header->r_inc[GMT_Y] - small)) - 1 + one_or_zero;
	}

	actual_col = GMT_memory (GMT, NULL, *width, unsigned int);
	if (geo) {
		small = 0.1 * header->inc[GMT_X];
		for (col = 0; col < (*width); col++) {
			x = GMT_col_to_x (GMT, col, wesn[XLO], wesn[XHI], header->inc[GMT_X], half_or_zero, *width);
			if (header->wesn[XLO] - x > small)
				x += 360.0;
			else if (x - header->wesn[XHI] > small)
				x -= 360.0;
			actual_col[col] = (unsigned int)GMT_grd_x_to_col (GMT, x, header);
		}
	}
	else {	/* Normal ordering */
		for (col = 0; col < (*width); col++) actual_col[col] = (*first_col) + col;
	}

	*index = actual_col;
	//GMT_Report (GMT->parent, GMT_MSG_DEBUG, "-> region: %g %g, grid: %g %g\n", wesn[XLO], wesn[XHI], header->wesn[XLO], header->wesn[XHI]);
	//printf("-> region: %g %g, grid: %g %g\n", wesn[XLO], wesn[XHI], header->wesn[XLO], header->wesn[XHI]);
	//GMT_Report (GMT->parent, GMT_MSG_DEBUG, "row: %d %d, col: %d %d\n", *first_row, *last_row, *first_col, *last_col);
	//printf("row: %d %d, col: %d %d\n", *first_row, *last_row, *first_col, *last_col);

	return (GMT_NOERROR);
}

double * GMT_grd_coord (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *header, int dir)
{	/* Allocate, compute, and return the x- or y-coordinates for a grid */
	unsigned int k;
	double *coord = NULL;
	assert (dir == GMT_X || dir == GMT_Y);
	if (dir == GMT_X) {
		coord = GMT_memory (GMT, NULL, header->nx, double);
		
		for (k = 0; k < header->nx; k++) coord[k] = GMT_grd_col_to_x (GMT, k, header);
	}
	else if (dir == GMT_Y) {
		coord = GMT_memory (GMT, NULL, header->ny, double);
		for (k = 0; k < header->ny; k++) coord[k] = GMT_grd_row_to_y (GMT, k, header);
	}
	
	return (coord);
}

int gmt_padspace (struct GMT_CTRL *GMT, struct GMT_GRID_HEADER *header, double *wesn, unsigned int *pad, struct GRD_PAD *P)
{	/* When padding is requested it is usually used to set boundary conditions based on
	 * two extra rows/columns around the domain of interest.  BCs like natural or periodic
	 * can then be used to fill in the pad.  However, if the domain is taken from a grid
	 * whose full domain exceeds the region of interest we are better off using the extra
	 * data to fill those pad rows/columns.  Thus, this function tries to determine if the
	 * input grid has the extra data we need to fill the BC pad with observations. */
	bool wrap;
	unsigned int side, n_sides = 0;
	double wesn2[4];

	/* First copy over original settings to the Pad structure */
	GMT_memset (P, 1, struct GRD_PAD);					/* Initialize to zero */
	GMT_memcpy (P->pad, pad, 4, int);					/* Duplicate the pad */
	if (!wesn) return (false);						/* No subset requested */
	if (wesn[XLO] == wesn[XHI] && wesn[YLO] == wesn[YHI]) return (false);	/* Subset not set */
	if (wesn[XLO] == header->wesn[XLO] && wesn[XHI] == header->wesn[XHI] && wesn[YLO] == header->wesn[YLO] && wesn[YHI] == header->wesn[YHI]) 
		return (false);	/* Subset equals whole area */
	GMT_memcpy (P->wesn, wesn, 4, double);					/* Copy the subset boundaries */
	if (pad[XLO] == 0 && pad[XHI] == 0 && pad[YLO] == 0 && pad[YHI] == 0) return (false);	/* No padding requested */

	/* Determine if data exist for a pad on all four sides.  If not we give up */
	wrap = GMT_grd_is_global (GMT, header);	/* If global wrap then we cannot be outside */
	if ((wesn2[XLO] = wesn[XLO] - pad[XLO] * header->inc[GMT_X]) < header->wesn[XLO] && !wrap)	/* Cannot extend west/xmin */
		{ n_sides++; wesn2[XLO] = wesn[XLO]; }
	else	/* OK to load left pad with data */
		P->pad[XLO] = 0;
	if ((wesn2[XHI] = wesn[XHI] + pad[XHI] * header->inc[GMT_X]) > header->wesn[XHI] && !wrap)	/* Cannot extend east/xmax */
		{ n_sides++; wesn2[XHI] = wesn[XHI]; }
	else	/* OK to load right pad with data */
		P->pad[XHI] = 0;
	if ((wesn2[YLO] = wesn[YLO] - pad[YLO] * header->inc[GMT_Y]) < header->wesn[YLO])	/* Cannot extend south/ymin */
		{ n_sides++; wesn2[YLO] = wesn[YLO]; }
	else	/* OK to load bottom pad with data */
		P->pad[YLO] = 0;
	if ((wesn2[YHI] = wesn[YHI] + pad[YHI] * header->inc[GMT_Y]) > header->wesn[YHI])	/* Cannot extend north/ymax */
		{ n_sides++; wesn2[YHI] = wesn[YHI]; }
	else	/* OK to load top pad with data */
		P->pad[YHI] = 0;
	if (n_sides == 4) return (false);	/* No can do */

	/* Here we know that there is enough input data to fill some or all of the BC pad with actual data values */
	/* We have temporarily set padding to zero (since the pad is now part of the region) for those sides we can extend */

	/* Temporarily enlarge the region so it now includes the padding we need */
	GMT_memcpy (P->wesn, wesn2, 4, double);

	/* Set BC */
	for (side = 0; side < 4; side++) {
		if (P->pad[side] == 0)
			header->BC[side] = GMT_BC_IS_DATA;
	}

	return (true);	/* Return true so the calling function can take appropriate action */
}

int GMT_read_grd (struct GMT_CTRL *GMT, char *file, struct GMT_GRID_HEADER *header, float *grid, double *wesn, unsigned int *pad, int complex_mode)
{	/* file:	- IGNORED -
	 * header:	grid structure header
	 * grid:	array with final grid
	 * wesn:	Sub-region to extract  [Use entire file if NULL or contains 0,0,0,0]
	 * padding:	# of empty rows/columns to add on w, e, s, n of grid, respectively
	 * complex_mode:	&1 | &2 if complex array is to hold real (1) and imaginary (2) parts (otherwise read as real only)
	 *		Note: The file has only real values, we simply allow space in the array
	 *		for imaginary parts when processed by grdfft etc.
	 */

	bool expand;		/* true or false */
	int err;		/* Implied by GMT_err_trap */
	struct GRD_PAD P;

	complex_mode &= GMT_GRID_IS_COMPLEX_MASK;	/* Remove any non-complex flags */
	/* If we are reading a 2nd grid (e.g., real, then imag) we must update info about the file since it will be a different file */
	if (header->complex_mode && (header->complex_mode & complex_mode) == 0) 
		GMT_err_trap (GMT_grd_get_format (GMT, file, header, true));
	
	expand = gmt_padspace (GMT, header, wesn, pad, &P);	/* true if we can extend the region by the pad-size to obtain real data for BC */

	gmt_grd_layout (GMT, header, grid, complex_mode & GMT_GRID_IS_COMPLEX_MASK, GMT_IN);	/* Deal with complex layout */

	
	GMT_err_trap ((*GMT->session.readgrd[header->type]) (GMT, header, grid, P.wesn, P.pad, complex_mode));

	if (expand) /* Must undo the region extension and reset nx, ny using original pad  */
		GMT_memcpy (header->wesn, wesn, 4, double);
	header->grdtype = gmt_get_grdtype (GMT, header);	/* Since may change if a subset */
	GMT_grd_setpad (GMT, header, pad);	/* Copy the pad to the header */
	GMT_set_grddim (GMT, header);		/* Update all dimensions */
	if (expand) GMT_grd_zminmax (GMT, header, grid);	/* Reset min/max since current extrema includes the padded region */
	GMT_pack_grid (GMT, header, grid, k_grd_unpack); /* revert scale and offset */
	GMT_BC_init (GMT, header);	/* Initialize grid interpolation and boundary condition parameters */

	return (GMT_NOERROR);
}

int GMT_grd_BC_set (struct GMT_CTRL *GMT, struct GMT_GRID *G, unsigned int direction)
{
	/* Set two rows of padding (pad[] can be larger) around data according
	   to desired boundary condition info in that header.
	   Returns -1 on problem, 0 on success.
	   If either x or y is periodic, the padding is entirely set.
	   However, if neither is true (this rules out geographical also)
	   then all but three corner-most points in each corner are set.

	   As written, not ready to use with "surface" for GMT 5, because
	   assumes left/right is +/- 1 and down/up is +/- mx.  In "surface"
	   the amount to move depends on the current mesh size, a parameter
	   not used here.

	   This is the revised, two-rows version (WHFS 6 May 1998).
	*/

	uint64_t mx;		/* Width of padded array; width as malloc'ed  */
	uint64_t mxnyp;		/* distance to periodic constraint in j direction  */
	uint64_t i, jmx;		/* Current i, j * mx  */
	uint64_t nxp2;		/* 1/2 the xg period (180 degrees) in cells  */
	uint64_t i180;		/* index to 180 degree phase shift  */
	uint64_t iw, iwo1, iwo2, iwi1, ie, ieo1, ieo2, iei1;  /* see below  */
	uint64_t jn, jno1, jno2, jni1, js, jso1, jso2, jsi1;  /* see below  */
	uint64_t jno1k, jno2k, jso1k, jso2k, iwo1k, iwo2k, ieo1k, ieo2k;
	uint64_t j1p, j2p;	/* j_o1 and j_o2 pole constraint rows  */
	unsigned int n_skip, n_set;
	unsigned int bok;		/* bok used to test that things are OK  */
	bool set[4] = {true, true, true, true};

	char *kind[5] = {"not set", "natural", "periodic", "geographic", "extended data"};
	char *edge[4] = {"left  ", "right ", "bottom", "top   "};

	if (G->header->complex_mode & GMT_GRID_IS_COMPLEX_MASK) return (GMT_NOERROR);	/* Only set up for real arrays */
	if (G->header->no_BC) return (GMT_NOERROR);	/* Told not to deal with BC stuff */

	for (i = n_skip = 0; i < 4; i++) {
		if (G->header->BC[i] == GMT_BC_IS_DATA) {set[i] = false; n_skip++;}	/* No need to set since there is data in the pad area */
	}
	if (n_skip == 4) {	/* No need to set anything since there is data in the pad area on all sides */
		//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "All boundaries set via extended data.\n");
		printf("All boundaries set via extended data.\n");
		return (GMT_NOERROR);
	}

	/* Check minimum size:  */
	if (G->header->nx < 1 || G->header->ny < 1) {
		//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "requires nx,ny at least 1.\n");
		printf("requires nx,ny at least 1.\n");
		return (GMT_NOERROR);
	}

	/* Check that pad is at least 2 */
	for (i = bok = 0; i < 4; i++) if (G->header->pad[i] < 2) bok++;
	if (bok > 0) {
		if (direction == GMT_IN) 
			{
			//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "called with a pad < 2; skipped.\n");
			printf("called with a pad < 2; skipped.\n");
			}
		return (GMT_NOERROR);
	}

	/* Initialize stuff:  */

	mx = G->header->mx;
	nxp2 = G->header->nxp / 2;	/* Used for 180 phase shift at poles  */

	iw = G->header->pad[XLO];	/* i for west-most data column */
	iwo1 = iw - 1;		/* 1st column outside west  */
	iwo2 = iwo1 - 1;	/* 2nd column outside west  */
	iwi1 = iw + 1;		/* 1st column  inside west  */

	ie = G->header->pad[XLO] + G->header->nx - 1;	/* i for east-most data column */
	ieo1 = ie + 1;		/* 1st column outside east  */
	ieo2 = ieo1 + 1;	/* 2nd column outside east  */
	iei1 = ie - 1;		/* 1st column  inside east  */

	jn = mx * G->header->pad[YHI];	/* j*mx for north-most data row  */
	jno1 = jn - mx;		/* 1st row outside north  */
	jno2 = jno1 - mx;	/* 2nd row outside north  */
	jni1 = jn + mx;		/* 1st row  inside north  */

	js = mx * (G->header->pad[YHI] + G->header->ny - 1);	/* j*mx for south-most data row  */
	jso1 = js + mx;		/* 1st row outside south  */
	jso2 = jso1 + mx;	/* 2nd row outside south  */
	jsi1 = js - mx;		/* 1st row  inside south  */

	mxnyp = mx * G->header->nyp;

	jno1k = jno1 + mxnyp;	/* data rows periodic to boundary rows  */
	jno2k = jno2 + mxnyp;
	jso1k = jso1 - mxnyp;
	jso2k = jso2 - mxnyp;

	iwo1k = iwo1 + G->header->nxp;	/* data cols periodic to bndry cols  */
	iwo2k = iwo2 + G->header->nxp;
	ieo1k = ieo1 - G->header->nxp;
	ieo2k = ieo2 - G->header->nxp;

	/* Duplicate rows and columns if nx or ny equals 1 */

	if (G->header->nx == 1) for (i = jn+iw; i <= js+iw; i += mx) G->data[i-1] = G->data[i+1] = G->data[i];
	if (G->header->ny == 1) for (i = jn+iw; i <= jn+ie; i++) G->data[i-mx] = G->data[i+mx] = G->data[i];

	/* Check poles for grid case.  It would be nice to have done this
		in GMT_boundcond_param_prep() but at that point the data
		array isn't passed into that routine, and may not have been
		read yet.  Also, as coded here, this bombs with error if
		the pole data is wrong.  But there could be an option to
		to change the condition to Natural in that case, with warning.  */

	if (G->header->registration == GMT_GRID_NODE_REG) {	/* A pole can only be a grid node with gridline registration */
		if (G->header->gn) {	/* North pole case */
			bok = 0;
			if (GMT_is_fnan (G->data[jn + iw])) {	/* First is NaN so all should be NaN */
				for (i = iw+1; i <= ie; i++) if (!GMT_is_fnan (G->data[jn + i])) bok++;
			}
			else {	/* First is not NaN so all should be identical */
				for (i = iw+1; i <= ie; i++) if (G->data[jn + i] != G->data[jn + iw]) bok++;
			}
			if (bok > 0)
				//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "Warning: %d (of %d) inconsistent grid values at North pole.\n", bok, G->header->nx);
				printf("Warning: %d (of %d) inconsistent grid values at North pole.\n", bok, G->header->nx);
		}

		if (G->header->gs) {	/* South pole case */
			bok = 0;
			if (GMT_is_fnan (G->data[js + iw])) {	/* First is NaN so all should be NaN */
				for (i = iw+1; i <= ie; i++) if (!GMT_is_fnan (G->data[js + i])) bok++;
			}
			else {	/* First is not NaN so all should be identical */
				for (i = iw+1; i <= ie; i++) if (G->data[js + i] != G->data[js + iw]) bok++;
			}
			if (bok > 0)
				//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "Warning: %d (of %d) inconsistent grid values at South pole.\n", bok, G->header->nx);
				printf( "Warning: %d (of %d) inconsistent grid values at South pole.\n", bok, G->header->nx);
		}
	}

	/* Start with the case that x is not periodic, because in that case we also know that y cannot be polar.  */

	if (G->header->nxp <= 0) {	/* x is not periodic  */

		if (G->header->nyp > 0) {	/* y is periodic  */

			for (i = iw, bok = 0; i <= ie; ++i) {
				if (G->header->registration == GMT_GRID_NODE_REG && !doubleAlmostEqualZero (G->data[jn+i], G->data[js+i]))
					++bok;
				if (set[YHI]) {
					G->data[jno1 + i] = G->data[jno1k + i];
					G->data[jno2 + i] = G->data[jno2k + i];
				}
				if (set[YLO]) {
					G->data[jso1 + i] = G->data[jso1k + i];
					G->data[jso2 + i] = G->data[jso2k + i];
				}
			}
			if (bok > 0)
				//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "Warning: %d (of %d) inconsistent grid values at South and North boundaries for repeated nodes.\n", bok, G->header->nx);
				printf("Warning: %d (of %d) inconsistent grid values at South and North boundaries for repeated nodes.\n", bok, G->header->nx);
			/* periodic Y rows copied.  Now do X naturals.
				This is easy since y's are done; no corner problems.
				Begin with Laplacian = 0, and include 1st outside rows
				in loop, since y's already loaded to 2nd outside.  */

			for (jmx = jno1; jmx <= jso1; jmx += mx) {
				if (set[XLO]) G->data[jmx + iwo1] = (float)(4.0 * G->data[jmx + iw]) - (G->data[jmx + iw + mx] + G->data[jmx + iw - mx] + G->data[jmx + iwi1]);
				if (set[XHI]) G->data[jmx + ieo1] = (float)(4.0 * G->data[jmx + ie]) - (G->data[jmx + ie + mx] + G->data[jmx + ie - mx] + G->data[jmx + iei1]);
			}

			/* Copy that result to 2nd outside row using periodicity.  */
			if (set[XLO]) {
				G->data[jno2 + iwo1] = G->data[jno2k + iwo1];
				G->data[jso2 + iwo1] = G->data[jso2k + iwo1];
			}
			if (set[XHI]) {
				G->data[jno2 + ieo1] = G->data[jno2k + ieo1];
				G->data[jso2 + ieo1] = G->data[jso2k + ieo1];
			}

			/* Now set d[laplacian]/dx = 0 on 2nd outside column.  Include 1st outside rows in loop.  */
			for (jmx = jno1; jmx <= jso1; jmx += mx) {
				if (set[XLO]) G->data[jmx + iwo2] = (G->data[jmx + iw - mx] + G->data[jmx + iw + mx] + G->data[jmx + iwi1])
					- (G->data[jmx + iwo1 - mx] + G->data[jmx + iwo1 + mx]) + (float)(5.0 * (G->data[jmx + iwo1] - G->data[jmx + iw]));
				if (set[XHI]) G->data[jmx + ieo2] = (G->data[jmx + ie - mx] + G->data[jmx + ie + mx] + G->data[jmx + iei1])
					- (G->data[jmx + ieo1 - mx] + G->data[jmx + ieo1 + mx]) + (float)(5.0 * (G->data[jmx + ieo1] - G->data[jmx + ie]));
			}

			/* Now copy that result also, for complete periodicity's sake  */
			if (set[XLO]) {
				G->data[jno2 + iwo2] = G->data[jno2k + iwo2];
				G->data[jso2 + iwo2] = G->data[jso2k + iwo2];
				G->header->BC[XLO] = GMT_BC_IS_NATURAL;
			}
			if (set[XHI]) {
				G->data[jno2 + ieo2] = G->data[jno2k + ieo2];
				G->data[jso2 + ieo2] = G->data[jso2k + ieo2];
				G->header->BC[XHI] = GMT_BC_IS_NATURAL;
			}

			/* DONE with X not periodic, Y periodic case.  Fully loaded.  */
			if (set[YLO] && set[YHI]) {
				G->header->BC[YLO] = G->header->BC[YHI] = GMT_BC_IS_PERIODIC;
				//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "Set boundary condition for bottom and top edges: %s\n", kind[G->header->BC[YLO]]);
				printf( "Set boundary condition for bottom and top edges: %s\n", kind[G->header->BC[YLO]]);
			}
			else if (set[YLO]) {
				G->header->BC[YLO] = GMT_BC_IS_PERIODIC;
				//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "Set boundary condition for %s edge: %s\n", edge[YLO], kind[G->header->BC[YLO]]);
			}
			else if (set[YHI]) {
				G->header->BC[YHI] = GMT_BC_IS_PERIODIC;
				//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "Set boundary condition for %s edge: %s\n", edge[YHI], kind[G->header->BC[YHI]]);
			}

			return (GMT_NOERROR);
		}
		else {	/* Here begins the X not periodic, Y not periodic case  */

			/* First, set corner points.  Need not merely Laplacian(f) = 0
				but explicitly that d2f/dx2 = 0 and d2f/dy2 = 0.
				Also set d2f/dxdy = 0.  Then can set remaining points.  */

	/* d2/dx2 */	if (set[XLO]) G->data[jn + iwo1]   = (float)(2.0 * G->data[jn + iw] - G->data[jn + iwi1]);
	/* d2/dy2 */	if (set[YHI]) G->data[jno1 + iw]   = (float)(2.0 * G->data[jn + iw] - G->data[jni1 + iw]);
	/* d2/dxdy */	if (set[XLO] && set[YHI]) G->data[jno1 + iwo1] = -(G->data[jno1 + iwi1] - G->data[jni1 + iwi1] + G->data[jni1 + iwo1]);

	/* d2/dx2 */	if (set[XHI]) G->data[jn + ieo1]   = (float)(2.0 * G->data[jn + ie] - G->data[jn + iei1]);
	/* d2/dy2 */	if (set[YHI]) G->data[jno1 + ie]   = (float)(2.0 * G->data[jn + ie] - G->data[jni1 + ie]);
	/* d2/dxdy */	if (set[XHI] && set[YHI]) G->data[jno1 + ieo1] = -(G->data[jno1 + iei1] - G->data[jni1 + iei1] + G->data[jni1 + ieo1]);

	/* d2/dx2 */	if (set[XLO]) G->data[js + iwo1]   = (float)(2.0 * G->data[js + iw] - G->data[js + iwi1]);
	/* d2/dy2 */	if (set[YLO]) G->data[jso1 + iw]   = (float)(2.0 * G->data[js + iw] - G->data[jsi1 + iw]);
	/* d2/dxdy */	if (set[XLO] && set[YLO]) G->data[jso1 + iwo1] = -(G->data[jso1 + iwi1] - G->data[jsi1 + iwi1] + G->data[jsi1 + iwo1]);

	/* d2/dx2 */	if (set[XHI]) G->data[js + ieo1]   = (float)(2.0 * G->data[js + ie] - G->data[js + iei1]);
	/* d2/dy2 */	if (set[YLO]) G->data[jso1 + ie]   = (float)(2.0 * G->data[js + ie] - G->data[jsi1 + ie]);
	/* d2/dxdy */	if (set[XHI] && set[YLO]) G->data[jso1 + ieo1] = -(G->data[jso1 + iei1] - G->data[jsi1 + iei1] + G->data[jsi1 + ieo1]);

			/* Now set Laplacian = 0 on interior edge points, skipping corners:  */
			for (i = iwi1; i <= iei1; i++) {
				if (set[YHI]) G->data[jno1 + i] = (float)(4.0 * G->data[jn + i]) - (G->data[jn + i - 1] + G->data[jn + i + 1] + G->data[jni1 + i]);
				if (set[YLO]) G->data[jso1 + i] = (float)(4.0 * G->data[js + i]) - (G->data[js + i - 1] + G->data[js + i + 1] + G->data[jsi1 + i]);
			}
			for (jmx = jni1; jmx <= jsi1; jmx += mx) {
				if (set[XLO]) G->data[iwo1 + jmx] = (float)(4.0 * G->data[iw + jmx]) - (G->data[iw + jmx + mx] + G->data[iw + jmx - mx] + G->data[iwi1 + jmx]);
				if (set[XHI]) G->data[ieo1 + jmx] = (float)(4.0 * G->data[ie + jmx]) - (G->data[ie + jmx + mx] + G->data[ie + jmx - mx] + G->data[iei1 + jmx]);
			}

			/* Now set d[Laplacian]/dn = 0 on all edge pts, including
				corners, since the points needed in this are now set.  */
			for (i = iw; i <= ie; i++) {
				if (set[YHI]) G->data[jno2 + i] = G->data[jni1 + i] + (float)(5.0 * (G->data[jno1 + i] - G->data[jn + i]))
					+ (G->data[jn + i - 1] - G->data[jno1 + i - 1]) + (G->data[jn + i + 1] - G->data[jno1 + i + 1]);
				if (set[YLO]) G->data[jso2 + i] = G->data[jsi1 + i] + (float)(5.0 * (G->data[jso1 + i] - G->data[js + i]))
					+ (G->data[js + i - 1] - G->data[jso1 + i - 1]) + (G->data[js + i + 1] - G->data[jso1 + i + 1]);
			}
			for (jmx = jn; jmx <= js; jmx += mx) {
				if (set[XLO]) G->data[iwo2 + jmx] = G->data[iwi1 + jmx] + (float)(5.0 * (G->data[iwo1 + jmx] - G->data[iw + jmx]))
					+ (G->data[iw + jmx - mx] - G->data[iwo1 + jmx - mx]) + (G->data[iw + jmx + mx] - G->data[iwo1 + jmx + mx]);
				if (set[XHI]) G->data[ieo2 + jmx] = G->data[iei1 + jmx] + (float)(5.0 * (G->data[ieo1 + jmx] - G->data[ie + jmx]))
					+ (G->data[ie + jmx - mx] - G->data[ieo1 + jmx - mx]) + (G->data[ie + jmx + mx] - G->data[ieo1 + jmx + mx]);
			}
			/* DONE with X not periodic, Y not periodic case.  Loaded all but three cornermost points at each corner.  */

			for (i = n_set = 0; i < 4; i++) if (set[i]) {
				n_set++;
				G->header->BC[i] = GMT_BC_IS_NATURAL;
			}
			if (n_set == 4) {
				//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "Set boundary condition for all edges: %s\n", kind[G->header->BC[XLO]]);
			}
			for (i = 0; i < 4; i++) if (set[i]) {
				//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "Set boundary condition for %s edge: %s\n", edge[i], kind[G->header->BC[i]]);
				printf("Set boundary condition for %s edge: %s\n", edge[i], kind[G->header->BC[i]]);
			}
			return (GMT_NOERROR);
		}
		/* DONE with all X not periodic cases  */
	}
	else {	/* X is periodic.  Load x cols first, then do Y cases.  */
		if (set[XLO]) G->header->BC[XLO] = GMT_BC_IS_PERIODIC;
		if (set[XHI]) G->header->BC[XHI] = GMT_BC_IS_PERIODIC;
		for (jmx = jn, bok = 0; jmx <= js; jmx += mx) {
			if (G->header->registration == GMT_GRID_NODE_REG && !doubleAlmostEqualZero (G->data[jmx+iw], G->data[jmx+ie]))
				++bok;
			if (set[XLO]) {
				G->data[iwo1 + jmx] = G->data[iwo1k + jmx];
				G->data[iwo2 + jmx] = G->data[iwo2k + jmx];
			}
			if (set[XHI]) {
				G->data[ieo1 + jmx] = G->data[ieo1k + jmx];
				G->data[ieo2 + jmx] = G->data[ieo2k + jmx];
			}
		}
		if (bok > 0)
			//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "Warning: %d (of %d) inconsistent grid values at West and East boundaries for repeated nodes.\n", bok, G->header->ny);
			printf("Warning: %d (of %d) inconsistent grid values at West and East boundaries for repeated nodes.\n", bok, G->header->ny);
		if (G->header->nyp > 0) {	/* Y is periodic.  copy all, including boundary cols:  */
			for (i = iwo2, bok = 0; i <= ieo2; ++i) {
				if (G->header->registration == GMT_GRID_NODE_REG && !doubleAlmostEqualZero (G->data[jn+i], G->data[js+i]))
					++bok;
				if (set[YHI]) {
					G->data[jno1 + i] = G->data[jno1k + i];
					G->data[jno2 + i] = G->data[jno2k + i];
				}
				if (set[YLO]) {
					G->data[jso1 + i] = G->data[jso1k + i];
					G->data[jso2 + i] = G->data[jso2k + i];
				}
			}
			if (bok > 0)
				//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "Warning: %d (of %d) inconsistent grid values at South and North boundaries for repeated nodes.\n", bok, G->header->nx);
			/* DONE with X and Y both periodic.  Fully loaded.  */
			printf("Warning: %d (of %d) inconsistent grid values at South and North boundaries for repeated nodes.\n", bok, G->header->nx);

			if (set[YLO] && set[YHI]) {
				G->header->BC[YLO] = G->header->BC[YHI] = GMT_BC_IS_PERIODIC;
				//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "Set boundary condition for bottom and top edges: %s\n", kind[G->header->BC[YLO]]);
			}
			else if (set[YLO]) {
				G->header->BC[YLO] = GMT_BC_IS_PERIODIC;
				//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "Set boundary condition for %s edge: %s\n", edge[YLO], kind[G->header->BC[YLO]]);
			}
			else if (set[YHI]) {
				G->header->BC[YHI] = GMT_BC_IS_PERIODIC;
				//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "Set boundary condition for %s edge: %s\n", edge[YHI], kind[G->header->BC[YHI]]);
			}
			return (GMT_NOERROR);
		}

		/* Do north (top) boundary:  */

		if (G->header->gn) {	/* Y is at north pole.  Phase-shift all, incl. bndry cols. */
			if (G->header->registration == GMT_GRID_PIXEL_REG) {
				j1p = jn;	/* constraint for jno1  */
				j2p = jni1;	/* constraint for jno2  */
			}
			else {
				j1p = jni1;		/* constraint for jno1  */
				j2p = jni1 + mx;	/* constraint for jno2  */
			}
			for (i = iwo2; set[YHI] && i <= ieo2; i++) {
				i180 = G->header->pad[XLO] + ((i + nxp2)%G->header->nxp);
				G->data[jno1 + i] = G->data[j1p + i180];
				G->data[jno2 + i] = G->data[j2p + i180];
			}
			if (set[YHI]) {
				G->header->BC[YHI] = GMT_BC_IS_GEO;
				//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "Set boundary condition for %s edge: %s\n", edge[YHI], kind[G->header->BC[YHI]]);
			}
		}
		else {
			/* Y needs natural conditions.  x bndry cols periodic.
				First do Laplacian.  Start/end loop 1 col outside,
				then use periodicity to set 2nd col outside.  */

			for (i = iwo1; set[YHI] && i <= ieo1; i++) {
				G->data[jno1 + i] = (float)(4.0 * G->data[jn + i]) - (G->data[jn + i - 1] + G->data[jn + i + 1] + G->data[jni1 + i]);
			}
			if (set[XLO] && set[YHI]) G->data[jno1 + iwo2] = G->data[jno1 + iwo2 + G->header->nxp];
			if (set[XHI] && set[YHI]) G->data[jno1 + ieo2] = G->data[jno1 + ieo2 - G->header->nxp];


			/* Now set d[Laplacian]/dn = 0, start/end loop 1 col out,
				use periodicity to set 2nd out col after loop.  */

			for (i = iwo1; set[YHI] && i <= ieo1; i++) {
				G->data[jno2 + i] = G->data[jni1 + i] + (float)(5.0 * (G->data[jno1 + i] - G->data[jn + i]))
					+ (G->data[jn + i - 1] - G->data[jno1 + i - 1]) + (G->data[jn + i + 1] - G->data[jno1 + i + 1]);
			}
			if (set[XLO] && set[YHI]) G->data[jno2 + iwo2] = G->data[jno2 + iwo2 + G->header->nxp];
			if (set[XHI] && set[YHI]) G->data[jno2 + ieo2] = G->data[jno2 + ieo2 - G->header->nxp];

			/* End of X is periodic, north (top) is Natural.  */
			if (set[YHI]) {
				G->header->BC[YHI] = GMT_BC_IS_NATURAL;
				//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "Set boundary condition for %s edge: %s\n", edge[YHI], kind[G->header->BC[YHI]]);
			}
		}

		/* Done with north (top) BC in X is periodic case.  Do south (bottom)  */

		if (G->header->gs) {	/* Y is at south pole.  Phase-shift all, incl. bndry cols. */
			if (G->header->registration == GMT_GRID_PIXEL_REG) {
				j1p = js;	/* constraint for jso1  */
				j2p = jsi1;	/* constraint for jso2  */
			}
			else {
				j1p = jsi1;		/* constraint for jso1  */
				j2p = jsi1 - mx;	/* constraint for jso2  */
			}
			for (i = iwo2; set[YLO] && i <= ieo2; i++) {
				i180 = G->header->pad[XLO] + ((i + nxp2)%G->header->nxp);
				G->data[jso1 + i] = G->data[j1p + i180];
				G->data[jso2 + i] = G->data[j2p + i180];
			}
			if (set[YLO]) {
				G->header->BC[YLO] = GMT_BC_IS_GEO;
				//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "Set boundary condition for %s edge: %s\n", edge[YLO], kind[G->header->BC[YLO]]);
			}
		}
		else {
			/* Y needs natural conditions.  x bndry cols periodic.
				First do Laplacian.  Start/end loop 1 col outside,
				then use periodicity to set 2nd col outside.  */

			for (i = iwo1; set[YLO] && i <= ieo1; i++) {
				G->data[jso1 + i] = (float)(4.0 * G->data[js + i]) - (G->data[js + i - 1] + G->data[js + i + 1] + G->data[jsi1 + i]);
			}
			if (set[XLO] && set[YLO]) G->data[jso1 + iwo2] = G->data[jso1 + iwo2 + G->header->nxp];
			if (set[XHI] && set[YHI]) G->data[jso1 + ieo2] = G->data[jso1 + ieo2 - G->header->nxp];


			/* Now set d[Laplacian]/dn = 0, start/end loop 1 col out,
				use periodicity to set 2nd out col after loop.  */

			for (i = iwo1; set[YLO] && i <= ieo1; i++) {
				G->data[jso2 + i] = G->data[jsi1 + i] + (float)(5.0 * (G->data[jso1 + i] - G->data[js + i]))
					+ (G->data[js + i - 1] - G->data[jso1 + i - 1]) + (G->data[js + i + 1] - G->data[jso1 + i + 1]);
			}
			if (set[XLO] && set[YLO]) G->data[jso2 + iwo2] = G->data[jso2 + iwo2 + G->header->nxp];
			if (set[XHI] && set[YHI]) G->data[jso2 + ieo2] = G->data[jso2 + ieo2 - G->header->nxp];

			/* End of X is periodic, south (bottom) is Natural.  */
			if (set[YLO]) {
				G->header->BC[YLO] = GMT_BC_IS_NATURAL;
				//GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "Set boundary condition for %s edge: %s\n", edge[YLO], kind[G->header->BC[YLO]]);
			}
		}

		/* Done with X is periodic cases.  */

		return (GMT_NOERROR);
	}
}

