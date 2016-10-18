/*
 * gmt_parce.cpp
 *
 *  Created on: Feb 22, 2015
 *      Author: nishita
 */
#include <stddef.h>
#include "memory.h"
#include "gmt_constants.h"
#include "gmt_defaults.h"
#include "surface.h"
#include "gmt_macros.h"
static inline struct GMTAPI_CTRL * gmt_get_api_ptr (struct GMTAPI_CTRL *ptr) {return (ptr);} //nishita

char ** GMT_Create_Args (void *V_API, int *argc, struct GMT_OPTION *head)
{	/* This function creates a character array with the command line options that
	 * correspond to the linked options provided.  It is the inverse of GMT_Create_Options.
	 * The number of array strings is returned via *argc.
	 */

	char **txt = NULL, buffer[GMT_BUFSIZ] = {""};
	unsigned int arg = 0;
	size_t n_alloc = GMT_SMALL_CHUNK;
	struct GMT_OPTION *opt = NULL;
	struct GMT_CTRL *G = NULL;
	struct GMTAPI_CTRL *API = NULL;

	if (V_API == NULL)
		return_null (V_API, GMT_NOT_A_SESSION);	/* GMT_Create_Session has not been called */
		//return(GMT_NOT_A_SESSION); //nishita
	if (head == NULL)
		return_null (V_API, GMT_OPTION_LIST_NULL);	/* No list of options was given */
		//return (GMT_OPTION_LIST_NULL); //nishita
	API = gmt_get_api_ptr ((struct GMTAPI_CTRL *)V_API);	/* Cast void pointer to a GMTAPI_CTRL pointer */

	*argc = 0;	/* Start off with no arguments */

	G = API->GMT;	/* GMT control structure */
	txt = GMT_memory (G, NULL, n_alloc, char *);

	for (opt = head; opt; opt = opt->next) {	/* Loop over all options in the linked list */
		if (!opt->option) continue;			/* Skip all empty options */
		if (opt->option == GMT_OPT_SYNOPSIS)		/* Produce special - option for synopsis */
			sprintf (buffer, "-");
		else if (opt->option == GMT_OPT_INFILE)		/* Option for input filename [or numbers] */
			sprintf (buffer, "%s", opt->arg);
		else if (opt->arg && opt->arg[0])			/* Regular -?arg commandline option with argument for some ? */
			sprintf (buffer, "-%c%s", opt->option, opt->arg);
		else							/* Regular -? commandline argument without argument */
			sprintf (buffer, "-%c", opt->option);

		txt[arg] = GMT_memory (G, NULL, strlen (buffer)+1, char);	/* Get memory for this item */

		/* Copy over the buffer contents */
		strcpy (txt[arg], buffer);
		arg++;	/* One more option added */
		if (arg == n_alloc) {	/* Need more space for our growing list */
			n_alloc += GMT_SMALL_CHUNK;
			txt = GMT_memory (G, txt, n_alloc, char *);
		}
	}
	/* OK, done processing all options */
	if (arg == 0) {	/* Found no options, so delete the list we allocated */
		GMT_free (G, txt);
	}
	else if (arg < n_alloc) {	/* Trim back on the list to fit what we want */
		txt = GMT_memory (G, txt, arg, char *);
	}

	*argc = arg;	/* Pass back the number of items created */
	return (txt);	/* Pass back the char* array to the calling module */
}

int GMT_Destroy_Args (void *V_API, int argc, char **args[])
{	/* Delete all text arguments, perhaps those created by GMT_Create_Args
	 * Note that a pointer to the args[] array is expected so that we can
	 * set it to NULL afterwards. */

	struct GMTAPI_CTRL *API = NULL;
	if (V_API == NULL)
		return_error (V_API, GMT_NOT_A_SESSION);		/* GMT_Create_Session has not been called */
		//return (GMT_NOT_A_SESSION); //nishita
		if (argc == 0 || !args)
			return_error (V_API, GMT_ARGV_LIST_NULL);	/* We were given no args to destroy, so there! */
			//return (GMT_ARGV_LIST_NULL); //nishita
	if (argc < 0) return_error (V_API, GMT_COUNTER_IS_NEGATIVE);		/* We were given a negative count! */
	//return (GMT_COUNTER_IS_NEGATIVE); //nishita
	API = gmt_get_api_ptr (V_API);	/* Cast void pointer to a GMTAPI_CTRL pointer */
	/* Just deallocate the space taken by the list of arguments */
	while (argc--) GMT_free (API->GMT, (*args)[argc]);
	GMT_free (API->GMT, *args);	/* Free the array itself */
	return (GMT_OK);		/* No error encountered */
}

char * GMT_Create_Cmd (void *V_API, struct GMT_OPTION *head)
{	/* This function creates a single character string with the command line options that
	 * correspond to the linked options provided.
	 */

	char *txt = NULL, buffer[GMT_BUFSIZ] = {""};
	bool first = true;
	size_t length = 0, inc, n_alloc = GMT_BUFSIZ;
	struct GMT_OPTION *opt = NULL;
	struct GMT_CTRL *G = NULL;
	struct GMTAPI_CTRL *API = NULL;

	if (V_API == NULL) return_null (V_API, GMT_NOT_A_SESSION);	/* GMT_Create_Session has not been called */
	if (head == NULL) return_null (V_API, GMT_OPTION_LIST_NULL);	/* No list of options was given */
	API = gmt_get_api_ptr (V_API);	/* Cast void pointer to a GMTAPI_CTRL pointer */

	G = API->GMT;	/* GMT control structure */
	txt = GMT_memory (G, NULL, n_alloc, char);

	for (opt = head; opt; opt = opt->next) {	/* Loop over all options in the linked list */
		if (!opt->option) continue;			/* Skip all empty options */
		if (opt->option == GMT_OPT_SYNOPSIS)		/* Produce special - option for synopsis */
			sprintf (buffer, "-");
		else if (opt->option == GMT_OPT_INFILE)		/* Option for input filename [or numbers] */
			sprintf (buffer, "%s", opt->arg);
		else if (opt->arg && opt->arg[0])			/* Regular -?arg commandline option with argument for some ? */
			sprintf (buffer, "-%c%s", opt->option, opt->arg);
		else							/* Regular -? commandline argument without argument */
			sprintf (buffer, "-%c", opt->option);

		inc = strlen (buffer);
		if (!first) inc++;	/* Count the space between args */
		if ((length + inc) >= n_alloc) {	/* Will need more memory */
			n_alloc <<= 1;
			txt = GMT_memory (G, txt, n_alloc, char);
		}
		if (!first) strcat (txt, " ");	/* Add space betwen args */
		strcat (txt, buffer);
		length += inc;
		first = false;
	}
	length++;	/* Need space for trailing \0 */
	/* OK, done processing all options */
	if (length == 1)	/* Found no options, so delete the string we allocated */
		GMT_free (G, txt);
	else if (length < n_alloc)	/* Trim back on the list to fit what we want */
		txt = GMT_memory (G, txt, length, char);

	return (txt);		/* Pass back the results to the calling module */
}




