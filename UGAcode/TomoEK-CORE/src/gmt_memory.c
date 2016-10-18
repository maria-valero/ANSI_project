/*
 * gmt_memory.cpp
 *
 *  Created on: Feb 24, 2015
 *      Author: nishita
 */
#include <stddef.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include "memory.h"
#include "gmt_type.h"
void GMT_free_func (struct GMT_CTRL *GMT,void *addr, bool align, const char *where)
{


//#ifdef MEMDEBUG
	//if (GMT->hidden.mem_keeper->active) {
		//bool is_safe_to_free = gmt_memtrack_sub (GMT, where, addr);
		//if (is_safe_to_free == false)
			//return; /* Address addr was not allocated by GMT_memory_func before */
	//}
//#endif

	if (align) {	/* Must free aligned memory */
#ifdef HAVE_FFTW3F
		fftwf_free (addr);
#else
		free (addr);
#endif
	}
	else
		free (addr);
}

	void *GMT_memory_func (struct GMT_CTRL *GMT, void *prev_addr, size_t nelem, size_t size, bool align, const char *where)
	{
		/* Multi-functional memory allocation subroutine.
		   If prev_addr is NULL, allocate new memory of nelem elements of size bytes.
			Ignore when nelem == 0.
		   If prev_addr exists, reallocate the memory to a larger or smaller chunk of nelem elements of size bytes.
			When nelem = 0, free the memory.
		   If align is true we seek to get aligned memory.
		*/

		void *tmp = NULL;

		if (nelem == SIZE_MAX) {	/* Probably 32-bit overflow */
		//	GMT_report_func (GMT, GMT_MSG_NORMAL, where, "Error: Requesting SIZE_MAX number of items (%" PRIuS ") - exceeding 32-bit counting?\n", nelem);
//#ifdef DEBUG
//			GMT_report_func (GMT, GMT_MSG_NORMAL, where, "GMT_memory called\n");
//#endif
			GMT_exit (GMT, EXIT_FAILURE); return NULL;
		}


		if (prev_addr) {
			if (nelem == 0) { /* Take care of n == 0 */
				GMT_free (GMT, prev_addr);
				return (NULL);
			}
			if (align) {
//#ifdef HAVE_FFTW3F
				tmp = NULL; /* currently defunct */


			}
			else
				tmp = realloc ( prev_addr, nelem * size);
			if (tmp == NULL)
				//die_if_memfail (GMT, nelem, size, where);
				printf("die fail\n");
		}
		else {
			if (nelem == 0) return (NULL); /* Take care of n == 0 */
			if (align) {
//#ifdef HAVE_FFTW3F
//				tmp = fftwf_malloc (nelem * size);
//#elif defined(WIN32) || defined(USE_MEM_ALIGNED)
	//			tmp = _aligned_malloc (nelem * size, 16U);
//#elif defined(HAVE_POSIX_MEMALIGN)
	//			(void)posix_memalign (&tmp, 16U, nelem * size);
//#elif defined(HAVE_MEMALIGN)
				tmp = memalign (16U, nelem * size);

				if (tmp != NULL)
					tmp = memset (tmp, 0, nelem * size);
			}
			else
				tmp = calloc (nelem, size);
			//if (tmp == NULL)
				//die_if_memfail (GMT, nelem, size, where);
		}

//#ifdef MEMDEBUG
	//	if (GMT->hidden.mem_keeper->active)
		//	gmt_memtrack_add (GMT, where, tmp, prev_addr, nelem * size);
//#endif
		return (tmp);
	}




void GMT_free_tmp_arrays (struct GMT_CTRL *GMT)
{
	/* Free temporary coordinate memory used by this session */
	size_t col;

	//if (GMT->hidden.mem_cols) GMT_Report (GMT->parent, GMT_MSG_DEBUG, "GMT memory: Free %" PRIuS " temporary column arrays, each of length : %" PRIuS "\n", GMT->hidden.mem_cols, GMT->hidden.mem_rows);
	for (col = 0; col < GMT->hidden.mem_cols; col++) {	/* For each column, free an array */
		if (GMT->hidden.mem_coord[col]) GMT_free (GMT, GMT->hidden.mem_coord[col]);
	}
	if (GMT->hidden.mem_coord) GMT_free (GMT, GMT->hidden.mem_coord);
	GMT->hidden.mem_rows = GMT->hidden.mem_cols = 0;
}
