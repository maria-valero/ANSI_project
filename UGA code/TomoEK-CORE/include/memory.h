#ifndef MEMORY_H_
#define MEMORY_H_
#include <stdbool.h>
enum GMT_enum_mem_alloc {	/* Initial memory for 2 double columns is 32 Mb */
	GMT_INITIAL_MEM_COL_ALLOC	= 2U,
	GMT_INITIAL_MEM_ROW_ALLOC	= 2097152U	/* 2^21 */	
};

//void GMT_free_func (struct GMT_CTRL *GMT,void *addr, bool align, const char *where);
//void *GMT_memory_func (struct GMT_CTRL *GMT, void *prev_addr, size_t nelem, size_t size, bool align, const char *where);
/* Set __func__ identifier */
#ifndef HAVE___FUNC__
#	ifdef HAVE___FUNCTION__
#		define __func__ __FUNCTION__
#	else
#		define __func__ "<unknown>"
#	endif
#endif

/* Convenience macro for GMT_memory_func */
#if defined (DEBUG) || defined (MEMDEBUG)
#define GMT_memory(C,ptr,n,type) GMT_memory_func(C,ptr,n,sizeof(type),false,__SOURCE_LINE_FUNC)
#define GMT_memory_aligned(C,ptr,n,type) GMT_memory_func(C,ptr,n,sizeof(type),true,__SOURCE_LINE_FUNC)
#else
#define GMT_memory(C, ptr,n,type) GMT_memory_func(C, ptr,n,sizeof(type),false,__func__)
#define GMT_memory_aligned(C,ptr,n,type) GMT_memory_func(C,ptr,n,sizeof(type),true,__func__)
#endif


/* Convenience macro for GMT_free_func */
#if defined (DEBUG) || defined (MEMDEBUG)
#define GMT_free(C,ptr) (GMT_free_func(C,ptr,false,__SOURCE_LINE_FUNC),(ptr)=NULL)
#define GMT_free_aligned(C,ptr) (GMT_free_func(C,ptr,true,__SOURCE_LINE_FUNC),(ptr)=NULL)
#else
#define GMT_free(C,ptr) (GMT_free_func(C,ptr,false,__func__),(ptr)=NULL)
#define GMT_free_aligned(C,ptr) (GMT_free_func(C,ptr,true,__func__),(ptr)=NULL)
#endif
#endif
