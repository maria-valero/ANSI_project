/*
 * gmt_api.h
 *
 *  Created on: Feb 23, 2015
 *      Author: nishita
 */

#ifndef GMT_API_H_
#define GMT_API_H_

/* Three different i/o status: unused, actively using, or used */
enum GMT_enum_status {
	GMT_IS_UNUSED = 0,	/* We have not yet read from/written to this resource */
	GMT_IS_USING,		/* Means we have started reading from/writing to this file */
	GMT_IS_USED};

void * GMT_Read_Data (void *V_API, unsigned int family, unsigned int method, unsigned int geometry, unsigned int mode, double wesn[], char *input, void *data);
void * GMT_Get_Record (void *V_API, unsigned int mode, int *retval);
int GMT_Destroy_Data (void *V_API, void *object);
void * GMT_Create_Data (void *V_API, unsigned int family, unsigned int geometry, unsigned int mode, uint64_t dim[], double *range, double *inc, unsigned int registration, int pad, void *data);
int GMT_Register_IO (void *V_API, unsigned int family, unsigned int method, unsigned int geometry, unsigned int direction, double wesn[], void *resource);
int GMTAPI_Unregister_IO (struct GMTAPI_CTRL *API, int object_ID, unsigned int direction);
void * GMT_Retrieve_Data (void *V_API, int object_ID);
int GMT_Put_Record (void *V_API, unsigned int mode, void *record);

#endif /* GMT_API_H_ */
