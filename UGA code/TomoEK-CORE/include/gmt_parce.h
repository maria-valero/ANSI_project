/*
 * gmt_parce.h
 *
 *  Created on: Feb 22, 2015
 *      Author: nishita
 */

#ifndef GMT_PARCE_H_
#define GMT_PARCE_H_


char ** GMT_Create_Args (void *V_API, int *argc, struct GMT_OPTION *head);
int GMT_Destroy_Args (void *V_API, int argc, char **args[]);
char * GMT_Create_Cmd (void *V_API, struct GMT_OPTION *head);


#endif /* GMT_PARCE_H_ */
