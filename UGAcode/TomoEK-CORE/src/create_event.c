//============================================================================
// Name        : test.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <glob.h>

//#include "../include/Tomography.h"
//#include "createArrivalTimeDiffFile.h"
//#define ARRIVALTIME_FILE "../data/arrival.dat"
#define DAT_FILE_LIST "list.dat"
int TotalEvent;

void createArrivalTime(char *arrivalTimefile) {
	char buffer[300];
	char aux[30];

	FILE * fpRead, *fpWrite;
	char * line = NULL;
	size_t len = 0;
	ssize_t read;
	int eventCounter = 1;
	//int eventFound = 0;
	char arrpath[100] ;//= "/home/nishita/buildGMTUsingMake/data";
	char dstPath[100];
	char path[100];
	if (getcwd(path, sizeof(path)) == NULL)
	{
		printf("Error \n\n\n");
		//return -1;
	}
	strcpy(path,"/home/labm/.core/TomoEK/data/");
	strcpy(arrpath,path);
	strcat(arrpath,arrivalTimefile);
	strcpy(dstPath,path);
	printf("dstPath :%s\n",arrpath);
	fpRead = fopen(arrpath, "r");
	if (fpRead == NULL)
	{
		printf("File %s Not found!!\n",arrivalTimefile);
	}
	while ((read = getline(&line, &len,fpRead)) != -1) {
		sprintf(aux,"Event %d\n",eventCounter);
		if(strcmp(line,aux) == 0)
		{
			if(eventCounter > 1)
				fclose(fpWrite);
			sprintf(buffer,"%sEvent_%d.dat",dstPath,eventCounter);
			fpWrite = fopen(buffer, "w");
			eventCounter++;
			continue;
		}
		else
		{
			fprintf(fpWrite, "%s", line);
		}
	}
	fclose(fpRead);
	fclose(fpWrite);
	if (line)
		free(line);
	TotalEvent = eventCounter - 1;
	//printf("TotalEvent %d\n", TotalEvent);
}


void calculateTiimeDiff(char* stationFile)
{

	double Arr_Vecs[50][2];
	double stationCor[50][3];
	char buffer[300];
	char outbuff[100], listbuff[100];
	FILE *diff, *fpRead, *staionfp, *listfp;
	double stationNum;
	double arrivalTime, lon,lat;
	int i,j,k;
	int stationCounter = 0;
	int totalStaCnt = 0;
	char stapath[100] ;
	char dstPath[100];
	char path[100];
		if (getcwd(path, sizeof(path)) == NULL)
		{
			printf("Error \n\n\n");
			//return -1;
		}
		strcpy(path,"/home/labm/.core/TomoEK/data/");
		strcpy(stapath,path);
		strcat(stapath,stationFile);
		strcpy(dstPath,path);
		//printf("dstPath :%s\n",dstPath);


	staionfp = fopen(stapath,"r");
	if (staionfp == NULL)
		printf("File %s Not found!!\n",stationFile);

	fscanf (staionfp, "%lf %lf %lf", &lon, &lat, &stationNum);
	while (!feof (staionfp))
	{
		//printf ("%lf %lf %lf\n", lon, lat, stationNum);
		stationCor[totalStaCnt][0] = lon;
		stationCor[totalStaCnt][1] = lat;
		stationCor[totalStaCnt][2] = stationNum;
		totalStaCnt++;
		fscanf (staionfp, "%lf %lf %lf", &lon, &lat, &stationNum);

	}
	fclose(staionfp);
	strcpy(listbuff,path);
	strcat(listbuff,DAT_FILE_LIST);
	listfp = fopen(listbuff,"w");
	for(i = 1; i <= TotalEvent;i++ )
	{
		sprintf(buffer,"%sEvent_%d.dat",dstPath,i);
		printf("File name %s\n", buffer);

		fpRead = fopen(buffer, "r");
		if (fpRead == NULL)
		{
			printf("File %s Not found!!\n",buffer);
		}
		fscanf (fpRead, "%lf %lf ", &stationNum, &arrivalTime);
		while (!feof (fpRead))
		{
			Arr_Vecs[stationCounter][0] = stationNum;
			Arr_Vecs[stationCounter][1] = arrivalTime;
			stationCounter++;
			fscanf (fpRead, "%lf %lf", &stationNum, &arrivalTime);
		}
		fclose(fpRead);
		sprintf(outbuff,"%s%.0f_Event_%d.dat",dstPath,i,Arr_Vecs[0][0]);

		fprintf(listfp,"%.0f_Event_%d.dat\n",i,Arr_Vecs[0][0]);

		diff = fopen(outbuff,"w");
		//printf("stationCounter %d \n",totalStationCounter);


		for(j = 0; j < stationCounter ;j++ ) //maximum number of stations
		{
			for(k = 0; k < totalStaCnt ;k++)
			{
				//printf ("%lf %lf\n", Arr_Vecs[j][0], stationCor[k][2]);
				if(Arr_Vecs[j][0] == stationCor[k][2])
				{
					fprintf(diff,"%lf %lf %lf\n", stationCor[k][0],stationCor[k][1],(Arr_Vecs[j][1] - Arr_Vecs[0][1]));
					break;
				}
			}
		}
		stationCounter = 0;
		fclose(diff);
	}
	fclose(listfp);
}
/*

int main(int argc,char *argv[]) {
	 if(argc != 3)
	    {
	      printf("usage:arrivalTime.dat station.dat %s %s",argv[1],argv[2]);
	      return 0;
	    }
	createArrivalTime(argv[1]);
	calculateTiimeDiff(argv[2]);
	return 0;
}


int main(int argc, char **argv) {
  const char *pattern = "/home/nishita/workspace/arrivalDiff/data/*.txt";
  glob_t pglob;

  glob(pattern, GLOB_ERR, NULL, &pglob);

  printf("Found %d matches\n", pglob.gl_pathc);
  printf("First match: %s\n", pglob.gl_pathv[0]);

  globfree(&pglob);


  return 0;
}*/
