
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#define DAT_FILE_LIST "list.dat"

double get_dist(double lat1,double lon1,double lat2,double lon2)
{
	double theta,pi,temp;
	double radius=6371;
	pi=4.0*atan(1.0);

	lat1=atan(0.993277*tan(lat1/180*pi))*180/pi;
	lat2=atan(0.993277*tan(lat2/180*pi))*180/pi;

	temp = sin((90 - lat1) / 180 * pi) * cos(lon1 / 180 * pi)
											* sin((90 - lat2) / 180 * pi) * cos(lon2 / 180 * pi)
											+ sin((90 - lat1) / 180 * pi) * sin(lon1 / 180 * pi)
											* sin((90 - lat2) / 180 * pi) * sin(lon2 / 180 * pi)
											+ cos((90 - lat1) / 180 * pi) * cos((90 - lat2) / 180 * pi);
	if (temp > 1) {
		printf("warning cos(theta)>1 and correct to 1!!\n",temp);
		temp=1;
	}
	if (temp < -1) {
		printf("warning cos(theta)<-1 and correct to -1!!\n",temp);
		temp=-1;
	}
	theta=fabs(acos(temp));
	return theta*radius;
}

int createSlowMap(double period, char *pa,  char *datapath) {

	FILE *fin, *fin2, *fin3, *fout, *file1;
	int i,j;
	int npts_x,npts_y;
	char buff1[1000]={""},stalist[1000]={""},path[1000]={""};
	double lat, lon, lat2, lon2,radius, pi;

	int marker_nn,marker_EN[2][2],marker_E,marker_N;
	char staName[100]={""};
	double dist;


	radius=6371.1391285;
	pi=4.0*atan(1.0);
	double dx,dy,x0,y0,x1,y1,temp,temp1,temp2,lat_temp;
	npts_x=176;
	npts_y=126;
	//fprintf(stderr, "Memory check!!\n");
	double tr_t[npts_x][npts_y];
	double dx_km[npts_y],dy_km;
	//fprintf(stderr, "Memory enough!!\n");

	dx=0.2;//degree
	dy=0.2;//degree
	x0=235.0;
	y0=25.0;
	x1=x0+(npts_x-1)*dx;
	y1=y0+(npts_y-1)*dy;
	for (j = 1; j < npts_y - 1; j++) {
		lat_temp=y0+j*dy;
		lat_temp=atan(0.993277*tan(lat_temp/180.0*pi))*180.0/pi;
		dx_km[j]=radius*sin((90.0-lat_temp)/180.0*pi)*dx/180.0*pi;
	}
	dy_km=radius*dy/180.0*pi;

	//if (getcwd(path, sizeof(path)) == NULL)
	//	{
	//		printf("Error \n\n\n");
	//		return -1;		
	//	}
	
	strcpy(stalist,pa);
	strcpy(path,datapath);
	//	strcpy(path,"/home/labm/.core/TomoEK/data/");
	//	strcpy(stalist,path);
	//	strcat(stalist,DAT_FILE_LIST);
	
	file1=fopen(stalist,"r");
	for(;;)
	{
		if(fscanf(file1,"%s",&staName)==EOF)
			break;
		sprintf(buff1,"%s%s.HD_0.0",path,staName);
printf("CREATING SLOWNESS FILE %s\npath: %s\nstaName %s\n\n",buff1,path,staName);
		if((fin=fopen(buff1,"r"))==NULL)
		{
			printf( "%s.HD_FILE not exists!!!\n",staName);
			return 1;
		}
		sprintf(buff1,"%s%s.HD_0.2",path,staName);
		if ((fin2 = fopen(buff1, "r")) == NULL) {
			printf("%s.HD_02_FILE not exists!!!\n",staName);
			return 1;
		}
		sprintf(buff1,"%sslow_azi_%s.HD",path,staName);

		fout=fopen(buff1,"w");
		for (i = 0; i < npts_x; i++) {
			for (j = 0; j < npts_y; j++) {
				tr_t[i][j]=0;
			}
		}
		sprintf(buff1, "%s%s",path, staName);
		if ((fin3 = fopen(buff1, "r")) == NULL) {
			printf("%s FILE not exists!!!",staName);
			return 1;
		}
		fclose(fin3);
		for (;;) {
			if (fscanf(fin, "%lf %lf %lf", &lon, &lat, &temp) == EOF)
				break;
			if (fscanf(fin2, "%lf %lf %lf", &lon2, &lat2, &temp2) == EOF)
				break;
			if (lon != lon2 || lat != lat2) {
				printf("HD and HD_0.2 files not compatiable!!!");
				return 0;
			}
			if (lon > x1 + 0.01 || lon < x0 - 0.01 | lat > y1 + 0.01
					|| lat < y0 - 0.01)
				continue;
			i=(int)((lon-x0)/dx+0.1);
			j=(int)((lat-y0)/dy+0.1);
			if (temp < temp2 - 1 || temp > temp2 + 1 || temp < 2.0 * period) {
				tr_t[i][j]=0;
				continue;
			}
			marker_nn=3;
			marker_EN[0][0]=0;
			marker_EN[0][1]=0;
			marker_EN[1][0]=0;
			marker_EN[1][1]=0;
			fin3=fopen(buff1,"r");
			for (;;) {
				if (fscanf(fin3, "%lf %lf %lf", &lon2, &lat2, &temp2) == EOF) {
					temp=0;
					fclose(fin3);
					break;
				}
				if(lon2-lon<0)
					marker_E=0;
				else
					marker_E=1;
				if(lat2-lat<0)
					marker_N=0;
				else
					marker_N=1;
				if(marker_EN[marker_E][marker_N]!=0)
					continue;
				dist=get_dist(lat,lon,lat2,lon2);
				if (dist < 150) {
					marker_nn--;
					if (marker_nn == 0) {
						fclose(fin3);
						break;
					}
					marker_EN[marker_E][marker_N]++;
				}
			}
			tr_t[i][j]=temp;
			
		}
		fclose(fin);
		fclose(fin2);

		//      double temp1,temp2;
		for (i = 1; i < npts_x - 1; i++) {
			for (j = 1; j < npts_y - 1; j++) {
				temp1=(tr_t[i+1][j]-tr_t[i-1][j])/2.0/dx_km[j];
				temp2=(tr_t[i][j+1]-tr_t[i][j-1])/2.0/dy_km;
				if (temp2 == 0) {
					temp2=0.00001;
				}
				temp=sqrt(temp1*temp1+temp2*temp2);

				if (temp > 0.5 || temp < 0.2 || tr_t[i + 1][j] == 0
						|| tr_t[i - 1][j] == 0 || tr_t[i][j + 1] == 0
						|| tr_t[i][j - 1] == 0) {
					fprintf(fout,"%lf %lf 0 999\n",x0+i*dx,y0+j*dy);
				} else {
					if(temp1>0&&temp2>0)
						fprintf(fout, "%lf %lf %lf %lf\n", x0 + i * dx, y0 + j * dy,
								temp, atan(temp2 / temp1) / pi * 180);
					if(temp1>0&&temp2<=0)
						fprintf(fout, "%lf %lf %lf %lf\n", x0 + i * dx, y0 + j * dy,
								temp, 360 + atan(temp2 / temp1) / pi * 180);
					if(temp1<=0)
						fprintf(fout, "%lf %lf %lf %lf\n", x0 + i * dx, y0 + j * dy,
								temp, 180 + atan(temp2 / temp1) / pi * 180);
				}
			}
		}
		fclose(fout);
	}
	fclose(file1);
	return 0;
}
