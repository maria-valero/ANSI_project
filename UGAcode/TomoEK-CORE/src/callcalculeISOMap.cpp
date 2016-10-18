#include <iostream>
#include <string.h>

using namespace std;
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#define DAT_FILE_LIST "list.dat"

//extern struct ISO_AXIS CalculateISOMap(double min, double max, int N_bin, char* fileName) ;
//extern int CalculateISOMap(double min, double max, int N_bin, char* fileName) ;
int CalculateISOMap(double min, double max, int N_bin, char* fileName);

int main()
{
    cout << "Hello world!" << endl;
    char filename[10];
    strcpy(filename,"map");
    CalculateISOMap(0, 360, 18, filename) ;
    return 0;
}

int CalculateISOMap(double min, double max, int N_bin, char* fileName) {
	FILE *fin,*file1,*fout, *file_iso, *file_ani;
	int i, j, k;
	int npts_x, npts_y;
	char buff1[100], name_iso[100], name_ani[100], name_ani_n[100];
	double t_lat, t_lon, radius, pi, sta1_lon, sta1_lat;
	int t_i, t_j, nsta;
	int ii, jj, kk, kkk, min_n;
	double d_bin;
	//  cout<<min<<" "<<max<<" "<<N_bin<<endl;
	double hist[N_bin];
	double slow_sum1[N_bin];
	double slow_un[N_bin];

	//struct ISO_AXIS coordinats;

	d_bin = (max - min) / N_bin;



	radius = 6371.1391285;
	pi = 4.0 * atan(1.0);
	double dx, dy, x0, y0, x1, y1, temp, lat_temp, trash1, trash2;
	double lon,lat,temp2;

	npts_x = 176;
	npts_y = 126;

	//coordinats.x_coor = npts_x;
	//coordinats.y_coor = npts_y;

	fprintf(stderr, "Memory check!!\n");
	static double slow[176][126][2000];
	static double azi[176][126][2000];
	double weight[2000];
	double weight_sum;
	int n[npts_x][npts_y];
	double slow_sum[npts_x][npts_y],slow_std[npts_x][npts_y];
	// double dx_km[npts_y],dy_km;

	fprintf(stderr,"Memory enough!!\n");

	dx=0.2;//degree
	dy=0.2;//degree
	x0=235.0;
	y0=25.0;
	x1=x0+(npts_x-1.0)*dx;
	y1=y0+(npts_y-1.0)*dy;
	//for(j=1;j<npts_y-1;j++)
	//  {
	//  lat_temp=y0+j*dy;
	//  lat_temp=atan(0.993277*tan(lat_temp/180*pi))*180/pi;
	//  dx_km[j]=radius*sin((90-lat_temp)/180*pi)*dx/180*pi;
	//}
	//dy_km=radius*dy/180*pi;
	for(i=0;i<npts_x;i++)
	{
		for(j=0;j<npts_y;j++)
		{
			slow_sum[i][j] = 0.0;
			n[i][j] = 0.0;
		}
	}
	char event_name[100];
	char stalist[100],path[100];
	if (getcwd(path, sizeof(path)) == NULL)
	{
		printf("Error \n\n\n");
		return -1.0;
	}

	strcpy(path,"/home/labm/.core/TomoEK/data/");
	strcpy(stalist,path);
	strcat(stalist,DAT_FILE_LIST);
	file1=fopen(stalist,"r");
	nsta=0.0;
    cout<<stalist<<endl;

	for(;;)
	{

		if(fscanf(file1,"%s",&event_name)==EOF)   //I change here &event_name for event_name
		 	break;

		nsta++;
		sprintf(buff1,"%sslow_azi_%s.HD",path,event_name);
printf("NOMBRE:..... %s",buff1);
		if ((fin = fopen(buff1, "r")) == NULL) {
			printf("NODESLOWNESSFILE no exist!!\n");
			return 1;
		}
		for (;;)
		{
			if (fscanf(fin, "%lf %lf %lf %lf", &lon, &lat, &temp, &temp2) == EOF)
				break;
			if (lon > x1 + 0.01 || lon < x0 - 0.01 | lat > y1 + 0.01
					|| lat < y0 - 0.01)
				continue;
			i = (int)((lon - x0) / dx + 0.1);
			j = (int)((lat - y0) / dy + 0.1);
			if (temp < 0.5 && temp > 0.2)
			{
				slow[i][j][n[i][j]] = temp;
				azi[i][j][n[i][j]] = temp2;
				n[i][j]++;
				//	     slow_sum[i][j]+=temp;
			}
		}
		fclose(fin);
	}
	fclose(file1);
	sprintf(name_iso, "%s%s.iso",path,fileName);
	//sprintf(name_iso,"%sslow_azi_%s.txt.HD",path,event_name);
	file_iso = fopen(name_iso, "w");
	nsta = 100.0;
	double w2;
	double temp_slow_sum;
	int temp_n;
	for (i = 0; i < npts_x; i++)
	{
		for (j = 0; j < npts_y; j++)
		{
			if (n[i][j] < 0.5 * nsta)
			{
				fprintf(file_iso, "%lf %lf 0 9999 %d\n", x0 + i * dx,
						y0 + j * dy, n[i][j]);
//printf("........................%d...............................\n",n[i][j]);
				continue;
			}
			w2 = 0.0;
			weight_sum = 0.0;
			for (k = 0; k < n[i][j]; k++)
			{
				weight[k] = 0.0;
				for (kk = 0; kk < n[i][j]; kk++)
				{
					if (fabs(azi[i][j][kk] - azi[i][j][k]) < 25.0)
						weight[k]++;
				}
				weight[k] = 1.0 / weight[k];
				weight_sum += weight[k];
			}
			for (k = 0; k < n[i][j]; k++)
			{
				weight[k] = weight[k] / weight_sum;
				slow_sum[i][j] += weight[k] * slow[i][j][k];
				w2 += weight[k] * weight[k];
			}
			//	  if(n[i][j]>=0.5*nsta)
				//{
			//
			//  slow_sum[i][j]=slow_sum[i][j]/n[i][j];
			temp_slow_sum = slow_sum[i][j];
			temp = 0.0;
			for (k = 0; k < n[i][j]; k++)
			{
				temp += weight[k] * (slow[i][j][k] - slow_sum[i][j])
																		* (slow[i][j][k] - slow_sum[i][j]);
			}
			slow_std[i][j] = sqrt(temp / (1.0 - w2));
			w2 = 0.0;
			weight_sum = 0.0;
			slow_sum[i][j] = 0.0;
			temp_n = 0.0;
			for (k = 0; k < n[i][j]; k++)
			{
				if (fabs(slow[i][j][k] - temp_slow_sum) > 2.0 * slow_std[i][j])
					continue;
				weight_sum += weight[k];
				temp_n++;
			}
			//	   if(i==marker_i&&j==marker_j)
			// fprintf(stderr,"%d %d\n",n[i][j],temp_n);
			for (k = 0; k < n[i][j]; k++)
			{
				if (fabs(slow[i][j][k] - temp_slow_sum) > 2.0 * slow_std[i][j])
					continue;
				weight[k] = weight[k] / weight_sum;
				slow_sum[i][j] += weight[k] * slow[i][j][k];
				w2 += weight[k] * weight[k];
			}
			temp = 0.0;
			for (k = 0; k < n[i][j]; k++)
			{
				if (fabs(slow[i][j][k] - temp_slow_sum) > 2.0 * slow_std[i][j])
					continue;
				temp += weight[k] * (slow[i][j][k] - slow_sum[i][j])
																		* (slow[i][j][k] - slow_sum[i][j]);
			}
			slow_std[i][j] = sqrt(temp / (1.0 - w2));
			temp = slow_std[i][j] * sqrt(w2) / slow_sum[i][j] / slow_sum[i][j];
			//	  temp=sqrt(temp/(n[i][j]-1)/n[i][j])/slow_sum[i][j]/slow_sum[i][j];
			// cout<<x0+i*dx<<" "<<y0+j*dy<<" "<<1/slow_sum[i][j]<<" "<<temp<<" "<<n[i][j]<<endl;

fprintf(file_iso, "%lf %lf %lf %lf %d\n", x0 + i * dx, y0 + j * dy,
					1.0 / slow_sum[i][j], temp, temp_n);
			//}
			//else
			//{
			//  //    cout<<x0+i*dx<<" "<<y0+j*dy<<" 0 "<<" 9999 "<<n[i][j]<<endl;
			//  fprintf(file_iso,"%lf %lf 0 9999 %d\n",x0+i*dx,y0+j*dy,n[i][j]);
			//}
		}
	}
	fclose(file_iso);

	sprintf(name_ani, "%s%s.ani", path,fileName);
	sprintf(name_ani_n, "%s%s_ani_n", path,fileName);
	file_ani = fopen(name_ani, "w");
	fout = fopen(name_ani_n, "w");

	//@maria: Printing to check
	printf("\n\nFILE NAME_ANI: %s",name_ani);
	printf("\n\nFILE NAME_ANI_N (OUTPUT): %s",name_ani_n);

	double test_lon, test_lat;
	//test_lon=240.6;
	//test_lat=36.6;

	for (i = 0; i < npts_x; i++)
	{
		for (j = 0; j < npts_y; j++)
		{
			//  i=int((test_lon-x0)/dx+0.1);
			//j=int((test_lat-y0)/dy+0.1);
			//cout<<i<<" "<<j<<endl;
			//fprintf(file_ani,"\n>\n%lf %lf\n",x0+i*dx,y0+j*dy);
			for (k = 0; k < N_bin; k++)
			{
				hist[k] = 0.0;
				slow_sum1[k] = 0.0;
				slow_un[k] = 0.0;
			}
			if (i - 3.0 < 0 || i + 3.0 >= npts_x || j - 3.0 < 0 || j + 3.0 >= npts_y)
				continue;
			kkk = 0;
			//	  min_n=999999999;
			for (ii = i - 3; ii <= i + 3; ii += 3)
			{
				for (jj = j - 3; jj <= j + 3; jj += 3)
				{
					if (n[ii][jj] < nsta * 0.5)
						continue;
					kkk += n[ii][jj];
					//  if(n[ii][jj]<min_n)
					//min_n=n[ii][jj];
				}
			}
			fprintf(fout, "%lf %lf %d\n", x0 + i * dx, y0 + j * dy, kkk);
			if (kkk < 9.0 * nsta * 0.5 || n[i][j] < nsta * 0.5)
				continue;
			for (ii = i - 3; ii <= i + 3; ii += 3)
			{
				for (jj = j - 3; jj <= j + 3; jj += 3)
				{
					if (n[ii][jj] < nsta * 0.5)
						continue;
					for (k = 0; k < n[ii][jj]; k++)
					{
						if (azi[ii][jj][k] > max || azi[ii][jj][k] < min)
						{
							fprintf(stderr, "out of range!!");
							return 1;
						}
						hist[(int)((azi[ii][jj][k] - min) / d_bin)]++;
						slow_sum1[(int)((azi[ii][jj][k] - min) / d_bin)] +=
								slow[ii][jj][k] - slow_sum[ii][jj];
					}
				}
			}
			kk = 0.0;

			for (k = 0; k < N_bin; k++)
			{
				if (hist[k] >= 10.0)
				{
					kk++;
				}
			}
			fprintf(file_ani, "%lf %lf %d\n", x0 + i * dx, y0 + j * dy, kk);
			for (k = 0; k < N_bin; k++)
			{
				if (hist[k] >= 10.0)
				{
					slow_sum1[k]=slow_sum1[k]/hist[k];
					//		  slow_un[k]=slow_un[k]-slow_sum1[k]*slow_sum1[k]*hist[k];
					//slow_un[k]=sqrt(slow_un[k]/(hist[k]-1)/hist[k])/(slow_sum[i][j]+slow_sum1[k])/(slow_sum[i][j]+slow_sum1[k]); //uncertainty of vel not slow
					slow_un[k]=slow_std[i][j]/sqrt((double)(hist[k]));
					slow_un[k]=slow_un[k]/(slow_sum[i][j]+slow_sum1[k])/(slow_sum[i][j]+slow_sum1[k]);//uncertainty of vel not slow
					//slow_un[k]=slow_std[i][j];
					//cout<<min+(0.5+k)*d_bin<<" "<<1/slow_sum1[k]<<" "<<slow_un[k]<<endl;
					fprintf(file_ani,"%lf %lf %lf\n",min+(0.5+k)*d_bin,1.0/(slow_sum[i][j]+slow_sum1[k]),slow_un[k]);
				}
			}
			//	  fprintf(stderr,"done!!\n");
			//return 0;
		}
	}
	fclose(file_ani);
	fclose(fout);

	sprintf(buff1, "%s%s.tmo", path,fileName);
	fout = fopen(buff1, "w");
	sprintf(name_iso, "%s%s.iso",path,fileName);

	//@maria: Printing to check
	printf("\n\nFINAL READING FROM: %s",name_iso);

	//sprintf(name_iso,"%sslow_azi_%s.txt.HD",path,event_name);
	printf("\nNombre del archivo: %s\n",buff1);
	file_iso = fopen(name_iso, "r");
	fprintf(fout, "0.000000\n0.000000\n0.000000\n0.000000\n0.000000\n0.000000\n0.000000\n0.000000\n0.000000\n0.000000\n0.000000\n0.000000\n%d\n%d\n1\n",npts_x, npts_y);
	for(;;)
		{
	if (fscanf(file_iso, "%lf %lf %lf %lf %d", &lon, &lat, &temp, &temp2, &nsta) == EOF)
				break;
	fprintf(fout, "%lf\n", temp);
		}
	fclose(fout);
	fclose(file_iso);
	return 0;
}
