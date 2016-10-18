#define MAIN
#include "./mysac64.h"
#include "./seis_header.h"

#include        <sys/types.h>
#include        <dirent.h>
#include        <fcntl.h>
#include        <sys/stat.h>
#include        <pthread.h>
#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
using namespace std;


int main(int na, char *arg[])
{
        int nmax = 1800000;
        FILE *fsac,*fp;

        int i,j,k;
        char fname[20];
        char station_name[8],knetwk[8];
        int station_number;
        char temp[8],st[8];
        char rw_filename[50];

        int hours,npts;
        SAC_HD *SHD = &SAC_HEADER;

        struct seis_header seis;
        struct seis_header *SEIS = &seis;
        float delta;
        int year,day,hour,min,sec;
        int year2,day2,hour2;
        int month,day3;
        
        int skip;

        if((fsac=fopen(arg[1], "rb"))==NULL){
                printf("Can't read or find the SAC file!\n");
                exit(0);
        }
        else{
                fread(SHD,sizeof(SAC_HD),1,fsac);
                npts = SHD->npts;
                hours=int((SHD->npts)/nmax);
                float sig[hours][nmax];
                float sigtem[nmax];
                printf("npts is %d\n",npts);
                printf("hours is %d\n",hours);

                delta = SHD->delta;
                for(i=0;i<5;i++) station_name[i] = SHD->kstnm[i];
                for(i=0;i<3;i++){
                        if(SHD->knetwk[i] == ' '){
                                knetwk[i] = '\0';
                                break;
                        }
                        else{
                                knetwk[i]=SHD->knetwk[i];
                        }
                }
                printf("station name is %s.\n",station_name);
                printf("delta is %lf.\n",delta);
                printf("network is %s.\n",knetwk);
                
                for(i=0;i<5;i++) temp[i] = st[i+2];
                station_number = atof(temp);
                year = SHD->nzyear;
                day = SHD->nzjday;
                hour = SHD->nzhour;
                min = SHD->nzmin;
                sec = SHD->nzsec;

                skip = 0;
                hour2 = hour;
                day2 = day;
                year2 = year;
                if(min != 0 || sec != 0){
                        hour2 = hour + 1;
                        if(hour2 >23){
                                hour2 = 0;
                                day2 = day + 1;
                        }
                        skip = (3600 - ((min*60)+sec))/(SHD->delta);
                }
                fread(sigtem,sizeof(float),skip,fsac);
                for(i=0;i<hours;i++){
			printf("%d\n",hours);
			fread(sig[i],sizeof(float),nmax,fsac);
			//for(j=0;j<3;i++){
			//	printf("%d %lf\n",j,sig[i][j]);
			//}
		}
		fclose(fsac);
                fsac = NULL;
                SHD->iftype = (int)ITIME;
                SHD->leven = (int)TRUE;
                SHD->lovrok = (int)TRUE;
		SHD->internal4 = 6L;
                SHD->npts = nmax;
		
		for(k=0;k<hours;k++){
                        if(day2<32) {month = 1;day3 = day2;}
                        else if(day2<60) {month = 2;day3 = day2-31;}
                        else if(day2<91) {month =3; day3 = day2 - 59;}
                        else if(day2<121) {month =4;day3 = day2 - 90;}
                        else if(day2<152) {month =5;day3 = day2 -120;}
                        else if(day2<182) {month =6;day3 = day2 - 151;}
                        else if(day2<213) {month =7;day3 = day2 - 181;}
                        else if(day2<244) {month =8;day3 = day2 - 212;}
                        else if(day2<274) {month =9;day3 = day2 - 243;}
                        else if(day2<305) {month =10;day3 = day2 - 274;}
                        else if(day2<335) {month =11;day3 = day2 - 304;}
                        else if(day2<366) {month =12;day3 = day2 - 334;}
                        sprintf(rw_filename, "%04d%02d%02d%02d0000.%s-0000.seis", year2,month,day3,hour2,station_name);
                        fsac=fopen(rw_filename,"wb");
	
			SHD->nzjday = day2;
                        SHD->nzhour = hour2;
                        SHD->nzmin = 0;
                        SHD->nzsec = 0;
                        SHD->depmin=sig[k][0];
                        SHD->depmax=sig[k][0];
                        printf("\n END of hour %d\n",hour2);
			
			for(j=0;j<nmax;j++){
                                if ( SHD->depmin > sig[k][j] ) SHD->depmin = sig[k][j];
                                if ( SHD->depmax < sig[k][j] ) SHD->depmax = sig[k][j];
                        }
                        SEIS->nt = SHD->npts;
                        SEIS->id = 0;
                        SEIS->rline = station_number;
                        SEIS->rpt = 0;
                        for(i=0;i<5;i++) SEIS->name[i] = station_name[i];
                        SEIS->x = SHD->stlo;
                        SEIS->y = SHD->stla;
                        SEIS->z = SHD->stel;
                        SEIS->t0 = 0.0;
                        SEIS->dt = delta;
			fwrite(SEIS,SEIS_HEADER_SZ,1,fsac);
                        fwrite(sig[k],sizeof(float),nmax,fsac);
                        fclose(fsac);
                        fsac=NULL;

                        hour2 = hour2 + 1;
                        if(hour2 >23){
                                hour2 = 0;
                                day2 = day2 + 1;
                        }

                }

        }
        return 0;
}
