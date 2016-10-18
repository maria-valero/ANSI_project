/* Perform correlations on LB data
 * R. Clayton, Caltech, 2011/08/31
 */        

/* Parameters: */
float	tdagcwind	= 1.0;	/*(sec)*/
float	fdagcwind	= 1.0;	/*(hz)*/
//float	fmax		=60.0;	/*(Hz)*/
int	lcorr		=200;	/* samples */
float	dtcorr		=0.01;	/* sec */
#define MAIN
#include	<stdio.h>
#include	<sys/types.h>
#include	<dirent.h>
#include	<stdlib.h>
#include	<unistd.h>
#include	<fcntl.h>
#include	<sys/stat.h>
#include	<pthread.h>
#include	"seis_header.h"
#include        "mysac64.h"
#define numsts 2300	        // This was set as 1025, changed to 1900 Oct 19,2015 Oner
#define numcents 2300          // This was set as 1025, changed to 1900 Oct 19, 2015 Oner
#define LAG 4000               // This was set as 4000, changed to  2000 December 2, 2015 Oner

//int c_rline,c_rpt, citr;
int c_b,c_e;
//static int NSMAX = 5306;
//static int NSMAX2 = 501;
int	ns	=1800; //number of center stations
int	csize	= 1800000;        // This was 900000. I changed it to 1800001   !!!!!
int	ftsize	=1048576;
int	lwind	= 250;
int     lag     = LAG;
int     nwin    = 5000;
char   g_line[numsts][20],g_pt[numsts][20];
int ng;
int ngns[numsts];
int ns_max=-1;
float cor_stack[numcents][numsts][2*LAG+1];
int nstack[numcents][numsts];
float gx[numsts],gy[numsts];

struct complex
   {
   	float re, im;
   };
struct trace
{
  int	rline;
  int	rpt;
  int	index;
  int	nt;
  float	*data;
  //float data[4000000];
  float     x,y;
  float t0;
};

struct trace tr[numsts];
struct trace tr_low[0];
struct trace tr_down[0];


int	nthreads= 1; //number of thread
struct arginfo
   {
   	int	itr;
	int	ntr;
	int	lwind;
	int	id;
   };
struct arginfo arg[16];

struct arginfo2
   {
   	int	itr;
	int	lfreq;
	int	id;
   };
struct arginfo2 arg2[16];

/*c/////////////////////////////////////////////////////////////////////////*/
/*--------------------------------------------------------------------------*/
        void write_sac (char *fname, float *sig, SAC_HD *SHD)
/*----------------------------------------------------------------------------
----------------------------------------------------------------------------*/
{
 FILE *fsac;
 int i;
/*..........................................................................*/
        fsac = fopen(fname, "wb");

        if ( !SHD ) SHD = &SAC_HEADER;


        SHD->iftype = (int)ITIME;
        SHD->leven = (int)TRUE;

        SHD->lovrok = (int)TRUE;
        SHD->internal4 = 6L;


  /*+++++++++++++++++++++++++++++++++++++++++*/
     SHD->depmin = sig[0];
     SHD->depmax = sig[0];
 
   for ( i = 0; i < SHD->npts ; i++ ) {
    if ( SHD->depmin > sig[i] ) SHD->depmin = sig[i];
    if ( SHD->depmax < sig[i] ) SHD->depmax = sig[i];
   }
   fwrite(SHD,sizeof(SAC_HD),1,fsac);

   fwrite(sig,sizeof(float),(int)(SHD->npts),fsac);


   fclose (fsac);
}

load_data(char *dir, char *dir2)
   {
	DIR *dirp, *opendir();
	struct dirent *dp, *readdir();
	FILE *fd, *fopen();
	struct seis_header seis;
	char fullname[128];
	int i,ii;

	if( (dirp= opendir(dir2)) == NULL )
	   {
	   	fprintf(stderr,"cannot open directory= %s\n",dir);
		exit(-1);
	   }
	fprintf(stdout,"%s\n",dir2);

	ns= 0;
	for(i=0;i<ng;i++)
	  {
	    ngns[i]=-1;
	    sprintf(fullname,"%s/%s-%s/%s.%s-%s.seis",dir2,g_line[i],g_pt[i],dir,g_line[i],g_pt[i]);
	    if( (fd= fopen(fullname,"r")) == NULL )
		   {
		     fprintf(stderr,"cannot open file= %s\n",fullname);
			continue;
		   }
	    fprintf(stderr,"read %s\n",fullname);
		fread(&seis,1,SEIS_HEADER_SZ,fd);
		if(seis.nt < csize) continue;

		fread(tr[ns].data,csize,4,fd);
		fclose(fd);
		tr[ns].rline= seis.rline;
		tr[ns].rpt= seis.rpt;
		tr[ns].nt= seis.nt;
		tr[ns].x=seis.x;
		tr[ns].y=seis.y;
		tr[ns].t0=seis.t0;
		tr[ns].index=i;
		for(ii=tr[ns].nt;ii<ftsize;ii++)      // use ii<4*ftsize  December 19, 2015 Oner
		  tr[ns].data[ii]=0;
		ngns[i]=ns;
		ns++;
		//		fprintf(stderr,"ns= %4d\n",ns);
		if(ns > numsts) break;
		
	   }
	closedir(dirp);
	return(ns);
   }

SAC_HD shd_cor2;
main(int ac, char **av)
   {
     shd_cor2=sac_null;

     FILE *ff;
     int i, j,k,temp1,temp2,ierr,lfreq;
     float val,x,y,cx,cy,cor1[2*LAG+1],dis_cri;
     char dir[100],file_name[100],buff[300],buff1[300];
       char dirtime[100];
     int prepare(), correlate();

     if(ac != 7)
       {
	 fprintf(stderr,"Usage: correlate time_list c_b c_e sta_line_pt_x_y.lst dis_cri seis_dir\n");
	 exit(-1);
       }
     //int c_b,c_e;
     c_b=atoi(av[2]);
     c_e=atoi(av[3]);
     if(c_e-c_b>numcents-1)
       {
	 fprintf(stderr,"too many center station!!\n");
	 return 0;
       }
    
  for(i=c_b;i<=c_e;i++)
    {
      //      fprintf(stderr,"%d ",i);
      for(j=0;j<numsts;j++)
	{
	  for(k=0;k<2*LAG+1;k++)
	       {
		 cor_stack[i-c_b][j][k]=0;
		 nstack[i-c_b][j]=0;
	       }
	   }
       }
	dis_cri=atof(av[5]);

	if((ff=fopen(av[4],"r"))==NULL)
	  {
	    fprintf(stderr,"problem open %s!!\n",av[4]);
	    return 0;
	  }
	ng=0;
	for(;;)
	  {
	    if(fscanf(ff,"%s %s %f %f",g_line[ng],g_pt[ng],&x,&y)==EOF)
	      {
		break;
	      }
	    gx[ng]=x;
	    gy[ng]=y;
	    ng++;
	  }
	fclose(ff);
	ftsize= nextpower2(csize);
	for(i=0;i<numsts;i++)
	  {
	    if( (tr[i].data = (float *)(malloc(4*ftsize))) == NULL )
	      {
		fprintf(stderr,"cannot alloc memory ns=%d\n",i);
		exit(-1);
	      }
	  }
       if( (tr_low[0].data = (float *)(malloc(4*ftsize))) == NULL )
	      {
              fprintf(stderr,"cannot alloc memory tr_low \n");
              exit(-1);
          }
       if( (tr_down[0].data = (float *)(malloc(4*ftsize))) == NULL )
	      {
              fprintf(stderr,"cannot alloc memory tr_low \n");
              exit(-1);
          }
       ff=fopen(av[1],"r");
       fscanf(ff,"%s",dirtime);
       fprintf(stderr,"timelist=%s\n",dirtime);
	for(i=c_b;i<c_e;i++)
        {
	    sprintf(buff,"if [ ! -d %s_%s_%s ]; then mkdir %s_%s_%s; fi",g_line[i],g_pt[i],dirtime,g_line[i],g_pt[i],dirtime);
	    system(buff);
	  }
       fclose(ff);
	if((ff=fopen(av[1],"r"))==NULL)
	  {
	    fprintf(stderr,"problem open %s!!\n",av[1]);
	    return 0;
	  }
	for(;;)
	  {
	    if(fscanf(ff,"%s",dir)==EOF)
	      {
		break;
	      }

	    fprintf(stdout,"load data %s\n",dir);
	    ns= load_data(dir,av[6]);
	    fprintf(stdout,"ns= %d\n",ns);
	    fprintf(stdout,"prepare data\n");
          
          for(i=0;i<csize;i++)
          {
              if(i%5==0)
              {
              //fprintf(stderr,"Frequency = %d,Original point= %f\n",i,tr[0].data[i]);
              }
          }
              fprintf(stderr,"---------------------------------11111--------------------------------\n");
	    ierr=launch(prepare,nthreads);
	    
          for(i=0;i<ftsize/2;i++)
          {
             //fprintf(stderr,"Frequency = %d,Original Real1 = %f Image1= %f\n",i,tr[0].data[i*2],tr[0].data[i*2+1]);
          }
          for(i=0;i<ftsize/2;i++)
          {
             //fprintf(stderr,"Frequency = %d,low_pass: Real2 = %f Image2= %f\n",i,tr_low[0].data[i*2],tr_low[0].data[i*2+1]);
          }
           fprintf(stderr,"------------------------------------22222-----------------------------\n");
	    if(ierr==-1)
	      continue;
	    fprintf(stdout,"correlate data ns=%d\n",ns);
	    launch2(correlate,nthreads);
          	  }
	fclose(ff);

	int nt2,nt2_ori;
	float ori_delta=0.002;    // Added by Oner Sufri, April 27 2016.

//	printf("ftsize depends on csize, so csize is %d\n",csize); /* Added by Oner Sufri December 2, 2015 */
//	printf("ftsize in the main is %d\n",ftsize);     /* Added by Oner Sufri December 2, 2015 */

	//lfreq= (ftsize/2)* 60 / 125; Changed to below. April 27, 2016.
	lfreq=(ftsize/2);
	nt2= nextpower2(csize);
	nt2_ori=ftsize;
	fprintf(stdout,"outputing CC\n");
       struct seis_header seis;
       i=c_b;
       sprintf(file_name,"%s.%s-%s.seis",dirtime,g_line[i],g_pt[i]);
       FILE *pFile;
       pFile=fopen(file_name,"w");
       seis.rline=tr[0].rline;
       seis.rpt=tr[0].rpt;
       seis.nt=nt2/5;
       seis.x=tr[0].x;
       seis.y=tr[0].y;
       seis.t0=tr[0].t0;
       fwrite(&seis,1,SEIS_HEADER_SZ,pFile);
       i=nt2/5;
       fwrite(tr_down[0].data,i,4,pFile);
       fclose(pFile);
	/*for(i=c_b;i<c_e;i++)
	  {
	    for(j=0;j<ng;j++)
	      {  
		if(nstack[i-c_b][j]==0)
		  continue;

		sprintf(file_name,"%s_%s_%s/COR_%s_%s_%s_%s.SAC",g_line[i],g_pt[i],dirtime,g_line[i],g_pt[i],g_line[j],g_pt[j]);
		shd_cor2.npts=2*lag+1;
		//shd_cor2.delta=1/250.0*nt2_ori/nt2; //below is the correct one, April 27 2016. From Yadong's email.
		shd_cor2.delta=ori_delta;
		shd_cor2.b=-1*shd_cor2.delta*lag;
		shd_cor2.user1=atoi(g_line[i]);
		shd_cor2.user2=atoi(g_pt[i]);
		shd_cor2.evla=gx[i];
		shd_cor2.evlo=gy[i]; 
		shd_cor2.user5=atoi(g_line[j]);
		shd_cor2.user6=atoi(g_pt[j]);
		shd_cor2.stla=gx[j];
		shd_cor2.stlo=gy[j];
		shd_cor2.user9=nstack[i-c_b][j];
              //fprintf(stdout,"filename:%s",file_name);
		write_sac (file_name, cor_stack[i-c_b][j], &shd_cor2 );
	      }
	  }*/
   }



prepare(struct arginfo *arg)
   {
     //     fprintf(stderr,"gill aaaa load %f\n",tr[0].data[257144]);
	struct trace *ptr, *lowtr;
	int ntr, itr, i, k, nt, nt2;
	int ii;
	float dist,temp_max;

	itr= arg->itr;
	ntr= arg->ntr;
   	nt= tr[itr].nt;
	nt2= nextpower2(nt);
	float *temp;

	/* double used for agc */
	if( (temp= (double *)(malloc(4*nt2/2))) == NULL)
	 {
	 	fprintf(stderr,"cannot alloc temp memory nt= %d\n",nt);
		exit(-1);
	 }
	for(k= 0; k< ntr; k++)
	   {
	     //ptr= &tr[itr+k];
           ptr=&tr[0];
		/* Time doamin agc to equivalize time windows */
	     taper(ptr->data,nt,lwind);
	     forfft(ptr->data,nt2,1);

	     //frequency whitening

	     /*for(i=0;i<nt2/2;i++)
	       {
		 temp[i]=sqrt(ptr->data[i*2]*ptr->data[i*2]+ptr->data[i*2+1]*ptr->data[i*2+1]);
	       }	
	     for(i=0;i<nt2/2;i++)
	       {
		 temp_max=0;
		 for(ii=i-100;ii<i+100;ii++)
		   {
		     if(ii<0||ii>=nt2/2)
		       continue;
		     if(temp[ii]>temp_max)
		       temp_max=temp[ii];
		   }
		 ptr->data[i*2]=ptr->data[i*2]/temp_max;
		 ptr->data[i*2+1]=ptr->data[i*2+1]/temp_max;
	       }
	     
	   }*/
       }
	free(temp);
           
           // low-pass filter
           lowtr=&tr_low[0];
          int lf;
           lf=nt2/2;
       
       fprintf(stderr,"lf = %d ______________________________________\n",lf);
       ptr=&tr[0];
           for(i=0;i<lf;i++)
           {
               if(i>=1800 && i<=5400)
               {
                   lowtr->data[i*2]=ptr->data[i*2];
                   lowtr->data[i*2+1]=ptr->data[i*2+1];
                  // fprintf(stderr,"do_low_pass_~~~~~~~ Frequency = %d,real = %f, image = %f \n",i,lowtr->data[i*2],lowtr->data[i*2+1]);
               }
               else
               {
                   lowtr->data[i*2]=0;
                   lowtr->data[i*2+1]=0;
               }
           }
           //for(i=0;i<lf;i++)
           //{
            //   fprintf(stderr,"Frequency = %d, Real = %d Image= %d\n",i,tr[0].data[i*2],tr[0].data[i*2+1]);
           //}
           //for(i=0;i<lf;i++)
           //{
           //    fprintf(stderr,"Frequency = %d, Real = %d Image= %d\n",i,lowtr->data[i*2],lowtr->data[i*2+1]);
           //}
   }

#ifdef OLDAGC
agc_scale(	float	*data,
		int	nt,
		float	*gain,
		int	lwind)
   {
	int lwind2, i, j, j1, j2;
	double sum, fabs();

	lwind2= lwind/2;

	for(i=0; i<nt; i++)
	   {
		j1= i-lwind2;
		if(j1 < 0) j1= 0;
		j2= i+lwind2;
		if(j2 >= nt) j2= nt-1;
		sum= 0.0;
		for(j=j1; j<= j2; j++) sum += fabs(data[j]);
		if(sum < 1.e-15) sum= 1.0e-15;
		gain[i]= 1.0/sum;
	   }
	for(i=0; i<nt; i++) data[i] *= gain[i];
   }
#endif
agc_scale(	float	*data,
		int	nt,
		double	*gain,
		int	lwind)
   {
	int lwind2, i, j, j1, j2;
	double sum, fabs();

	lwind2= lwind/2;


	gain[0]= fabs(data[0]);
	for(i=1; i<nt; i++)
		gain[i]= gain[i-1] + fabs(data[i]);
	for(i=0; i<nt; i++)
	   {
		j1= i-lwind2;
		if(j1 < 0) j1= 0;
		j2= i+lwind2;
		if(j2 >= nt) j2= nt-1;
		sum= gain[j2] - gain[j1];
		if(sum < 1.e-15) sum= 1.0e-15;
		gain[i]= 1.0/sum;
	   }
	for(i=0; i<nt; i++) data[i] *= gain[i];
   }

taper(		float	*data,
		int	nt,
		int	lramp)
   {
	double dwt, wt;
	int i;

   	dwt= 1.0/(lramp);
	wt= 0.0;
	for(i=0; i<lramp; i++)
	   {
	   	data[i] *= wt;
		data[nt-i-1] *= wt;
		wt += dwt;
	   }
   }
	   	
pthread_t	thr[16];

int launch(void(*update)(void *), int nthread)
   {
	/* this routine divides the work, and launches a thread for
	   each part.  It then waits for all of them to finish
	 */
	int load, ithr;
	load= (ns + nthread -1)/nthread;
	int i;
	for(ithr= 0; ithr < nthread; ithr++)
	   {
		/* launch threads */
		arg[ithr].itr = ithr*load;
		arg[ithr].ntr = load;
		if(arg[ithr].itr +arg[ithr].ntr > ns) arg[ithr].ntr= ns-arg[ithr].itr;
		arg[ithr].id  = ithr;
		thr[ithr]= ithr;
		if(pthread_create(&thr[ithr],NULL,(void *)update,&arg[ithr]))
		   {
		   	fprintf(stderr,"failure to create thread\n");
			exit(-1);
		   }
	   }
	for(ithr= 0; ithr < nthread; ithr++)
	   {
		/* wait for thread to finish */
	   	if(pthread_join(thr[ithr],NULL))
		   {
		   	fprintf(stderr,"failure to rejoin thread\n");
			exit(-1);
		   }
	   }
	return 0;
   }

launch2(void(*update)(void *), int nthread)
   {
	/* this routine divides the work, and launches a thread for
	   each part.  It then waits for all of them to finish
	 */
	int i, ithr,lfreq, nt, nt2;
	nt= tr[0].nt;
	nt2= nextpower2(nt);
	lfreq= (nt2/2);    // April 27, 2016 changed from (nt/2)*60/125
	
//	printf("\n from launch2 lfreq is %d\n",lfreq);                          /* Oner  December 2, 2015 */
//	printf("from launch2 nt is %d\n",nt);				   /* Oner  December 2, 2015 */
//	printf("from launch2 nt2, which is nexpower2(nt), is %d\n",nt2);        /* Oner  December 2, 2015 */ 

	for(i=c_b; i<c_e; )
	{
		/* launch threads */
		for(ithr= 0; ithr < nthread; ithr++)
		   {
		   	if(i+ithr >=c_e) break;
			if(ngns[i+ithr]==-1)
			  {
			    i++;
			    ithr--;
			    continue;
			  }
			arg2[ithr].itr = i+ithr;
			arg2[ithr].lfreq = lfreq;
			arg2[ithr].id  = ithr;
			thr[ithr]= ithr;
			if(pthread_create(&thr[ithr],NULL,(void *)update,&arg2[ithr]))
			   {
			   	fprintf(stderr,"failure to create thread\n");
				exit(-1);
			   }
		   }
		for(ithr= 0; ithr < nthread; ithr++)
		   {
			/* wait for thread to finish */
		   	if(i+ithr >= c_e) break;
		   	if(pthread_join(thr[ithr],NULL))
			   {
			   	fprintf(stderr,"failure to rejoin thread\n");
				exit(-1);
			   }
		   }
		i += nthread;
	   }
   }

correlate(struct arginfo2 *arg2)
   {
     
     struct complex *ctr1, *ctr2, *ccor;
     char file_name[100];
     //SAC_HD shd_cor1;
     float *cor,cor1[2*LAG+1];
     float max;
     int i, j, itr, lfreq;
     int ii;
     itr= ngns[arg2->itr];
     

     lfreq= arg2->lfreq;
     int  nt2,nt2_ori;
     nt2= nextpower2(lfreq);
     nt2_ori=nextpower2(tr[itr].nt);
     //    fprintf(stderr,"lfreq= %d nt2 %d nt2 %d\n",lfreq,nt2,nextpower2(tr[itr].nt));
     if( (cor= (float *)(malloc(4*lfreq*2))) == NULL )
       {
	 fprintf(stderr,"cannot alloc memory for cor, lfreq=%d\n",lfreq);
	 exit(-1);
       }
     ccor= (struct complex *)(cor);
     
     ctr1= (struct complex *)(tr_low[0].data);


    // for(i=0; i<ns; i++)
    //   {
	 ctr2= (struct complex *)(tr_down[0].data);
       fprintf(stderr,"lfreq= %d~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",lfreq);
	 for(j=0; j<lfreq; j++)
	   {
	     /*ccor[j].re=  ctr1[j].re * ctr2[j].re
	       + ctr1[j].im * ctr2[j].im;
	     ccor[j].im=  ctr1[j].im * ctr2[j].re
	       - ctr1[j].re * ctr2[j].im;*/
           ccor[j].re=ctr1[j].re;
           ccor[j].im=ctr1[j].im;
           
	  // }
       }
       fprintf(stderr,"before invfft~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	 invfft(cor,nt2,-1);
       fprintf(stderr,"after invfft~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
       for(i=0;i<nt2;i++)
       {
           //fprintf(stderr,"Frequency = %d,downsampling point= %f\n",i,cor[i]);
       }
           //downsampling
           j=0;
           for(ii=0;ii<nt2;ii++)
           {
               if(ii%5==0)
               {
                   //fprintf(stderr,"Frequency = %d,low_pass: number = %f \n",ii,cor[ii]);
                   tr_down[0].data[j]=cor[ii];
                                      j=j+1;
               }
           }
       fprintf(stderr,"last j= %d~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",j);
	 /*if(lag>nt2/2)
	   lag=nt2/2;
	 max=0;
	 for(ii=nt2-lag;ii<nt2;ii++)
	   {
	     if(fabs(cor[ii])>max)
	       {
		 max=fabs(cor[ii]);
	       }
	     cor1[ii-nt2+lag]=cor[ii];
	   }
	 for(ii=0;ii<=lag;ii++)
	   {
	     if(fabs(cor[ii])>max)
	       max=fabs(cor[ii]);
	     cor1[ii+lag]=cor[ii];
	    
	   }
	 for(ii=0;ii<=2*lag;ii++)
	   {
	     cor_stack[arg2->itr-c_b][tr[i].index][2*lag-ii]+=cor1[ii]/max;
	   }
	 nstack[arg2->itr-c_b][tr[i].index]++;
       }*/
     free(cor);
       fprintf(stderr,"end invfft~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
   }
