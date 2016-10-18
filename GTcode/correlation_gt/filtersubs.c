
#include	<stdio.h>
float sin_table[] =
   {
	1.0000000e+00,	/* sin(pi/2) */
	7.0710678e-01,	/* sin(pi/4) */
	3.8268343e-01,	/* sin(pi/8) */
	1.9509032e-01,	/* sin(pi/16) */
	9.8017140e-02,	/* sin(pi/32) */
	4.9067674e-02,	/* sin(pi/64) */
	2.4541228e-02,	/* sin(pi/128) */
	1.2271538e-02,	/* sin(pi/256) */
	6.1358846e-03,	/* sin(pi/512) */
	3.0679568e-03,	/* sin(pi/1024) */
	1.5339802e-03,	/* sin(pi/2048) */
	7.6699032e-04,	/* sin(pi/4096) */
	3.8349519e-04,	/* sin(pi/8192) */
	1.9174760e-04,	/* sin(pi/16384) */
	9.5873799e-05, /* sin(pi/32768) */
	4.7936899e-05, /* sin(pi/65536) */
	2.3968499e-05, /* sin(pi/131072) */
	1.1984224e-05, /* sin(pi/262144) */
	5.9921125e-06, /* sin(pi/524288) */
	2.9960562e-06, /* sin(pi/1048576) */
	1.4980281e-06, /* sin(pi/2097152) */

	7.49014068e-07, /* sin(pi/4194304) */
	3.74507034e-07, /* sin(pi/8388608) */
	1.87253517e-07, /* sin(pi/16777216) */
	9.36267585e-08, /* sin(pi/33554432) */
	4.68133792e-08, /* sin(pi/67108864) */
	2.34066896e-08  /* sin(pi/134217728) */
   };
/*
 * cfft - radix 2 FFT for complex data
 *
 * n is the number of complex points.
 * x is the single precision (not double!!) complex vector.
 * isign is the sign of the transform.
 *
 * the routine does no normalization
 */
struct complex { float re; float im; };
cfft(x,n,isign)
struct complex *x;
int n,isign;
   {
	register struct complex *px, *qx, *rx;
	struct complex *limit, *qlimit, dtemp;
	float cn, sn, cd, sd, temp, real, imag;
	int m, j, istep;
	float *psintab;
	extern float sin_table[];

	limit= x + n;
	j= 0;
	for(px=x; px<limit; px++)
	   {
		if(px < (qx= x+j))
		   {	dtemp= *qx; *qx= *px; *px= dtemp;	}
		m = n>>1;
		while( m>=1 && j>=m )
		   { j-= m; m>>= 1;    }
		j+= m;
	   }
	rx= x+1;
	for(px=x; px<limit; px+= 2, rx+= 2)
	   {
		temp= rx->re;
		rx->re= px->re -temp;
		px->re += temp;
		temp= rx->im;
		rx->im= px->im -temp;
		px->im += temp;
	   }
	j=2;
	psintab= sin_table;
	while( j < n )
	   {
		istep= j<<1;
		sd= *psintab++;
		temp= *psintab;
		cd= 2.0 * temp * temp;
		cn= 1.0;
		sn= 0.0;
		if( isign < 0 ) sd= -sd;
		qlimit= x+j;
		for(qx=x; qx< qlimit; qx++)
		   {
			for(px=qx; px<limit; px+= istep)
			   {
				rx= px + j;
				real= cn * rx->re - sn * rx->im;
				imag= sn * rx->re + cn * rx->im;
				rx->re = px->re - real;
				rx->im = px->im - imag;
				px->re += real;
				px->im += imag;
			   }
			temp= cd * cn + sd * sn;
			sn += (sd * cn - cd * sn);
			cn -= temp;
		   }
		j= istep;
	   }
	return 0;
   }

forfft(x,n,isign)
register struct complex *x;
int n, isign;
   {
	register struct complex *px, *rx;
	float cn, sn, cd, sd, arg;
	float are, aim, bre, bim, real, imag;
	float *psintab;
	extern float sin_table[];
	int k;

	cfft(x,n/2,isign);

	/* do DC and Nyquist */
	real= x[0].re;
	imag= x[0].im;
	x[0].re= real + imag;
	x[0].im= real - imag;

	/* set up for sine recurrsion */
	psintab= sin_table;
	for(k=4; k<n; k <<= 1) psintab++;
	sd= *psintab++;
	real= *psintab;
	cd= 2.0 * real * real;

	sn= 0.0;
	cn= 1.0;
	if(isign < 0) sd= -sd;
	px= x + 1;
	rx= x + n/2 -1;
	while( px <= rx )
	   {
		real = cd*cn + sd*sn;
		imag = sd*cn - cd*sn;
		cn -= real;
		sn += imag;

		are= 0.5*(px->re + rx->re);
		aim= 0.5*(px->im - rx->im);
		bre= 0.5*(px->im + rx->im);
		bim= 0.5*(rx->re - px->re);

		real= bre*cn - bim*sn;
		imag= bre*sn + bim*cn;

		px->re = are + real;
		px->im = aim + imag;
		rx->re = are - real;
		rx->im = imag - aim;

		px++;
		rx--;
	   }
	if(abs(isign) > 1)
	   {
		x[n/2].re= x[0].im;
		x[0].im= x[n/2].im = 0.0;
	   }
   }

invfft(x,n,isign)
register struct complex *x;
int n, isign;
   {
	register struct complex *px, *rx;
	float cn, sn, cd, sd;
	float are, aim, bre, bim, real, imag;
	float *psintab;
	extern float sin_table[];
	int k;

	if(abs(isign) > 1) x[0].im= x[n/2].re;

	/* do DC and Nyquist */
	real= x[0].re;
	imag= x[0].im;
	x[0].re= real + imag;
	x[0].im= real - imag;

	/* set up for sine recurrsion */
	psintab= sin_table;
	for(k=4; k<n; k <<= 1) psintab++;
	sd= *psintab++;
	real= *psintab;
	cd= 2.0 * real * real;

	sn= 0.0;
	cn= 1.0;
	if(isign < 0) sd= -sd;
	px= x + 1;
	rx= x + n/2 -1;
	while( px <= rx )
	   {
		real = cd*cn + sd*sn;
		imag = sd*cn - cd*sn;
		cn -= real;
		sn += imag;

		are= (px->re + rx->re);
		aim= (px->im - rx->im);
		bre= (px->re - rx->re);
		bim= (px->im + rx->im);

		real= bre*cn - bim*sn;
		imag= bre*sn + bim*cn;

		px->re = are - imag;
		px->im = aim + real;
		rx->re = are + imag;
		rx->im = real - aim;

		px++;
		rx--;
	   }
	cfft(x,n/2,isign);
   }


filter(f,dt,x,lx,y,ly,notch)
float *f, dt, *x, *y;
int lx, ly, notch;
   {
	int i, j, nyq, ifilt[4];
	float tmp, wt, dwt;
	struct complex *cy;

	if(x != y)	/* out of place */
	   {
		if(ly > lx)
		   {
			for(i=0;  i<lx; i++) y[i]= x[i];
			for(i=lx; i<ly; i++) y[i] = 0.0;
		   }
		 else	for(i=0; i<ly; i++) y[i]= x[i];
	   }
	 else		/* in place */
	   {
		if(ly > lx)
			for(i=lx; i<ly; i++) y[i] = 0.0;
	   }
	
	for(j=0; j<3; j++)
		for(i=j+1; i<4; i++)
			if(f[j] > f[i])
			   {
				tmp= f[i]; f[i]= f[j]; f[j]= tmp;
			   }
	for(i=0; i<4; i++)
	   {
		ifilt[i] = (int)(dt *ly * f[i]);
		if(ifilt[i] < 0) ifilt[i]= 0;
		if(ifilt[i] > ly/2 +1) ifilt[i]= ly/2 +1;
	   }
	
/*
fprintf(stderr,"filter: lx=%d ly=%d notch=%d dt=%6.4f ifilt=%d,%d,%d,%d\n",
	lx,ly,notch,dt,ifilt[0],ifilt[1],ifilt[2],ifilt[3]);
*/
	forfft(y,ly,1);
	cy= (struct complex *)(y);

	if(notch == 0)
	   {
		for(i=0; i< ifilt[0]; i++) cy[i].re= cy[i].im= 0.0;

		if(ifilt[1] > ifilt[0])
		   {
			dwt= 1.0/(float)(ifilt[1] - ifilt[0]);
			wt= 0.0;
			for(i= ifilt[0]; i< ifilt[1]; i++, wt += dwt)
			   {
				cy[i].re *= wt;
				cy[i].im *= wt;
			   }
		   }

		if(ifilt[3] > ifilt[2])
		   {
			dwt= 1.0/(float)(ifilt[3] - ifilt[2]);
			wt= 1.0;
			for(i= ifilt[2]; i< ifilt[3]; i++, wt -= dwt)
			   {
				cy[i].re *= wt;
				cy[i].im *= wt;
			   }
		   }
		nyq= ly/2+1;
		for(i=ifilt[3]; i<nyq; i++) cy[i].re= cy[i].im= 0.0;
	   }
	 else
	   {
		if(ifilt[1] > ifilt[0])
		   {
			dwt= 1.0/(float)(ifilt[1] - ifilt[0]);
			wt= 1.0;
			for(i= ifilt[0]; i< ifilt[1]; i++, wt -= dwt)
			   {
				cy[i].re *= wt;
				cy[i].im *= wt;
			   }
		   }

		for(i=ifilt[1]; i< ifilt[2]; i++) cy[i].re= cy[i].im= 0.0;

		if(ifilt[3] > ifilt[2])
		   {
			dwt= 1.0/(float)(ifilt[3] - ifilt[2]);
			wt= 0.0;
			for(i= ifilt[2]; i< ifilt[3]; i++, wt += dwt)
			   {
				cy[i].re *= wt;
				cy[i].im *= wt;
			   }
		   }
	   }

	invfft(cy,ly,-1);

	wt= 1.0/(float)(ly);
	for(i=0; i<ly; i++) y[i] *= wt;

	for(i=lx; i<ly; i++) y[i] = 0.0;
   }

/* determine power spectrun of 'x'. Note that x must be of length lfft.
   'lfft is the length of fft to use.  The result 'amp' is of length
   lfft/2 +1.  Composite=1 means add result to existing result.
 */
spectrum(x,lx,amp,lfft,composite)
float *x, *amp;
int lx, lfft, composite;
   {
	int i, nyq;
	struct complex *cx;

	cx= (struct complex *)(x);
	nyq= lx/2 + 1;
	for(i=lx; i<lfft; i++) x[i]= 0.0;
	if(!composite) for(i=0; i<nyq; i++) amp[i]= 0.0;

	forfft(x,lx,1);

	for(i=0; i<nyq; i++)
		amp[i] += cx[i].re * cx[i].re + cx[i].im * cx[i].im;
   }

nextpower2(int n)
   {
	int n2;
	n2= 1;
	while(n2 <= n) n2 *= 2;
	return(n2);
   }
