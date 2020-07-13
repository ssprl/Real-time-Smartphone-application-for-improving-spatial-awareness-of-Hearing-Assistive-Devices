#include "twoDOA.h"

/*
 twoDOA.c
 Created by Anshuman Ganguly, Yiya Hao, and Abdullah Kucuk on 4/28/2017

 DOA
 Frame Based
 with VAD

Inputs are:
 1- 2 channel input 
 
 Outputs are:
 1- corrtheta_est which using VAD

*/
double max1(double a, double b) /* Finds max of inputs*/
{
	if (a > b)
		return a;
	else
		return b;
}

twomicDOA*
newDOA(int NFFT)  /*Initializations for twomicDOA function*/
{

	
	twomicDOA* newDOA= (twomicDOA*)malloc(sizeof(twomicDOA));

	newDOA->prevmag_fft_Framex= (double*)malloc((NFFT/2+1)*sizeof(double));
	newDOA->prevmag_fft_Frameh= (double*)malloc((NFFT / 2 + 1) * sizeof(double));


	for (int i = 0; i < NFFT/2+1; i++)
	{
		newDOA->prevmag_fft_Framex[i] = 0;
		newDOA->prevmag_fft_Frameh[i] = 0;
	}
	newDOA->prevcorrtheta_est = 0;
	newDOA->SFxmax = 0;
	newDOA->SFxavg = 0;
	newDOA->flagSFx = 0;
	
	newDOA->doIT = *twomicDOAv_00;

	
	return newDOA;
}

double *fftshift(double *x) /*This function for fftshift operations*/
{
	int npoint = 1024;
	double *output = (double*)malloc((npoint * 2) * sizeof(double));
	int k;
	for (k = 0; k < 1024; k++) {
		output[k + npoint] = x[k];
		output[k] = x[k + npoint];
	}
	return output;
}

int findmax(double *arr3) { // finds max index of array index that bigger than sca

	
	double tmp_max = -999999;
	int tmp;
	for (int ii = 0; ii < 1599; ii++) {


		if (tmp_max < arr3[ii]  ) {
			tmp_max = arr3[ii];
			tmp = ii;
		}

	}

	return tmp;

}

void twomicDOAv_00(twomicDOA *twoDOA,double *x, double *h, double *win, int i, double *prevmag_fft_Framex, double *prevmag_fft_Frameh,double prevcorrtheta_est,double SFxmax, double SFxavg,double flagSFx, int fs, int l,int NFFT ) { /*main function that finds angle of speaker location*/

	
	
	int n, N;
	n = l;
	N = n;

	int Ncorr = 1599;  

	double * lags = (double*)malloc(l * sizeof(double));

	for (int i = 0; i < l-1; i++)
	{
		lags[i] = -799 + i;
		lags[i] = lags[i] / fs;
	}
	



	int n_F = 1;
	twoDOA->corrtheta_est = 0;

	double deltau = 0;

	int Trainingframe = 10; /*%Training frames for threshold calculation*/
	int	DurationThresh = 1; /*% Duration thresold for smooth transition*/

	double *X_real= (double*)malloc(NFFT * sizeof(double));  /* Define memory space for fft of real part of x*/
	double *X_img = (double*)malloc(NFFT * sizeof(double));  /* Define memory space for fft of imaginaty part of x */

	double *H_real = (double*)malloc(NFFT * sizeof(double)); /* Define memory space for fft of real part of h*/
	double *H_img = (double*)malloc(NFFT * sizeof(double));  /* Define memory space for fft of imaginary part of h */

	for (int j = 0; j < NFFT; j++)
	{
		X_real[j] = 0;
		X_img[j] = 0;
		H_real[j] = 0;
		H_img[j] = 0;

	}


	double *mag_fft_Framex= (double*)malloc((NFFT/2+1) * sizeof(double));  /* Define memory space for magnitude of fft of first channel input*/
	double *sfx = (double*)malloc((NFFT / 2 + 1) * sizeof(double));        /* Define memory space for spectral flux of first channel input*/

	double *mag_fft_Frameh = (double*)malloc((NFFT / 2 + 1) * sizeof(double)); /* Define memory space for magnitude of fft of  second channel input*/
	double *sfh = (double*)malloc((NFFT / 2 + 1) * sizeof(double));        /* Define memory space for spectral flux of  second channel input*/

	int idx = 0;

	for (int i = 0; i < NFFT/2+1; i++)
	{
		mag_fft_Framex[i] = 0;		
		sfx[i] = 0;

		mag_fft_Frameh[i] = 0;
		sfh[i] = 0;



	}

	double *x_frame_2 = (double*)malloc((NFFT) * sizeof(double)); /* Define memory space for windowed first channel input*/
	double *h_frame_2 = (double*)malloc((NFFT) * sizeof(double)); /* Define memory space for windowed second channel input*/

	for (int i = 0; i < n; i++)
	{

		x_frame_2[i] = win[i] * x[i];
		h_frame_2[i] = win[i] * h[i];
	

	}

	/*We need this for loop for FFT and NFFT=2048 */
	for (int i = n; i < NFFT; i++)
	{
		x_frame_2[i] = 0;
		h_frame_2[i] = 0;
	}

	Transform *trans1; /*Initializaton for FFT*/
	trans1 = newTransform(NFFT); /*Initializaton for FFT*/
	trans1->doTransform(trans1, x_frame_2); /*This fuction takes FFT of input*/

	Transform *trans2; /*Initializaton for FFT*/
	trans2 = newTransform(NFFT); /*Initializaton for FFT*/
	trans2->doTransform(trans2, h_frame_2); /*This fuction takes FFT of input*/



	double *fft_Framex_real = (double*)malloc((NFFT) * sizeof(double)); /* Define memory space for real part of fft of x*/
	double *fft_Framex_img = (double*)malloc((NFFT) * sizeof(double)); /* Define memory space for imaginary part of fft of x*/

	double *fft_Frameh_real = (double*)malloc((NFFT) * sizeof(double)); /* Define memory space for real part of fft of h*/
	double *fft_Frameh_img = (double*)malloc((NFFT) * sizeof(double)); /* Define memory space for imaginary part of fft of h*/

	for (int j = 0; j < NFFT; j++)
	{
		X_real[j] = trans1->real[j];
		X_img[j] = trans1->imaginary[j];
		H_real[j] = trans2->real[j];
		H_img[j] = trans2->imaginary[j];


		fft_Framex_real[j] = X_real[j] / n;
		fft_Framex_img[j] = X_img[j] / n;

		fft_Frameh_real[j] = H_real[j] / n;
		fft_Frameh_img[j] = H_img[j] / n;

	}

	/*%% Implementing VAD
     % Calculation of Spectral flux(SF)*/

	double *currmag_fft_Framex = (double*)malloc((NFFT / 2 + 1) * sizeof(double));
	double *currmag_fft_Frameh = (double*)malloc((NFFT / 2 + 1) * sizeof(double));

	double sum_sfx=0;
	double sum_sfh=0;

	for (int i = 0; i < NFFT/2+1; i++)
	{
		
		mag_fft_Framex[i] = sqrt(pow(fft_Framex_real[i],2)+ pow(fft_Framex_img[i], 2));
		mag_fft_Frameh[i] = sqrt(pow(fft_Frameh_real[i], 2) + pow(fft_Frameh_img[i], 2));
		
		currmag_fft_Framex[i] = mag_fft_Framex[i];
		currmag_fft_Frameh[i] = mag_fft_Frameh[i];

		sfx[i] = currmag_fft_Framex[i] - prevmag_fft_Framex[i];
		sfh[i] = currmag_fft_Frameh[i] - prevmag_fft_Frameh[i];

		sum_sfx += pow(sfx[i],2);
		sum_sfh += pow(sfh[i], 2);
		
	
		twoDOA->prevmag_fft_Framex[i] = currmag_fft_Framex[i];
		twoDOA->prevmag_fft_Frameh[i] = currmag_fft_Frameh[i];

	}

	double SFx, SFh;
	SFx = sum_sfx / (n_F*NFFT);
	SFh = sum_sfh / (n_F*NFFT);

	/*END of Spectral Flux Calculation*/

	/* Cross Correlation calculation*/
	double *R12_real = (double*)malloc((NFFT ) * sizeof(double)); /* Define memory space for real part of fft of crosscorrelation*/
	double *R12_img = (double*)malloc((NFFT) * sizeof(double)); /* Define memory space for imaginary part of fft of crosscorrelation*/

	for (int i = 0; i < NFFT; i++)
	{
		R12_real[i] = X_real[i] * H_real[i] - ((-1 * H_img[i]) * X_img[i]); 
		R12_img[i] = X_real[i] * (-1 * H_img[i]) + (X_img[i] * H_real[i]);
		
	}
	Transform *trans3; /*Initializaton for IFFT*/
	trans3 = newTransform(NFFT); /*Initializaton for IFFT*/
	trans3->invTransform(trans3, R12_real, R12_img); /*This fuction takes FFT of input*/



	double *r12_temp = (double*)malloc((NFFT) * sizeof(double));
	r12_temp = fftshift(trans3->real);

	

	double *r12= (double*)malloc((n) * sizeof(double));
	double *r12_abs = (double*)malloc((n) * sizeof(double));

	for (int i = (NFFT / 2 + 1 - (Ncorr - 1) / 2)-1; i < NFFT / 2 + 1 + (Ncorr - 1) / 2; i++)
	{
		r12[i-(NFFT / 2 + 1 - (Ncorr - 1) / 2-1)] = r12_temp[i];
		r12_abs[i - (NFFT / 2 + 1 - (Ncorr - 1) / 2 - 1)] = fabs(r12[i - (NFFT / 2 + 1 - (Ncorr - 1) / 2 - 1)]);

	}
	

	 idx = findmax(r12_abs);
	 
	 /* END of Cross Correlation calculation*/

	 deltau = lags[idx];

	 double  ret, val, degree, theta_est;

	
	 val = 180.0 / PI;

	 degree = deltau * 343 / (0.13);
	// printf("i=%1.9g\n", degree);
	 if (degree>1)
		 theta_est = 0;
	 else if (degree<-1)
		 theta_est = PI * val;            /*% Estimated angle for GCC*/
	 else
		 theta_est = acos(degree) * val; /*% Estimated angle for GCC*/


	double  corrtheta_est = 0;
	
	/*%% Finding direction of speech only based on SF*/
	/*%% Setting SF Threshold*/


	if (i==0) {
		SFxavg = SFx;
		SFxmax = SFx;
	}
	else {

		if (i<=Trainingframe) {
			if (SFx > SFxmax) {
			
				SFxmax = SFx;
			
			}
			SFxavg = (SFxavg*(i-1 ) + SFx) / i; /*% Spectral Flux Calculation for VAD*/
		}
	}

	if ((i >= Trainingframe) && (SFx > SFxavg) && (SFx > SFxmax)) {
	
		flagSFx +=   1;								
	}
	else {
		flagSFx = 0;
	}


	if (i >= Trainingframe) {
		if (flagSFx > DurationThresh)
			corrtheta_est = theta_est;   /*% Estimated angle with VAD*/
		
		else
			corrtheta_est = prevcorrtheta_est; /*% Estimated angle with VAD*/
	
	
	}
	else {
	
		corrtheta_est = 0;
	}

	
	twoDOA->corrtheta_est = corrtheta_est;
	twoDOA->prevcorrtheta_est = corrtheta_est;

	twoDOA->SFxmax = SFxmax;

	twoDOA->SFxavg = SFxavg;

	twoDOA->flagSFx = flagSFx;


}

/*Transform are used for FFT and IFFT operations*/
Transform*
newTransform(int points)
{
	Transform* newTransform = (Transform*)malloc(sizeof(Transform));

	newTransform->points = points;
	newTransform->real = (double*)malloc(points * sizeof(double));
	newTransform->imaginary = (double*)malloc(points * sizeof(double));
	newTransform->sine = NULL;
	newTransform->cosine = NULL;
	newTransform->doTransform = *FFT;
	newTransform->invTransform = *IFFT;
	newTransform->sine = (double*)malloc((points / 2) * sizeof(double));
	newTransform->cosine = (double*)malloc((points / 2) * sizeof(double));
	//precompute twiddle factors
	double arg;
	int i;
	for (i = 0; i<points / 2; i++)
	{
		arg = -2 * M_PI*i / points;
		newTransform->cosine[i] = cos(arg);
		newTransform->sine[i] = sin(arg);
	}
	return newTransform;
}

void
FFT(Transform* fft, double* input)
{
	int i, j, k, L, m, n, o, p, q, r;
	float tempReal, tempImaginary, cos, sin, xt, yt;
	k = fft->points;
	for (i = 0; i<k; i++)
	{
		fft->real[i] = input[i];
		fft->imaginary[i] = 0;
	}

	j = 0;
	m = k / 2;
	//bit reversal
	for (i = 1; i<(k - 1); i++)
	{
		L = m;
		while (j >= L)
		{
			j = j - L;
			L = L / 2;
		}
		j = j + L;
		if (i<j)
		{
			tempReal = fft->real[i];
			tempImaginary = fft->imaginary[i];
			fft->real[i] = fft->real[j];
			fft->imaginary[i] = fft->imaginary[j];
			fft->real[j] = tempReal;
			fft->imaginary[j] = tempImaginary;
		}
	}
	L = 0;
	m = 1;
	n = k / 2;
	//computation
	for (i = k; i>1; i = (i >> 1))
	{
		L = m;
		m = 2 * m;
		o = 0;
		for (j = 0; j<L; j++)
		{
			cos = fft->cosine[o];
			sin = fft->sine[o];
			o = o + n;
			for (p = j; p<k; p = p + m)
			{
				q = p + L;
				xt = cos*fft->real[q] - sin*fft->imaginary[q];
				yt = sin*fft->real[q] + cos*fft->imaginary[q];
				fft->real[q] = (fft->real[p] - xt);
				fft->imaginary[q] = (fft->imaginary[p] - yt);
				fft->real[p] = (fft->real[p] + xt);
				fft->imaginary[p] = (fft->imaginary[p] + yt);
			}
		}
		n = n >> 1;
	}
}
void
IFFT(Transform* fft, double* inputreal, double* inputimaginary)
{

	int i, j, k, L, m, n, o, p, q, r;
	double tempReal, tempImaginary, cos, sin, xt, yt;
	k = fft->points;
	for (i = 0; i<k; i++)
	{
		fft->real[i] = inputreal[i];
		fft->imaginary[i] = (-1)*inputimaginary[i];
	}

	j = 0;
	m = k / 2;
	//bit reversal
	for (i = 1; i<(k - 1); i++)
	{
		L = m;
		while (j >= L)
		{
			j = j - L;
			L = L / 2;
		}
		j = j + L;
		if (i<j)
		{
			tempReal = fft->real[i];
			tempImaginary = fft->imaginary[i];
			fft->real[i] = fft->real[j];
			fft->imaginary[i] = fft->imaginary[j];
			fft->real[j] = tempReal;
			fft->imaginary[j] = tempImaginary;
		}
	}
	L = 0;
	m = 1;
	n = k / 2;
	//computation
	for (i = k; i>1; i = (i >> 1))
	{
		L = m;
		m = 2 * m;
		o = 0;
		for (j = 0; j<L; j++)
		{
			cos = fft->cosine[o];
			sin = fft->sine[o];
			o = o + n;
			for (p = j; p<k; p = p + m)
			{
				q = p + L;
				xt = cos*fft->real[q] - sin*fft->imaginary[q];
				yt = sin*fft->real[q] + cos*fft->imaginary[q];
				fft->real[q] = (fft->real[p] - xt);
				fft->imaginary[q] = (fft->imaginary[p] - yt);
				fft->real[p] = (fft->real[p] + xt);
				fft->imaginary[p] = (fft->imaginary[p] + yt);
			}
		}
		n = n >> 1;
	}



	for (i = 0; i<k; i++)
	{
		fft->real[i] = fft->real[i] / k;
		fft->imaginary[i] = fft->imaginary[i] / k;
	}

}
void
transformMagnitude(Transform* transform, double* output)
{
	int n;
	for (n = 0; n<transform->points; n++)
	{
		output[n] = sqrt(transform->real[n] * transform->real[n] + transform->imaginary[n] * transform->imaginary[n]);
	}
}

void
invtranMagnitude(Transform* transform, double* output)
{
	int n;
	float a;
	a = 1.0 / transform->points;
	for (n = 0; n < transform->points; n++)
	{
		output[n] = a * sqrt(transform->real[n] * transform->real[n] + transform->imaginary[n] * transform->imaginary[n]);
	}
}

void
destroyTransform(Transform** transform)
{
	if (*transform != NULL) {
		if ((*transform)->cosine != NULL) {
			free((*transform)->cosine);
			(*transform)->cosine = NULL;
		}
		if ((*transform)->sine != NULL) {
			free((*transform)->sine);
			(*transform)->sine = NULL;
		}
		if ((*transform)->real != NULL) {
			free((*transform)->real);
			(*transform)->real = NULL;
		}
		if ((*transform)->imaginary != NULL) {
			free((*transform)->imaginary);
			(*transform)->imaginary = NULL;
		}
		free(*transform);
		*transform = NULL;
	}
}
