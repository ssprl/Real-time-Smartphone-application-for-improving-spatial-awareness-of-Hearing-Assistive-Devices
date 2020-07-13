#ifndef TWODOA_H
#define TWODOA_H
#define _USE_MATH_DEFINES
#define PI 3.14159265
#include <stdlib.h>
#include <stdio.h>
#include <math.h>




typedef struct twomicDOA {
	double  *prevmag_fft_Framex, *prevmag_fft_Frameh;
	double corrtheta_est,prevcorrtheta_est, SFxmax, SFxavg, flagSFx;

	void(*doIT)(struct twomicDOA *twoDOA,double *x, double *h, double *win, int i, double *prevmag_fft_Framex, double *prevmag_fft_Frameh,double prevcorrtheta_est,double SFxmax, double SFxavg,double flagSFx, int fs, int l,int NFFT ) ;

} twomicDOA;




twomicDOA * newDOA(int NFFT);
double *fftshift(double *x);

void twomicDOAv_00(twomicDOA *twoDOA, double *x, double *h, double *win, int i, double *prevmag_fft_Framex, double *prevmag_fft_Frameh, double prevcorrtheta_est, double SFxmax, double SFxavg, double flagSFx, int fs, int l, int NFFT);

typedef struct Transform {
	int points;
	double* sine;
	double* cosine;
	double* real;
	double* imaginary;
	void(*doTransform)(struct Transform* transform, double* input);
	void(*invTransform)(struct Transform* transform, double* inputreal, double* inputimaginary);
} Transform;

Transform* newTransform(int points);
void transformMagnitude(Transform* transform, double* output);
void invtranMagnitude(Transform* transform, double* output);
void destroyTransform(Transform** transform);
void FFT(Transform* fft, double* input);
void IFFT(Transform* fft, double* inputreal, double* inputimaginary);





#endif