//
// Created by axk166230 on 9/11/2017.
//

#ifndef TWOCHANDOA_V4_TWODOA_H
#define TWOCHANDOA_V4_TWODOA_H
float max1(float a, float b);


#define _USE_MATH_DEFINES
#define PI 3.14159265
#define MAXIT 1000
#define EULER 0.5772156649
#define FPMIN 1.0e-30
#define EPSS  1.0e-7f
#define pi 3.1415926535897932384626433832795
#define vall  180.0 / PI;
//#define fss 16000
//#define NFFTT  2048
//#define nn = 1600
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <android/log.h>
#include <time.h>
//#include "audio_main.h"





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
void destroyTransform(Transform* transform);
void FFT(Transform* fft, double* input);
void IFFT(Transform* fft, double* inputreal, double* inputimaginary);

typedef struct logMMSE_parameters {
    double  *ksi, *Xk_prev, *noise_mean, *sig1, *spec_real, *spec_imag, *spec_real2, *spec_imag2, *output, *output2, *wind, *sys_wind, *windowed_signal, *windowed_signal2, *x_old, *x_final, *y_old, *y_final;

    double *noise_mu2, *prev_frame, *prev_frame2, *last_output, *up_noise_mu2;
    double vad_decisionn, avg_timing;
    int NFFT, n, PERC, len1, len2;

    ////////* Memory Update *////////////

    double *insign_x,*insign_h,*prev_x,*prev_h,*sig2, *sigforgama, *gammak, *maxgammak, *log_sigma_k;
    //double *sig2, *sigforgama, *gammak, *maxgammak, *log_sigma_k;
    double *A, *vk, * ei_vk, *hw, *spec_hw_real, *spec_hw_imag, *spec_hw_real2, *spec_hw_imag2;
    Transform *trans33, *trans, *trans_y, *trans22, *trans22_y;

    ///////////////////////

    void(*enhance_signal)(struct logMMSE_parameters *logMMSE, double *insign, double *insign2, int frame_counter, int Trainingframe);

} logMMSE_parameters;

void enhance_signal(logMMSE_parameters *logMMSE, double *insign, double *insign2, int frame_counter, int Trainingframe);
logMMSE_parameters* newparameters(int NFFT, int n, int PERC);
double *Hamming(int np);


typedef struct twomicDOA {
    double  *prevmag_fft_Framex, *prevmag_fft_Frameh;
    /* Update for memory issue */
    double  *lags, *X_real,*X_img,*H_real,*H_img,*mag_fft_Framex,*sfx,*mag_fft_Frameh,*sfh,*x_frame_2,*h_frame_2,*fft_Framex_real,*fft_Framex_img,*fft_Frameh_real,*fft_Frameh_img;
    double  *currmag_fft_Framex,*currmag_fft_Frameh,*R12_real,*R12_img,*r12,*r12_abs,*r12_temp,*windd ;
    int Trainingframe, DurationThresh,n,Ncorr,idx;
    double deltau,SFx, degree, theta_est,corrtheta_est1,sum_sfx;
    Transform *trans1,*trans2,*trans3;

    ///////////////////////////
    double corrtheta_est,prevcorrtheta_est, SFxmax, SFxavg, flagSFx,mic_dist;

    void(*doIT)(struct twomicDOA *twoDOA,double *x, double *h, int i, double *prevmag_fft_Framex, double *prevmag_fft_Frameh,double prevcorrtheta_est,double SFxmax, double SFxavg,double flagSFx, int fs, int l,int NFFT, float ThreadTime, float DurationTime, logMMSE_parameters *logMMSE, bool isEnchanced ) ;

} twomicDOA;



twomicDOA * newDOA(int NFFT, int n, int fs,float ThreadTime,float DurationTime);
//double *fftshift(double *x);

void twomicDOAv_00(twomicDOA *twoDOA, double *x, double *h, int i, double *prevmag_fft_Framex, double *prevmag_fft_Frameh, double prevcorrtheta_est, double SFxmax, double SFxavg, double flagSFx, int fs, int l, int NFFT,  float ThreadTime, float DurationTime, logMMSE_parameters *logMMSE, bool isEnchanced );



#endif //TWOCHANDOA_V4_TWODOA_H
