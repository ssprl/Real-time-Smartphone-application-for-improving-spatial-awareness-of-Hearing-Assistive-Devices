//
// Created by Abdullah Kucuk on 9/11/2017.
//



#include "twoDOA.h"



float max1(float a, float b)
{
    if(a>b)
        return a;
    else
        return b;
}

twomicDOA*
newDOA(int NFFT, int n, int fs,float ThreadTime,float DurationTime)
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

    newDOA->doIT = twomicDOAv_00;
    newDOA->mic_dist=0.13; // mic distance in meters

    ////////////* Update for memory issue *///////////////////////

    newDOA->lags = (double*)malloc(n * sizeof(double));

    for (int i = 0; i < n-1; i++)
    {
        newDOA->lags[i] = -1*(n/2-1) + i;
        newDOA->lags[i] = newDOA->lags[i] / fs;
        //printf("lags=%1.7g\t i=%d\n", lags[i],i);
    }


    newDOA->X_real= (double*)calloc(NFFT , sizeof(double));
    newDOA->X_img = (double*)calloc(NFFT , sizeof(double));

    newDOA->H_real = (double*)calloc(NFFT , sizeof(double));
    newDOA->H_img = (double*)calloc(NFFT , sizeof(double));


    newDOA->mag_fft_Framex= (double*)calloc((NFFT/2+1) , sizeof(double));
    newDOA->sfx = (double*)calloc((NFFT / 2 + 1) , sizeof(double));
    newDOA->mag_fft_Frameh = (double*)calloc((NFFT / 2 + 1) , sizeof(double));
    newDOA->sfh = (double*)calloc((NFFT / 2 + 1) , sizeof(double));


    newDOA->x_frame_2 = (double*)calloc((NFFT) , sizeof(double));
    newDOA->h_frame_2 = (double*)calloc((NFFT) , sizeof(double));


    newDOA->fft_Framex_real = (double*)calloc((NFFT) , sizeof(double));
    newDOA->fft_Framex_img = (double*)calloc((NFFT) , sizeof(double));
    newDOA->fft_Frameh_real = (double*)calloc((NFFT) , sizeof(double));
    newDOA->fft_Frameh_img = (double*)calloc((NFFT) , sizeof(double));

    newDOA->currmag_fft_Framex = (double*)calloc((NFFT / 2 + 1) , sizeof(double));
    newDOA->currmag_fft_Frameh = (double*)calloc((NFFT / 2 + 1) , sizeof(double));

    newDOA->R12_real = (double*)malloc((NFFT ) * sizeof(double));
    newDOA->R12_img = (double*)malloc((NFFT) * sizeof(double));

    newDOA->r12_temp = (double*)malloc((NFFT) * sizeof(double));
    newDOA->r12= (double*)malloc((n) * sizeof(double));
    newDOA->r12_abs = (double*)malloc((n) * sizeof(double));

    newDOA->windd = (double*)malloc((n) * sizeof(double));
    newDOA->windd = Hamming(n);

    newDOA->Trainingframe =  (int) ThreadTime;
    newDOA->DurationThresh = (int) DurationTime;

    newDOA->n=n;
    newDOA->Ncorr = n-1;

    newDOA->deltau = 0;
    newDOA->idx=0;

    newDOA->SFx=0;
    newDOA->sum_sfx=0;

    newDOA->corrtheta_est1 = 0;


    newDOA->trans1 = newTransform(NFFT);
    newDOA->trans2 = newTransform(NFFT);
    newDOA->trans3 = newTransform(NFFT);
    ////////////////////////////////////////

    return newDOA;
}


void fftshift(double *x,int NFFT,double *output)
{
    //int npoint = NFFT/2;
    //double *output = (double*)malloc((NFFT) * sizeof(double));
    int k;
    for (k = 0; k < NFFT/2; k++) {
        output[k + NFFT/2] = x[k];
        output[k] = x[k + NFFT/2];
    }

}//*/


int findmax(double *arr3, int n) { // finds max index of array index that bigger than sca


    double tmp_max = -999999;
    int tmp;
    for (int ii = 0; ii < n-1; ii++) {


        if (tmp_max < arr3[ii]  ) {
            tmp_max = arr3[ii];
            tmp = ii;
            //tmp_size++;
            //printf("max = %f,   arr = %f,    index = %d", tmp_maxtmp);
        }





    }
    int temp_idx = tmp;
    return temp_idx;

}

void twomicDOAv_00(twomicDOA *twoDOA,double *x, double *h, int framecounter, double *prevmag_fft_Framex, double *prevmag_fft_Frameh,double prevcorrtheta_est,double SFxmax, double SFxavg,double flagSFx, int fs, int l,int NFFT, float ThreadTime , float DurationTime, logMMSE_parameters *logMMSE , bool isEnchanced) {


  


    if(isEnchanced) {
        



        if (framecounter < twoDOA->Trainingframe) {

            for (int kk = 0; kk < logMMSE->n; kk++) {
                logMMSE->insign_x[kk] = twoDOA->windd[kk] * x[kk];
                logMMSE->insign_h[kk] = twoDOA->windd[kk] * h[kk];


            }
            /*for (int ii = logMMSE->n; ii < logMMSE->NFFT; ii++) {
                insign_x[ii] = 0;
                insign_h[ii] = 0;

            }*/
            logMMSE->enhance_signal(logMMSE, logMMSE->insign_x,logMMSE->insign_h, framecounter, twoDOA->Trainingframe);



        }
        else {

            if (framecounter == twoDOA->Trainingframe) {
                for (int ii = 0; ii < logMMSE->n; ii++) {
                    logMMSE->insign_x[ii] = twoDOA->windd[ii]* x[ii];
                    logMMSE->prev_frame[ii] = x[ii];

                    logMMSE->insign_h[ii] =  twoDOA->windd[ii] * h[ii];
                    logMMSE->prev_frame2[ii] = h[ii];

                }
                /* for (int ii = logMMSE->n; ii < logMMSE->NFFT; ii++) {
                     insign_x[ii] = 0;
                     insign_h[ii] = 0;

                 }*/
                logMMSE->enhance_signal(logMMSE, logMMSE->insign_x, logMMSE->insign_h, framecounter, twoDOA->Trainingframe);
                for (int ii = 0; ii < logMMSE->len1; ii++) {
                    x[ii] = logMMSE->x_final[ii];
                    h[ii] = logMMSE->y_final[ii];
                }




            }
            else if (framecounter == twoDOA->Trainingframe + 1) {


                for (int ii = 0; ii < logMMSE->n; ii++) {
                    if (ii < logMMSE->len2) {
                        logMMSE->insign_x[ii] =
                                twoDOA->windd[ii] * logMMSE->prev_frame[ii + logMMSE->len2];
                        logMMSE->insign_h[ii] =
                                twoDOA->windd[ii] * logMMSE->prev_frame2[ii + logMMSE->len2];

                    }
                    else {
                        logMMSE->insign_x[ii] = twoDOA->windd[ii] * x[ii - logMMSE->len2];
                        logMMSE->insign_h[ii] = twoDOA->windd[ii] * h[ii - logMMSE->len2];

                    }
                    logMMSE->prev_frame[ii] = x[ii];
                    logMMSE->prev_frame2[ii] = h[ii];

                }


                logMMSE->enhance_signal(logMMSE, logMMSE->insign_x, logMMSE->insign_h, framecounter, twoDOA->Trainingframe);
                for (int ii = 0; ii < logMMSE->len1; ii++) {
                    x[ii] = logMMSE->x_final[ii];
                    h[ii] = logMMSE->y_final[ii];
                }

            }
            else {

                /* if (framecounter==25) {
                     int ake = 0;
                     ake=+1;
                 }*/


                for (int ii = 0; ii < logMMSE->n; ii++) {
                    logMMSE->insign_x[ii] =  twoDOA->windd[ii] * logMMSE->prev_frame[ii];
                    logMMSE->insign_h[ii] =  twoDOA->windd[ii] * logMMSE->prev_frame2[ii];
                    logMMSE->prev_x[ii]=x[ii];
                    logMMSE->prev_h[ii]=h[ii];

                }





                logMMSE->enhance_signal(logMMSE, logMMSE->insign_x, logMMSE->insign_h, framecounter, twoDOA->Trainingframe);




                for (int ii = 0; ii < logMMSE->len1; ii++) {
                    x[ii] = logMMSE->x_final[ii];
                    h[ii] = logMMSE->y_final[ii];

                }




                for (int ii = 0; ii < logMMSE->n; ii++) {
                    if (ii < logMMSE->len2) {
                        logMMSE->insign_x[ii] =
                                twoDOA->windd[ii] * logMMSE->prev_frame[ii + logMMSE->len2];
                        logMMSE->insign_h[ii] =
                                twoDOA->windd[ii] * logMMSE->prev_frame2[ii + logMMSE->len2];
                    }
                    else {
                        logMMSE->insign_x[ii] = twoDOA->windd[ii] * x[ii - logMMSE->len2];
                        logMMSE->insign_h[ii] = twoDOA->windd[ii] * h[ii - logMMSE->len2];

                    }
                    logMMSE->prev_frame[ii] = logMMSE->prev_x[ii];
                    logMMSE->prev_frame2[ii] = logMMSE->prev_h[ii];

                }
               // clock_t t;
                //t = clock();

                logMMSE->enhance_signal(logMMSE, logMMSE->insign_x, logMMSE->insign_h, framecounter, twoDOA->Trainingframe);

               /* t = clock() - t;
                double ake_timing=(((double)t)/ CLOCKS_PER_SEC)+(logMMSE->avg_timing*(framecounter-Trainingframe-2));

                logMMSE->avg_timing=ake_timing/(framecounter-Trainingframe-1);
                __android_log_print(ANDROID_LOG_INFO,"Time", "Time1= %1.9g",((float)logMMSE->avg_timing));// */

                for (int ii = logMMSE->len1; ii < 2*logMMSE->len1; ii++) {
                    x[ii] = logMMSE->x_final[ii-logMMSE->len1];
                    h[ii] = logMMSE->y_final[ii-logMMSE->len1];

                }



            }


        }


    }

   
    twoDOA->corrtheta_est = 0;


    for (int i = 0; i < twoDOA->n; i++)
    {

        twoDOA->x_frame_2[i] = twoDOA->windd[i] * x[i];
        twoDOA->h_frame_2[i] = twoDOA->windd[i] * h[i];
        //printf("h_frame_2=%1.7g\t i=%d\n", h_frame_2[i], i);

    }

   

    twoDOA->trans2->doTransform(twoDOA->trans2, twoDOA->h_frame_2);
    twoDOA->trans1->doTransform(twoDOA->trans1, twoDOA->x_frame_2);



    for (int j = 0; j < NFFT; j++)
    {
        twoDOA->X_real[j] = twoDOA->trans1->real[j];
        twoDOA->X_img[j] = twoDOA->trans1->imaginary[j];
        twoDOA->H_real[j] = twoDOA->trans2->real[j];
        twoDOA->H_img[j] = twoDOA->trans2->imaginary[j];


        twoDOA->fft_Framex_real[j] = twoDOA->X_real[j] / twoDOA->n;
        twoDOA->fft_Framex_img[j] = twoDOA->X_img[j] / twoDOA->n;

        twoDOA->fft_Frameh_real[j] = twoDOA->H_real[j] / twoDOA->n;
        twoDOA->fft_Frameh_img[j] = twoDOA->H_img[j] / twoDOA->n;

        //printf("x_frame_2_real=%1.7g\t i=%d\n", trans2->real[j], j);
    }


   
    twoDOA->sum_sfx=0;
    //double sum_sfh=0;

    for (int i = 0; i < NFFT/2+1; i++)
    {

        twoDOA->mag_fft_Framex[i] = sqrt(pow(twoDOA->fft_Framex_real[i],2)+ pow(twoDOA->fft_Framex_img[i], 2));
        twoDOA->mag_fft_Frameh[i] = sqrt(pow(twoDOA->fft_Frameh_real[i], 2) + pow(twoDOA->fft_Frameh_img[i], 2));

        twoDOA->currmag_fft_Framex[i] =  twoDOA->mag_fft_Framex[i];
        twoDOA->currmag_fft_Frameh[i] =  twoDOA->mag_fft_Frameh[i];

        twoDOA->sfx[i] = twoDOA->currmag_fft_Framex[i] - prevmag_fft_Framex[i];
        twoDOA->sfh[i] = twoDOA->currmag_fft_Frameh[i] - prevmag_fft_Frameh[i];

        twoDOA->sum_sfx += pow( twoDOA->sfx[i],2);
        //sum_sfh += pow( twoDOA->sfh[i], 2);

        //printf("mag_fft_Framex_real=%1.7g\t i=%d\n", mag_fft_Frameh[i], i);

        twoDOA->prevmag_fft_Framex[i] = twoDOA->currmag_fft_Framex[i];
        twoDOA->prevmag_fft_Frameh[i] = twoDOA->currmag_fft_Frameh[i];

    }

   // double SFx;// SFh;
    twoDOA->SFx = twoDOA->sum_sfx / (NFFT);
    //SFh = sum_sfh / (n_F*NFFT);

 
    for (int i = 0; i < NFFT; i++)
    {
        twoDOA->R12_real[i] = twoDOA->X_real[i] * twoDOA->H_real[i] - ((-1 * twoDOA->H_img[i]) * twoDOA->X_img[i]);
        twoDOA->R12_img[i] = twoDOA->X_real[i] * (-1 * twoDOA->H_img[i]) + (twoDOA->X_img[i] * twoDOA->H_real[i]);
        //printf("R12_img=%1.7g\t i=%d\n", R12_img[i], i);

    }
    /*Transform *trans3;
    trans3 = newTransform(NFFT);//*/
    twoDOA->trans3->invTransform(twoDOA->trans3, twoDOA->R12_real, twoDOA->R12_img);

    fftshift(twoDOA->trans3->real,NFFT,twoDOA->r12_temp);




    for (int i = (NFFT / 2 + 1 - (twoDOA->Ncorr - 1) / 2)-1; i < NFFT / 2 + 1 + (twoDOA->Ncorr - 1) / 2; i++)
    {
        twoDOA->r12[i-(NFFT / 2 + 1 - (twoDOA->Ncorr - 1) / 2-1)] = twoDOA->r12_temp[i];
        twoDOA->r12_abs[i - (NFFT / 2 + 1 - (twoDOA->Ncorr - 1) / 2 - 1)] = fabs(twoDOA->r12[i - (NFFT / 2 + 1 - (twoDOA->Ncorr - 1) / 2 - 1)]);

    }
      twoDOA->idx = findmax(twoDOA->r12_abs,twoDOA->n);

    //int temp_idx = idx;
    //
    twoDOA->deltau = twoDOA->lags[twoDOA->idx];


    twoDOA->degree = twoDOA->deltau * 343 / (twoDOA->mic_dist);
    // printf("i=%1.9g\n", degree);
    if (twoDOA->degree>1){
        twoDOA->theta_est = 0;}
    else if (twoDOA->degree<-1){
        twoDOA->theta_est = PI * vall;}
    else{
        twoDOA->theta_est = acos(twoDOA->degree) * vall;}



    /*%% Setting SF Threshold*/


    if (framecounter==0) {
        SFxavg = twoDOA->SFx;
        SFxmax = twoDOA->SFx;
    }
    else {

        if (framecounter<=twoDOA->Trainingframe) {
            if (twoDOA->SFx > SFxmax) {

                SFxmax = twoDOA->SFx;

            }
            SFxavg = (SFxavg*(framecounter-1 ) + twoDOA->SFx) / framecounter;
        }
    }

    if ((framecounter >= twoDOA->Trainingframe) && (twoDOA->SFx > SFxavg) && (twoDOA->SFx > SFxmax)) {

        flagSFx +=   1;										/*   % Duration*/
    }
    else {
        flagSFx = 0;
    }


    if (framecounter >= twoDOA->Trainingframe) {
        if (flagSFx > twoDOA->DurationThresh)
            twoDOA->corrtheta_est1 = twoDOA->theta_est;

        else
            twoDOA->corrtheta_est1 = prevcorrtheta_est;


    }
    else {

        twoDOA->corrtheta_est1 = 0;
    }

    //printf("corrtheta_est=%1.7g\t \n", corrtheta_est);
    twoDOA->corrtheta_est = twoDOA->corrtheta_est1;
    twoDOA->prevcorrtheta_est = twoDOA->corrtheta_est1;

    twoDOA->SFxmax = SFxmax;

    twoDOA->SFxavg = SFxavg;

    twoDOA->flagSFx = flagSFx;

  

}

Transform*
newTransform(int points)
{
    Transform* newTransform = (Transform*)malloc(sizeof(Transform));

    newTransform->points = points;
    newTransform->real = (double*)malloc(points * sizeof(double));
    newTransform->imaginary = (double*)malloc(points * sizeof(double));
    newTransform->sine = NULL;
    newTransform->cosine = NULL;
    newTransform->doTransform = FFT;
    newTransform->invTransform = IFFT;
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
    int i, j, k, L, m, n, o, p, q;
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

    int i, j, k, L, m, n, o, p, q;
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
destroyTransform(Transform* transform)
{
    if (transform != NULL) {
        if ((transform)->cosine != NULL) {
            free((transform)->cosine);
            (transform)->cosine = NULL;
        }
        if ((transform)->sine != NULL) {
            free((transform)->sine);
            (transform)->sine = NULL;
        }
        if ((transform)->real != NULL) {
            free((transform)->real);
            (transform)->real = NULL;
        }
        if ((transform)->imaginary != NULL) {
            free((transform)->imaginary);
            (transform)->imaginary = NULL;
        }
        free(transform);
        transform = NULL;
    }
}

//////////////////////*mod_logMMSE*////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
//  mod_logMMSE.c
//
//
//  Created by Abdullah Kucuk on 08/30/2017.
//  Copyright  2017 Abdullah Kucuk. All rights reserved.
*/

#include "twoDOA.h"
//double *Hamming( int np);
double aa, mu, eta, ksi_min;
//int k;


logMMSE_parameters*
newparameters(int NFFT, int n, int PERC)
{
    //int nn = sizeof(ip);
    logMMSE_parameters* newparameters = (logMMSE_parameters*)malloc(sizeof(logMMSE_parameters));
    newparameters->PERC = PERC; /*% window overlap in percent of frame size*/
    newparameters->len1 = floor(n*PERC / 100);
    newparameters->len2 = n - newparameters->len1;

    newparameters->NFFT = NFFT;
    newparameters->n = n;
    newparameters->ksi = (double*)malloc((NFFT) * sizeof(double));

    newparameters->Xk_prev = (double*)malloc((NFFT) * sizeof(double));

    newparameters->noise_mu2 = (double*)malloc((NFFT) * sizeof(double));


    newparameters->noise_mean = (double*)malloc((NFFT) * sizeof(double));
    newparameters->sig1 = (double*)malloc((NFFT) * sizeof(double));
    newparameters->spec_real = (double*)malloc((NFFT) * sizeof(double));
    newparameters->spec_imag = (double*)malloc((NFFT) * sizeof(double));
    newparameters->spec_real2 = (double*)malloc((NFFT) * sizeof(double));
    newparameters->spec_imag2 = (double*)malloc((NFFT) * sizeof(double));

    newparameters->windowed_signal = (double*)malloc((NFFT) * sizeof(double));
    newparameters->windowed_signal2 = (double*)malloc((NFFT) * sizeof(double));

    newparameters->wind = (double*)malloc((n) * sizeof(double));

    newparameters->wind=Hamming(n);
    /*for (int a = 0; a < n; a++) {

        printf("newparameters->wind[%d]=%.17f\n", a, newparameters->wind[a]);

    }*/
    newparameters->sys_wind = (double*)malloc((newparameters->len1) * sizeof(double));
    newparameters->x_old = (double*)malloc((newparameters->len1) * sizeof(double));
    newparameters->x_final = (double*)malloc((newparameters->len1) * sizeof(double));
    newparameters->y_old = (double*)malloc((newparameters->len1) * sizeof(double));
    newparameters->y_final = (double*)malloc((newparameters->len1) * sizeof(double));

    newparameters->prev_frame = (double*)malloc((newparameters->n) * sizeof(double));
    newparameters->prev_frame2 = (double*)malloc((newparameters->n) * sizeof(double));

    newparameters->output = (double*)malloc((NFFT) * sizeof(double));
    newparameters->output2 = (double*)malloc((NFFT) * sizeof(double));
    newparameters->avg_timing=0;
    for (int ii = 0; ii < newparameters->len1; ii++)
    {
        newparameters->sys_wind[ii] = 1 / (newparameters->wind[ii] + newparameters->wind[ii + newparameters->len1]);
        newparameters->x_old[ii] = 0;
        newparameters->x_final[ii] = 0;
        newparameters->y_old[ii] = 0;
        newparameters->y_final[ii] = 0;

    }

    for (int k = 0; k <newparameters->NFFT; k++)
    {
        newparameters->noise_mean[k] = 0;


    }

    newparameters->enhance_signal = enhance_signal;

    ////////////////* Memory Update *///////////////////

    newparameters->insign_x = (double *) calloc(NFFT, sizeof(double));
    newparameters->insign_h = (double *) calloc(NFFT, sizeof(double));
    newparameters->prev_x = (double *) calloc(NFFT, sizeof(double));
    newparameters->prev_h = (double *) calloc(NFFT, sizeof(double));

    newparameters->sig2 = (double*)malloc((NFFT) * sizeof(double));
    newparameters->sigforgama = (double*)malloc((NFFT) * sizeof(double));
    newparameters->gammak = (double*)malloc((NFFT) * sizeof(double));
    newparameters->maxgammak = (double*)malloc((NFFT) * sizeof(double));
    newparameters->log_sigma_k = (double*)malloc((NFFT) * sizeof(double));//*/

    newparameters->A = (double*)malloc((NFFT) * sizeof(double));
    newparameters->vk = (double*)calloc((NFFT) , sizeof(double));
    newparameters->ei_vk = (double*)calloc((NFFT) , sizeof(double));
    newparameters->hw = (double*)calloc((NFFT), sizeof(double));
    newparameters->spec_hw_real = (double*)calloc((NFFT) , sizeof(double));
    newparameters->spec_hw_imag = (double*)calloc((NFFT) , sizeof(double));
    newparameters->spec_hw_real2 = (double*)calloc((NFFT) , sizeof(double));
    newparameters->spec_hw_imag2 = (double*)calloc((NFFT) , sizeof(double));
    newparameters->trans33 = newTransform(NFFT);
    newparameters->trans = newTransform(NFFT);
    newparameters->trans_y = newTransform(NFFT);
    newparameters->trans22 = newTransform(NFFT);
    newparameters->trans22_y = newTransform(NFFT);


    //*/


    ////////////////////////////////////////////////////


    return newparameters;
}


void  minn(double *x, int size, double *y, double min) { /*If array has bigger value than 40, this function set this value as 40*/
   // double min = 40;

    for (int i = 0; i < size; i++)
    {
        if (min < x[i])
            y[i] = min;
        else
            y[i] = x[i];


    }

}

void maxx(double *x, int size, double maxvalue, double value, double *y ) { /*This function sets array values to 0 if given value less than 0 */


    for (int i = 0; i < size; i++)
    {
        if (maxvalue >(x[i] - value))
            y[i] = maxvalue;
        else
            y[i] = x[i] - value;

    }


}

double summ(double *x, int size) {

    double sum = 0;
    for (int i = 0; i < size; i++)
    {
        sum += x[i];
    }

    return sum;

}
double expint_new(int n, double x)
{
    int i, ii, nm1;
    double a, b, c, d, del, fact, h, psi;
    double ans=0;

    nm1 = n - 1;
    if (n < 0 || x < 0.0 || (x == 0.0 && (n == 0 || n == 1)))
        ans = 0;
    else {
        if (n == 0)
            ans = exp(-x) / x;
        else {
            if (x == 0.0)
                ans = 1.0 / nm1;

            else {
                if (x > 1.0) {
                    b = x + n;
                    c = 1.0 / FPMIN;
                    d = 1.0 / b;
                    h = d;
                    for (i = 1; i <= MAXIT; i++) {
                        a = -i*(nm1 + i);
                        b += 2.0;
                        d = 1.0 / (a*d + b);
                        c = b + a / c;
                        del = c*d;
                        h *= del;
                        if (fabs(del - 1.0) < EPSS) {
                            ans = h*exp(-x);
                            return ans;
                        }
                    }
                }
                else {
                    ans = (nm1 != 0 ? 1.0 / nm1 : -log(x) - EULER);
                    fact = 1.0;
                    for (i = 1; i <= MAXIT; i++) {
                        fact *= -x / i;
                        if (i != nm1) del = -fact / (i - nm1);
                        else {
                            psi = -EULER;
                            for (ii = 1; ii <= nm1; ii++) psi += 1.0 / ii;
                            del = fact*(-log(x) + psi);
                        }
                        ans += del;
                        if (fabs(del) < fabs(ans)*EPSS)
                            return ans;
                    }
                }
            }
        }
    }
    return ans;
}
double *Hamming( int np) {// Hamming definition is changed according to Hamming.m file
    double *y = (double*)calloc(np, sizeof(double));
    int a;
    for (a = 0; a < np; a++) {
        y[a] = (0.54 - 0.46*cos((2 * pi*a) / (np - 1)));
        //printf("y[%d]=%.17f\n", a, y[a]);

    }

    return y;
    free(y);
    y=NULL;
}

void enhance_signal(logMMSE_parameters *logMMSE, double *insign, double *insign2, int frame_counter, int Trainingframe) {



    /*for (int kk = 0; kk <logMMSE->NFFT - 800; kk++)
    {
    //logMMSE->windowed_signal[kk] = logMMSE->wind[kk] * sig[kk];
    printf("noise_mu2[%d]=%.17g\n", kk + 1, logMMSE->noise_mu2[kk]);

    }*/



    if (frame_counter<Trainingframe) {
       // Transform* trans33;
       // trans33 = newTransform(logMMSE->NFFT);
        logMMSE->trans33->doTransform(logMMSE->trans33, insign);
        for (int k = 0; k <logMMSE->NFFT; k++)
        {
            logMMSE->noise_mean[k] += sqrt(pow(logMMSE->trans33->real[k], 2) + pow(logMMSE->trans33->imaginary[k], 2)); /*UPDATE MAGNITUDE LEFT AND RIGHT*/
            logMMSE->noise_mu2[k] = pow(logMMSE->noise_mean[k] / Trainingframe, 2);


        }

       // destroyTransform(trans33);
        //trans33=NULL;

    }
    else {



        // k = 0;
        aa = 0.98f;
        mu = 0.98f;
        eta = 0.15f;
        ksi_min = pow(10, -25 / 10);



        logMMSE->trans->doTransform(logMMSE->trans, insign);

        //Transform *trans_y;
        //logMMSE->trans_y = newTransform(logMMSE->NFFT);
        logMMSE->trans_y->doTransform(logMMSE->trans_y, insign2); //spectrum of y

        for (int kk = 0; kk < logMMSE->NFFT; kk++)
        {
            logMMSE->sig1[kk] = sqrt(pow(logMMSE->trans->real[kk], 2) + pow(logMMSE->trans->imaginary[kk], 2));

            logMMSE->spec_real[kk] = logMMSE->trans->real[kk];
            logMMSE->spec_imag[kk] = logMMSE->trans->imaginary[kk];

            logMMSE->spec_real2[kk] = logMMSE->trans_y->real[kk];
            logMMSE->spec_imag2[kk] = logMMSE->trans_y->imaginary[kk];

            logMMSE->sig2[kk] = pow(logMMSE->sig1[kk], 2);
            logMMSE->sigforgama[kk] = logMMSE->sig2[kk] / logMMSE->noise_mu2[kk];

        }

        // UPDATED as follows//logMMSE->gammak = minn(logMMSE->sigforgama, logMMSE->NFFT);/*% limit post SNR to avoid overflows*/

        minn(logMMSE->sigforgama, logMMSE->NFFT, logMMSE->gammak ,40);/*% limit post SNR to avoid overflows*/

        //double * maxgammak = (double*)malloc((logMMSE->NFFT) * sizeof(double));



        // UPDATED as follows//logMMSE->maxgammak = maxx(logMMSE->gammak, logMMSE->NFFT, 0, 1);
        maxx(logMMSE->gammak, logMMSE->NFFT, 0, 1,logMMSE->maxgammak );

        if (frame_counter == Trainingframe) {
            for (int l = 0; l < logMMSE->NFFT; l++)
            {
                logMMSE->ksi[l] = aa + (1 - aa)*logMMSE->maxgammak[l];

            }
        }
        else {

            for (int lx = 0; lx < logMMSE->NFFT; lx++)
            {
                logMMSE->ksi[lx] = aa*logMMSE->Xk_prev[lx] / logMMSE->noise_mu2[lx] + (1 - aa)*logMMSE->maxgammak[lx];/* % a priori SNR*/

            }

            // UPDATED as follows//logMMSE->ksi = maxx(logMMSE->ksi, logMMSE->NFFT, ksi_min, 0); /*% limit ksi to -25 dB*/
                maxx(logMMSE->ksi, logMMSE->NFFT, ksi_min, 0,logMMSE->ksi ); /*% limit ksi to -25 dB*/
        }


        //double * log_sigma_k = (double*)malloc((logMMSE->NFFT) * sizeof(double));

        for (int ll = 0; ll < logMMSE->NFFT; ll++)
        {
            logMMSE->log_sigma_k[ll] = logMMSE->gammak[ll] * logMMSE->ksi[ll] / (1 + logMMSE->ksi[ll]) - log(1 + logMMSE->ksi[ll]);

        }


        logMMSE->vad_decisionn = summ(logMMSE->log_sigma_k, logMMSE->NFFT) / logMMSE->n;


        if (logMMSE->vad_decisionn < eta)
        {
            // double *up_noise_mu2 = (double*)malloc((logMMSE->NFFT) * sizeof(double));
            for (int l = 0; l < logMMSE->NFFT; l++)
            {
                /* up_noise_mu2[l] = mu*logMMSE->noise_mu2[l] + (1 - mu)*sig2[l];
                 logMMSE->noise_mu2[l] = up_noise_mu2[l];*/
                logMMSE->noise_mu2[l]=mu*logMMSE->noise_mu2[l] + (1 - mu)*logMMSE->sig2[l];

            }
            // free(up_noise_mu2);
        }




        for (int ll = 0; ll < logMMSE->NFFT; ll++)
        {
            logMMSE->A[ll] = logMMSE->ksi[ll] / (1 + logMMSE->ksi[ll]);
            logMMSE->vk[ll] = logMMSE->A[ll] * logMMSE->gammak[ll];
            logMMSE->ei_vk[ll] = 0.5*expint_new(1, logMMSE->vk[ll]);
            logMMSE->hw[ll] = logMMSE->A[ll] * exp(logMMSE->ei_vk[ll]);
            logMMSE->sig1[ll] = logMMSE->sig1[ll] *logMMSE->hw[ll];
            logMMSE->Xk_prev[ll] = pow(logMMSE->sig1[ll], 2);

            logMMSE->spec_hw_real[ll] = logMMSE->spec_real[ll] * logMMSE->hw[ll];
            logMMSE->spec_hw_imag[ll] = logMMSE->spec_imag[ll] * logMMSE->hw[ll];

            logMMSE->spec_hw_real2[ll] = logMMSE->spec_real2[ll] * logMMSE->hw[ll];
            logMMSE->spec_hw_imag2[ll] = logMMSE->spec_imag2[ll] * logMMSE->hw[ll];
        }
        //clock_t t;
        //t = clock();


        logMMSE->trans22->invTransform(logMMSE->trans22, logMMSE->spec_hw_real, logMMSE->spec_hw_imag);


        logMMSE->trans22_y->invTransform(logMMSE->trans22_y, logMMSE->spec_hw_real2, logMMSE->spec_hw_imag2);
        /*t = clock() - t;
        ake_timing=(((double)t)/ CLOCKS_PER_SEC)+ake_timing;
        __android_log_print(ANDROID_LOG_INFO,"Time", "Time1= %1.9g",ake_timing/(frame_counter-10));// */


        for (int ii = 0; ii < logMMSE->n; ii++)
        {
            logMMSE->output[ii] = logMMSE->trans22->real[ii];
            logMMSE->output2[ii] = logMMSE->trans22_y->real[ii];

        }


        for (int kk = 0; kk < logMMSE->len1; kk++)
        {

            logMMSE->x_final[kk] = (logMMSE->x_old[kk] + logMMSE->output[kk])*logMMSE->sys_wind[kk];
            logMMSE->x_old[kk] = logMMSE->output[kk + logMMSE->len1];

            logMMSE->y_final[kk] = (logMMSE->y_old[kk] + logMMSE->output2[kk])*logMMSE->sys_wind[kk];
            logMMSE->y_old[kk] = logMMSE->output2[kk + logMMSE->len1];

            //printf("logMMSE->x_final[%d]=%1.7g\n", kk+1,logMMSE->x_final[kk]);
        }


    

    }

}
