#include "main.h"




int main() {

	/*Initialization for current speech*/
	int l = 206721;
	int n = 1600;
	int fs = 16000;
	int n_F = round(l / n);

	int NFFT = 2048;

	
	/*TwomicDOA function estimates source location*/
	twomicDOA * twoDOA;
	twoDOA = newDOA(NFFT);

	double * x_frame = (double*)malloc(n * sizeof(double)); /* Define memory space for Frame_x */
	double * h_frame = (double*)malloc(n * sizeof(double)); /* Define memory space for Frame_h */
	double * output = (double*)malloc(n * sizeof(double));  /* Define memory space for output */
	int j = 0;
	
	for (int i = 0; i < n_F; i++) {  /*This for loop makes input frame based*/

		for (int k = i*n; k < (i + 1)* n; k++)
		{

			x_frame[k - (i*n)] = x[k];
			h_frame[k - (i*n)] = h[k];

		}
		
		clock_t t;
		t = clock();  /*For timing*/
		twoDOA->doIT(twoDOA, x_frame, h_frame, win, i, twoDOA->prevmag_fft_Framex, twoDOA->prevmag_fft_Frameh, twoDOA->prevcorrtheta_est, twoDOA->SFxmax, twoDOA->SFxavg, twoDOA->flagSFx, fs, n, NFFT);/*TwomicDOA function estimates source location*/
		t = clock() - t;  /*For timing*/

		//printf("It took me %d clicks (%1.9g  seconds)\n", t, (float)t / CLOCKS_PER_SEC);
	output[i] = twoDOA->corrtheta_est;
	}

	


	_getch();


}