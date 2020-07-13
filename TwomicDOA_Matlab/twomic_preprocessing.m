
%%
% twomic_preprocessing.m
% Created by Anshuman Ganguly and Abdullah Kucuk on 4/28/2017
% 
% DOA
% non realtime
% with VAD
%
%Inputs are:
% 1- MergedDOA_HINT.wav for clean speech
% 2- NoisyBabble15dB.wav for noisy speeech
%
% Outputs are:
% 1- RMSEVAD;
% 2- RMSENOVAD;
% 3- RMSEGCC_PHAT;
% 4- RMSEGCC_SCOT;
% 5- RMSEGCC_ECKART;
% 6- RMSEGCC_ML;
%%
clear all;
clc;
close all;

[sig,fs]=audioread('MergedDOA_HINT.wav'); % Read clean speech data



%% Upsampling to 48kHz
 newfs = 48000;
  sig = resample(sig,newfs,fs); % Upsample the clean speech from 16kHz to 48kHz
%%

 [noisysig,fs]=audioread('NoisyBabble15dB.wav');% Read noisy speech data

 noisysig = resample(noisysig,newfs,fs);% Upsample the noisy speech from 16kHz to 48kHz
%% 
 
      x = noisysig(:,1); h = noisysig(:,2); % x is first channel input; h is second channel input

x=x(0.3*newfs:end);
h=h(0.3*newfs:end);

l=length(x);  % length of input signal


%%

n = round(0.1*newfs); % Frame length is 100 ms 

N = n;
Ncorr = round((N))-1;
NFFT = 2^nextpow2(n);  % finds NFFT point according to Frame Length
lags = (-(Ncorr-1)/2:(Ncorr-1)/2).';
lags = lags/newfs;

n_F = round(l/n);    % Number of frame 
win = hamming(n);    % Creates hamming window

x_frame = zeros(n,n_F); % Initialization for Frame_x
h_frame = zeros(n,n_F); % Initialization for Frame_h

SFx = zeros(1,n_F);SFh = zeros(1,n_F);  % Initialization for Spectral Flux 
SFxcnt = zeros(1,n_F);                  % Initialization for Spectral Flux count

theta_est = zeros(1,n_F);                % Initialization for estimated angle without VAD

corrtheta_est = zeros(1,n_F);            % Initialization for estimated angle with VAD
deltau = zeros(1,n_F);
meanSFx = 0; maxSFx =0;
meanSFh = 0; maxSFh =0;
Trainingframe = 5;               % Training frames for threshold calculation
DurationThresh = 1;              % Duration thresold for smooth transition
cnt = 0;


[GT_F,time_ind] = GTspeech(sig,-40,n,newfs,3);  % This function determines the ground truth according to the clean speech
% time_ind=round(time_ind*newfs/n);               % Time indexes shows where the speech is.

for i = 1:n_F                                   % This for loop makes input frame based
    i;
    tic;
    x_frame = win.*x((i-1)*n+1:i*n);            % Frame based input signal multiply by window 
    h_frame = win.*h((i-1)*n+1:i*n);            % Frame based input signal multiply by window 
   
   
    
    
    X = fft(x_frame,NFFT);                      % Takes FFT of windowed signal
    H = fft(h_frame,NFFT);                      % Takes FFT of windowed signal
    
    fft_Framex = fft(x_frame,NFFT)/n;
    mag_fft_Framex(:,i) = abs(fft_Framex(1:NFFT/2+1));  % Magnitude of spectrum of input
    fft_Frameh = fft(h_frame,NFFT)/n;
    mag_fft_Frameh(:,i) = abs(fft_Frameh(1:NFFT/2+1)); % Magnitude of spectrum of input
%% Implementing VAD
   % Calculation of Spectral flux(SF)
    if i == 1
        sfx = 0;            % For first frame, Spectral Flux is zero.
        sfh = 0;            % For first frame, Spectral Flux is zero.
        %mag_fft_Frame(:,i) = 0;
        %sf =  mag_fft_Frame(:,i) ;
    else
        sfx = (mag_fft_Framex(:,i))-(mag_fft_Framex(:,i-1));  % Spectral Flux Calculation
        sfh = (mag_fft_Frameh(:,i))-(mag_fft_Frameh(:,i-1));  % Spectral Flux Calculation
    end
    SFx(i) = sum(sfx.^2)/(n_F*NFFT);
    SFh(i) = sum(sfh.^2)/(n_F*NFFT);

%% GCC
    
 R12 = X.*conj(H);               % Cross Correlation calculation

R12_hat = R12;
r12_temp = fftshift(ifft(R12_hat),1);
r12 = r12_temp(NFFT/2+1-(Ncorr-1)/2:NFFT/2+1+(Ncorr-1)/2,:);

[~,idx(i)] = max(abs(r12));
deltau(i) = lags(idx(i));
theta_est(i) = real((acos(deltau(i)*343/(0.13)))*180/pi); % Estimated angle for GCC

%% PHAT processor
    R12_PHAT = R12./(abs(R12)+eps);     % Calculation of PHAT processor
 R12_hat_PHAT = R12_PHAT;
 r12_temp_PHAT = fftshift(ifft(R12_hat_PHAT),1);
 r12_PHAT = r12_temp_PHAT(NFFT/2+1-(Ncorr-1)/2:NFFT/2+1+(Ncorr-1)/2,:);

 [~,idx_PHAT(i)] = max(abs(r12_PHAT));
 deltau_PHAT(i) = lags(idx_PHAT(i));
 theta_est_PHAT(i) = real((acos(deltau_PHAT(i)*343/(0.13)))*180/pi);    % Estimated angle for GCC_PHAT processor


    
%% SCOT(Smooth COherence Transform) processor
GS = 1./sqrt((X.*conj(X)).*(H.*conj(H)));
R12_SCOT = R12.*GS;                       % Calculation of SCOT processor

R12_hat_SCOT = R12_SCOT;
r12_temp_SCOT = fftshift(ifft(R12_hat_SCOT),1);
r12_SCOT = r12_temp_SCOT(NFFT/2+1-(Ncorr-1)/2:NFFT/2+1+(Ncorr-1)/2,:);

[~,idx_SCOT(i)] = max(abs(r12_SCOT));
deltau_SCOT(i) = lags(idx_SCOT(i));
theta_est_SCOT(i) = real((acos(deltau_SCOT(i)*343/(0.13)))*180/pi);  % Estimated angle for GCC_SCOT processor

%% Eckart filter

% if ec == 1
    GX = (X.*conj(X)); GH = (H.*conj(H));
    GE = abs(R12).*(GX - abs(R12)).*(GH-abs(R12));
    R12_hat_Eckart = GE.*R12;               % Calculation of SCOT processor
r12_temp_Eckart = fftshift(ifft(R12_hat_Eckart),1);

r12_Eckart = r12_temp_Eckart(NFFT/2+1-(Ncorr-1)/2:NFFT/2+1+(Ncorr-1)/2,:);

[~,idx_Eckart(i)] = max(abs(r12_Eckart));
deltau_Eckart(i) = lags(idx_Eckart(i));
theta_est_ECKART(i) = real((acos(deltau_Eckart(i)*343/(0.13)))*180/pi); % Estimated angle for GCC_ECKART processor
%% Hannan Thomson(HT) filter
GXX = (X.*conj(X)); GHH = (H.*conj(H));
gamma = ((R12))./(eps+sqrt((GXX).*(GHH)));
GHT = (((abs(gamma).^2)./(eps+(1-(abs(gamma).^2)))))./(eps+abs(R12));
R12_hat_HT = GHT.*R12;                     % Calculation of ML processor

r12_temp_HT = fftshift(ifft(R12_hat_HT),1);
r12_HT = r12_temp_HT(NFFT/2+1-(Ncorr-1)/2:NFFT/2+1+(Ncorr-1)/2,:);

[~,idx_HT(i)] = max(abs(r12_HT));
deltau_HT(i) = lags(idx_HT(i));
theta_est_ML(i) = real((acos(deltau_HT(i)*343/(0.13)))*180/pi);    % Estimated angle for GCC_ML processor


%% Finding direction of speech only based on SF


    % Setting SF Threshold
    if i == 1
        SFxavg(i) = SFx(i);
        SFxmax= SFx(i);
    else
        if i <= Trainingframe
            if SFx(i) > SFxmax
                SFxmax = SFx(i);
            end
            
            SFxavg = (SFxavg*(i-1)+SFx(i))/i;   % Spectral Flux Calculation for VAD
        end
        
    end
    
    if (i>=Trainingframe) && ((SFx(i) > SFxavg)&& (SFx(i)> SFxmax))
      
        flagSFx = flagSFx+1;                            % Duration
    else
        flagSFx = 0;
    end
    
    if (i>=Trainingframe) 
        if (flagSFx > DurationThresh)
            corrtheta_est(i) = theta_est(i);
            SFxcnt(i) = 1; 
        else 
            corrtheta_est(i) = corrtheta_est(i-1);  % Estimated angle with VAD
        end
    else
        corrtheta_est(i) = 0;                       % Estimated angle with VAD
    end
    


timer(i) = toc;

end
SF(1,:) =SFx; SF(2,:) = SFh; 

RMSENOVAD1=sqrt(mean((GT_F(time_ind(1):time_ind(2))-theta_est(time_ind(1):time_ind(2))).^2));   % Calculation of RMSE for 0 degree without VAD
RMSENOVAD2=sqrt(mean((GT_F(time_ind(3):time_ind(4))-theta_est(time_ind(3):time_ind(4))).^2));   % Calculation of RMSE for 45 degree without VAD
RMSENOVAD3=sqrt(mean((GT_F(time_ind(5):time_ind(6))-theta_est(time_ind(5):time_ind(6))).^2));   % Calculation of RMSE for 90 degree without VAD
RMSENOVAD4=sqrt(mean((GT_F(time_ind(7):time_ind(8))-theta_est(time_ind(7):time_ind(8))).^2));   % Calculation of RMSE for 135 degree without VAD
RMSENOVAD5=sqrt(mean((GT_F(time_ind(9):time_ind(10))-theta_est(time_ind(9):time_ind(10))).^2)); % Calculation of RMSE for 180 degree without VAD
 RMSENOVAD=(RMSENOVAD1+RMSENOVAD2+RMSENOVAD3+RMSENOVAD4+RMSENOVAD5)/5;                          % Taking average of RMSE

RMSEVAD1=sqrt(mean((GT_F(time_ind(1):time_ind(2))-corrtheta_est(time_ind(1):time_ind(2))).^2));  % Calculation of RMSE for 0 degree with VAD
RMSEVAD2=sqrt(mean((GT_F(time_ind(3):time_ind(4))-corrtheta_est(time_ind(3):time_ind(4))).^2));  % Calculation of RMSE for 45 degree with VAD
RMSEVAD3=sqrt(mean((GT_F(time_ind(5):time_ind(6))-corrtheta_est(time_ind(5):time_ind(6))).^2));  % Calculation of RMSE for 90 degree with VAD
RMSEVAD4=sqrt(mean((GT_F(time_ind(7):time_ind(8))-corrtheta_est(time_ind(7):time_ind(8))).^2));  % Calculation of RMSE for 135 degree with VAD
RMSEVAD5=sqrt(mean((GT_F(time_ind(9):time_ind(10))-corrtheta_est(time_ind(9):time_ind(10))).^2)); % Calculation of RMSE for 180 degree with VAD

RMSEVAD=(RMSEVAD1+RMSEVAD2+RMSEVAD3+RMSEVAD4+RMSEVAD5)/5;                                           % Taking average of RMSE






figure;subplot(6,1,1);plot([1:length(x)]/newfs,x);title('Noisy Signal'); xlabel('Time(sec)'); ylabel('');grid on;     %Plot the noisy signal
subplot(6,1,2);plot([1:length(theta_est)]*n/newfs,[SF] ); title('Spectral Flux'); xlabel('Time(sec)'); ylabel('');grid on;  %Plot Spectral Flux
subplot(6,1,3);plot([1:length(theta_est)]*n/newfs,[SFx-SFh] ); title(' Diff. between Spectral Flux bet. Two microphones '); xlabel('Time(sec)'); ylabel('');grid on; %Plot Spectral Flux Difference between Two microphones
subplot(6,1,4);plot([1:length(theta_est)]*n/newfs, [theta_est]); title('Estimated Angle without VAD '); xlabel('Time(sec)'); ylabel('Angle(Degree)');grid on;        %Plot Estimated angle without VAD 
subplot(6,1,5);plot([1:length(theta_est)]*n/newfs, [GT_F(1:length(theta_est))]);title('Ground Truth '); xlabel('Time(sec)'); ylabel('Angle(Degree)');grid on;        %Plot Ground Truth
subplot(6,1,6);plot([1:length(theta_est)]*n/newfs, [corrtheta_est]);title('Corrected Angle with VAD '); xlabel('Time(sec)'); ylabel('Angle(Degree)');grid on;        %Plot Estimated angle with VAD



RMSEGCC_PHAT1=sqrt(mean((GT_F(time_ind(1):time_ind(2))-theta_est_PHAT(time_ind(1):time_ind(2))).^2));  % Calculation of RMSE for 0 degree for PHAT processor
RMSEGCC_PHAT2=sqrt(mean((GT_F(time_ind(3):time_ind(4))-theta_est_PHAT(time_ind(3):time_ind(4))).^2));  % Calculation of RMSE for 45 degree for PHAT processor
RMSEGCC_PHAT3=sqrt(mean((GT_F(time_ind(5):time_ind(6))-theta_est_PHAT(time_ind(5):time_ind(6))).^2));  % Calculation of RMSE for 90 degree for PHAT processor
RMSEGCC_PHAT4=sqrt(mean((GT_F(time_ind(7):time_ind(8))-theta_est_PHAT(time_ind(7):time_ind(8))).^2));  % Calculation of RMSE for 135 degree for PHAT processor
RMSEGCC_PHAT5=sqrt(mean((GT_F(time_ind(9):time_ind(10))-theta_est_PHAT(time_ind(9):time_ind(10))).^2)); % Calculation of RMSE for 180 degree for PHAT processor
 RMSEGCC_PHAT=(RMSEGCC_PHAT1+RMSEGCC_PHAT2+RMSEGCC_PHAT3+RMSEGCC_PHAT4+RMSEGCC_PHAT5)/5;                % Taking average of RMSE
 
 RMSEGCC_SCOTT1=sqrt(mean((GT_F(time_ind(1):time_ind(2))-theta_est_SCOT(time_ind(1):time_ind(2))).^2)); % Calculation of RMSE for 0 degree for SCOT processor
RMSEGCC_SCOTT2=sqrt(mean((GT_F(time_ind(3):time_ind(4))-theta_est_SCOT(time_ind(3):time_ind(4))).^2));  % Calculation of RMSE for 45 degree for SCOT processor
RMSEGCC_SCOTT3=sqrt(mean((GT_F(time_ind(5):time_ind(6))-theta_est_SCOT(time_ind(5):time_ind(6))).^2));  % Calculation of RMSE for 90 degree for SCOT processor
RMSEGCC_SCOTT4=sqrt(mean((GT_F(time_ind(7):time_ind(8))-theta_est_SCOT(time_ind(7):time_ind(8))).^2));  % Calculation of RMSE for 135 degree for SCOT processor
RMSEGCC_SCOTT5=sqrt(mean((GT_F(time_ind(9):time_ind(10))-theta_est_SCOT(time_ind(9):time_ind(10))).^2)); % Calculation of RMSE for 180 degree for SCOT processor
 RMSEGCC_SCOT=(RMSEGCC_SCOTT1+RMSEGCC_SCOTT2+RMSEGCC_SCOTT3+RMSEGCC_SCOTT4+RMSEGCC_SCOTT5)/5;             % Taking average of RMSE

 RMSEGCC_ECKART1=sqrt(mean((GT_F(time_ind(1):time_ind(2))-theta_est_ECKART(time_ind(1):time_ind(2))).^2)); % Calculation of RMSE for 0 degree for ECKART processor
RMSEGCC_ECKART2=sqrt(mean((GT_F(time_ind(3):time_ind(4))-theta_est_ECKART(time_ind(3):time_ind(4))).^2));  % Calculation of RMSE for 45 degree for ECKART processor
RMSEGCC_ECKART3=sqrt(mean((GT_F(time_ind(5):time_ind(6))-theta_est_ECKART(time_ind(5):time_ind(6))).^2));  % Calculation of RMSE for 90 degree for ECKART processor
RMSEGCC_ECKART4=sqrt(mean((GT_F(time_ind(7):time_ind(8))-theta_est_ECKART(time_ind(7):time_ind(8))).^2));  % Calculation of RMSE for 135 degree for ECKART processor
RMSEGCC_ECKART5=sqrt(mean((GT_F(time_ind(9):time_ind(10))-theta_est_ECKART(time_ind(9):time_ind(10))).^2)); % Calculation of RMSE for 180 degree for ECKART processor
 RMSEGCC_ECKART=(RMSEGCC_ECKART1+RMSEGCC_ECKART2+RMSEGCC_ECKART3+RMSEGCC_ECKART4+RMSEGCC_ECKART5)/5;       % Taking average of RMSE
 
 
  RMSEGCC_ML1=sqrt(mean((GT_F(time_ind(1):time_ind(2))-theta_est_ML(time_ind(1):time_ind(2))).^2));  % Calculation of RMSE for 0 degree for ML processor
RMSEGCC_ML2=sqrt(mean((GT_F(time_ind(3):time_ind(4))-theta_est_ML(time_ind(3):time_ind(4))).^2));    % Calculation of RMSE for 45 degree for ML processor
RMSEGCC_ML3=sqrt(mean((GT_F(time_ind(5):time_ind(6))-theta_est_ML(time_ind(5):time_ind(6))).^2));    % Calculation of RMSE for 90 degree for ML processor
RMSEGCC_ML4=sqrt(mean((GT_F(time_ind(7):time_ind(8))-theta_est_ML(time_ind(7):time_ind(8))).^2));    % Calculation of RMSE for 135 degree for ML processor
RMSEGCC_ML5=sqrt(mean((GT_F(time_ind(9):time_ind(10))-theta_est_ML(time_ind(9):time_ind(10))).^2));  % Calculation of RMSE for 180 degree for ML processor
 RMSEGCC_ML=(RMSEGCC_ML1+RMSEGCC_ML2+RMSEGCC_ML3+RMSEGCC_ML4+RMSEGCC_ML5)/5;                         % Taking average of RMSE






RMSEVAD;
RMSENOVAD;
RMSEGCC_PHAT;
RMSEGCC_SCOT;
 RMSEGCC_ECKART;
 RMSEGCC_ML;

