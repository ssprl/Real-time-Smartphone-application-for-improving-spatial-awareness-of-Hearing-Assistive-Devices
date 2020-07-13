% This function determines the ground truth according to the clean speech
function [GT,time_ind] = GTSpeech(s,Thres,n,fs,SpeechDuration)





s_hat = s;


thres = 10^(Thres/10);

s_hat(abs(s_hat)<10^-(thres))=eps;
s_sil(abs(s_hat)>10^-(thres))=eps;



l = length(s);

n_F = round(l/n);
GT=zeros(1,n_F);
for i = 1 : n_F
    F_start(i) = (i-1)*n+1;
    F_end(i) = min(i*n,l);
    VADmeasures(i,1) = F_start(i);
    VADmeasures(i,2) = F_end(i);
 
    Frame = s(F_start(i):F_end(i));
    RMS(i) = rms(Frame);
    if ((RMS(i) > (thres))&&(0.68*fs<F_start(i))&&(2.1*fs>F_end(i)))
        GT(i) = 0;
         idx_0(i)=i;
    elseif ((RMS(i) > (thres))&&(4.2*fs<F_start(i))&&(5.6*fs>F_end(i)))
        GT(i) = 45;
         idx_45(i)=i;
    elseif ((RMS(i) > (thres))&&(7.7*fs<F_start(i))&&(9.09*fs>F_end(i)))
        GT(i) = 90;
         idx_90(i)=i;
    elseif ((RMS(i) > (thres))&&(11.2*fs<F_start(i))&&(12.5*fs>F_end(i)))
        GT(i) = 135;
         idx_135(i)=i;
        
    elseif ((RMS(i) > (thres))&&(14.75*fs<F_start(i))&&(16.08*fs>F_end(i)))
        GT(i) = 180;  
         idx_180(i)=i;
    end
    
    
    
end

time_ind(1)=idx_0(min(find(idx_0>0)));
time_ind(2)=max(idx_0);

time_ind(3)=idx_45(min(find(idx_45>0)));
time_ind(4)=max(idx_45);

time_ind(5)=idx_90(min(find(idx_90>0)));
time_ind(6)=max(idx_90);

time_ind(7)=idx_135(min(find(idx_135>0)));
time_ind(8)=max(idx_135);

time_ind(9)=idx_180(min(find(idx_180>0)));
time_ind(10)=max(idx_180);
end
