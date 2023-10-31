clear;close all;clc;
%%Read audio
SNR = 5;
[speech, fs] = audioread('clean speech.wav'); %we have this file in the assignment
[noise, fs2] = audioread('babble noise.wav');%we have this file in the assignment
% [noise, fs2] = audioread('stationary speech-shaped noise.wav');%we have this file in the assignment

N = length(speech);
noise = noise(1:N);%equalize the number of samples

noise = noise./db2mag(SNR)*norm(speech)/norm(noise);%we put the noise level 5 dB below the signal
SNR2 = snr(speech,noise);%Matlab command to check SNR has been aplied (SNR2 has to be equal to SNR)

speech_with_noise = speech+noise; % additive noise



audiowrite('speech_with_noise.wav',speech_with_noise,fs); % creates audio file with original clean signal with noise

%% CHOOSE WIENER OR SPECTRAL SUBTRACTION

method = 'wiener'; %change to 'subtraction' in case you want to try spectral subtraction

%% PSD NOISE

psd_avrg_noise=0;
L_noise = 1024;
overlap_noise = L_noise/2;
K_noise = floor((N - overlap_noise) / (L_noise - overlap_noise));

win = ones(L_noise,1); %rectanguler
% win = hamming(L_noise); %hamming
% win = bartlett(L_noise); %bartlett
% win = blackman(L_noise); %blackman

for i = 0:(K_noise -1)
    
    index_start = (L_noise-overlap_noise)*i+1;
    index_end = index_start+L_noise-1;
    
    frame_noise = noise(index_start:index_end);
    frame_noise_fft = fft(frame_noise).*win;
    
    psd_frame_noise = abs(frame_noise_fft).^2/L_noise/K_noise;
    
    
    psd_avrg_noise = (psd_avrg_noise +  psd_frame_noise);

    
    
end


plot(linspace(0,pi,L_noise/2),10*log10(psd_avrg_noise(1:L_noise/2)));
set(gca,'XTick',0:pi/4:pi) 
set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'})
title("Noise PSD")
hold on


%% NOISE ESTIMATION 1.2 seconds

psd_avrg_noise_es=0;
L_noise_es = 1024;
overlap_noise_es = L_noise_es/2;
K_noise_es = floor((fs*1.2 - overlap_noise_es) / (L_noise_es - overlap_noise_es));

for i = 0:(K_noise_es -1)
    
    index_start = (L_noise_es-overlap_noise_es)*i+1;
    index_end = index_start+L_noise_es-1;
    
    frame_noise_es = speech_with_noise(index_start:index_end);
    frame_noise_fft_es = fft(frame_noise_es).*win;
    
    psd_frame_noise_es = abs(frame_noise_fft_es).^2/L_noise_es/K_noise_es;
    
    
    psd_avrg_noise_es = (psd_avrg_noise_es +  psd_frame_noise_es);
    
    
    
end


plot(linspace(0,pi,L_noise_es/2),10*log10(psd_avrg_noise_es(1:L_noise_es/2)));
set(gca,'XTick',0:pi/4:pi) 
set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'})
title("Noise Estimation")
legend('Noise PSD','Noise Estimation PSD')
hold off

%% SPEECH with noise PSD

L_speech = 1024; %Empirical
%120 ms L = 192

overlap_speech = L_speech/2;
K_speech = floor((N - overlap_speech) / (L_speech - overlap_speech));

psd_avrg_speech = 0;
enhancedSpeech = zeros(N,1);
psd_avrg_noise_wiener = zeros(L_speech,1);
factor = 0.03;

count = 0;
% psd_avrg_noise_es = 0.25.*psd_avrg_noise_es;
for i = 0:(K_speech-1)
    
    index_start = (L_speech-overlap_speech)*i+1;
    index_end = index_start+L_speech-1;

    frame_speech = speech_with_noise(index_start:index_end).*win;

    frame_fft_speech = fft(frame_speech);
    phase = angle(frame_fft_speech);
  
   
    psd_frame_speech = abs(frame_fft_speech).^2./L_speech/K_speech;
    
    psd_avrg_speech = (psd_avrg_speech + psd_frame_speech);
   
    switch method
        case 'wiener'
            gain = psd_frame_speech ./ (psd_frame_speech +psd_avrg_noise_es);
            psd_enhanced = psd_frame_speech .*  gain;
    
        case 'subtraction'
            psd_enhanced =   max(0,psd_frame_speech - factor*psd_avrg_noise_es);

    end
   
   
  
    enhanced_frame =ifft(sqrt(psd_enhanced.*L_speech*K_speech).*exp(phase*1i)) ;
    
  
    
    enhancedSpeech(index_start:index_end) = enhancedSpeech(index_start:index_end) + enhanced_frame;
    
end


%Created audio file with enhanced speech
audiowrite('clean.wav',real(enhancedSpeech),fs);



figure;
subplot(2,2,1)
plot(real(enhancedSpeech));
title("Enhanced Speech")
axis([0 N -max(speech_with_noise) max(speech_with_noise)])

subplot(2,2,4)
plot(speech_with_noise);
title("Speech with noise")
axis([0 N -max(speech_with_noise) max(speech_with_noise)])

subplot(2,2,3)
plot(real(enhancedSpeech)-speech_with_noise);
title("Removed noise")

subplot(2,2,2)
plot(speech);
title("Original Speech")

%% SNR Comparation


power_clean = sum(abs(speech).^2);
power_noise = sum(abs(noise).^2);
power_noise_after = sum(abs(enhancedSpeech-speech).^2);
power_enhanced = sum(abs(enhancedSpeech).^2);


snr_before = 10*log10(power_clean/power_noise);
snr_after = 10*log10(power_enhanced/power_noise_after);


snr_ip = snr_after - snr_before;





