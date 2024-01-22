clc; clear;

% Loading the sound file
[y, samplingFrequency] = audioread('eric.wav');
StartTime = 0; 
EndTime = length(y)/samplingFrequency;
duration = EndTime-StartTime;

% Plotting the sound file in time domain
t = linspace(StartTime,EndTime,length(y));
figure (1);
plot(t,y);
title('Signal in time domain ');
xlabel('Time '); 
ylabel('Amplitude ');

% Playing the sound file
fprintf('Original Sound is play :\n');
sound (y, samplingFrequency);
pause (duration) ;

% Converting the sound signal to the frequency domain
YFrequency = fftshift(fft(y));

% Plotting the spectrum of the sound signal
f = linspace(-samplingFrequency/2, samplingFrequency/2, length(YFrequency));
figure (2);
subplot (2,1,1);
plot (f, abs(YFrequency)); 
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Original Frequency Spectrum');

subplot (2,1,2);
plot (f, angle(YFrequency)); 
title ('Phase spectrum ');
xlabel ('Frequency '); 
ylabel('Phase ');

% Eliminating frequencies above 4KHz using a low pass filter

% Low Pass Filter (Butterworth) %
nOnes  = ceil((8000 * length(YFrequency)) / samplingFrequency);
nZeros = floor((length(YFrequency) - nOnes) / 2); 
LPF = [ zeros(nZeros,1) ; ones(nOnes,1) ; zeros(length(YFrequency)-nOnes-nZeros,1)];
YFiltered = YFrequency.*LPF;

% Plotting the spectrum of filtered signal
F = linspace(-samplingFrequency/2, samplingFrequency/2, length(YFiltered));
figure (3);
subplot (2,1,1);
plot (F, abs(YFiltered));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Filtered Frequency Spectrum');

subplot (2,1,2);
plot (F, angle(YFiltered));
title ('Phase spectrum of filtered signal ');
xlabel ('Frequency '); 
ylabel('Phase ');

% Converting the filtered spectrum to time domain
y_filtered = real(ifft(ifftshift(YFiltered)));

% Playing the filtered sound signal
duration = length(y_filtered)/samplingFrequency;
fprintf('Filtered Sound is play :\n');
sound(y_filtered, samplingFrequency);
pause(duration) ;

%& Plotting the filtered signal in time domain
tl = linspace (0, length(y_filtered)/samplingFrequency, length(y_filtered));
figure(4); 
plot(tl, y_filtered);
title('Filtered signal in time domain ');
xlabel('Time '); 
ylabel('Amplitude ');


% Generating the message
CarrierFrequency = 100000; 
samplingFrequency = 48000;
newSamplingFreq = 5*CarrierFrequency;
m = resample(y_filtered, newSamplingFreq, samplingFrequency);
T = linspace (0, length(m)/newSamplingFreq, length(m)) ;
% Carrier signal
CarrierSignal = cos(2*pi*CarrierFrequency*T) ;
% Modulating the message using DSB-TC
Amplitude = 2*max(m) ;
xDSBTC = (Amplitude+m)' .* CarrierSignal;
% DSB-TC signal in frequency domain
sDSBTC = fftshift(fft(xDSBTC));
% Plotting the spectrum of DSB-TC signal
F_1 = linspace(-newSamplingFreq/2, newSamplingFreq/2, length(sDSBTC)) ;
figure (5); 
plot(F_1,abs(sDSBTC));
title('DSB-TC Modulated Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
% Modulating the message using DSB-SC
xDSBSC = m' .* CarrierSignal;
% DSB-SC in frequency domain
sDSBSC = fftshift(fft(xDSBSC));
% Plotting the spectrum of DSB-SC signal
F_2 = linspace(-newSamplingFreq/2, newSamplingFreq/2, length(sDSBSC));
figure (6);
plot(F_2, abs(sDSBSC));
title('DSB-SC Modulated Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Detection of DSB-TC signal using envelope detector
envelopDSBTC = abs(hilbert(xDSBTC));
% Plotting the detected envelope in time domain
figure (7); 
plot(T,envelopDSBTC) ;
title('Recieved signal after DSB-TC modulation ');
xlabel('Time' ); 
ylabel('Amplitude ');
% Decreasing the sampling frequency
Fs=0.2*newSamplingFreq;
% Generating the recieved sound and playing it
recSound = resample(envelopDSBTC, Fs, newSamplingFreq) ;
duration = length(recSound)/Fs;
fprintf('Recived Sound after DSB-TC Envelope is play :\n');
sound(recSound, Fs);
pause(duration) ;
% Detection of DSB-SC signal using envelope detector
envelopDSBSC = abs(hilbert(xDSBSC));
% Plotting the recieved envelope in time domain
figure(8); 
plot(T,envelopDSBSC) ;
title('Recieved signal after DSB-SC modulation ');
xlabel('Time '); 
ylabel('Amplitude ');
% Generating the recieved sound and playing it
recSound_2 = resample (envelopDSBSC,Fs,newSamplingFreq) ;
duration = length(recSound_2) /Fs;
fprintf('Recived Sound after DSB-SC Envelope is play :\n');
sound(recSound_2, Fs);
pause(duration) ;

% Use coherent detection for DSB-SC with different SNR values
SNR_values = [0, 10, 30]; % Specify SNR values in dB

for snr = SNR_values
	noise = awgn(xDSBTC, snr);
	r1 = 2*noise.*CarrierSignal;
	
	% Signal after coherent detection in frequency domain
	R = fftshift(fft(r1));
	
	% Low pass filter with cutoff frequency=4KHz
	Fc = 4000;
	No_1 = ceil(length(R)*2*Fc/newSamplingFreq);
	No_0 = floor((length(R)-No_1)/2);
	LPF = [zeros(No_0,1); ones(No_1, 1); zeros(No_0,1)];
	
	% Spectrum of recieved signal after the LPF
	R_1 = R' .* LPF;
	F_1 = linspace (-newSamplingFreq/2, newSamplingFreq/2, length(R));
	figure;
	subplot(2,1,1);
	plot(F_1, abs(R_1)) ;
	title(['Recieved signal for SNR = ' num2str(snr) ' dB in frequency domain']);
	xlabel('Frequency '); 
	ylabel('Magnitude ');
  
	% Recieved signal in time domain
	R_s1 = real(ifft(ifftshift(R_1)));
	R_p1 = real(ifft(ifftshift(R)));
	t_s1 = linspace (0, length(R_s1)/newSamplingFreq,length(R_s1));
	subplot(2,1,2);
	plot(t_s1, R_s1);
	title(['Recieved signal for SNR = ' num2str(snr) ' dB in time domain']);
	xlabel('Time '); 
	ylabel('Amplitude ');
	
	% Generating the recieved sound and playing it
	recSound = resample(R_p1, Fs, newSamplingFreq);
	duration = length(recSound)/Fs;
	fprintf('Recived Sound after Noise Addition is play :\n');
	sound(recSound, Fs);
	pause(duration) ;
end


% Coherent detection with a frequency error of 100Hz (Beet effect)
C_coh = cos(2*pi*(CarrierFrequency+100)*T);
noise = awgn(xDSBTC, 30);
r = noise .* C_coh;

% Spectrum of recieved signal after coherent detection
R_4 = fftshift(fft(r));

% Low pass filter with cutoff frequency=4KHz
Fc = 4000;
No_1 = ceil(length(R_4)*2*Fc/newSamplingFreq);
No_0 = floor((length(R_4)-No_1)/2);
LPF = [zeros(No_0,1); ones(No_1, 1); zeros(No_0,1)];
	
% Recieved signal after LPF
R4 = R_4' .* LPF;

f4 = linspace(-newSamplingFreq/2, newSamplingFreq/2, length(R4));
figure (12);
subplot (2,1,1);
plot (f4, abs(R4));
title('Recieved signal with frequency error in frequency domain ');
xlabel('Frequency ');
ylabel(' Magnitude ');

% Recieved signal in time domain
R_s4 = real(ifft(ifftshift(R4)));
R_p4 = real(ifft(ifftshift(R_4)));
t_s4 = linspace (0, length(R_s4)/newSamplingFreq, length(R_s4));
subplot (2,1,2);
plot (t_s4,R_s4);
title ('Recieved signal with frequency error in time domain ');
xlabel('Time '); 
ylabel('Amplitude ');

% Generating the recieved sound and playing it
recSound = resample (R_p4, Fs, newSamplingFreq);
duration = length (recSound) /Fs;
fprintf('Recived Sound after 100MHz Frequency Error is play :\n');
sound(recSound, Fs);
pause(duration) ;


% Coherent detection with phase error of 20 degrees (Attenuation)
C_coh = cos(2*pi*CarrierFrequency*T+20*pi/180);
noise = awgn(xDSBTC, 30);
r0 = noise .* C_coh;

% Spectrum of recieved signal after coherent detection
R_5=fftshift(fft(r0));


% Low pass filter with cutoff frequency=4KHz
Fc = 4000;
No_1 = ceil(length(R_5)*2*Fc/newSamplingFreq);
No_0 = floor((length(R_5)-No_1)/2);
LPF = [zeros(No_0,1); ones(No_1, 1); zeros(No_0,1)];

% Recieved signal after the LPF
R5 = R_5' .* LPF;

f5 = linspace (-newSamplingFreq/2, newSamplingFreq/2, length(R5));
figure (13);
subplot (2,1,1);
plot(f5, abs(R5));
title('Recieved signal with phase error in frequency domain ');
xlabel('Frequency '); 
ylabel(' Magnitude ');

% Recieved signal in time domain
R_s5 = real(ifft(ifftshift(R5)));
R_p5 = real(ifft(ifftshift(R_5)));
t_s5 = linspace (0, length (R_s5)/newSamplingFreq, length(R_s5));
subplot (2,1,2);
plot(t_s5, R_s5);
title('Recieved signal with phase error in time domain ');
xlabel('Time '); 
ylabel('Amplitude ');

% Generating the recieved sound and playing it
recSound = resample (R_p5, Fs,newSamplingFreq) ;
duration = length (recSound) /Fs;
fprintf('Recived Sound after 20Degree Phase Error is play :\n');
sound(recSound, Fs);
pause(duration) ;