clear;clc;

% Importing sound wave %
FS = 48000 ;
[Y, FS] = audioread('eric.wav');

% Plot time domain signal %
T = [(1 : length(Y)) / FS];
figure(1);
subplot(1, 2, 1);
plot(T, Y);
grid on ;
title('Time Domain Signal');

% Play sound wave %
fprintf('Original Sound is play :\n');
sound(Y, FS);

% Wait for the wave %
pause(8)

% Convert the Signal into Frequency Domain %
YF = fftshift(fft(Y));
YF_mag = abs(YF);

% Plot frequency domain signal %
F = linspace(-FS/2 , FS/2 , length(Y));
subplot(1, 2, 2);
plot(F, YF_mag);
grid on ;
title('Frequency Domain Signal');

% Low Pass Filter (Butterworth) %
nOnes  = ceil((8000 * length(Y)) / FS);
nZeros = floor((length(YF) - nOnes) / 2); 
LPF = [ zeros(nZeros,1) ; ones(nOnes,1) ; zeros(length(YF)-nOnes-nZeros,1)];
signal_frequency = YF.*LPF;
Y_3 = real(ifft(ifftshift(signal_frequency)));

% Play the Sound after the Filter %
fprintf('Sound after Filter is play :\n');
sound(Y_3, FS); 

% Plot the Filter output in Time Domain %
T = linspace(0, length(Y_3)/FS, length(Y_3));
figure(2)
plot(T, Y_3);
xlabel('Time'); ylabel('Amplitude'); 
title('The Resived Signal in Time Domain');

% Plot the Filter output in Frequency Domain %
amplitude_signal_frequency = abs(signal_frequency);
phase_signal_frequency = angle(signal_frequency);
ResFreq = linspace(-FS/2, FS/2, length(Y_3));
figure(3)

% Wait for the wave %
pause(8);

% Plot the Magnitude Spectrum %
subplot(1,2,1) 
plot(ResFreq, amplitude_signal_frequency)
xlabel('Frequency'); 
ylabel('Amplitude'); 
title('Magnitude Spectrum');

% Plot the Phase Spectrum %
subplot(1,2,2) 
plot(ResFreq, phase_signal_frequency)
xlabel('Frequency'); 
ylabel('Phase'); 
title('Phase Spectrum');

% Narrow Band FM Generation %
CarrierFrequency=100000;
newFS=5*CarrierFrequency;
Amplitude=15;
KF = 0.1;
newRES = resample(Y_3,500,48);
newT = 0:1/newFS:(1/newFS*length(newRES)-1/newFS);
theta=(2*pi*CarrierFrequency*newT)+(KF.*cumsum(newRES'));
S = Amplitude*cos(theta); 
SF = fftshift(fft(S));
SF_Real = linspace(-newFS/2,newFS/2,length(SF));

% Plotting Narrow Band FM %
figure(4);
subplot(2,1,1); 
plot(newT,S); 
title('NBFM Signal - Time Domain');
subplot(2,1,2); 
plot(SF_Real,abs(SF));
title('NBFM Signal - Frequency Domain');

% Signal Demodulation %
S_diff= diff(S);
Envelop = abs(hilbert(S_diff));
Envelop = Envelop - mean(Envelop);
DeMod_Time = linspace(0,(length(Envelop)/newFS),length(Envelop));
R_env = fftshift(fft(Envelop))./length(Envelop);
n = length(R_env)/newFS;
DC_F = (newFS/2)- 1;
R_env( floor(DC_F*n) : end-floor(DC_F*n) ) = 0;
f = linspace(-newFS/2,newFS/2,length(R_env));
r = real(ifft(ifftshift(R_env))).*length(R_env);

% Plotting Demodulated Signal %
figure(5)
subplot(2,1,1)
plot(DeMod_Time,r)
title('Received Signal - Time Domain')
subplot(2,1,2)
plot(f,abs(R_env)*length(R_env))
title('Received signal - Frequency Domain')
out = resample(r,48,500);
fprintf('Demodulated Sound is play :\n');
sound(out,FS);

% Wait for the wave %
pause(8);

