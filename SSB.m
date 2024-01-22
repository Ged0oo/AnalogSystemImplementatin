clc; clear;

% constants
PI = 3.14;

% STEP 1
% -------- read in the audio file and get its spectrum -------------------------
[audio, Fs] = audioread("Eric.wav");

%player = audioplayer(audio, Fs);
%play(player);
%pause(8);

audio_spectrum = fftshift(fft(audio));

t = linspace(0, length(audio)./Fs, length(audio));
fvec = linspace(-Fs/2, Fs/2, length(audio_spectrum));

figure;
plot(t, audio);
title("audio in time domain");

figure;
plot(fvec, abs(audio_spectrum));
title("audio in frequency domain");

% STEP 2
% -------- apply LPF at Fcutoff = 4KHz -----------------------------------------
Fcutoff = 4e3;
audio_spectrum(abs(fvec) > Fcutoff) = 0;

figure;
plot(fvec, abs(audio_spectrum));
title("audio spectrum after applying LPF @ Fcutoff=4KHz");

% STEP 3
% -------- obtain audio in time domain after LPF -------------------------------
audio = ifft(ifftshift(audio_spectrum));

%player = audioplayer(audio, Fs);
%play(player);
%pause(8);

% STEP 4
% -------- generate a DSB-SC signal with Fs = Fcarrier * 5, Fc = 100KHz --------
Fcarrier = 1e5;
Fs_new = Fcarrier * 5;

audio = resample(audio, Fs_new, Fs);

t = linspace(0, length(audio)./Fs_new, length(audio));

carrier = cos(2 * PI * Fcarrier * t)';

DSB_SC_signal = carrier .* audio;
DSB_SC_spectrum = fftshift(fft(DSB_SC_signal));

fvec = linspace(-Fs_new/2, Fs_new/2, length(DSB_SC_spectrum));

figure;
plot(t, real(DSB_SC_signal));
title("DSB-SC signal");

figure;
plot(fvec, abs(DSB_SC_spectrum));
title("DSB-SC spectrum");

% STEP 5 (first time - using ideal LPF)
% -------- obtain SSB (the LSB of the DSB-SC signal) ---------------------------
SSB_SC_spectrum = DSB_SC_spectrum;
SSB_SC_spectrum(abs(fvec) > Fcarrier) = 0;

figure;
plot(fvec, abs(SSB_SC_spectrum));
title("SSB-SC (LSB) spectrum");

% STEP 6 (first time - using ideal LPF)
% -------- use coherent detection to retrieve the msg from the SSB-SC msg ------
SSB_SC_signal = real(ifft(ifftshift(SSB_SC_spectrum)));
received_signal = carrier .* SSB_SC_signal;
%received_signal = resample(received_signal, Fs, Fs_new);

%player = audioplayer(received_signal, Fs);
%play(player);
pause(8);
received_spectrum = fftshift(fft(received_signal));

t = linspace(0, length(received_signal)./Fs, length(received_signal));

figure;
plot(t, received_signal);
title("received signal");

figure;
plot(fvec, abs(received_spectrum));
title("received spectrum");

% STEP 7
% STEP 5 (second time - using 4th order Butterworth filter as LPF) -------------
[b, a] = butter(2, [96000/(5*Fcarrier/2) 100000/(5*Fcarrier/2)], 'bandpass');
ssb_signal = filter(b, a, DSB_SC_signal);

% % Plot the spectrum of the SSB signal
Butter_LSB_spectrum = fftshift(fft(ssb_signal));
figure;
plot(fvec, abs(Butter_LSB_spectrum));
title('Butterworth LSB spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

received_signal = carrier .* ssb_signal;

% Plot the spectrum of the SSB modulated signal
spectrum_ssb_modulated = fftshift(fft(received_signal));
figure;
plot(fvec, abs(spectrum_ssb_modulated));
title('Spectrum of SSB Modulated Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Low-pass filtering to extract the baseband signal
cutoff_frequency = 4000;
N = length(spectrum_ssb_modulated);
n = N / Fs;
num_zeros  = floor((Fs / 2 - cutoff_frequency) * n);

frequencies = linspace(-Fs_new / 2, Fs_new / 2, N);
spectrum_ssb_modulated( [1:num_zeros    (N - num_zeros + 1):N]) = 0;

% Plot the spectrum of the filtered SSB signal
figure;
subplot(2,1,2);
plot(frequencies, abs(spectrum_ssb_modulated));
title('Spectrum of Filtered SSB Modulated Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Time Domain Plot of Filtered SSB Signal
ssb_modulated = real(ifft(ifftshift(spectrum_ssb_modulated)));

subplot(2,1,1);
plot(t, real(ssb_modulated));
title('Filtered SSB Signal in Time Domain');
xlabel('Time (seconds)');
ylabel('Amplitude');

% Resample the filtered SSB signal to play the sound
ssb_resampled = resample(ssb_modulated, Fs, Fs_new);

% Play the resampled SSB signal
sound(ssb_resampled, Fs);
pause(8);

% Fnq = Fs_new / 2;
% Wn = [Fcarrier - Fcutoff, Fcarrier] ./ Fnq;
% [b, a] = butter(4, Wn, 'bandpass');
% Butter_LSB_spectrum = filter(b, a, DSB_SC_spectrum);
% 
% fvec = linspace(-Fs_new/2, Fs_new/2, length(Butter_LSB_spectrum));
% 
% figure;
% plot(fvec, abs(Butter_LSB_spectrum));
% title("Butterworth LSB spectrum");
% 
% % STEP 6 (second time - using coherent detection after Butterworth filter) -----
% Butter_LSB_signal = real(ifft(ifftshift(Butter_LSB_spectrum)));
% received_signal = carrier .* Butter_LSB_signal;
% received_signal = resample(received_signal, Fs, Fs_new);
% 
% t = linspace(0, length(received_signal)./Fs, length(received_signal));
% 
% figure;
% plot(t, received_signal);
% title("Butterworth signal after coherent detection");

% STEP 8
% -------- receive the SSB-SC signal but with noise (SNR=0, 10, 30) ------------
SNRs = [0 10 30];

t = linspace(0, length(SSB_SC_signal)/Fs_new, length(SSB_SC_signal));
carrier = 2 * cos(2 * PI * Fcarrier * t)';

for i=1:1:length(SNRs)
  SNR = SNRs(i);
  sig_with_noise = awgn(SSB_SC_signal, SNR, 'measured');
  received_signal = sig_with_noise .* carrier;
  received_signal = resample(received_signal, Fs, Fs_new);
  received_spectrum = fftshift(fft(received_signal));

  t = linspace(0, length(received_signal)./Fs, length(received_signal));
  fvec = linspace(-Fs/2, Fs/2, length(received_spectrum));

  title1 = sprintf("Received Message At SNR=%d (time domain)", SNR);
  title2 = sprintf("Received Message At SNR=%d (frequency domain)", SNR);

  player = audioplayer(received_signal, Fs);
  play(player);
  pause(8); 
  
  figure;
  subplot(1, 2, 1);
  plot(t, received_signal);
  title(title1);
  subplot(1, 2, 2);
  plot(fvec, abs(received_spectrum));
  title(title2);
end

% STEP 9
% -------- generate an SSB-TC signal and demodulate it -------------------------
% ##### generation #####
[audio Fs] = audioread("Eric.wav");

Fs_new = 5 * Fcarrier;

audio = resample(audio, Fs_new, Fs);

t = linspace(0, length(audio)./Fs_new, length(audio));

carrier = cos(2 * PI * Fcarrier * t)';

A = 2 * max(audio);

SSB_TC_signal = (A + audio) .* carrier;
SSB_TC_spectrum = fftshift(fft(SSB_TC_signal));
fvec = linspace(-Fs_new/2, Fs_new/2, length(SSB_TC_spectrum));

SSB_TC_spectrum(abs(fvec) > Fcarrier) = 0;
SSB_TC_signal = real(ifft(ifftshift(SSB_TC_spectrum)));

% ##### demodulation #####
received_signal = SSB_TC_signal .* carrier .* 2;
received_signal = resample(received_signal, Fs, Fs_new);
received_spectrum = fftshift(fft(received_signal));

player = audioplayer(received_signal, Fs);
play(player);
pause(8);

t = linspace(0, length(received_signal)./Fs, length(received_signal));
fvec = linspace(-Fs/2, Fs/2, length(received_spectrum));

figure;
plot(t, received_signal);
title("received signal after SSB-TC");

figure;
plot(fvec, abs(received_spectrum));
title("received spectrum after SSB-TC");
