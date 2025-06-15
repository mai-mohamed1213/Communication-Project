clc;
clear;
close all;
%%%%%%%%%%%%%%%%%%%%%%Time and frequency parametres %%%%%%%%%%%%%%%%%%
fs = 100;
T = 100;
df = 0.01;

dt = 1/fs;
N = ceil(T/dt);
t = -(N*dt/2) : dt : ((N*dt/2) - dt);
t = t(1:N);

f = -(0.5*fs) : df : (0.5*fs - df);
f = f(1:N);

%%%%%%%%%%%%%%%%%%%%%%%%% 1- Plot X(t) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_time = zeros(size(t));
for i = 1:length(t)
    if t(i) >= -4 && t(i) < 0
        x_time(i) = t(i) + 5;
    elseif t(i) >= 0 && t(i) <= 4
        x_time(i) = -t(i) + 5;
    end
end

figure(1);
plot(t, x_time, 'r', 'LineWidth', 1.5);
title('Function X(t)');
xlabel('Time (s)');
ylabel('X(t)');
grid on;
xlim([-10 10]);

%%%%%%%%%%%%%%%%%%%%%%%%% 2- Analytical expression for X(f) %%%%%%%%%%%%%%%%%
x_analytical_freq = 16 * (sinc(4 * f)).^2 + 8 * sinc(8 * f);

%%%%%%%%%%%%%%%%%%%%%%%%% 3- Plot X(f) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_freq = fftshift(fft(x_time) * dt);

figure(2);
plot(f, abs(x_freq), 'r', 'LineWidth', 1.5);
hold on;
plot(f, abs(x_analytical_freq), 'b', 'LineWidth', 1.5);
title('Magnitude Spectrum |X(f)|');
legend('X', 'X (Analytical)');
xlabel('Frequency (Hz)');
ylabel('|X(f)|');
xlim([-10 10]);
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%% 4- Estimate Bandwidth of X %%%%%%%%%%%%%%%%%%%%
psd_x = abs(x_freq).^2; % Power Spectral Density

figure(3);
plot(f, psd_x, 'r');
title('Power Spectrum of X(t)');
xlabel('Frequency (Hz)');
ylabel('Power');
grid on;

threshold = 0.05 * max(psd_x);
indices = find(psd_x >= threshold);
BW_x = max(abs(f(indices)));

fprintf('X Estimated Bandwidth = %1.3f Hz\n', BW_x);

% Analytical bandwidth
psd_analytical_x = abs(x_analytical_freq).^2;
threshold = 0.05 * max(psd_analytical_x);
indices = find(psd_analytical_x >= threshold);
BW_analytical_x = max(abs(f(indices)));

fprintf('X Analytical Estimated Bandwidth = %1.3f Hz\n', BW_analytical_x);

%%%%%%%%%%%%%%%%%%%%%%%%% 5- Low Pass Filter (BW = 1 Hz) %%%%%%%%%%%%%%%%%%%
LPF_BW1 = 1;
H = abs(f) < LPF_BW1;
x_filtered1 = ifft(ifftshift(H .* x_freq)) / dt;

figure(4);
plot(t, x_time, 'r', 'LineWidth', 1.5);
hold on;
plot(t, real(x_filtered1), 'b', 'LineWidth', 1.5);
title('X(t) Before and After LPF (BW = 1 Hz)');
legend('Original X(t)', 'Filtered X(t)');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([-10 10]);
ylim([0 5]);
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%% 6- Low Pass Filter (BW = 0.3 Hz) %%%%%%%%%%%%%%%%%
LPF_BW2 = 0.3;
H = abs(f) < LPF_BW2;
x_filtered2 = ifft(ifftshift(H .* x_freq)) / dt;

figure(5);
plot(t, x_time, 'r', 'LineWidth', 1.5);
hold on;
plot(t, real(x_filtered2), 'b', 'LineWidth', 1.5);
title('X(t) Before and After LPF (BW = 0.3 Hz)');
legend('Original X(t)', 'Filtered X(t)');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([-10 10]);
ylim([0 5]);
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%% 7- Repeat for m(t) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a- Define m(t)
m_time = cos(2 * pi * 0.5 * t);
m_time(t > 4) = 0;
m_time(t < 0) = 0;

figure(6);
plot(t, m_time);
title('Function M(t)');
xlabel('Time (s)');
ylabel('M(t)');
grid on;
xlim([-10 10]);

% b- Analytical expression for M(f)
m_analytical_freq = 2 * sinc(4 * (f - 0.5)) .* exp(-1i * 2 * pi * 2 * (f - 0.5))
+ 2 * sinc(4 * (f + 0.5)) .* exp(-1i * 2 * pi * 2 * (f + 0.5));

% c- Plot M(f)
m_freq = fftshift(fft(m_time) * dt);

figure(7);
plot(f, abs(m_freq), 'r', 'LineWidth', 1.5);
hold on;
plot(f, abs(m_analytical_freq), 'b', 'LineWidth', 1.5);
legend('M (Numerical)', 'M (Analytical)');
xlabel('Frequency (Hz)');
ylabel('|M(f)|');
xlim([-10 10]);
grid on;

% d- Calculate bandwidth of m(t)
psd_m = abs(m_freq).^2;

figure(8);
plot(f, psd_m, 'r');
title('Power Spectrum of M(t)');
xlabel('Frequency (Hz)');
ylabel('Power');
grid on;

threshold = 0.05 * max(psd_m);
indices = find(psd_m >= threshold);
BW_m = max(abs(f(indices)));

fprintf('M Estimated Bandwidth = %1.3f Hz\n', BW_m);

% e- Analytical bandwidth
psd_analytical_m = abs(m_analytical_freq).^2;
threshold = 0.05 * max(psd_analytical_m);
indices = find(psd_analytical_m >= threshold);
BW_analytical_m = max(abs(f(indices)));

fprintf('M Analytical Estimated Bandwidth = %1.3f Hz\n', BW_analytical_m);

%%%%%%%%%%%%%%%%%%%%%%%%% 8- Frequency Division Multiplexing (FDM) %%%%%%%%%%%%%
c1 = cos(2 * pi * 20 * t);
x_modulated = x_filtered1 .* c1;  % s1(t)


%%%%%%%%%%%%%%%%%%%%%%%%% 9- Second carrier and modulated m(t) %%%%%%%%%%%%%%%%%
c2 = cos(2 * pi * 24 * t);
m_modulated = m_time .* c2;  % s2(t)
m_modulated_freq = fftshift(fft(m_modulated) * dt);

%%%%%%%%%%%%%%%%%%%%%%%%% 10- Apply Bandpass Filter for LSB %%%%%%%%%%%%%%%%%%%%
BBF = zeros(size(f));
BBF(f >= 22 & f <= 24) = 1;
BBF(f <= -22 & f >= -24) = 1;

m_modulated_freq2 = m_modulated_freq .* BBF;
m_modulated2 = ifft(ifftshift(m_modulated_freq2)) / dt;

%%%%%%%%%%%%%%%%%%%%%%%%% 11- Combine and plot FDM signals %%%%%%%%%%%%%%%%%%%%%
s1 = x_modulated;
s2 = m_modulated2;
s_time = s1 + s2;
s_freq = fftshift(fft(s_time) * dt);

figure(12);
plot(t, s_time, 'r', 'LineWidth', 1.5);
title('Combined FDM Signal s(t)');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([-4 4]);
ylim([-5 5]);
grid on;

figure(13);
plot(f, abs(s_freq), 'r', 'LineWidth', 1.5);
title('Spectrum of FDM Signal S(f)');
xlabel('Frequency (Hz)');
ylabel('|S(f)|');
xlim([-30 30]);
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%% 12- Coherent Demodulation for x(t) %%%%%%%%%%%%%%%%%%%
x_demodulated= s_time .* cos(2 * pi * 20 * t);
x_demodulated_freq = fftshift(fft(x_demodulated)) * dt;

LPF_BW = 1;
H = abs(f) < LPF_BW;
x_rec = 2 * real(ifft(ifftshift(H .* x_demodulated_freq)) / dt);

figure(14);
plot(t, x_time, 'r', 'LineWidth', 1.5);
hold on;
plot(t, x_rec, 'b', 'LineWidth', 1.5);
legend('Original X(t)', 'Demodulated X(t)');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([-10 10]);
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%% 13- Coherent Demodulation for m(t) %%%%%%%%%%%%%%%%%%%
m_demodulated = s_time .* cos(2 * pi * 24 * t);
m_demodulated_freq = fftshift(fft(m_demodulated)) * dt;

H = abs(f) < LPF_BW;
m_rec = 4 * real(ifft(ifftshift(H .* m_demodulated_freq)) / dt);

figure(15);
plot(t, m_time, 'r', 'LineWidth', 1.5);
hold on;
plot(t, m_rec, 'b', 'LineWidth', 1.5);
legend('Original M(t)', 'Demodulated M(t)');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([-10 10]);
grid on;

