% ---------------------Parameters-------------------%
bit_rate = 1000;                 % bits per second
fc = 5000;                       % Carrier frequency (Hz)
num_bits = 64;                   % Number of bits
samples_per_bit = 100;           % Oversampling
fs = bit_rate * samples_per_bit; % Sampling frequency
T = 1 / bit_rate;
t = 0:1/fs:num_bits*T - 1/fs;
% -------------------random binary data--------------------%
data_bits = randi([0 1], 1, num_bits)
data_upsampled = repelem(data_bits, samples_per_bit);
%  -------------------ASK Modulation-------------------%
carrier = cos(2*pi*fc*t);
ask_signal = data_upsampled .* carrier;
% -------------Time-Domain Plot (Transmitter Output)-------------%
figure;
plot(t(1:1000), ask_signal(1:1000), 'LineWidth', 1.5);
title('ASK Signal at Transmitter (Time Domain)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
% -----------Frequency-Domain Plot (Transmitter Output)------------%
n = length(ask_signal);
df=fs/n;
if (rem(n,2)==0)
    f = (-fs/2):df:(fs/2-df);
else
    f = -(fs/2-df/2):df:(fs/2-df/2);
end
ASK_fft = fftshift(abs(fft(ask_signal))/n);
figure;
plot(f, ASK_fft, 'LineWidth', 1.5);
title('Spectrum of ASK Signal at Transmitter');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 2*fc]);
grid on;
% ---------Coherent Receiver: Demodulation with Phase Offsets-------------%
phases = [30, 60, 90];
for i = 1:length(phases)
    phase_deg = phases(i);
    phase_rad = deg2rad(phase_deg);
    % receiver_carrier with phase offset
    receiver_carrier = cos(2*pi*fc*t + phase_rad);
    % Coherent detection
    received = ask_signal .* receiver_carrier;
    % Integrate and decision
    demod_bits = zeros(1, num_bits);
    for i = 1:num_bits
        idx_start = (i-1)*samples_per_bit + 1;
        idx_end = i*samples_per_bit;
        segment = received(idx_start:idx_end);
        avg_val = mean(segment);
        demod_bits(i) = avg_val > 0.15;  % simple thresholding
    end
   % Plot demodulated signal (partial for visualization)
    figure;
    plot(t(1:1000), received(1:1000), 'LineWidth', 1.5);
    title(['Demodulated Signal (Phase = ', num2str(phase_deg), '°)']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
    % -----------Frequency-Domain Plot (Receiver Output)------------%
    received_fft = fftshift(abs(fft(received)) / n);
    figure;
    plot(f, received_fft, 'LineWidth', 1.5);
    title(['Spectrum of Received Signal - Phase = ', num2str(phase_deg), '°']);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    grid on;
   end

