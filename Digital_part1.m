% ------------Parameters---------------%
num_bits = 64;
bit_rate = 100;               % bits per second
Fs = 1000;                    % Sampling frequency
Tb = 1/bit_rate;              % Bit period
t = 0:1/Fs:num_bits*Tb - 1/Fs;

% ----------Random bit stream----------%
bits = randi([0 1], 1, num_bits);

% ---------Polar NRZ Encoding----------%
polar_nrz = zeros(1, length(t));
for i_1 = 1:num_bits
    idx_p = (i_1-1)*Fs*Tb + 1:i_1*Fs*Tb;
    if bits(i_1) == 1
        polar_nrz(idx_p) = 1;
    else
        polar_nrz(idx_p) = -1;
    end
end
% -------------AMI Encoding--------------%
AMI =zeros(1, length(t));
no_of_ones=0;
for i_2 = 1:num_bits
    idx_A = (i_2-1)*Fs*Tb + 1:i_2*Fs*Tb;
    if bits(i_2) == 0
        AMI(idx_A) = 0;
    else
     no_of_ones=no_of_ones+1;
     if(rem(no_of_ones,2)==0)
         AMI(idx_A) = -1;
     else
         AMI(idx_A) = 1;

     end

    end
end

% -----------Plot  AMI Time Domain-----------%
figure (1);
plot(t, AMI, 'LineWidth', 2,'Color','R');
title('AMI Line Code - Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');
ylim([-1.5 1.5]);
grid on;
% -----------Plot Polar-NRZ Time Domain-----------%
figure(2)
plot(t,polar_nrz, 'LineWidth', 2,'color','b');
title('Polar-NRZ Line Code - Time Domain');
xlabel('Time (s)');

ylabel('Amplitude');
ylim([-1.5 1.5]);
grid on;

% -----------Plot AMI Frequency Domain------------%
N = length(AMI);
f = (-N/2:N/2-1)*(Fs/N);
AMI_fft = fftshift(abs(fft(AMI)/length(AMI)));
figure(3);
plot(f, AMI_fft, 'LineWidth', 1.5);
title('AMI Line Code - Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% ----------Plot Polar-NRZ Frequency Domain--------%
N = length(polar_nrz);
f_2 = (-N/2:N/2-1)*(Fs/N);
polar_nrz_fft = fftshift(abs(fft(polar_nrz)/length(polar_nrz)));
figure(4);
plot(f_2, polar_nrz_fft, 'LineWidth', 1.5);
title('Polar-NRZ Code - Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
