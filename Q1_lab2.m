clear all;
close all;

[serial_1, time_1, signal_1] = textread('mono2048.txt', '%f %f %f');
[serial_2, time_2, signal_2] = textread('monowavefil2.txt', '%f %f %f'); % change it with the DUT
%[serial_2, time_2, signal_2] = textread('monowaveguide.txt', '%f %f %f');
%[serial_2, time_2, signal_2] = textread('monotube.txt', '%f %f %f');


N1 = size(time_1, 1);
N2 = size(time_2, 1);



time_i_1 = zeros(1, N1);

for i = 1:N1 - 1
    time_i_1(i) = time_1(i + 1) - time_1(i);
end

time_d_avg_1 = sum(time_i_1.')./ N1;

time_i_2 = zeros(1, N2);

for i = 1:N2 - 1
    time_i_2(i) = time_2(i + 1) - time_2(i);
end

time_d_avg_2 = sum(time_i_1.')./ N2;


hann_ = hanning(N1);
signal_hann_1 =  detrend(signal_1);
%signal_hann_1 = hann_ .* detrend(signal_1);

signal_freq_1 = fft((signal_hann_1), N1);

Pyy_1 = signal_freq_1 .* conj(signal_freq_1) ./ N1;

f_max_1 = 1 ./ (2.* time_d_avg_1);

fs_1 = 2 .* f_max_1;

signal_hann_2 = detrend(signal_2);

%signal_hann_2 = hann_ .* detrend(signal_2);

signal_freq_2 = fft((signal_hann_2), N2);

Pyy_2 = signal_freq_2 .* conj(signal_freq_2) ./ N2;

f_max_2 = 1 ./ (2.* time_d_avg_2);

fs_2 = 2 .* f_max_2;

frequencies_1 = fs_1 .* (0:N1/2) ./ N1;

plot(frequencies_1(1:end), 0.5 * db(Pyy_1(1:N1/2+1)), 'LineWidth', 2);
hold on;

frequencies_2 = fs_2 .* (0:N2/2) ./ N2;

plot(frequencies_2(1:end), 0.5 * db(Pyy_2(1:N2/2+1)), 'LineWidth', 2);


xlabel('Frequencies(GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('P_{xx}, P_{yy} (dB scale)', 'FontSize', 12, 'FontWeight', 'bold');
title('Power Spectrum of mono pulse and response from in DUT', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Pulse', 'DUT'}, 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0.1 12])
%print('Power_Spectrum_db_wfil', '-depsc');
%print('Power_Spectrum_db_wg', '-depsc');
%print('Power_Spectrum_db_wtube', '-depsc');

figure(2);
plot(frequencies_1(1:end), (Pyy_1(1:N1/2+1)), 'LineWidth', 2);
hold on;

plot(frequencies_2(1:end), (Pyy_2(1:N2/2+1)), 'LineWidth', 2);


xlabel('Frequencies(GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('P_{xx}, P_{yy} (Linear Scale)', 'FontSize', 12, 'FontWeight', 'bold');
title('Power Spectrum of mono pulse and response from in DUT', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Pulse', 'DUT'}, 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0.1 12])
%print('Power_Spectrum_linear_wfil', '-depsc');
%print('Power_Spectrum_linear_wg', '-depsc');
%print('Power_Spectrum_linear_tube', '-depsc');


%% Transfer Function: 

h = (signal_freq_2)./(signal_freq_1);

%h1 = (signal_freq_2 .* hann_)./(signal_freq_1 .* hann_);


figure(3);
plot(frequencies_2(1:end), db(abs(h(1:N1/2+1))), 'LineWidth', 2, 'color', [0.6350, 0.0780, 0.1840]);

xlabel('Frequencies(GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('V_{hh} (dB scale)', 'FontSize', 12, 'FontWeight', 'bold');
title('Absolute value Spectrum of the transfer function', 'FontSize', 12, 'FontWeight', 'bold');

grid on;
xlim([0.1 13])
%print('Transfer_dB_wfil_mag', '-depsc');
%print('Transfer_dB_wg_mag', '-depsc');
%print('Transfer_dB_tube_mag', '-depsc');

figure(4);
plot(frequencies_2(1:end), abs(h(1:N1/2+1)), 'LineWidth', 2, 'color', [0.6350, 0.0780, 0.1840]);

xlabel('Frequencies(GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('V_{hh} (Linear Scale (mVolt))', 'FontSize', 12, 'FontWeight', 'bold');
title('Absolute value Spectrum of the transfer function', 'FontSize', 12, 'FontWeight', 'bold');

grid on;
xlim([0.1 13])

%print('Transfer_linear_wfil_mag', '-depsc');
%print('Transfer_linear_wg_mag', '-depsc');
%print('Transfer_linear_tube_mag', '-depsc');
% Windowing to remove noise

win_ = zeros(size(h, 1), size(h, 2));
win_a = zeros(size(h, 1), size(h, 2));
freq = fs_1 .* (0:N1-1) ./ N1;

% Abrupt window

% for i = 1:size(freq, 2)
%     if freq(i) > 0.1 && freq(i) < 10
%         win_a(i) = 1;
%         win_a(N1 - i) = 1;
%     end
% end

% Manual Hann window

count = 0;
for i = 1:size(freq, 2)
    if freq(i) > 0.1 && freq(i) < 13 %change it to 13 for wavefilter 
        count = count + 1;
    end
end

for i = 1:size(freq, 2)
    if freq(i) > 0.1 && freq(i) < 13 %change it to 13 for wavefilter
        win_(i) = sin((i - 0.1) .* pi ./ count).^2;
        win_(N1 - i) = win_(i);
    end
end




% Time domain impulse response

h_ = h .* win_;
%h_new = h1 .* win_;
h_a = h .* win_a;

h_t = ifft(h_);
h_t_1 = ifft(h_a);




figure(6);
plot(freq(1:end), (win_(1:end)), 'LineWidth', 2); % plot if you want to see the hann window

figure(5);

plot(time_1(1:end), real(h_t(1:end)), 'LineWidth', 2, 'color',  [0.6350, 0.0780, 0.1840]);
% hold on;
% plot(time_1(1:end), real(h_t_1(1:end)), '--', 'LineWidth', 2,  'color', [0.25, 0.25, 0.25]);

xlabel('time(nS)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('h(t) [mVolt]', 'FontSize', 12, 'FontWeight', 'bold');
title('Transfer function in time domain', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

%print('Transfer_time_wfil', '-depsc');
%print('Transfer_time_wg', '-depsc');
%print('Transfer_time_tube', '-depsc');

% figure(7);
% plot(time_1(1:end), real(h_t_1(1:end)), 'LineWidth', 2,  'color', [0.6350, 0.0780, 0.1840]);
% 
% xlabel('time(nS)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('h(t) [mVolt]', 'FontSize', 12, 'FontWeight', 'bold');
% title('Transfer function in time domain after applying Hanning window', 'FontSize', 12, 'FontWeight', 'bold');
% grid on;

%print('Transfer_time_wfil', '-depsc');
%print('Transfer_time_wg', '-depsc');
%print('Transfer_time_tube', '-depsc');

% figure(5);
% 
% s = detrend(signal_1);
% plot(time_1, signal_1);
% hold on;
% plot(time_1, s);

%-------------------------------

% Bonus Thing

N_noise = length(frequencies_1(frequencies_1>20));

Noise = sum((signal_freq_1(frequencies_1 > 20))) ./ N_noise;

SNR = max(abs(signal_freq_1))./abs(Noise);

sigma = max(signal_1)./SNR;

P_n = N1 .* sigma.^2;

beta = P_n ./ (mean(Pyy_1));

H = (signal_freq_2 .* conj(signal_freq_1)) ./ (signal_freq_1 .* conj(signal_freq_1) + beta);

H_new = H .* win_;

H_t = ifft(H_new);

figure(9);
plot(frequencies_2(1:end), abs(H(1:N1/2+1)), 'LineWidth', 2, 'color', [0.6350, 0.0780, 0.1840]);

xlabel('Frequencies(GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('V_{hh} (Linear Scale (mVolt))', 'FontSize', 12, 'FontWeight', 'bold');
title('Absolute value Spectrum of the transfer function', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0.1 13]);
print('bonus_wavefil_l', '-depsc');


figure(10);

plot(frequencies_2(1:end), 0.5 .* db(abs(H(1:N1/2+1))), 'LineWidth', 2, 'color', [0.25, 0.25, 0.25]);

xlim([0.1 6]);

xlabel('Frequencies(GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('V_{hh} (dB scale)', 'FontSize', 12, 'FontWeight', 'bold');
title('Absolute value Spectrum of the transfer function', 'FontSize', 12, 'FontWeight', 'bold');

grid on;
%print('bonus_wg_db', '-depsc');

figure(11);

plot(time_1(1:end), real(H_t(1:end)), 'LineWidth', 2, 'color', [0.6350, 0.0780, 0.1840]);

xlabel('time(nS)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('h(t) [mVolt]', 'FontSize', 12, 'FontWeight', 'bold');
title('Transfer function in time domain', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

%print('bonus_wg_t', '-depsc');
