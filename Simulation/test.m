%% 1. 参数定义
clear; clc; close all;

% 系统参数
Fs = 25e6;          % 基带采样率 (Hz)
PRT = 232.76e-6;    % 脉冲重复周期 (s)
c = 3e8;            % 光速 (m/s)

% 窄脉冲参数
tau_narrow = 0.16e-6;
gap1_duration = 11.4e-6;

% 中脉冲参数
tau_medium = 8e-6;
B_medium = 20e6;
k_medium = -B_medium / tau_medium; % 负调频率
gap2_duration = 31.8e-6;

% 长脉冲参数
tau_long = 28e-6;
B_long = 20e6;
k_long = B_long / tau_long; % 正调频率

% 目标参数
R = 0.3e3;           % 目标距离 (m)
SNR_dB = 20;        % 信噪比 (dB)

%% 2. 时间和样本数计算
% PRT总样本数
N_total = round(PRT * Fs);
t_total = (0:N_total-1) / Fs; % PRT内的总时间轴

% 各段样本数
N_narrow = round(tau_narrow * Fs);
N_gap1 = round(gap1_duration * Fs);
N_medium = round(tau_medium * Fs);
N_gap2 = round(gap2_duration * Fs);
N_long = round(tau_long * Fs);

% 各脉冲起始样本点
start_idx_narrow = 1; % MATLAB索引从1开始
end_idx_narrow = start_idx_narrow + N_narrow - 1;

start_idx_medium = end_idx_narrow + N_gap1 + 1;
end_idx_medium = start_idx_medium + N_medium - 1;

start_idx_long = end_idx_medium + N_gap2 + 1;
end_idx_long = start_idx_long + N_long - 1;

%% 3. 生成发射信号 (一个PRT)
Tx_signal = zeros(1, N_total);

% (1) 生成窄脉冲 (简单脉冲，基带为常数)
Tx_signal(start_idx_narrow : end_idx_narrow) = 1.0 + 0i;

% (2) 生成中脉冲 (负线性调频)
t_medium_pulse = ((0:N_medium-1) / Fs) - tau_medium / 2;
medium_pulse = exp(1i * pi * k_medium * t_medium_pulse.^2);
Tx_signal(start_idx_medium : end_idx_medium) = medium_pulse;

% (3) 生成长脉冲 (正线性调频)
t_long_pulse = ((0:N_long-1) / Fs) - tau_long / 2;
long_pulse = exp(1i * pi * k_long * t_long_pulse.^2);
Tx_signal(start_idx_long : end_idx_long) = long_pulse;

%% 4. 生成回波信号
% 计算目标回波延迟
time_delay = 2 * R / c;
sample_delay = round(time_delay * Fs);

% 通过延时和补零生成无噪声的回波信号
Rx_signal_noiseless = circshift(Tx_signal, sample_delay);

% 计算信号功率并添加噪声
signal_power = mean(abs(Tx_signal(Tx_signal ~= 0)).^2);
noise_power = signal_power / (10^(SNR_dB / 10));

% 生成复高斯白噪声
noise = (randn(1, N_total) + 1i * randn(1, N_total)) * sqrt(noise_power / 2);
Rx_signal = Rx_signal_noiseless + noise;

%% 5. 可视化发射信号
figure('Name', '发射信号分析');
subplot(2, 1, 1);
plot(t_total * 1e6, real(Tx_signal));
title('发射信号 (Tx) - 一个PRT');
ylabel('实部'); grid on;
subplot(2, 1, 2);
plot(t_total * 1e6, imag(Tx_signal));
xlabel('时间 (us)');
ylabel('虚部'); grid on;

figure('Name', '回波信号分析');
plot(t_total * 1e6, real(Rx_signal));
title(sprintf('接收回波信号 (Rx) - 目标距离 %.1f km, SNR %.0f dB', R/1e3, SNR_dB));
xlabel('时间 (us)');
ylabel('幅度'); grid on;



%% 6. 脉冲压缩处理
% (1) 生成三个脉冲的匹配滤波器
MF_narrow = conj(fliplr(Tx_signal(start_idx_narrow : end_idx_narrow)));
MF_medium = conj(fliplr(Tx_signal(start_idx_medium : end_idx_medium)));
MF_long = conj(fliplr(Tx_signal(start_idx_long : end_idx_long)));

% (2) 截取各段用于处理的回波数据 (从发射结束时刻开始)
Rx_segment_narrow = Rx_signal(end_idx_narrow + 1 : end);
Rx_segment_medium = Rx_signal(end_idx_medium + 1 : end);
Rx_segment_long = Rx_signal(end_idx_long + 1 : end);

% (3) 执行匹配滤波 (卷积)
PC_narrow_output = conv(Rx_segment_narrow, MF_narrow, 'full');
PC_medium_output = conv(Rx_segment_medium, MF_medium, 'full');
PC_long_output = conv(Rx_segment_long, MF_long, 'full');


%% 7. CORRECTED: 按采样波门拼接处理结果
% (1) 定义各段采样波门点数
N_gate_narrow = 228;
N_gate_medium = 723;
N_gate_long = 2453;
N_total_gate = N_gate_narrow + N_gate_medium + N_gate_long;

% (2) --- 核心逻辑 ---
% 从每段处理结果中，截取它所负责的距离范围的数据
% 定义拼接的起始和结束点
idx1_start = 1;
idx1_end = N_gate_narrow;

idx2_start = idx1_end + 1;
idx2_end = idx1_end + N_gate_medium;

idx3_start = idx2_end + 1;
idx3_end = idx2_end + N_gate_long;

% 从各自的处理输出中截取正确的距离段
piece1 = PC_narrow_output(idx1_start:idx1_end);
piece2 = PC_medium_output(idx2_start:idx2_end);
piece3 = PC_long_output(idx3_start:idx3_end);

% 拼接成最终的剖面
final_profile = [piece1, piece2, piece3];

% (3) 创建拼接后的统一距离轴 (这部分不变)
final_range_axis = (0:N_total_gate-1) * (c / (2 * Fs));

% (4) 计算分界点的距离值，用于绘图 (这部分不变)
range_boundary1 = (N_gate_narrow) * (c / (2 * Fs));
range_boundary2 = (N_gate_narrow + N_gate_medium) * (c / (2 * Fs));

%% 8. NEW: 可视化拼接后的距离剖面
figure('Name', '拼接后的全距离剖面');
plot(final_range_axis / 1000, abs(final_profile));
title('窄/中/长脉冲处理拼接后的全距离剖面');
xlabel('距离 (km)');
ylabel('幅度');
grid on;

figure('Name', '拼接后的全距离剖面');
plot(final_range_axis / 1000, 20*log10(abs(final_profile)));
title('窄/中/长脉冲处理拼接后的全距离剖面');
xlabel('距离 (km)');
ylabel('幅度(dB)');
grid on;

% 绘制目标真实位置
line([R/1000 R/1000], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);

% 绘制分界线
line([range_boundary1/1000 range_boundary1/1000], ylim, 'Color', 'g', 'LineStyle', ':');
line([range_boundary2/1000 range_boundary2/1000], ylim, 'Color', 'm', 'LineStyle', ':');

% 添加文字说明
text(range_boundary1/2000, max(abs(final_profile))*0.9, '窄脉冲处理段', 'Color', 'k', 'HorizontalAlignment', 'center');
text(range_boundary1/1000 + (range_boundary2-range_boundary1)/2000, max(abs(final_profile))*0.9, '中脉冲处理段', 'Color', 'k', 'HorizontalAlignment', 'center');
text(range_boundary2/1000 + (max(final_range_axis)-range_boundary2)/2000, max(abs(final_profile))*0.9, '长脉冲处理段', 'Color', 'k', 'HorizontalAlignment', 'center');

legend('拼接后剖面', '目标真实位置', '窄/中 分界', '中/长 分界');

