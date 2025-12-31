%% 1. 仿真参数定义
clc; clear; close all;

% --- 雷达与阵列参数 ---
fc = 9450e6;                     % 雷达工作中心频率 (Hz)
c = 2.99792458e8;                % 光速 (m/s)
wavelength = c / fc;             % 波长
N_ch = 16;                       % 接收通道数
d = 0.0138;                      % 阵元间距 (m) (与您代码一致)

% --- 目标参数 ---
target_range = 800;              % 目标距离 (m)
target_velocity = 30;            % 目标速度 (m/s)
target_angle_deg = 2.0;          % <--- 目标真实角度 (偏离中心2度)

%% 2. 生成和波束(Σ)与差波束(Δ)的DBF系数
% 这是一个教科书式的标准方法，用于为均匀线性阵列生成系数

% --- 和波束系数 (Σ) ---
% 使用汉明窗对所有通道进行加权，以降低旁瓣
win_sum = hamming(N_ch)';
coeffs_sum = win_sum; 

% --- 差波束系数 (Δ) ---
% 在和波束的幅度加权基础上，叠加一个反对称的线性相位
win_diff = hamming(N_ch)';
phase_taper_diff = -pi/2 : pi/(N_ch-1) : pi/2; % 线性相位斜坡
coeffs_diff = win_diff .* exp(1j * phase_taper_diff);

%% 3. 可视化和/差波束方向图
figure('Name', '和波束与差波束方向图');
angles = -90:0.1:90;
steering_vectors = exp(1j * 2 * pi * d * (0:N_ch-1)' * sind(angles) / wavelength);

pattern_sum = abs(coeffs_sum * steering_vectors);
pattern_diff = abs(coeffs_diff * steering_vectors);

plot(angles, 20*log10(pattern_sum / max(pattern_sum)), 'b-', 'LineWidth', 2);
hold on;
plot(angles, 20*log10(pattern_diff / max(pattern_diff)), 'r-', 'LineWidth', 2);
grid on;
ylim([-40, 5]);
title('和波束(Σ)与差波束(Δ)方向图');
xlabel('角度 (度)');
ylabel('归一化增益 (dB)');
legend('和波束 (Σ)', '差波束 (Δ)');
xline(0, '--k', '轴上 (Boresight)');

%% 4. 模拟目标回波并进行波束形成
% 只需模拟单个脉冲在目标峰值处的数据即可
% (为简化，忽略距离和多普勒，只关注角度)

% --- 接收到的16通道信号 (在目标峰值处) ---
target_angle_rad = deg2rad(target_angle_deg);
% 目标的阵列流形矢量
target_steering_vec = exp(1j * 2 * pi * d * (0:N_ch-1)' * sin(target_angle_rad) / wavelength);
% 假设信号幅度为1，接收到的信号就是其流形矢量
signal_16ch_peak = target_steering_vec.';

% --- 应用和/差波束形成 ---
output_sum = signal_16ch_peak * coeffs_sum';
output_diff = signal_16ch_peak * coeffs_diff';

%% 5. 执行单脉冲测角
% --- 5.1 计算单脉冲斜率 k (系统标定) ---
% 单脉冲斜率是天线固有的，需要预先标定。
% 我们这里通过仿真一个偏离轴上很小的角度(如0.1度)来估算它。
calib_angle_deg = 0.1;
calib_angle_rad = deg2rad(calib_angle_deg);
calib_steering_vec = exp(1j * 2 * pi * d * (0:N_ch-1)' * sin(calib_angle_rad) / wavelength).';
calib_ratio = (calib_steering_vec * coeffs_diff') / (calib_steering_vec * coeffs_sum');
% 斜率 k ≈ 角度 / (差/和比值的实部)
k_monopulse = calib_angle_deg / real(calib_ratio);
fprintf('系统标定：单脉冲斜率 k ≈ %.4f (度/单位比值)\n\n', k_monopulse);

% --- 5.2 计算目标的单脉冲比 ---
monopulse_ratio = output_diff / output_sum;

% --- 5.3 估算目标角度 ---
% 估算出的角度 ≈ k * (差/和比值的实部)
estimated_angle_deg = k_monopulse * real(monopulse_ratio);

%% 6. 结果分析
fprintf('--- 单脉冲测角结果 ---\n');
fprintf('目标真实角度: %.4f 度\n', target_angle_deg);
fprintf('和波束(Σ)输出幅度: %.2f\n', abs(output_sum));
fprintf('差波束(Δ)输出幅度: %.2f\n', abs(output_diff));
fprintf('单脉冲比 (Δ/Σ): %.4f + %.4fi\n', real(monopulse_ratio), imag(monopulse_ratio));
fprintf('----------------------------------------\n');
fprintf('>>> 估算出的目标角度: %.4f 度 <<<\n', estimated_angle_deg);
fprintf('----------------------------------------\n');