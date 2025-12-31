% run_music_algorithm.m
% 描述: 演示MUSIC超分辨率测向算法如何分辨两个角度相近的信号。
clc; clear; close all;

%% 1. 仿真参数定义
fprintf('--- 1. 配置仿真所需参数 ---\n');
fc = 9450e6;
c = 2.99792458e8;
channel_num = 16;
element_spacing = 0.0138;
wavelength = c / fc;

% --- 关键：定义两个角度相近的信号，模拟多径 ---
source_angles = [2.0, -1.5]; % [直达波角度, 反射波角度] (度)
source_amplitudes = [1, 0.7]; % 反射波幅度稍弱
num_sources = length(source_angles);

% --- 快拍数与信噪比 ---
num_snapshots = 256; % 采集的快拍数 (类似于CPI内的脉冲数)
SNR_dB = 15;         % 信噪比

%% 2. 生成16通道的接收信号
fprintf('--- 2. 生成16通道接收信号 ---\n');
n_indices = (0:channel_num-1)';

% --- 生成两个信源的理论导向矢量 ---
A = exp(1j * 2 * pi * element_spacing * n_indices * sind(source_angles) / wavelength);

% --- 生成随机的基带信号 (多快拍) ---
S = (randn(num_sources, num_snapshots) + 1j * randn(num_sources, num_snapshots)) / sqrt(2);
S(1,:) = S(1,:) * source_amplitudes(1); % 施加幅度
S(2,:) = S(2,:) * source_amplitudes(2);

% --- 生成噪声 ---
noise_power = 1 / (10^(SNR_dB/10));
N = sqrt(noise_power/2) * (randn(channel_num, num_snapshots) + 1j * randn(channel_num, num_snapshots));

% --- 最终接收信号 = 信号 + 噪声 ---
X = A * S + N;

%% 3. 执行MUSIC算法
fprintf('--- 3. 执行MUSIC算法 ---\n');

% --- 3.1 计算协方差矩阵 ---
Rxx = (X * X') / num_snapshots;

% --- 3.2 特征值分解 ---
[E, D] = eig(Rxx);
eigenvalues = real(diag(D)); % 特征值为实数
[eigenvalues, sort_idx] = sort(eigenvalues, 'descend'); % 降序排列
E = E(:, sort_idx); % 对应地调整特征向量的顺序

% --- 3.3 划分子空间 ---
% 信号子空间 (对应最大的 num_sources 个特征值)
Es = E(:, 1:num_sources);
% 噪声子空间 (对应剩下的 N-num_sources 个特征值)
En = E(:, num_sources+1:end);

% --- 3.4 谱峰搜索 ---
scan_angles = -20:0.1:20;
music_spectrum = zeros(size(scan_angles));
for i = 1:length(scan_angles)
    % 获取当前扫描角度的导向矢量
    a = exp(1j * 2 * pi * element_spacing * n_indices * sind(scan_angles(i)) / wavelength);
    
    % 核心公式: 计算空间谱值
    % a' * En * En' * a 计算了导向矢量a在噪声子空间上的投影能量
    % 当a是真实信号方向时，它与噪声子空间正交，分母接近0，谱值出现峰值
    music_spectrum(i) = 1 / (a' * En * En' * a);
end

%% 4. 结果可视化
fprintf('--- 4. 生成结果图表 ---\n');
figure('Name', 'MUSIC算法 vs. 传统波束形成');

% --- 绘制MUSIC空间谱 ---
plot(scan_angles, 10*log10(abs(music_spectrum) / max(abs(music_spectrum))), 'r-', 'LineWidth', 2, 'DisplayName', 'MUSIC 空间谱');
hold on;

% --- 作为对比，绘制传统DBF的响应 ---
% 假设一个传统波束指向0度
w_conventional = hamming(channel_num);
A_scan = exp(1j * 2 * pi * element_spacing * n_indices * sind(scan_angles) / wavelength);
pattern_conventional = w_conventional.' * A_scan;
plot(scan_angles, 20*log10(abs(pattern_conventional) / max(abs(pattern_conventional))), 'b--', 'LineWidth', 1.5, 'DisplayName', '传统波束 (指向0°)');

grid on;
title('MUSIC超分辨率算法 vs. 传统波束形成');
xlabel('角度 (度)');
ylabel('归一化空间谱 (dB)');
legend;
ylim([-50, 5]);

% 标记真实信号位置
xline(source_angles(1), '--k', sprintf('真实目标 1: %.1f°', source_angles(1)));
xline(source_angles(2), '--k', sprintf('真实目标 2: %.1f°', source_angles(2)));