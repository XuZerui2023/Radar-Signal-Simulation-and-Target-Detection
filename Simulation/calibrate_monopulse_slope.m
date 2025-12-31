% calibrate_monopulse_slope.m
% 描述: 这是一个用于仿真标定用于和差比幅单脉冲测角中单脉冲斜率(K值)的专用脚本。

%  修改记录
%  date       by      version   modify
%  25/10/10   XZR      v1.0     创建

clc; clear; close all;

%% 1. 参数配置 (与您的主程序和方向图一致)
fprintf('--- 1. 配置标定所需参数 ---\n');
dbf_coef_path = 'C:\Users\a\Desktop\9.8 问题\Simulation\X8数据采集250522_DBFcoef.csv'; % <--- 修改我
fc = 9450e6;
c = 2.99792458e8;
channel_num = 16;
element_spacing = 0.0138;
wavelength = c / fc;

% --- 关键：从您的方向图图例中读取的精确波束指向 ---
beam_angles_deg = [-16, -9.6, -3.2, 3.2, 9.6, 16, 22.6, 29.2, 36.1, 43.3, 51, 59.6, 70.3];

%% 2. 加载并反转DBF系数
fprintf('--- 2. 加载并解析DBF系数 ---\n');
DBF_coeffs_data = readmatrix(dbf_coef_path);
DBF_coeffs_C = DBF_coeffs_data(:, 1:2:end) + 1j * DBF_coeffs_data(:, 2:2:end);
DBF_coeffs_C = fliplr(DBF_coeffs_C); % 保持与主程序一致

%% 3. 生成鉴角曲线并计算K值
fprintf('--- 3. 生成鉴角曲线并计算斜率 K ---\n');

% --- 选择一对相邻波束进行标定 (例如中心波束#7和#8) ---
beam_A_idx = 12;
beam_B_idx = 13;
w_A = DBF_coeffs_C(beam_A_idx, :);
w_B = DBF_coeffs_C(beam_B_idx, :);

% 确定这对波束的中心交点角度
crossover_angle = (beam_angles_deg(beam_A_idx) + beam_angles_deg(beam_B_idx)) / 2;

% 在交点附近设置一个精细的角度扫描范围
calib_angles = linspace(crossover_angle - 3, crossover_angle + 3, 501);

% --- 计算每个角度下的单脉冲比 ---
n_indices = (0:channel_num-1)';
steering_vectors = exp(1j * 2 * pi * element_spacing * n_indices * sind(calib_angles) / wavelength);

% 计算和、差波束在每个角度上的响应
response_A = w_A * steering_vectors;
response_B = w_B * steering_vectors;
response_sum = response_A + response_B;
response_diff = response_A - response_B;
monopulse_ratio = response_diff ./ response_sum;

% --- 计算斜率 K ---
% K值是鉴角曲线在零点附近的斜率
% 角度偏移 = K * real(单脉冲比) => K = 角度偏移 / real(单脉冲比)
% 我们在零点附近取一小段线性区域来估算斜率
[~, center_idx] = min(abs(calib_angles - crossover_angle));
delta_idx = 10; % 在中心点附近取10个点
angles_for_slope = calib_angles(center_idx-delta_idx : center_idx+delta_idx) - crossover_angle;
ratios_for_slope = real(monopulse_ratio(center_idx-delta_idx : center_idx+delta_idx));
% 使用线性拟合得到斜率
p = polyfit(ratios_for_slope, angles_for_slope, 1);
k_slope = p(1);

fprintf('\n标定完成 (波束对 #%d, #%d):\n', beam_A_idx, beam_B_idx);
fprintf('----------------------------------------\n');
fprintf('>>> 单脉冲斜率 (K值) ≈ %.4f (度/单位比值) <<<\n', k_slope);
fprintf('----------------------------------------\n');

%% 4. 可视化鉴角曲线
figure('Name', '单脉冲鉴角曲线');
plot(calib_angles, real(monopulse_ratio), 'b-', 'LineWidth', 2);
hold on;
plot(calib_angles, imag(monopulse_ratio), 'r--', 'LineWidth', 1.5);
grid on;
ax = gca;
ax.XAxisLocation = 'origin'; % 将X轴移到y=0处
ax.YAxisLocation = 'origin'; % 将Y轴移到x=0处
xlabel('目标角度 (度)');
ylabel('单脉冲比值');
title(sprintf('波束对 (#%d, #%d) 的鉴角曲线', beam_A_idx, beam_B_idx));
legend('比值实部 (用于测角)', '比值虚部 (通常为零)');
text(crossover_angle + 0.5, 0.1, sprintf('计算出的 K 值 ≈ %.2f', k_slope), 'FontSize', 12);