% plot_beam_patterns.m
% 描述: 这是一个用于加载DBF系数并绘制13个雷达波束方向图的专用脚本。

%  修改记录
%  date       by      version   modify
%  25/10/10   XZR      v1.0     创建


clc; clear; close all;

%% 1. 参数配置
% =========================================================================
fprintf('--- 1. 配置绘图所需参数 ---\n');

% --- 1.1 文件路径 ---
% !!! 注意: 请将DBF系数文件的路径修改为您的实际位置 !!!
dbf_coef_path = 'C:\Users\a\Desktop\9.8 问题\Simulation\X8数据采集250522_DBFcoef.csv'; % <--- 修改我

% --- 1.2 雷达与阵列参数 ---
fc = 9500e6;                     % 雷达工作中心频率 (Hz)
c = 2.99792458e8;                % 光速 (m/s)
channel_num = 16;                % 接收通道数
element_spacing = 0.0138;        % 阵元间距 (m)

% --- 1.3 派生参数 ---
wavelength = c / fc;             % 波长

%% 2. 加载并解析DBF系数
% =========================================================================
fprintf('--- 2. 加载并解析DBF系数 ---\n');

if ~exist(dbf_coef_path, 'file')
    error('找不到DBF系数文件，请检查路径: %s', dbf_coef_path);
end
DBF_coeffs_data = readmatrix(dbf_coef_path);
% 将系数转换为复数形式, 尺寸应为 [13, 16]
DBF_coeffs_C = DBF_coeffs_data(:, 1:2:end) + 1j * DBF_coeffs_data(:, 2:2:end);
DBF_coeffs_C = fliplr(DBF_coeffs_C);
num_beams = size(DBF_coeffs_C, 1);
fprintf('成功加载 %d 个波束的系数。\n', num_beams);

%% 3. 计算阵列流形 (Array Manifold) / 导向矢量 (Steering Vectors)
% =========================================================================
fprintf('--- 3. 计算阵列导向矢量 ---\n');
% 需要计算在每个可能的来波方向上，16个通道接收到的信号相位关系

% 定义感兴趣的角度范围
angles_deg = -90:0.1:100; % 从-90度到+90度，步进0.1度
angles_rad = deg2rad(angles_deg);

% 阵元索引 (从0到15)
n_indices = (1:channel_num)';

% 计算阵列流形矩阵, 每一列对应一个角度的导向矢量
% 矩阵尺寸为 [16, 1801]
steering_vectors = exp(1j * 2 * pi * element_spacing * n_indices * sin(angles_rad) / wavelength);

%% 4. 计算并绘制13个波束的方向图
% =========================================================================
fprintf('--- 4. 计算并绘制方向图 ---\n');

figure('Name', '13个波束的雷达方向图');
hold on;

% 存储每个波束的峰值方向
peak_angles = zeros(1, num_beams);
%b = 1;
for b = 1:num_beams
    % 提取当前波束的DBF加权系数 (一个 1x16 的行向量)
    w = DBF_coeffs_C(b, :);

    % 计算方向图: 将DBF系数与所有角度的导向矢量做内积
    % (1x16) * (16x1801) -> (1x1801)
    beam_pattern = w * steering_vectors;

    % 将方向图转换为dB，并进行归一化 (使主瓣峰值为0dB)
    beam_pattern_abs = abs(beam_pattern);
   % beam_pattern_db = 20 * log10(beam_pattern_abs / max(beam_pattern_abs));
    beam_pattern_db = 20 * log10(beam_pattern_abs);
    % 绘制当前波束的方向图
    plot(angles_deg, beam_pattern_db, 'LineWidth', 1.5);

    % 找到当前波束的峰值角度
    [~, max_idx] = max(beam_pattern_db);
    peak_angles(b) = angles_deg(max_idx);
end

hold off;
grid on;
title('雷达13波束方向图');
xlabel('角度 (度)');
ylabel('归一化增益 (dB)');
legend(arrayfun(@(b, ang) sprintf('波束 #%d (指向 %.1f°)', b, ang), 1:num_beams, peak_angles, 'UniformOutput', false));
ylim([-50, 5]); % 设置一个好的观察范围，以便看清旁瓣
ax = gca;
ax.FontSize = 12;

fprintf('方向图绘制完成。\n');