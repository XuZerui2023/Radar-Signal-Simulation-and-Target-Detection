% calibrate_all_monopulse_slopes.m
% 描述: 自动为所有12个相邻波束对进行仿真标定，生成单脉冲斜率(K值)查找表。

%  修改记录
%  date       by      version   modify
%  25/10/12   XZR      v1.0     创建

clc; clear; close all;

%% 1. 参数配置
fprintf('--- 1. 配置标定所需参数 ---\n');
dbf_coef_path = 'C:\Users\a\Desktop\9.8 问题\Simulation\X8数据采集250522_DBFcoef.csv'; % DBF系数文件路径
fc = 9450e6;                     % 雷达工作中心频率 (Hz)
c = 2.99792458e8;                % 光速 (m/s)
channel_num = 16;                % 接收通道数
element_spacing = 0.0138;        % 阵元间距 (m)
wavelength = c / fc;             % 计算波长

% --- 关键：从 plot_beam_patterns.m 程序绘制的方向图中得到13个波束的中心指向角度 ---
beam_angles_deg = [-16, -9.6, -3.2, 3.2, 9.6, 16, 22.6, 29.2, 36.1, 43.3, 51, 59.6, 70.3];

%% 2. 加载DBF系数
fprintf('--- 2. 加载并解析DBF系数 ---\n');
DBF_coeffs_data = readmatrix(dbf_coef_path);
DBF_coeffs_C = DBF_coeffs_data(:, 1:2:end) + 1j * DBF_coeffs_data(:, 2:2:end);
DBF_coeffs_C = fliplr(DBF_coeffs_C);  % 修正通道顺序，与主程序中的回波仿真模型保持一致

%% 3. 循环标定所有波束对
fprintf('--- 3. 循环标定所有12个波束对 ---\n');
k_slopes_LUT = zeros(1, 12); % 初始化K值查找表

figure('Name', '所有波束对的鉴角曲线');
hold on;

for pair_idx = 1:12
    % --- 准备当前波束对的数据 ---
    beam_A_idx = pair_idx;      % 左波束的索引
    beam_B_idx = pair_idx + 1;  % 右波束的索引
    w_A = DBF_coeffs_C(beam_A_idx, :); % 提取左波束的16个DBF通道系数
    w_B = DBF_coeffs_C(beam_B_idx, :); % 提取右波束的16个DBF通道系数
    
    % --- 定义仿真扫描的角度范围 ---
    % 计算两个波束指向的中心点，即鉴角曲线的理论过零点
    crossover_angle = (beam_angles_deg(beam_A_idx) + beam_angles_deg(beam_B_idx)) / 2;
    % 定义扫描宽度，以覆盖两个波束的重叠区域
    scan_width = abs(beam_angles_deg(beam_A_idx) - beam_angles_deg(beam_B_idx));
    calib_angles = linspace(crossover_angle - scan_width, crossover_angle + scan_width, 501);

    % --- 计算理论回波信号 ---
    n_indices = (0:channel_num-1)';
    % 核心: 为每个扫描角度，计算一个16x1的"阵列导向矢量"(Steering Vector)
    % 它描述了来自该角度的平面波到达16个阵元时的理论相位关系
    steering_vectors = exp(1j * 2 * pi * element_spacing * n_indices * sind(calib_angles) / wavelength);
    
    % --- 生成鉴角曲线 --- 
    response_A = w_A * steering_vectors; % 计算波束A对来自所有扫描角度的信号的响应 (复数)
    response_B = w_B * steering_vectors; % 计算波束B对来自所有扫描角度的信号的响应 (复数)
    monopulse_ratio = (response_A - response_B) ./ (response_A + response_B); % 核心: 对每个角度，计算和差比值，得到完整的鉴角曲线

    % 计算斜率 K
    % K值是鉴角曲线在中心点附近的斜率。我们通过线性拟合来精确计算它。
    % 找到离理论中心点最近的角度索引
    [~, center_idx] = min(abs(calib_angles - crossover_angle));
    delta_idx = 5; % 在中心点附近左右各取5个点，形成一个小窗口用于拟合 5+5+1
    
    % 提取这个小窗口对应的角度偏移量和比值实部
    angles_for_slope = calib_angles(center_idx-delta_idx : center_idx+delta_idx) - crossover_angle;
    ratios_for_slope = real(monopulse_ratio(center_idx-delta_idx : center_idx+delta_idx));
    % 核心: 对"比值(x)"和"角度(y)"进行一阶多项式拟合(即线性拟合 y=ax+b)
    % p(1)就是拟合出的直线斜率a，即我们需要的K值
    p = polyfit(ratios_for_slope, angles_for_slope, 1);
    k_slopes_LUT(pair_idx) = p(1);
    
    % 绘制当前波束对的鉴角曲线
    plot(calib_angles, real(monopulse_ratio), 'LineWidth', 1.5, 'DisplayName', sprintf('Pair %d-%d', beam_A_idx, beam_B_idx));
end

hold off; grid on;
title('12个波束对的鉴角曲线 (实部)');
xlabel('角度 (度)'); ylabel('单脉冲比值');
legend('show', 'Location', 'eastoutside');

% --- 打印最终的查找表 ---
fprintf('\n--- K值查找表 (k_slopes_LUT) ---\n');
fprintf('波束对\tK值\n');
for i = 1:12
    fprintf('%d-%d\t\t%.4f\n', i, i+1, k_slopes_LUT(i));
end
disp('以下为12组波束对的单脉冲测角（斜率）K值:');
disp(['k_slopes_LUT = [', num2str(k_slopes_LUT), '];']);