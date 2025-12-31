% main_simulate_and_process_v3.m
% 版本: v3.0
% 描述: 这是一个集成了高保真回波仿真、DBF、脉冲压缩（含偏移校正）、
%       MTD和CFAR检测的全流程雷达信号处理仿真与分析脚本。
%       脚本基于吴晓燕女士的论文思想，对关键偏移量进行了校正，
%       并在每个处理阶段都提供了详细的可视化输出。
% 日期: 2025年09月12日

clc; clear; close all;

%% 1. 用户配置区
% =========================================================================
fprintf('--- 1. 开始进行用户参数配置 ---\n');

% --- 1.1 定义仿真的目标 ---
targets = struct(...
    'Range', {6000}, ...         % 目标距离 (m)
    'Velocity', {30}, ...          % 目标速度 (m/s, 正为远离)
    'RCS', {50}, ...                % 雷达散射截面 (m^2)
    'ElevationAngle', {15} ...        % 目标俯仰角 (度)
);

% --- 1.2 文件路径配置 ---
% !!! 注意: 请将以下路径修改为您项目中文件的实际位置 !!!
base_path = 'C:\Users\a\Desktop\9.8 问题\Simulation'; % <--- 修改我
dbf_coef_path = fullfile(base_path, 'X8数据采集250522_DBFcoef.csv');
% fir_filter_path = fullfile(base_path, 'Simulation', 'FIR.mat'); % FIR系数保存在此

% --- 1.3 CFAR 参数配置 ---
cfar_params.refCells_V = 4;      % 速度维参考单元
cfar_params.guardCells_V = 6;   % 速度维保护单元
cfar_params.refCells_R = 4;      % 距离维参考单元
cfar_params.guardCells_R = 6;    % 距离维保护单元
cfar_params.T_CFAR = 6;          % 检测门限因子
cfar_params.method = 'GOCA';       % 'CA' (平均), 'GOCA' (选大), 'SOCA' (选小)

%% 2. 雷达系统参数配置
% =========================================================================
fprintf('--- 2. 配置雷达系统参数 ---\n');

% --- 2.1 基础参数 (根据您的文档) ---
config.Sig_Config.c = 2.99792458e8;             % 光速 (m/s)
config.Sig_Config.fs = 25e6;                    % 基带采样率 (Hz)
config.Sig_Config.fc = 9450e6;                  % 雷达工作中心频率 (Hz)
config.Sig_Config.prtNum = 332;                 % 一帧 PRT 个数
config.Sig_Config.point_PRT = 3404;             % 每PRT总采样点数
config.Sig_Config.channel_num = 16;             % 接收通道数
config.Sig_Config.beam_num = 13;                % DBF后波束数量
config.Sig_Config.prt = 232.76e-6;              % PRT周期 (s)
config.Sig_Config.B = 20e6;                     % 带宽 (Hz)
config.Sig_Config.tao = [0.16e-6, 8e-6, 28e-6]; % 脉宽 [窄, 中, 长] (s)
config.Sig_Config.point_prt_segments = [228, 723, 2453]; % 各脉冲段对应的采样波门点数

% --- 2.2 天线阵列参数 ---
config.Array.element_spacing = 0.0138;          % 阵元间距 (m)

% --- 2.3 派生参数 ---
config.Sig_Config.prf = 1/config.Sig_Config.prt;
config.Sig_Config.wavelength = config.Sig_Config.c / config.Sig_Config.fc;
config.Sig_Config.K2 = -config.Sig_Config.B / config.Sig_Config.tao(2); % 中脉冲调频斜率
config.Sig_Config.K3 = config.Sig_Config.B / config.Sig_Config.tao(3);  % 长脉冲调频斜率
ts = 1 / config.Sig_Config.fs;

%% 3. 生成发射参考波形 (物理真实脉冲)
% =========================================================================
fprintf('--- 3. 生成发射参考波形 ---\n');
% (此部分代码已根据之前的讨论修正)
num_all_prt = round(config.Sig_Config.prt * config.Sig_Config.fs );        % 一个完整PRT信号的总长度
num_samples_1 = round(config.Sig_Config.tao(1) * config.Sig_Config.fs);
num_samples_2 = round(config.Sig_Config.tao(2) * config.Sig_Config.fs);
num_samples_3 = round(config.Sig_Config.tao(3) * config.Sig_Config.fs);
t1 = linspace(-config.Sig_Config.tao(1)/2, config.Sig_Config.tao(1)/2, num_samples_1);
t2 = linspace(-config.Sig_Config.tao(2)/2, config.Sig_Config.tao(2)/2, num_samples_2);
t3 = linspace(-config.Sig_Config.tao(3)/2, config.Sig_Config.tao(3)/2, num_samples_3);

pulse1 = sin(2*pi*t1 + pi/2);
pulse2 = exp(1j*2*pi*(0.5*config.Sig_Config.K2*(t2.^2)));
pulse3 = exp(1j*2*pi*(0.5*config.Sig_Config.K3*(t3.^2)));


% pulse1_prt_s = 
tx_pulse = zeros(1, config.Sig_Config.point_PRT);
tx_pulse(1:num_samples_1) = pulse1;
offset1 = config.Sig_Config.point_prt_segments(1);
tx_pulse(offset1+1 : offset1+num_samples_2) = pulse2;
offset2 = offset1 + config.Sig_Config.point_prt_segments(2);
tx_pulse(offset2+1 : offset2+num_samples_3) = pulse3;

figure;
plot(real(tx_pulse))
title('发射信号实部时域波形图');


%% 4. 模拟16通道目标回波
% =========================================================================
fprintf('--- 4. 模拟16通道目标回波 ---\n');
% (此部分代码保持不变，仅修改了噪声和增益)
raw_iq_data = complex(zeros(config.Sig_Config.prtNum, config.Sig_Config.point_PRT, config.Sig_Config.channel_num));
for m = 1:config.Sig_Config.prtNum
    pulse_echo_all_channels = complex(zeros(config.Sig_Config.point_PRT, config.Sig_Config.channel_num));
    for k = 1:length(targets)
        range = targets(k).Range;
        delay = 2 * range / config.Sig_Config.c;
        delay_samples = round(delay / ts);
        doppler_freq = 2 * targets(k).Velocity / config.Sig_Config.wavelength;
        doppler_phase_shift = exp(1j * 2 * pi * doppler_freq * (m-1) * config.Sig_Config.prt);
        amplitude_gain = 1e11; 
        lambda_sq = config.Sig_Config.wavelength^2;
        % amplitude = amplitude_gain * sqrt(targets(k).RCS * lambda_sq) / (range^2 * (4*pi)^(3/2));
        amplitude = 1;
        target_echo_base = zeros(1, config.Sig_Config.point_PRT);
        if (delay_samples > 0) && (delay_samples < config.Sig_Config.point_PRT)
            len_echo = min(length(tx_pulse), config.Sig_Config.point_PRT - delay_samples);
            target_echo_base(delay_samples+1 : delay_samples+len_echo) = tx_pulse(1:len_echo);
        end
        target_echo_base = amplitude * target_echo_base * doppler_phase_shift;
        phase_shifts_rad = deg2rad(calculate_phase_shifts(targets(k).ElevationAngle, config.Array.element_spacing, config.Sig_Config.wavelength));
        channel_phasors = exp(1j * phase_shifts_rad);
        target_echo_multichannel = target_echo_base.' * channel_phasors;
        pulse_echo_all_channels = pulse_echo_all_channels + target_echo_multichannel;
    end
    raw_iq_data(m, :, :) = reshape(pulse_echo_all_channels, [1, config.Sig_Config.point_PRT, config.Sig_Config.channel_num]);
end
noise_power = 1e-3;
noise = sqrt(noise_power/2) * (randn(size(raw_iq_data)) + 1j * randn(size(raw_iq_data)));
raw_iq_data = raw_iq_data + noise;

figure;
plot(real(raw_iq_data(1,:,2)))
title('16通道信号实部时域波形图');


%% 5. 数字波束形成 (DBF)
% =========================================================================
fprintf('--- 5. 执行数字波束形成 (DBF) ---\n');
% --- 加载DBF系数 ---
DBF_coeffs_data = readmatrix(dbf_coef_path);
DBF_coeffs_data_C = DBF_coeffs_data(:, 1:2:end) + 1j * DBF_coeffs_data(:, 2:2:end);
% --- 逐脉冲DBF ---
iq_data_13beam = complex(zeros(config.Sig_Config.prtNum, config.Sig_Config.point_PRT, config.Sig_Config.beam_num));
for i = 1:config.Sig_Config.prtNum
    single_pulse_16ch = squeeze(raw_iq_data(i, :, :));
    single_pulse_13beam = single_pulse_16ch * DBF_coeffs_data_C';
    iq_data_13beam(i, :, :) = single_pulse_13beam;
end
% --- DBF结果可视化 ---
figure('Name', '阶段1: DBF后波束功率');
beam_to_check = 4; % 选择一个靠近目标的波束来显示
imagesc(db(abs(squeeze(iq_data_13beam(:,:,beam_to_check)))));
title(sprintf('DBF后波束 #%d 的功率 (dB)', beam_to_check));
xlabel('距离单元 (快时间)'); ylabel('脉冲数 (慢时间)'); colorbar;

%% 6. 脉冲压缩 (含偏移量校正)
% =========================================================================
fprintf('--- 6. 执行分段脉冲压缩 (含偏移校正) ---\n');
% --- 加载窄脉冲FIR滤波器系数 ---
fir_coeffs = [794,1403,2143,2672,2591,1711,-58,-2351,-4592,-5855,-5338,-2389,3005,10341,18410,25779,30907,32768,30907,25779,18410,10341,3005,-2389,-5338,-5855,-4592,-2351,-58,1711,2591,2672,2143,1403,794];
fir_coeffs = fir_coeffs/max(fir_coeffs); % 归一化

% --- 分段处理 ---
p1_len = config.Sig_Config.point_prt_segments(1);
p2_len = config.Sig_Config.point_prt_segments(2);
p3_len = config.Sig_Config.point_prt_segments(3);
pc_results_13beam = complex(zeros(size(iq_data_13beam)));

for b = 1:config.Sig_Config.beam_num
    beam_data = squeeze(iq_data_13beam(:,:,b));
    % 1. 分割信号
    signal_01 = beam_data(:, 1:p1_len);
    signal_02 = beam_data(:, p1_len+1 : p1_len+p2_len);
    signal_03 = beam_data(:, p1_len+p2_len+1 : end);
    
    % 2. 独立脉压与偏移校正
    % 窄脉冲 (FIR滤波)
    pc_01 = filter(fir_coeffs, 1, signal_01, [], 2);
    fir_delay = round(mean(grpdelay(fir_coeffs)));
    pc_01 = circshift(pc_01, -fir_delay, 2); % 校正FIR延迟
    
    % 中脉冲 (匹配滤波)
    pc_02 = ifft(fft(signal_02, [], 2) .* repmat(fft(conj(fliplr(pulse2)), p2_len), 332, 1), [], 2);
    pc_02 = circshift(pc_02, -(num_samples_2-1), 2); % 校正脉压偏移
    
    % 长脉冲 (匹配滤波)
    pc_03 = ifft(fft(signal_03, [], 2) .* repmat(fft(conj(fliplr(pulse3)), p3_len), 332, 1), [], 2);
    pc_03 = circshift(pc_03, -(num_samples_3-1), 2); % 校正脉压偏移
    
    % 3. 拼接结果
    pc_results_13beam(:,:,b) = [pc_01, pc_02, pc_03];
end
% --- 脉冲压缩结果可视化 ---
figure('Name', '阶段2: 脉冲压缩后距离剖面');
deltaR = config.Sig_Config.c / (2*config.Sig_Config.fs);
range_axis = (0:config.Sig_Config.point_PRT-1) * deltaR;
plot(range_axis, 20*log10(abs(pc_results_13beam(1,:,beam_to_check)))); % 显示第一个脉冲的脉压结果
grid on; hold on;
for i = 1:length(targets)
    xline(targets(i).Range, '--r', sprintf('预设目标 %d', i));
end
title(sprintf('波束 #%d 脉压后距离-幅度图', beam_to_check));
xlabel('距离 (m)'); ylabel('幅度');

%% 7. MTD 处理
% =========================================================================
fprintf('--- 7. 执行MTD处理 ---\n');
mtd_results_13beam = complex(zeros(size(pc_results_13beam)));
window_mtd = kaiser(config.Sig_Config.prtNum, 8);
for b = 1:config.Sig_Config.beam_num
    pc_beam_data = squeeze(pc_results_13beam(:,:,b));
    mtd_results_13beam(:,:,b) = fftshift(fft(pc_beam_data .* window_mtd, [], 1), 1);
end

% --- MTD结果可视化 (RDM图) ---
figure('Name', '阶段3: MTD后RDM图');
velocity_axis = linspace(-config.Sig_Config.prf/2, config.Sig_Config.prf/2, config.Sig_Config.prtNum) * config.Sig_Config.wavelength / 2;
rdm_to_plot = abs(squeeze(mtd_results_13beam(:,:,beam_to_check)));
imagesc(range_axis, velocity_axis, 20*log10(rdm_to_plot));
title(sprintf('波束 #%d 的距离-多普勒图 (RDM)', beam_to_check));
xlabel('距离 (m)'); ylabel('速度 (m/s)'); colorbar; axis xy;


% % 速度剖面图
% % 1. 创建一个新的图形窗口
% figure('Name', 'MTD后速度剖面图');
% % 2. 找到RDM中能量最强的点，以该点所在的距离单元进行切片
% [~, max_idx] = max(abs(mtd_results_13beam(:,:,beam_to_check))); % 找到能量最大点的线性索引
% [max_v_idx, max_r_idx] = ind2sub(size(mtd_results_13beam(:,:,beam_to_check)), max_idx); % 转换为二维(行,列)索引
% 
% % 3. 提取该距离单元（即RDM矩阵的一整列）的速度维数据
% velocity_profile = abs(mtd_results_13beam(:, max_r_idx, beam_to_check));
% % 4. 绘制速度剖面图
% plot(velocity_axis, 20*log10(velocity_profile)); % 将幅度转换为dB单位显示
% % 5. 美化图形
% grid on;
% title(sprintf('速度维剖面图 @ 距离 = %.2f m', range_axis(max_r_idx)));
% xlabel('速度 (m/s)');
% ylabel('幅度 (dB)');
% hold on;
% peak_velocity = velocity_axis(max_v_idx);
% xline(peak_velocity, '--r', sprintf('峰值速度 @ %.2f m/s', peak_velocity)); % 在峰值速度处画一条垂直虚线
% hold off;





%% 8. CFAR 检测
% =========================================================================
fprintf('--- 8. 执行CFAR检测 ---\n');
rdm_for_cfar = abs(squeeze(mtd_results_13beam(:,:,beam_to_check)));
[cfar_detections, cfar_threshold] = local_cfar_detector(rdm_for_cfar, cfar_params);
% --- CFAR结果可视化 ---
figure('Name', '阶段4: CFAR检测结果');
% 找到最强目标进行剖面分析
[~, max_idx] = max(rdm_for_cfar(:));
[~, max_r_idx] = ind2sub(size(rdm_for_cfar), max_idx);
subplot(2,1,1);
imagesc(range_axis, velocity_axis, cfar_detections);
title(sprintf('波束 #%d CFAR检测结果 (1=目标)', beam_to_check));
xlabel('距离 (m)'); ylabel('速度 (m/s)'); axis xy;
subplot(2,1,2);
plot(velocity_axis, rdm_for_cfar(:, max_r_idx), 'b-', 'DisplayName', '信号强度');
hold on;
plot(velocity_axis, cfar_threshold(:, max_r_idx), 'r--', 'DisplayName', 'CFAR检测门限');
grid on; legend;
title(sprintf('信号与门限对比 @ 距离 %.2f m', range_axis(max_r_idx)));
xlabel('速度 (m/s)'); ylabel('幅度');

%% 本地子函数
% =========================================================================
function channel_phase_shifts_deg = calculate_phase_shifts(arrival_angle_deg, element_spacing_m, wavelength_m)
    arrival_angle_rad = deg2rad(arrival_angle_deg);
    delta_phi_rad = (2 * pi * element_spacing_m * sin(arrival_angle_rad)) / wavelength_m;
    channel_indices = 0:15;
    channel_phase_shifts_rad = channel_indices * delta_phi_rad;
    channel_phase_shifts_deg = rad2deg(channel_phase_shifts_rad);
end

function [detection_map, threshold_map] = local_cfar_detector(rdm_abs, params)
    % 简化的2D-CFAR检测器
    [num_vel, num_range] = size(rdm_abs);
    detection_map = false(num_vel, num_range);
    threshold_map = zeros(num_vel, num_range);
    
    % 速度维CFAR
    for r = 1:num_range
        range_slice = rdm_abs(:, r);
        for v = 1:num_vel
            % 定义参考窗和保护窗
            g_start1 = max(1, v - params.guardCells_V/2);
            g_end1 = v-1;
            r_start1 = max(1, g_start1 - params.refCells_V);
            r_end1 = g_start1-1;
            
            g_start2 = v+1;
            g_end2 = min(num_vel, v + params.guardCells_V/2);
            r_start2 = g_end2 + 1;
            r_end2 = min(num_vel, r_start2 + params.refCells_V - 1);
            
            noise_cells = [range_slice(r_start1:r_end1); range_slice(r_start2:r_end2)];
            if isempty(noise_cells)
                noise_level = 0;
            else
                noise_level = mean(noise_cells);
            end
            
            threshold = noise_level * params.T_CFAR;
            threshold_map(v,r) = threshold;
            
            if rdm_abs(v,r) > threshold
                detection_map(v,r) = true;
            end
        end
    end
end