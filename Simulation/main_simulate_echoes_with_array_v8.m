% main_multiframe_simulation_v9.m
%
% 描述:
%   这是一个多帧雷达仿真主脚本。
%   它负责初始化系统参数、模拟目标随时间的运动，并循环调用
%   一个“处理核”函数 (fun_process_single_frame) 来处理每一帧。
%
% 架构:
%   %% 0. 多帧仿真配置 (要跑多少帧?)
%   %% 1-3. 一次性设置 (加载 v8 的所有参数和滤波器)
%   %% 4.   主循环 (for frame_idx = 1:N)
%           4.1. 状态演进 (目标移动)
%           4.2. 调用处理核 (处理当前帧)
%           4.3. 累积结果 (保存航迹)
%   %% 5.   最终可视化 (绘制航迹)
%
clc; clear; close all;
%% 0. 多帧仿真配置
% =========================================================================
total_frames_to_run = 3;  % <--- 在此设置要模拟的总帧数
cumulative_final_log = []; % 用于累积所有帧的最终目标

%% 1. & 2. 一次性设置 (从 v8 复制)
% =========================================================================
fprintf('--- 1. & 2. 正在加载所有配置参数 (来自 v8)... ---\n');

% --- 1.1 目标 *初始* 状态 (t=0) ---
clear targets;
targets(1).Range = 3000;
targets(1).Velocity = 15;
targets(1).ElevationAngle = 10;
targets(1).SNR_dB = 10;

targets(2).Range = 10000;
targets(2).Velocity = 20;
targets(2).ElevationAngle = 10;
targets(2).SNR_dB = 15;
P_noise_floor = 1;

% --- 1.2 文件路径 ---
base_path = 'C:\Users\a\Desktop\9.8 问题\Simulation'; % <--- !! 检查路径
dbf_coef_path = fullfile(base_path, 'X8数据采集250522_DBFcoef.csv');

% --- 1.3 CFAR 参数 ---
cfar_params.refCells_V = 5;      cfar_params.guardCells_V = 10;
cfar_params.refCells_R = 5;      cfar_params.guardCells_R = 10;
cfar_params.T_CFAR = 8;          cfar_params.method = 'GOCA';
% --- 1.4 聚类参数 ---
cluster_params.max_range_sep = 30;
cluster_params.max_vel_sep = 0.4;
cluster_params.max_angle_sep = 5.0; % (v8 - Stage 1)

% --- 2.1 & 2.3 基础与派生参数 ---
config.Sig_Config.c = 2.99792458e8;
config.Sig_Config.fs = 25e6;
config.Sig_Config.fc = 9450e6;
config.Sig_Config.prtNum = 332;
config.Sig_Config.prt = 232.76e-6;
config.Sig_Config.B = 20e6;
config.Sig_Config.tao = [0.16e-6, 8e-6, 28e-6];
config.Sig_Config.gap_duration = [11.4e-6, 31.8e-6, 153.4e-6];
config.Sig_Config.point_prt_segments = [228, 723, 2453];
config.Sig_Config.channel_num = 16;
config.Sig_Config.beam_num = 13;
config.Array.element_spacing = 0.0138;
config.Sig_Config.wavelength = config.Sig_Config.c / config.Sig_Config.fc;
ts = 1 / config.Sig_Config.fs;
num_all_prt = round(config.Sig_Config.prt * config.Sig_Config.fs ); % 5819
config.Sig_Config.point_PRT = num_all_prt; 
N_total_gate = sum(config.Sig_Config.point_prt_segments); % 3404

%% 3. 预计算 (从 v8 提取)
% =========================================================================
fprintf('--- 3. 正在预计算所有滤波器、波形和系数... ---\n');

% (新) 创建一个结构体来存放所有预计算的数据
precomputed_data = struct();

% --- 3.1 生成发射波形 (v8 %% 3) ---
tau1 = config.Sig_Config.tao(1); tau2 = config.Sig_Config.tao(2); tau3 = config.Sig_Config.tao(3);
gap_duration1 = config.Sig_Config.gap_duration(1); gap_duration2 = config.Sig_Config.gap_duration(2);
k2 = -config.Sig_Config.B/tau2; k3 = config.Sig_Config.B/tau3;
num_samples_1 = round(tau1 * config.Sig_Config.fs);
num_samples_2 = round(tau2 * config.Sig_Config.fs);
num_samples_3 = round(tau3 * config.Sig_Config.fs);
t2 = linspace(-tau2/2, tau2/2, num_samples_2);
t3 = linspace(-tau3/2, tau3/2, num_samples_3);
pulse1 = ones(1, num_samples_1);
pulse2 = exp(1j*2*pi*(0.5*k2*(t2.^2)));
pulse3 = exp(1j*2*pi*(0.5*k3*(t3.^2)));
tx_pulse = complex(zeros(1, num_all_prt));
tx_pulse(1:num_samples_1) = pulse1;
offset1 = round((tau1+gap_duration1) * config.Sig_Config.fs);
tx_pulse(offset1+1 : offset1+num_samples_2) = pulse2;
offset2 = offset1 + round((tau2+gap_duration2) * config.Sig_Config.fs);
tx_pulse(offset2+1 : offset2+num_samples_3) = pulse3;
precomputed_data.tx_pulse = tx_pulse;
precomputed_data.P_signal_unscaled = mean(abs(tx_pulse(tx_pulse ~= 0)).^2);

% --- 3.2 生成匹配滤波器 (v8 %% 6.1) ---
fir_coeffs = [794,1403,2143,2672,2591,1711,-58,-2351,-4592,-5855,-5338,-2389,3005,10341,18410,25779,30907,32768,30907,25779,18410,10341,3005,-2389,-5338,-5855,-4592,-2351,-58,1711,2591,2672,2143,1403,794];
fir_coeffs = 6 * fir_coeffs/max(fir_coeffs); 
precomputed_data.MF_narrow = fir_coeffs;
precomputed_data.fir_delay = round(mean(grpdelay(fir_coeffs)));

win_medium = kaiser(length(pulse2), 4.5); 
precomputed_data.MF_medium_win = fliplr(conj(pulse2 .* win_medium'));
win_long = kaiser(length(pulse3), 4.5); 
precomputed_data.MF_long_win = fliplr(conj(pulse3 .* win_long'));

% --- 3.3 (优化) 预计算频域滤波器 ---
L_h_med = length(precomputed_data.MF_medium_win);
L_h_long = length(precomputed_data.MF_long_win);
gap_duration1_num = gap_duration1 * config.Sig_Config.fs;
gap_duration2_num = gap_duration2 * config.Sig_Config.fs;
seg_start_medium = num_samples_1 + gap_duration1_num + num_samples_2 + 1;
seg_start_long = num_samples_1 + gap_duration1_num + num_samples_2 + gap_duration2_num + num_samples_3 + 1;
L_s_med = num_all_prt - seg_start_medium + 1;
L_s_long = num_all_prt - seg_start_long + 1;
precomputed_data.N_fft_med = 2^nextpow2(L_s_med + L_h_med - 1);
precomputed_data.N_fft_long = 2^nextpow2(L_s_long + L_h_long - 1);
precomputed_data.MF_medium_fft = fft(precomputed_data.MF_medium_win, precomputed_data.N_fft_med, 2);
precomputed_data.MF_long_fft = fft(precomputed_data.MF_long_win, precomputed_data.N_fft_long, 2);

% --- 3.4 预计算拼接参数 (v8 %% 6.2) ---
precomputed_data.N_gate_narrow = config.Sig_Config.point_prt_segments(1);
precomputed_data.N_gate_medium = config.Sig_Config.point_prt_segments(2);
precomputed_data.N_gate_long = config.Sig_Config.point_prt_segments(3);
precomputed_data.N_total_gate = N_total_gate;
precomputed_data.seg_start_narrow = num_samples_1 + 1;
precomputed_data.seg_start_medium = seg_start_medium;
precomputed_data.seg_start_long = seg_start_long;

% --- 3.5 预计算 MTD 窗 (v8 %% 7) ---
precomputed_data.MTD_win = kaiser(config.Sig_Config.prtNum, 4.5);

% --- 3.6 预计算坐标轴 (v8 %% 9) ---
k_slopes_LUT = [-4.6391,-4.6888,-4.7578,-4.7891,-4.7214,-4.7513,-5.2343,-5.4529,-5.7323,-6.1685,-7.0256,-8.7612]; % 测角K值
v_max = config.Sig_Config.wavelength / (2 * config.Sig_Config.prt);
precomputed_data.velocity_axis = linspace(-v_max/2, v_max/2, config.Sig_Config.prtNum);
precomputed_data.range_axis = (0:N_total_gate-1) * (config.Sig_Config.c / (2 * config.Sig_Config.fs));
precomputed_data.deltaR = config.Sig_Config.c*ts/2; 
precomputed_data.deltaV = v_max / config.Sig_Config.prtNum;
precomputed_data.beam_angles_deg = [-16, -9.6, -3.2, 3.2, 9.6, 16, 22.6, 29.2, 36.1, 43.3, 51, 59.6, 70.3];
precomputed_data.k_slopes_LUT = k_slopes_LUT; 

% --- 3.7 加载 DBF 系数 (v8 %% 5) ---
try
    DBF_coeffs_data = readmatrix(dbf_coef_path);
    precomputed_data.DBF_coeffs_data_C = DBF_coeffs_data(:, 1:2:end) + 1j * DBF_coeffs_data(:, 2:2:end);
catch E
    fprintf('!!! 致命错误: 无法加载 DBF 系数文件 !!!\n');
    fprintf('请检查路径: %s\n', dbf_coef_path);
    rethrow(E);
end

fprintf('--- 预计算完成。进入主仿真循环。 ---\n');

%% 4. 主循环 (多帧仿真)
% =========================================================================
T_frame = config.Sig_Config.prtNum * config.Sig_Config.prt; % 一帧的持续时间 (约 0.077 秒)
tic;

for frame_idx = 1 : total_frames_to_run
    fprintf('\n--- 正在处理第 %d 帧 / 共 %d 帧 ---\n', frame_idx, total_frames_to_run);

    % --- 4.1. 状态演进 (State Evolution) ---
    % 在生成回波之前，更新每个目标的物理位置
    % (我们只更新距离，速度保持不变)
    for k = 1:length(targets)
        % (注意) v 为正 (靠近) 时，距离 R 应该 *减小*
        targets(k).Range = targets(k).Range - (targets(k).Velocity * T_frame);
    end
    
    % --- 4.2. 运行 v8 的处理核 ---
    % (此函数包含了 v8 的 %% 4 到 %% 11 的所有逻辑)
    [final_targets] = fun_process_single_frame(targets, config, cfar_params, ...
        cluster_params, precomputed_data, frame_idx);
    
    % --- 4.3. 累积结果 (注入 iFrame) ---
    if ~isempty(final_targets)
        % (使用我们为实测项目写的 "iFrame 注入" 逻辑)
        num_goals = length(final_targets);
        frame_num_cell = num2cell(repmat(frame_idx, 1, num_goals));
        [final_targets.iFrame] = deal(frame_num_cell{:});
        
        cumulative_final_log = [cumulative_final_log, final_targets];
    end
    
end % (结束主循环)
simulation_time = toc;
fprintf('\n--- 多帧仿真全部完成 ---\n');
fprintf('总共处理了 %d 帧, 累积了 %d 个检测点。\n', total_frames_to_run, length(cumulative_final_log));
fprintf('仿真总耗时: %.2f 秒 (平均每帧 %.2f 秒)。\n', simulation_time, simulation_time / total_frames_to_run);

%% 5. 最终可视化 (绘制航迹)
% =========================================================================
if ~isempty(cumulative_final_log)
    fprintf('--- 正在生成可视化航迹图 ---\n');
    
    % 提取所有目标的 R, V, A 和 帧号
    all_ranges = [cumulative_final_log.Range];
    all_velocities = [cumulative_final_log.Velocity];
    all_angles = [cumulative_final_log.Angle];
    all_frames = [cumulative_final_log.iFrame];
    
    % --- 绘制 距离-帧号 (时间) 航迹图 ---
    figure('Name', '多帧仿真航迹 (距离 vs. 时间)');
    scatter(all_frames, all_ranges, 30, all_velocities, 'filled');
    xlabel('帧号 (iFrame)');
    ylabel('距离 (m)');
    title('目标航迹 (距离 vs. 时间)');
    grid on;
    set(gca, 'XLim', [0, total_frames_to_run+1]);
    h_cb = colorbar;
    ylabel(h_cb, '速度 (m/s)');
    
    % --- 绘制 距离-角度 航迹图 ---
    figure('Name', '多帧仿真航迹 (距离 vs. 角度)');
    scatter(all_ranges, all_angles, 30, all_frames, 'filled');
    xlabel('距离 (m)');
    ylabel('俯仰角 (度)');
    title('目标航迹 (距离 vs. 角度)');
    grid on;
    h_cb = colorbar;
    ylabel(h_cb, '帧号 (iFrame)');
    
else
    fprintf('--- 仿真结束，未检测到任何目标。 ---\n');
end



