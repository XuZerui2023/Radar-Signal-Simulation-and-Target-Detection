% main_plot_snr_vs_angle_error.m
%
% 描述:
%   这是一个基于 v7.6 仿真代码的【完整】蒙特卡洛分析脚本。
%   它通过在不同的信噪比 (SNR) 下重复运行仿真，
%   来测量并绘制单脉冲测角精度随 SNR 变化的曲线。
%
%   (v-Sim 完整版: 包含所有子函数实现)
%
clc; clear; close all;
%% 0. 蒙特卡洛仿真配置
% =========================================================================
fprintf('--- 0. 配置蒙特卡洛仿真 ---\n');
% --- 在此定义要测试的SNR点 (dB) ---
snr_db_vector = -10:2:30;
% --- 在此定义每个SNR点要重复试验的次数 ---
% (100次较快，500次更平滑)
num_trials = 100; 

% --- 定义目标的"真值" (Ground Truth) ---
true_target.Range = 10000;
true_target.Velocity = 20;
true_target.ElevationAngle = 10; % <-- "真" 角度
true_target.RCS = 1;
% 找出这个"真"角度对应的K值 (来自 v7.6 %% 9)
% 角度10度, 位于波束对 [9.6, 16] (索引 5, 6) 之间
true_target.pair_idx = 5; 
k_slopes_LUT = [-4.6391,-4.6888,-4.7578,-4.7891,-4.7214,-4.7513,-5.2343,-5.4529,-5.7323,-6.1685,-7.0256,-8.7612];
true_target.k_slope = k_slopes_LUT(true_target.pair_idx);

% --- 初始化结果存储 ---
angle_error_std = zeros(size(snr_db_vector)); % 存储每个SNR下的角度标准差
detection_probability = zeros(size(snr_db_vector)); % 存储探测概率

%% 1. & 2. 复制 v7.6 的雷达参数配置
% =========================================================================
fprintf('--- 1. & 2. 加载雷达系统参数 (来自 v7.6) ---\n');
% --- 1.1 目标定义 (仅用于S4循环) ---
clear targets;
targets(1).Range = true_target.Range;
targets(1).Velocity = true_target.Velocity;
targets(1).RCS = true_target.RCS;
targets(1).ElevationAngle = true_target.ElevationAngle;
P_noise_floor = 1; 

% --- 1.2 文件路径 ---
% !!! 重要 !!!
% !!! 请确保此路径指向您本地的 DBF 系数文件 !!!
base_path = 'C:\Users\a\Desktop\9.8 问题\Simulation'; 
dbf_coef_path = fullfile(base_path, 'X8数据采集250522_DBFcoef.csv');

% --- 1.3 CFAR 参数 ---
cfar_params.refCells_V = 5;      cfar_params.guardCells_V = 10;
cfar_params.refCells_R = 5;      cfar_params.guardCells_R = 10;
cfar_params.T_CFAR = 8;          cfar_params.method = 'GOCA';
% --- 1.4 聚类参数 ---
cluster_params.max_range_sep = 30;
cluster_params.max_vel_sep = 0.4;
cluster_params.max_angle_sep = 5.0; % (v7.7 第一级聚类)

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
% (重要) 确保 num_all_prt 与 v7.6 匹配
num_all_prt = round(config.Sig_Config.prt * config.Sig_Config.fs ); % 5819
config.Sig_Config.point_PRT = num_all_prt; 
N_total_gate = sum(config.Sig_Config.point_prt_segments); % 3404

%% 3. 预计算: 波形、滤波器、坐标轴、DBF系数
% =========================================================================
fprintf('--- 3. 预计算所有不变的参数... ---\n');
% --- 3.1 生成发射波形 (v7.6 %% 3) ---
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
P_signal_unscaled = mean(abs(tx_pulse(tx_pulse ~= 0)).^2);

% --- 3.2 生成匹配滤波器 (v7.6 %% 6.1) ---
fir_coeffs = [794,1403,2143,2672,2591,1711,-58,-2351,-4592,-5855,-5338,-2389,3005,10341,18410,25779,30907,32768,30907,25779,18410,10341,3005,-2389,-5338,-5855,-4592,-2351,-58,1711,2591,2672,2143,1403,794];
fir_coeffs = 6 * fir_coeffs/max(fir_coeffs); MF_narrow = fir_coeffs;
fir_delay = round(mean(grpdelay(MF_narrow)));
win_medium = kaiser(length(pulse2), 4.5); MF_medium_win = fliplr(conj(pulse2 .* win_medium'));
win_long = kaiser(length(pulse3), 4.5); MF_long_win = fliplr(conj(pulse3 .* win_long'));

% (优化) 预计算频域滤波器，用于S6中的快速PC
L_h_med = length(MF_medium_win);
L_h_long = length(MF_long_win);
% (计算S6中 segment_medium/long 的长度)
gap_duration1_num = gap_duration1 * config.Sig_Config.fs;
gap_duration2_num = gap_duration2 * config.Sig_Config.fs;
seg_start_medium = num_samples_1 + gap_duration1_num + num_samples_2 + 1;
seg_start_long = num_samples_1 + gap_duration1_num + num_samples_2 + gap_duration2_num + num_samples_3 + 1;
L_s_med = num_all_prt - seg_start_medium + 1;
L_s_long = num_all_prt - seg_start_long + 1;
% (计算PC所需的FFT点数)
N_fft_med = 2^nextpow2(L_s_med + L_h_med - 1);
N_fft_long = 2^nextpow2(L_s_long + L_h_long - 1);
% (预计算FFT)
MF_medium_fft = fft(MF_medium_win, N_fft_med, 2);
MF_long_fft = fft(MF_long_win, N_fft_long, 2);

% --- 3.3 生成坐标轴 (v7.6 %% 7 & 9) ---
N_gate_narrow = config.Sig_Config.point_prt_segments(1);
N_gate_medium = config.Sig_Config.point_prt_segments(2);
N_gate_long = config.Sig_Config.point_prt_segments(3);
seg_start_narrow = num_samples_1 + 1;
% (MTD窗)
beta_MTD = 4.5; MTD_win = kaiser(config.Sig_Config.prtNum, beta_MTD);
% (坐标轴)
v_max = config.Sig_Config.wavelength / (2 * config.Sig_Config.prt);
velocity_axis = linspace(-v_max/2, v_max/2, config.Sig_Config.prtNum);
range_axis = (0:N_total_gate-1) * (config.Sig_Config.c / (2 * config.Sig_Config.fs));
deltaR = config.Sig_Config.c*ts/2; 
deltaV = v_max / config.Sig_Config.prtNum;
beam_angles_deg = [-16, -9.6, -3.2, 3.2, 9.6, 16, 22.6, 29.2, 36.1, 43.3, 51, 59.6, 70.3];

% --- 3.4 加载 DBF 系数 (v7.6 %% 5) ---
try
    DBF_coeffs_data = readmatrix(dbf_coef_path);
    DBF_coeffs_data_C = DBF_coeffs_data(:, 1:2:end) + 1j * DBF_coeffs_data(:, 2:2:end);
catch E
    fprintf('!!! 致命错误: 无法加载 DBF 系数文件 !!!\n');
    fprintf('请检查路径: %s\n', dbf_coef_path);
    rethrow(E);
end
fprintf('--- 所有参数和波形准备就绪 ---\n');

%% 4. 开始蒙特卡洛仿真
% =========================================================================
tic;
for i_snr = 1:length(snr_db_vector)
    
    current_snr_db = snr_db_vector(i_snr);
    fprintf('--- 正在执行: SNR = %.1f dB ( %d / %d )...\n', current_snr_db, i_snr, length(snr_db_vector));
    
    % --- 4.1 初始化内层循环的结果 ---
    errors_for_this_snr = NaN(num_trials, 1); % 存储单次角度误差
    detections_count = 0; % 统计探测次数
    
    % --- 4.2 内层循环: 重复试验 ---
    parfor i_trial = 1:num_trials % (使用 parfor 并行计算)
        
        % (S4) 生成回波 (不加噪)
        % (为适配 parfor，我们必须在循环内重新声明关键变量)
        local_targets = targets;
        local_config = config;
        
        raw_iq_data = complex(zeros(local_config.Sig_Config.prtNum, num_all_prt, local_config.Sig_Config.channel_num));
        for m = 1:local_config.Sig_Config.prtNum
            pulse_echo_all_channels = complex(zeros(num_all_prt, local_config.Sig_Config.channel_num));
            for k = 1:length(local_targets)
                delay = 2 * local_targets(k).Range / local_config.Sig_Config.c;
                delay_samples = round(delay / ts);
                doppler_freq = 2 * local_targets(k).Velocity / local_config.Sig_Config.wavelength;
                doppler_phase_shift = exp(1j * 2 * pi * doppler_freq * (m-1) * local_config.Sig_Config.prt);
                
                % (核心) 使用当前外层循环的 SNR
                SNR_k_lin = 10^(current_snr_db / 10);
                P_signal_k = SNR_k_lin * P_noise_floor;
                amplitude = sqrt(P_signal_k / P_signal_unscaled);
                
                target_echo_base = complex(zeros(1, num_all_prt));
                if (delay_samples > 0) && (delay_samples < num_all_prt)
                    len_echo = min(length(tx_pulse), num_all_prt - delay_samples);
                    target_echo_base(delay_samples+1 : delay_samples+len_echo) = tx_pulse(1:len_echo);
                end
                target_echo_base = amplitude * target_echo_base * doppler_phase_shift;
                phase_shifts_rad = deg2rad(calculate_phase_shifts(local_targets(k).ElevationAngle, local_config.Array.element_spacing, local_config.Sig_Config.wavelength));
                channel_phasors = exp(1j * phase_shifts_rad);
                target_echo_multichannel = target_echo_base.' * channel_phasors;
                pulse_echo_all_channels = pulse_echo_all_channels + target_echo_multichannel;
            end
            raw_iq_data(m, :, :) = reshape(pulse_echo_all_channels, [1, num_all_prt, local_config.Sig_Config.channel_num]);
        end
        raw_iq_data_reset0 = raw_iq_data;

        % (S4.1) 加噪 (核心: 每次试验都加 *新* 噪声)
        raw_iq_data_noise = complex(zeros(size(raw_iq_data_reset0)));
        for c = 1:16
            raw_iq_data_temp = squeeze(raw_iq_data_reset0(:,:,c));
            I_noise = randn(size((raw_iq_data_temp)));
            Q_noise = randn(size(raw_iq_data_temp));
            noise = (I_noise + 1j * Q_noise) .* sqrt(P_noise_floor/2); % 缩放到基底功率
            raw_iq_data_noise(:,:,c) = raw_iq_data_temp + noise;
        end
        
        % (S5) DBF
        iq_data_13beam = complex(zeros(local_config.Sig_Config.prtNum, num_all_prt, local_config.Sig_Config.beam_num));
        for i = 1:local_config.Sig_Config.prtNum
            single_pulse_16ch = squeeze(raw_iq_data_noise(i, :, :));
            single_pulse_13beam = single_pulse_16ch * DBF_coeffs_data_C';
            iq_data_13beam(i, :, :) = single_pulse_13beam;
        end
        
        % (S6) 脉冲压缩 (优化版: 使用频域)
        pc_results_13beam = complex(zeros(local_config.Sig_Config.prtNum, N_total_gate, local_config.Sig_Config.beam_num));
        for b = 1:local_config.Sig_Config.beam_num
            beam_data = squeeze(iq_data_13beam(:,:,b));
            segment_narrow = beam_data(:, seg_start_narrow:end);
            segment_medium = beam_data(:, seg_start_medium:end);
            segment_long = beam_data(:, seg_start_long:end);
            
            % 窄脉冲
            pc_out_narrow_full = filter(MF_narrow, 1, segment_narrow, [], 2);
            pc_out_narrow_full_shift = circshift(pc_out_narrow_full, -fir_delay, 2);
            
            % 中脉冲 (快速频域)
            segment_medium_fft = fft(segment_medium, N_fft_med, 2);
            pc_out_medium_full = ifft(segment_medium_fft .* MF_medium_fft, N_fft_med, 2);
            
            % 长脉冲 (快速频域)
            segment_long_fft = fft(segment_long, N_fft_long, 2);
            pc_out_long_full = ifft(segment_long_fft .* MF_long_fft, N_fft_long, 2);
            
            % 拼接
            piece1 = pc_out_narrow_full_shift(:, 1: N_gate_narrow);
            piece2 = pc_out_medium_full(:, N_gate_narrow + 1 : N_gate_narrow + N_gate_medium);
            piece3 = pc_out_long_full(:, N_gate_narrow + N_gate_medium + 1 : N_total_gate);
            pc_results_13beam(:,:,b) = [piece1, piece2, piece3];
        end
        
        % (S7) MTD
        rdm_13beam = complex(zeros(size(pc_results_13beam)));
        for b = 1:local_config.Sig_Config.beam_num
            pc_data_single_beam = squeeze(pc_results_13beam(:,:,b));
            pc_data_single_beam_win = pc_data_single_beam .* MTD_win;
            rdm_13beam(:,:,b) = fftshift(fft(pc_data_single_beam_win, [], 1), 1);
        end
        
        % (S8) CFAR
        [all_raw_detections, rdm_for_cfar_all] = fun_run_goca_cfar_8(rdm_13beam, cfar_params, local_config);
        
        % (S9) 参数估计 (v7.6 逻辑)
        parameterized_detections = fun_parameter_estimation_9(all_raw_detections, rdm_for_cfar_all, rdm_13beam, ...
            range_axis, velocity_axis, deltaR, deltaV, beam_angles_deg, k_slopes_LUT);
        
        % (S10) 波束内聚类 (v7.7 逻辑)
        intra_beam_targets = fun_cluster_stage1_10(parameterized_detections, cluster_params);
        
        % (S11) 波束间聚类 (v7.7 逻辑)
        final_targets = fun_cluster_stage2_11(intra_beam_targets, cluster_params);
        
        % --- 4.3 记录本次试验的结果 ---
        if ~isempty(final_targets)
            % (假设我们只关心检测到的第一个目标)
            if ~isempty(final_targets)
                % (一个更鲁棒的方法是循环 final_targets，找到 R/V/A 最接近真值的那个)
                measured_angle = final_targets(1).Angle;
                errors_for_this_snr(i_trial) = measured_angle - true_target.ElevationAngle;
                detections_count = detections_count + 1;
            end
        end
        
    end % (结束内层循环 i_trial)
    
    % --- 4.4 计算这个SNR点的统计数据 ---
    angle_error_std(i_snr) = std(errors_for_this_snr, 'omitnan');
    detection_probability(i_snr) = detections_count / num_trials;
    
    fprintf('  > SNR %.1f dB 完成: 探测率 = %.1f%%, 角度误差 STDEV = %.4f 度\n', ...
        current_snr_db, detection_probability(i_snr)*100, angle_error_std(i_snr));
    
end % (结束外层循环 i_snr)
simulation_time = toc;
fprintf('--- 蒙特卡洛仿真总耗时: %.2f 秒 ---\n', simulation_time);

%% 5. 绘制最终结果图
% =========================================================================
fprintf('--- 5. 仿真完成，正在绘图 ---\n');

% --- 5.1 绘制测角误差 vs. SNR ---
figure('Name', '测角误差 vs. 信噪比 (蒙特卡洛仿真)', 'Position', [100, 100, 800, 1000]);
subplot(2, 1, 1);
plot(snr_db_vector, angle_error_std, 'bo-', 'LineWidth', 2, 'DisplayName', 'v7.6 仿真结果 (样条插值 R/V, 整数索引 A)');
hold on;

% --- (可选) 绘制理论曲线进行对比 ---
% 理论公式: Error_Std = (BeamWidth / k_mono) / sqrt(SNR_linear)
% 简化版 (您 v7.6 使用的): Error_Std = (k_slope * sqrt(2)) / sqrt(SNR_linear)
snr_linear_vector = 10.^(snr_db_vector / 10);
% (注意) k_slope 在 v7.6 中是 (A-B)/(A+B)，理论公式需要 *2*
theoretical_error_std = (abs(true_target.k_slope) * sqrt(2)) ./ sqrt(snr_linear_vector); 
plot(snr_db_vector, theoretical_error_std, 'r--', 'LineWidth', 2, 'DisplayName', '理论极限 (k * sqrt(2) / sqrt(SNR))');

grid on;
xlabel('信噪比 (dB)');
ylabel('测角误差 (标准差, 度)');
title(sprintf('单脉冲测角精度 (蒙特卡洛 %d 次试验)', num_trials));
legend('show', 'Location', 'southwest');
set(gca, 'FontSize', 12);
ylim([0, max(max(angle_error_std), max(theoretical_error_std)) * 1.2]); % 动态 Y 轴

% --- 5.2 绘制探测概率 vs. SNR ---
subplot(2, 1, 2);
plot(snr_db_vector, detection_probability * 100, 'ms-', 'LineWidth', 2);
grid on;
xlabel('信噪比 (dB)');
ylabel('探测概率 (Pd, %)');
title(sprintf('目标探测概率 (CFAR T=%d)', cfar_params.T_CFAR));
ylim([-5, 105]);
set(gca, 'FontSize', 12);

fprintf('--- 绘图完成 ---\n');


%% 6. 本地子函数 (v7.6/v7.7 算法的完整实现)
% =========================================================================

% (S4 - 本地子函数)
function channel_phase_shifts_deg = calculate_phase_shifts(arrival_angle_deg, element_spacing_m, wavelength_m)
    arrival_angle_rad = deg2rad(arrival_angle_deg);
    delta_phi_rad = (2 * pi * element_spacing_m * sin(arrival_angle_rad)) / wavelength_m;
    channel_indices = 0:15;
    channel_phase_shifts_rad = channel_indices * delta_phi_rad;
    channel_phase_shifts_deg = rad2deg(channel_phase_shifts_rad);
end

% (S8 - v7.6 逻辑)
function [all_raw_detections, rdm_for_cfar_all] = fun_run_goca_cfar_8(rdm_13beam, cfar_params, config)
    
    all_raw_detections = []; 
    rdm_for_cfar_all = complex(zeros(config.Sig_Config.prtNum, config.Sig_Config.point_prt_segments(1)+config.Sig_Config.point_prt_segments(2)+config.Sig_Config.point_prt_segments(3), config.Sig_Config.beam_num-1));
    
    T_CFAR = cfar_params.T_CFAR;
    guard_R = cfar_params.guardCells_R; guard_V = cfar_params.guardCells_V;
    ref_R = cfar_params.refCells_R; ref_V = cfar_params.refCells_V;

    for pair_idx = 1:(config.Sig_Config.beam_num - 1)
        beam_idx_A = pair_idx;
        beam_idx_B = pair_idx + 1;
        rdm_beamA = abs(squeeze(rdm_13beam(:,:,beam_idx_A)));
        rdm_beamB = abs(squeeze(rdm_13beam(:,:,beam_idx_B)));
        rdm_for_cfar = rdm_beamA + rdm_beamB;
        rdm_for_cfar_all(:,:,pair_idx) = rdm_for_cfar;
        
        [num_V, num_R] = size(rdm_for_cfar);
        cfar_detections = zeros(num_V, num_R);
        
        for r = (ref_R + guard_R + 1) : (num_R - ref_R - guard_R)
            for v = (ref_V + guard_V + 1) : (num_V - ref_V - guard_V)
                cut_power = rdm_for_cfar(v, r);
                
                % 距离维 GOCA
                range_win_leading = rdm_for_cfar(v, r - guard_R - ref_R : r - guard_R - 1);
                range_win_trailing = rdm_for_cfar(v, r + guard_R + 1 : r + guard_R + ref_R);
                noise_R = max(mean(range_win_leading), mean(range_win_trailing));
                
                % 速度维 GOCA
                doppler_win_leading = rdm_for_cfar(v - guard_V - ref_V : v - guard_V - 1, r);
                doppler_win_trailing = rdm_for_cfar(v + guard_V + 1 : v + guard_V + ref_V, r);
                noise_V = max(mean(doppler_win_leading), mean(doppler_win_trailing));
                
                noise_estimate = max(noise_R, noise_V);
                threshold = T_CFAR * noise_estimate;
                
                if cut_power > threshold
                    cfar_detections(v, r) = 1;
                end
            end
        end
        
        [detected_v_idx, detected_r_idx] = find(cfar_detections);
        for i = 1:length(detected_v_idx)
            v_idx = detected_v_idx(i);
            r_idx = detected_r_idx(i);
            rdm_abs = rdm_for_cfar(v_idx, r_idx);
            all_raw_detections(end+1, :) = [v_idx, r_idx, pair_idx, rdm_abs];
        end
    end
end

% (S9 - v7.6 逻辑)
function parameterized_detections = fun_parameter_estimation_9(all_raw_detections, rdm_for_cfar_all, rdm_13beam, range_axis, velocity_axis, deltaR, deltaV, beam_angles_deg, k_slopes_LUT)
    
    num_total_raw_detections = size(all_raw_detections, 1);
    if num_total_raw_detections == 0
        parameterized_detections = [];
        return;
    end
    
    [num_V, num_R] = size(squeeze(rdm_13beam(:,:,1)));
    
    % --- 插值参数 (来自 v7.6 %% 9) ---
    extraDots = 2; rInterpTimes = 8; vInterpTimes = 4;
    
    parameterized_detections = repmat(struct('Range', 0, 'Velocity', 0, 'Angle', 0, 'Power', 0, 'PairIndex', 0), num_total_raw_detections, 1);
    
    for i = 1:num_total_raw_detections
        v_idx = all_raw_detections(i, 1);   % 整数索引
        r_idx = all_raw_detections(i, 2);   % 整数索引
        pair_idx = all_raw_detections(i, 3);
        power = all_raw_detections(i, 4);
        
        rdm_interp = rdm_for_cfar_all(:,:,pair_idx);
        
        % --- 距离维插值 (v7.6 %% 9 逻辑) ---
        rCellsFix = (r_idx - extraDots) : (r_idx + extraDots);
        rCellsFix = rCellsFix(rCellsFix >= 1 & rCellsFix <= num_R);
        if length(rCellsFix) < 3
            rCellMax = r_idx;
        else
            rCellsFix = sort(rCellsFix);
            mtdDataUsed_r = rdm_interp(v_idx, rCellsFix);
            rCellsFixQ = rCellsFix(1): 1/rInterpTimes :rCellsFix(end);
            mtdDataUsedQ_r = interp1(rCellsFix-rCellsFix(1),mtdDataUsed_r,rCellsFixQ-rCellsFix(1),'spline');
            [~,I1] = max(mtdDataUsedQ_r);
            rCellMax = rCellsFixQ(I1(1));
        end
        est_range = range_axis(r_idx) + (rCellMax-r_idx)*deltaR;

        % --- 速度维插值 (v7.6 %% 9 逻辑) ---
        vCellsFix = (v_idx - extraDots) : (v_idx + extraDots);
        vCellsFix = vCellsFix(vCellsFix >= 1 & vCellsFix <= num_V);
        if length(vCellsFix) < 3
            vCellMax = v_idx;
        else
            vCellsFix = sort(vCellsFix);
            mtdDataUsed_v = rdm_interp(vCellsFix,r_idx);
            vCellsFixQ = vCellsFix(1): 1/vInterpTimes :vCellsFix(end);
            mtdDataUsedQ_v = interp1(vCellsFix-vCellsFix(1),mtdDataUsed_v,vCellsFixQ-vCellsFix(1),'spline');
            [~,I2] = max(mtdDataUsedQ_v);
            vCellMax = vCellsFixQ(I2(1));
        end
        % (v7.6 修正后的速度计算)
        est_velocity = velocity_axis(v_idx) + (vCellMax-v_idx)*deltaV;
        
        % --- 角度估计 (v7.6 %% 9 逻辑 - *包含缺陷*) ---
        % (缺陷: 使用 v_idx, r_idx 整数索引，而不是 rCellMax, vCellMax)
        S_A = rdm_13beam(v_idx, r_idx, pair_idx);
        S_B = rdm_13beam(v_idx, r_idx, pair_idx + 1);
        
        monopulse_ratio = (S_A - S_B) / (S_A + S_B + eps); % (v7.6 使用了复数比)
        k_slope = k_slopes_LUT(pair_idx);
        angle_A = beam_angles_deg(pair_idx);
        angle_B = beam_angles_deg(pair_idx + 1);
        angle_offset = k_slope * real(monopulse_ratio);
        est_angle = (angle_A + angle_B)/2 + angle_offset;
        
        % --- 存储结果 ---
        parameterized_detections(i).Range = est_range;
        parameterized_detections(i).Velocity = est_velocity;
        parameterized_detections(i).Angle = est_angle;
        parameterized_detections(i).Power = power;
        parameterized_detections(i).PairIndex = pair_idx;
    end
end

% (S10 - v7.7 逻辑)
function intra_beam_targets = fun_cluster_stage1_10(parameterized_detections, cluster_params)
    
    num_total_raw_detections = length(parameterized_detections);
    if num_total_raw_detections == 0
        intra_beam_targets = [];
        return;
    end
    
    cluster_IDs_stage1 = zeros(num_total_raw_detections, 1);
    current_cluster_ID = 0;
    
    for i = 1:num_total_raw_detections
        if cluster_IDs_stage1(i) == 0
            current_cluster_ID = current_cluster_ID + 1;
            points_to_visit = i;
            while ~isempty(points_to_visit)
                current_idx = points_to_visit(1);
                points_to_visit(1) = [];
                if cluster_IDs_stage1(current_idx) == 0
                    cluster_IDs_stage1(current_idx) = current_cluster_ID;
                    for j = 1:num_total_raw_detections
                        if cluster_IDs_stage1(j) == 0
                            dist_r = abs(parameterized_detections(current_idx).Range - parameterized_detections(j).Range);
                            dist_v = abs(parameterized_detections(current_idx).Velocity - parameterized_detections(j).Velocity);
                            dist_a = abs(parameterized_detections(current_idx).Angle - parameterized_detections(j).Angle);
                            
                            if dist_r <= cluster_params.max_range_sep && dist_v <= cluster_params.max_vel_sep && dist_a <= cluster_params.max_angle_sep
                                points_to_visit(end+1) = j;
                            end
                        end
                    end
                end
            end
        end
    end
    num_clusters_stage1 = current_cluster_ID;
    
    % --- 合并: 功率加权平均 ---
    intra_beam_targets = repmat(struct('Range', 0, 'Velocity', 0, 'Angle', 0, 'Power', 0), num_clusters_stage1, 1);
    for i = 1:num_clusters_stage1
        cluster_mask = (cluster_IDs_stage1 == i);
        detections_in_cluster = parameterized_detections(cluster_mask);
        powers = [detections_in_cluster.Power]';
        total_power = sum(powers);
        
        intra_beam_targets(i).Range = sum([detections_in_cluster.Range]' .* powers) / total_power;
        intra_beam_targets(i).Velocity = sum([detections_in_cluster.Velocity]' .* powers) / total_power;
        intra_beam_targets(i).Angle = sum([detections_in_cluster.Angle]' .* powers) / total_power;
        intra_beam_targets(i).Power = total_power;
    end
end

% (S11 - v7.7 逻辑)
function final_targets = fun_cluster_stage2_11(intra_beam_targets, cluster_params)

    num_intra_beam_targets = length(intra_beam_targets);
    if num_intra_beam_targets == 0
        final_targets = [];
        return;
    end

    cluster_IDs_stage2 = zeros(num_intra_beam_targets, 1);
    current_cluster_ID = 0;
    
    for i = 1:num_intra_beam_targets
        if cluster_IDs_stage2(i) == 0
            current_cluster_ID = current_cluster_ID + 1;
            points_to_visit = i;
            while ~isempty(points_to_visit)
                current_idx = points_to_visit(1);
                points_to_visit(1) = [];
                if cluster_IDs_stage2(current_idx) == 0
                    cluster_IDs_stage2(current_idx) = current_cluster_ID;
                    for j = 1:num_intra_beam_targets
                        if cluster_IDs_stage2(j) == 0
                            % --- 核心：只基于 (R, V) 分组 ---
                            dist_r = abs(intra_beam_targets(current_idx).Range - intra_beam_targets(j).Range);
                            dist_v = abs(intra_beam_targets(current_idx).Velocity - intra_beam_targets(j).Velocity);
                            
                            if dist_r <= cluster_params.max_range_sep && dist_v <= cluster_params.max_vel_sep
                                points_to_visit(end+1) = j;
                            end
                        end
                    end
                end
            end
        end
    end
    num_clusters_stage2 = current_cluster_ID;
    
    % --- 合并: "赢家通吃" ---
    final_targets = repmat(struct('Range', 0, 'Velocity', 0, 'Angle', 0, 'Power', 0), num_clusters_stage2, 1);
    for i = 1:num_clusters_stage2
        cluster_mask = (cluster_IDs_stage2 == i);
        detections_in_cluster = intra_beam_targets(cluster_mask);
        
        powers = [detections_in_cluster.Power]';
        [max_power, idx_winner] = max(powers);
        winner_detection = detections_in_cluster(idx_winner);
        
        final_targets(i).Range = winner_detection.Range;
        final_targets(i).Velocity = winner_detection.Velocity;
        final_targets(i).Angle = winner_detection.Angle;
        final_targets(i).Power = winner_detection.Power;
    end
end