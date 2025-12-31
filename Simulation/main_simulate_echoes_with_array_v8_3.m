% main_multiframe_simulation_v9_Scan.m (v9.2 - 真实航迹版)
%
% 描述:
%   这是一个多帧雷达仿真主脚本。
%   (v9.1): 增加了伺服机构的模拟 (方位角)。
%   (v9.2): (根据您的反馈) 修正了 "状态演进" 逻辑。
%           现在模拟目标在 *固定高度* 和 *恒定地面速度* 下的
%           匀速直线运动，正确地计算 *瞬时* 俯仰角和径向速度。
%  修改记录
%  date       by      version   modify
%  25/11/13   XZR     V8.3      加入帧间聚类算法和匀速直线运动动目标模拟

clc; clear; close all;
%% 0. 多帧仿真配置
% =========================================================================
total_frames_to_run = 50;  % <--- 在此设置要模拟的总帧数
cumulative_final_log = []; % 用于累积所有帧的最终目标

%% 1. & 2. 一次性设置 (从 v8 复制)
% =========================================================================
fprintf('--- 1. & 2. 正在加载所有配置参数 (来自 v8)... ---\n');

% --- 1.0 (v9.1) 扫描配置 ---
config.scan.rpm = 6; % 天线转速 (RPM, 6 RPM = 36 deg/sec)
config.scan.start_azimuth = 0; % 扫描起始方位 (度)

% --- 1.1 目标 *初始* 状态 (t=0) ---
clear targets;
% (注意) 这些值现在被视为 t=0 时的 *初始测量值*
targets(1).Range = 3000;
targets(1).Velocity = 20;
targets(1).ElevationAngle = 10;
targets(1).SNR_dB = 10;
targets(2).Range = 10000;
targets(2).Velocity = 25;
targets(2).ElevationAngle = 10;
targets(2).SNR_dB = 15;
P_noise_floor = 1;

% --- 1.2 文件路径 ---
base_path = 'C:\Users\a\Desktop\9.8 问题\Simulation'; % <--- !! 检查路径
dbf_coef_path = fullfile(base_path, 'X8数据采集250522_DBFcoef.csv');

% --- 1.3 CFAR 参数 ---
cfar_params.refCells_V = 5;      
cfar_params.guardCells_V = 10;
cfar_params.refCells_R = 5;      
cfar_params.guardCells_R = 10;
cfar_params.T_CFAR = 8;          
cfar_params.method = 'GOCA';
% --- 1.4 聚类参数 (波束内/间) ---
cluster_params.max_range_sep = 30;
cluster_params.max_vel_sep = 0.4;
cluster_params.max_angle_sep = 5.0; % (v8 - Stage 1)

% --- 1.5 帧间聚类 (航迹关联) 参数 ---
config.inter_frame_cluster.enable = true;
% 1. 门限倍数（K值）
config.inter_frame_cluster.K = 1;
% 2. 关联门限
config.inter_frame_cluster.Gate_R = cluster_params.max_range_sep * config.inter_frame_cluster.K;  % (m) 帧内距离门限
config.inter_frame_cluster.Gate_V = cluster_params.max_vel_sep * config.inter_frame_cluster.K;    % (m/s) 帧内速度门限
config.inter_frame_cluster.Gate_El = cluster_params.max_angle_sep * config.inter_frame_cluster.K; % (度) 帧内俯仰角门限
config.inter_frame_cluster.Gate_Az = 10.0;
config.inter_frame_cluster.Max_Frame_Gap = 3; 

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

%% 3. 预计算 (v9.2 升级)
% =========================================================================
fprintf('--- 3. 正在预计算所有滤波器、波形和系数... ---\n');

precomputed_data = struct();

% --- 3.0 (v9.1) 帧时间 和 扫描增量 ---
T_frame = config.Sig_Config.prtNum * config.Sig_Config.prt;
deg_per_sec = config.scan.rpm * (360 / 60);
config.scan.deg_per_frame = deg_per_sec * T_frame;
precomputed_data.T_frame = T_frame;
fprintf('  > 仿真配置: %.1f RPM = %.2f deg/sec = %.3f deg/frame\n', ...
    config.scan.rpm, deg_per_sec, config.scan.deg_per_frame);

% --- 3.1 (v9.2 新增) 航迹物理模型初始化 ---
% (我们在这里向 "targets" 结构体中添加 *不变* 的物理参数)
fprintf('  > 正在初始化 %d 个目标的物理航迹...\n', length(targets));
for k = 1:length(targets)
    R_0 = targets(k).Range;
    El_0 = targets(k).ElevationAngle;
    V_rad_0 = targets(k).Velocity;
    
    % (计算并存储 *恒定* 的物理参数)
    targets(k).const_H = R_0 * sind(El_0);
    targets(k).const_V_ground = V_rad_0 / cosd(El_0);
    
    % (存储 *当前* 状态)
    targets(k).current_R_ground = R_0 * cosd(El_0);
    
    fprintf('    > 目标 %d: 初始 R=%.0fm, El=%.1f deg, V_rad=%.1f m/s\n', k, R_0, El_0, V_rad_0);
    fprintf('    > -> 已计算: 固定 H=%.1f m, 固定 V_ground=%.2f m/s\n', targets(k).const_H, targets(k).const_V_ground);
end

% (复制 %% 3.1 到 %% 3.7 的所有 v9.1 代码, 此处省略以保持简洁)
% ... (v9.1 %% 3.1)
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
% ... (v9.1 %% 3.2)
fir_coeffs = [794,1403,2143,2672,2591,1711,-58,-2351,-4592,-5855,-5338,-2389,3005,10341,18410,25779,30907,32768,30907,25779,18410,10341,3005,-2389,-5338,-5855,-4592,-2351,-58,1711,2591,2672,2143,1403,794];
fir_coeffs = 6 * fir_coeffs/max(fir_coeffs); 
precomputed_data.MF_narrow = fir_coeffs;
precomputed_data.fir_delay = round(mean(grpdelay(fir_coeffs)));
win_medium = kaiser(length(pulse2), 4.5); 
precomputed_data.MF_medium_win = fliplr(conj(pulse2 .* win_medium'));
win_long = kaiser(length(pulse3), 4.5); 
precomputed_data.MF_long_win = fliplr(conj(pulse3 .* win_long'));
% ... (v9.1 %% 3.3)
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
% ... (v9.1 %% 3.4)
precomputed_data.N_gate_narrow = config.Sig_Config.point_prt_segments(1);
precomputed_data.N_gate_medium = config.Sig_Config.point_prt_segments(2);
precomputed_data.N_gate_long = config.Sig_Config.point_prt_segments(3);
precomputed_data.N_total_gate = N_total_gate;
precomputed_data.seg_start_narrow = num_samples_1 + 1;
precomputed_data.seg_start_medium = seg_start_medium;
precomputed_data.seg_start_long = seg_start_long;
% ... (v9.1 %% 3.5)
precomputed_data.MTD_win = kaiser(config.Sig_Config.prtNum, 4.5);
% ... (v9.1 %% 3.6)
v_max = config.Sig_Config.wavelength / (2 * config.Sig_Config.prt);
precomputed_data.velocity_axis = linspace(-v_max/2, v_max/2, config.Sig_Config.prtNum);
precomputed_data.range_axis = (0:N_total_gate-1) * (config.Sig_Config.c / (2 * config.Sig_Config.fs));
precomputed_data.deltaR = config.Sig_Config.c*ts/2; 
precomputed_data.deltaV = v_max / config.Sig_Config.prtNum;
precomputed_data.beam_angles_deg = [-16, -9.6, -3.2, 3.2, 9.6, 16, 22.6, 29.2, 36.1, 43.3, 51, 59.6, 70.3];
precomputed_data.k_slopes_LUT = [-4.6391,-4.6888,-4.7578,-4.7891,-4.7214,-4.7513,-5.2343,-5.4529,-5.7323,-6.1685,-7.0256,-8.7612];
% ... (v9.1 %% 3.7)
try
    DBF_coeffs_data = readmatrix(dbf_coef_path);
    precomputed_data.DBF_coeffs_data_C = DBF_coeffs_data(:, 1:2:end) + 1j * DBF_coeffs_data(:, 2:2:end);
catch E
    fprintf('!!! 致命错误: 无法加载 DBF 系数文件 !!!\n');
    fprintf('请检查路径: %s\n', dbf_coef_path);
    rethrow(E);
end

fprintf('--- 预计算完成。进入主仿真循环。 ---\n');

%% 4. 主循环 (多帧仿真) (v9.2 升级)
% =========================================================================
current_azimuth = config.scan.start_azimuth;
tic;

% (新) 创建一个临时的 `targets` 副本，用于在循环中传递给处理函数
targets_for_processing = targets;

for frame_idx = 1 : total_frames_to_run
    fprintf('\n--- 正在处理第 %d 帧 / 共 %d 帧 ---\n', frame_idx, total_frames_to_run);

    % --- 4.1. 状态演进 (v9.2 真实航迹模型) ---
    fprintf('  (S4.1) 正在更新状态 (帧 %d): Az = %.2f deg\n', frame_idx, current_azimuth);
    
    % a) 演进伺服方位角
    current_azimuth = mod(current_azimuth + config.scan.deg_per_frame, 360);
    
    % b) 演进每个目标的状态
    for k = 1:length(targets)
        % 1. 演进 *基本* 状态 (地面距离)
        R_g_new = targets(k).current_R_ground - (targets(k).const_V_ground * T_frame);
        targets(k).current_R_ground = R_g_new; % 更新 "k" 循环外的状态
        
        % 2. 重新计算 *瞬时* 测量值
        R_new = sqrt(R_g_new^2 + targets(k).const_H^2);
        El_new = asind(targets(k).const_H / R_new);
        V_rad_new = targets(k).const_V_ground * cosd(El_new);
        
        % 3. 更新 *将要被处理* 的结构体
        targets_for_processing(k).Range = R_new;
        targets_for_processing(k).ElevationAngle = El_new;
        targets_for_processing(k).Velocity = V_rad_new;
        % (SNR 保持不变)
        
        fprintf('    > 目标 %d (H=%.0f m): R_g=%.0f m -> R=%.0f m, El=%.2f deg, V_rad=%.2f m/s\n', ...
            k, targets(k).const_H, R_g_new, R_new, El_new, V_rad_new);
    end
    
    % --- 4.2. 运行 v8 的处理核 ---
    % (将 *更新后* 的 targets_for_processing 传入)
    [final_targets] = fun_process_single_frame(targets_for_processing, config, cfar_params, ...
        cluster_params, precomputed_data, frame_idx);
    
    % --- 4.3. 累积结果 (注入 iFrame 和 iAntAngle) ---
    if ~isempty(final_targets)
        num_goals = length(final_targets);
        
        frame_num_cell = num2cell(repmat(frame_idx, 1, num_goals));
        ant_angle_cell = num2cell(repmat(current_azimuth, 1, num_goals));
        
        [final_targets.iFrame] = deal(frame_num_cell{:});
        [final_targets.iAntAngle] = deal(ant_angle_cell{:});
        
        cumulative_final_log = [cumulative_final_log, final_targets];
    end
    
end % (结束主循环)
simulation_time = toc;
fprintf('\n--- 多帧仿真全部完成 ---\n');
fprintf('总共处理了 %d 帧, 累积了 %d 个检测点。\n', total_frames_to_run, length(cumulative_final_log));

%% 5. 帧间聚类 (航迹关联) (v9.2 最终版)
% =========================================================================
if config.inter_frame_cluster.enable
    fprintf('--- 5. 正在执行帧间聚类 (5D 航迹关联)...\n');
    
    % --- 5.1 提取输入数据和门限 ---
    detection_log = cumulative_final_log;
    num_detections = length(detection_log);
    
    Gate_R = config.inter_frame_cluster.Gate_R;
    Gate_V = config.inter_frame_cluster.Gate_V;
    Gate_Az = config.inter_frame_cluster.Gate_Az;
    Gate_El = config.inter_frame_cluster.Gate_El;
    Max_Frame_Gap = config.inter_frame_cluster.Max_Frame_Gap;
    
    if num_detections == 0
        fprintf('  > 总日志为空，无需进行帧间聚类。\n');
        final_tracks_log = [];
    else
        % --- 5.2 聚类算法核心 (BFS) ---
        cluster_IDs_final = zeros(num_detections, 1);
        current_cluster_ID = 0;
        
        for i = 1:num_detections
            if cluster_IDs_final(i) == 0
                current_cluster_ID = current_cluster_ID + 1;
                points_to_visit = i;
                while ~isempty(points_to_visit)
                    current_idx = points_to_visit(1);
                    points_to_visit(1) = [];
                    if cluster_IDs_final(current_idx) == 0
                        cluster_IDs_final(current_idx) = current_cluster_ID;
                        for j = 1:num_detections
                            if cluster_IDs_final(j) == 0
                                % --- 核心: 5D 门限 (R, V, Az, El, Frame) ---
                                dist_r = abs(detection_log(current_idx).Range - detection_log(j).Range);
                                dist_v = abs(detection_log(current_idx).Velocity - detection_log(j).Velocity);
                                dist_az = abs(detection_log(current_idx).iAntAngle - detection_log(j).iAntAngle);
                                dist_el = abs(detection_log(current_idx).Angle - detection_log(j).Angle);
                                dist_frame = abs(detection_log(current_idx).iFrame - detection_log(j).iFrame);
                                
                                if (dist_r <= Gate_R) && (dist_v <= Gate_V) && ...
                                   (dist_az <= Gate_Az) && (dist_el <= Gate_El) && ...
                                   (dist_frame <= Max_Frame_Gap)
                                
                                    points_to_visit(end+1) = j;
                                end
                            end
                        end
                    end
                end
            end
        end
        num_final_clusters = current_cluster_ID;
        fprintf('  > 帧间聚类完成，将 %d 个总检测点合并为 %d 条唯一航迹。\n', num_detections, num_final_clusters);

        % --- 5.3 对每个簇进行合并 (混合归并策略) ---
        final_tracks_log = repmat(struct('Range', 0, 'Velocity', 0, 'Angle', 0, 'Azimuth', 0, 'Power', 0, 'FirstFrame', 0, 'LastFrame', 0, 'NumPoints', 0), num_final_clusters, 1);
        
        for i = 1:num_final_clusters
            cluster_mask = (cluster_IDs_final == i);
            detections_in_cluster = detection_log(cluster_mask);
            
            powers = [detections_in_cluster.Power]';
            total_power = sum(powers);
            
            [max_power, idx_winner] = max(powers);
            winner_detection = detections_in_cluster(idx_winner);
            
            angles_az = [detections_in_cluster.iAntAngle]';
            final_azimuth = sum(angles_az .* powers) / total_power;
            
            frames_in_cluster = [detections_in_cluster.iFrame];
            
            final_tracks_log(i).Range = winner_detection.Range;
            final_tracks_log(i).Velocity = winner_detection.Velocity;
            final_tracks_log(i).Angle = winner_detection.Angle;
            final_tracks_log(i).Azimuth = final_azimuth;
            final_tracks_log(i).Power = max_power;
            final_tracks_log(i).FirstFrame = min(frames_in_cluster);
            final_tracks_log(i).LastFrame = max(frames_in_cluster);
            final_tracks_log(i).NumPoints = length(frames_in_cluster);
        end
    end
else
    fprintf('--- 5. 帧间聚类被禁用，使用未聚类的日志进行绘图 ---\n');
    % (格式转换)
    num_detections = length(cumulative_final_log);
    final_tracks_log = repmat(struct('Range', 0, 'Velocity', 0, 'Angle', 0, 'Azimuth', 0, 'Power', 0, 'FirstFrame', 0, 'LastFrame', 0, 'NumPoints', 0), num_detections, 1);
    for i=1:num_detections
       final_tracks_log(i).Range = cumulative_final_log(i).Range;
       final_tracks_log(i).Velocity = cumulative_final_log(i).Velocity;
       final_tracks_log(i).Angle = cumulative_final_log(i).Angle;
       final_tracks_log(i).Azimuth = cumulative_final_log(i).iAntAngle;
       final_tracks_log(i).Power = cumulative_final_log(i).Power;
       final_tracks_log(i).FirstFrame = cumulative_final_log(i).iFrame;
       final_tracks_log(i).LastFrame = cumulative_final_log(i).iFrame;
       final_tracks_log(i).NumPoints = 1;
    end
end

%% 6. 最终可视化 (基于 "帧间聚类" 的结果)
% =========================================================================
if ~isempty(final_tracks_log)
    fprintf('--- 6. 正在生成最终航迹的可视化图表 ---\n');
    
    all_ranges = [final_tracks_log.Range];
    all_velocities = [final_tracks_log.Velocity];
    all_angles_el = [final_tracks_log.Angle];
    all_angles_az = [final_tracks_log.Azimuth];
    all_num_points = [final_tracks_log.NumPoints];
    
    % --- 绘制 距离-方位角 (PPI) 航迹图 ---
    figure('Name', '最终航迹图 (PPI 视图)');
    polarscatter(deg2rad(all_angles_az), all_ranges, all_num_points*10 + 20, all_velocities, 'filled');

    title(sprintf('最终 %d 条航迹 (距离 vs. 方位角)', length(final_tracks_log)));
    
    % --- 绘制 距离-俯仰角 (RHI) 航迹图 ---
    figure('Name', '最终航迹图 (RHI 视图)');
    scatter(all_ranges, all_angles_el, all_num_points*10 + 20, all_velocities, 'filled');
    xlabel('距离 (m)');
    ylabel('俯仰角 (度)');
    title(sprintf('最终 %d 条航迹 (距离 vs. 俯仰角)', length(final_tracks_log)));
    grid on;
    h_cb = colorbar;
    ylabel(h_cb, '速度 (m/s)');
    
    % --- (v9.2 新增) 绘制 R/El/V_rad vs. 时间 ---
    figure('Name', '目标状态演进 (vs. 帧号)');
    % (找到点数最多的航迹)
    [~, main_track_idx] = max([final_tracks_log.NumPoints]);
    if ~isempty(main_track_idx)
        main_track_id = find([cluster_IDs_final == main_track_idx]);
        main_track_points = cumulative_final_log(main_track_id);
        
        % 提取这条航迹的 *所有* 原始点 (聚类前)
        frames = [main_track_points.iFrame];
        ranges = [main_track_points.Range];
        angles = [main_track_points.Angle];
        vels = [main_track_points.Velocity];
        [frames, sort_idx] = sort(frames); % 确保按时间排序
        
        subplot(3,1,1);
        plot(frames, ranges(sort_idx), 'bo-'); title('主航迹 (距离 vs. 时间)'); ylabel('距离 (m)'); grid on;
        subplot(3,1,2);
        plot(frames, angles(sort_idx), 'ro-'); title('主航迹 (俯仰角 vs. 时间)'); ylabel('俯仰角 (度)'); grid on;
        subplot(3,1,3);
        plot(frames, vels(sort_idx), 'go-'); title('主航迹 (径向速度 vs. 时间)'); ylabel('速度 (m/s)'); grid on;
        xlabel('帧号');
    end
    
else
    fprintf('--- 仿真结束，未检测到任何最终航迹。 ---\n');
end

% --- (可选) 绘制 "聚类前" 和 "聚类后" 的对比图 ---
if config.inter_frame_cluster.enable && ~isempty(cumulative_final_log)
    figure('Name', '帧间聚类对比 (PPI)');
    % (聚类前)
    subplot(1, 2, 1);
    plot_az_pre = [cumulative_final_log.iAntAngle];
    plot_r_pre = [cumulative_final_log.Range];
    polarscatter(deg2rad(plot_az_pre), plot_r_pre, 20, 'r', 'filled', 'MarkerFaceAlpha', 0.5);

    title(sprintf('聚类前 (%d 个检测点)', length(cumulative_final_log)));
    % (聚类后)
    subplot(1, 2, 2);
    plot_az_post = [final_tracks_log.Azimuth];
    plot_r_post = [final_tracks_log.Range];
    plot_size_post = [final_tracks_log.NumPoints]*5 + 20;
    polarscatter(deg2rad(plot_az_post), plot_r_post, plot_size_post, 'b', 'filled');

    title(sprintf('聚类后 (%d 条航迹)', length(final_tracks_log)));
end