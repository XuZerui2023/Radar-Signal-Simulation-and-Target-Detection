% fun_process_single_frame.m
%
% 描述:
%   这是 v8 仿真的“处理核”函数。
%   它接收 *当前* 帧的目标状态，以及预计算好的所有参数/滤波器，
%   然后执行 v8 的完整处理链 (%% 4 到 %% 11)，
%   最终返回该帧的 `final_targets` 列表。
%
% 优化:
%   - %% 6 (脉冲压缩) 已被替换为频域 FFT 方法，以提高速度。
%   - 所有本地子函数已包含在内。
%
function final_targets = fun_process_single_frame(targets, config, cfar_params, cluster_params, precomputed_data, frame_idx)

% --- 从 config 和 precomputed_data 中解包所需变量 ---
P_noise_floor = 1; % (来自 v8 %% 1)
num_all_prt = config.Sig_Config.point_PRT;
ts = 1 / config.Sig_Config.fs;

% (来自 precomputed_data)
tx_pulse = precomputed_data.tx_pulse;
P_signal_unscaled = precomputed_data.P_signal_unscaled;
DBF_coeffs_data_C = precomputed_data.DBF_coeffs_data_C;
MF_narrow = precomputed_data.MF_narrow;
fir_delay = precomputed_data.fir_delay;
MF_medium_fft = precomputed_data.MF_medium_fft;
MF_long_fft = precomputed_data.MF_long_fft;
N_fft_med = precomputed_data.N_fft_med;
N_fft_long = precomputed_data.N_fft_long;
N_gate_narrow = precomputed_data.N_gate_narrow;
N_gate_medium = precomputed_data.N_gate_medium;
N_gate_long = precomputed_data.N_gate_long;
N_total_gate = precomputed_data.N_total_gate;
seg_start_narrow = precomputed_data.seg_start_narrow;
seg_start_medium = precomputed_data.seg_start_medium;
seg_start_long = precomputed_data.seg_start_long;
MTD_win = precomputed_data.MTD_win;
range_axis = precomputed_data.range_axis;
velocity_axis = precomputed_data.velocity_axis;
deltaR = precomputed_data.deltaR;
deltaV = precomputed_data.deltaV;
beam_angles_deg = precomputed_data.beam_angles_deg;
k_slopes_LUT = precomputed_data.k_slopes_LUT;

% --- (v8 - %% 4) 生成16通道目标回波信号 ---
fprintf('  (S4) 正在生成第 %d 帧的回波与噪声...\n', frame_idx);
raw_iq_data = complex(zeros(config.Sig_Config.prtNum, num_all_prt, config.Sig_Config.channel_num));
for m = 1:config.Sig_Config.prtNum 
    pulse_echo_all_channels = complex(zeros(num_all_prt, config.Sig_Config.channel_num));
    for k = 1:length(targets)
        range = targets(k).Range;
        velocity = targets(k).Velocity;
        elevation_angle = targets(k).ElevationAngle;
        
        delay = 2 * range / config.Sig_Config.c;
        delay_samples = round(delay / ts);
        doppler_freq = 2 * velocity / config.Sig_Config.wavelength;
        doppler_phase_shift = exp(1j * 2 * pi * doppler_freq * (m-1) * config.Sig_Config.prt);
        
        % (使用 v8 %% 4 的 SNR 逻辑)
        SNR_k_lin = 10^(targets(k).SNR_dB / 10);
        P_signal_k = SNR_k_lin * P_noise_floor;
        amplitude = sqrt(P_signal_k / P_signal_unscaled);
        
        target_echo_base = complex(zeros(1, num_all_prt));
        if (delay_samples > 0) && (delay_samples < num_all_prt)
            len_echo = min(length(tx_pulse), num_all_prt - delay_samples);
            target_echo_base(delay_samples+1 : delay_samples+len_echo) = tx_pulse(1:len_echo);
        end
        target_echo_base = amplitude * target_echo_base * doppler_phase_shift;
        phase_shifts_rad = deg2rad(calculate_phase_shifts(elevation_angle, config.Array.element_spacing, config.Sig_Config.wavelength));
        channel_phasors = exp(1j * phase_shifts_rad);
        target_echo_multichannel = target_echo_base.' * channel_phasors;
        pulse_echo_all_channels = pulse_echo_all_channels + target_echo_multichannel;
    end
    raw_iq_data(m, :, :) = reshape(pulse_echo_all_channels, [1, num_all_prt, config.Sig_Config.channel_num]);
end
raw_iq_data_reset0 = raw_iq_data;

% (v8 - %% 4.1) 加噪 (核心: 每次都加 *新* 噪声)
raw_iq_data_noise = complex(zeros(size(raw_iq_data_reset0)));
for c = 1:16
    raw_iq_data_temp = squeeze(raw_iq_data_reset0(:,:,c));
    I_noise = randn(size((raw_iq_data_temp)));
    Q_noise = randn(size(raw_iq_data_temp));
    noise = (I_noise + 1j * Q_noise) .* sqrt(P_noise_floor/2);
    raw_iq_data_noise(:,:,c) = raw_iq_data_temp + noise;
end

% --- (v8 - %% 5) 数字波束形成 (DBF) ---
fprintf('  (S5) 正在执行 DBF...\n');
iq_data_13beam = complex(zeros(config.Sig_Config.prtNum, num_all_prt, config.Sig_Config.beam_num));
for i = 1:config.Sig_Config.prtNum
    single_pulse_16ch = squeeze(raw_iq_data_noise(i, :, :));
    single_pulse_13beam = single_pulse_16ch * DBF_coeffs_data_C';
    iq_data_13beam(i, :, :) = single_pulse_13beam;
end

% --- (v8 - %% 6) 脉冲压缩 (已优化为频域 FFT) ---
fprintf('  (S6) 正在执行脉冲压缩 (FFT优化版)...\n');
pc_results_13beam = complex(zeros(config.Sig_Config.prtNum, N_total_gate, config.Sig_Config.beam_num));
for b = 1:config.Sig_Config.beam_num
    beam_data = squeeze(iq_data_13beam(:,:,b));
    
    % 提取段
    segment_narrow = beam_data(:, seg_start_narrow:end);
    segment_medium = beam_data(:, seg_start_medium:end);
    segment_long = beam_data(:, seg_start_long:end);
    
    % 窄脉冲 (Filter)
    pc_out_narrow_full = filter(MF_narrow, 1, segment_narrow, [], 2);
    pc_out_narrow_full_shift = circshift(pc_out_narrow_full, -fir_delay, 2);
    
    % 中脉冲 (FFT)
    segment_medium_fft = fft(segment_medium, N_fft_med, 2);
    pc_out_medium_full = ifft(segment_medium_fft .* MF_medium_fft, N_fft_med, 2);
    
    % 长脉冲 (FFT)
    segment_long_fft = fft(segment_long, N_fft_long, 2);
    pc_out_long_full = ifft(segment_long_fft .* MF_long_fft, N_fft_long, 2);
    
    % 拼接
    piece1 = pc_out_narrow_full_shift(:, 1: N_gate_narrow);
    piece2 = pc_out_medium_full(:, N_gate_narrow + 1 : N_gate_narrow + N_gate_medium);
    piece3 = pc_out_long_full(:, N_gate_narrow + N_gate_medium + 1 : N_total_gate);
    pc_results_13beam(:,:,b) = [piece1, piece2, piece3];
end

% --- (v8 - %% 7) MTD 处理 ---
fprintf('  (S7) 正在执行 MTD...\n');
rdm_13beam = complex(zeros(size(pc_results_13beam)));
for b = 1:config.Sig_Config.beam_num
    pc_data_single_beam = squeeze(pc_results_13beam(:,:,b));
    pc_data_single_beam_win = pc_data_single_beam .* MTD_win;
    rdm_13beam(:,:,b) = fftshift(fft(pc_data_single_beam_win, [], 1), 1);
end

% --- (v8 - %% 8) 多波束CFAR检测 ---
fprintf('  (S8) 正在执行 2D-GOCA-CFAR...\n');
[all_raw_detections, rdm_for_cfar_all] = fun_run_goca_cfar_8(rdm_13beam, cfar_params, config);
num_total_raw_detections = size(all_raw_detections, 1);

% --- (v8 - %% 9) 逐点参数估计 (样条插值) ---
fprintf('  (S9) 正在执行参数精测...\n');
parameterized_detections = fun_parameter_estimation_9(all_raw_detections, rdm_for_cfar_all, rdm_13beam, ...
    range_axis, velocity_axis, deltaR, deltaV, beam_angles_deg, k_slopes_LUT);

% --- (v8 - %% 10) 第一级: 波束【内】聚类 (R/V/A, 加权平均) ---
fprintf('  (S10) 正在执行波束内聚类...\n');
intra_beam_targets = fun_cluster_stage1_10(parameterized_detections, cluster_params);

% --- (v8 - %% 11) 第二级: 波束【间】聚类 (R/V, 赢家通吃) ---
fprintf('  (S11) 正在执行波束间聚类...\n');
final_targets = fun_cluster_stage2_11(intra_beam_targets, cluster_params);

fprintf('  > 第 %d 帧处理完毕，发现 %d 个唯一目标。\n', frame_idx, length(final_targets));

end % (主函数结束)

%% 12. 本地子函数 (从 v8 复制)
% =========================================================================
% (S4 - 本地子函数)
function channel_phase_shifts_deg = calculate_phase_shifts(arrival_angle_deg, element_spacing_m, wavelength_m)
    arrival_angle_rad = deg2rad(arrival_angle_deg);
    delta_phi_rad = (2 * pi * element_spacing_m * sin(arrival_angle_rad)) / wavelength_m;
    channel_indices = 0:15;
    channel_phase_shifts_rad = channel_indices * delta_phi_rad;
    channel_phase_shifts_deg = rad2deg(channel_phase_shifts_rad);
end

% (S8 - v8 %% 8 的逻辑)
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

% (S9 - v8 %% 9 的逻辑)
function parameterized_detections = fun_parameter_estimation_9(all_raw_detections, rdm_for_cfar_all, rdm_13beam, range_axis, velocity_axis, deltaR, deltaV, beam_angles_deg, k_slopes_LUT)
    
    num_total_raw_detections = size(all_raw_detections, 1);
    if num_total_raw_detections == 0
        parameterized_detections = [];
        return;
    end
    
    [num_V, num_R] = size(squeeze(rdm_13beam(:,:,1)));
    
    % --- 插值参数 (来自 v8 %% 9) ---
    extraDots = 2; rInterpTimes = 8; vInterpTimes = 4;
    
    parameterized_detections = repmat(struct('Range', 0, 'Velocity', 0, 'Angle', 0, 'Power', 0, 'PairIndex', 0), num_total_raw_detections, 1);
    
    for i = 1:num_total_raw_detections
        v_idx = all_raw_detections(i, 1);   % 整数索引
        r_idx = all_raw_detections(i, 2);   % 整数索引
        pair_idx = all_raw_detections(i, 3);
        power = all_raw_detections(i, 4);
        
        rdm_interp = rdm_for_cfar_all(:,:,pair_idx);
        
        % --- 距离维插值 (v8 %% 9 逻辑) ---
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

        % --- 速度维插值 (v8 %% 9 逻辑) ---
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
        % (v8 修正后的速度计算)
        est_velocity = velocity_axis(v_idx) + (vCellMax-v_idx)*deltaV;
        
        % --- 角度估计 (v8 %% 9 逻辑 - *包含缺陷*) ---
        % (缺陷: 使用 v_idx, r_idx 整数索引，而不是 rCellMax, vCellMax)
        S_A = abs(rdm_13beam(v_idx, r_idx, pair_idx));
        S_B = abs(rdm_13beam(v_idx, r_idx, pair_idx + 1));
        
        monopulse_ratio = (S_A - S_B) / (S_A + S_B + eps);
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

% (S10 - v8 %% 10 逻辑)
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

% (S11 - v8 %% 11 逻辑)
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
        [~, idx_winner] = max(powers); % (修复 v8 的一个小 bug, max_power 未使用)
        winner_detection = detections_in_cluster(idx_winner);
        
        final_targets(i).Range = winner_detection.Range;
        final_targets(i).Velocity = winner_detection.Velocity;
        final_targets(i).Angle = winner_detection.Angle;
        final_targets(i).Power = winner_detection.Power; % (修复 v8 的一个小 bug)
    end
end