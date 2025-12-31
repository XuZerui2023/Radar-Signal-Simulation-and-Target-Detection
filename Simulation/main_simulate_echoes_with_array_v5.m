% main_simulate_and_process_v3.m
% 版本: v4.0
% 描述: 这是一个集成了高保真回波仿真、DBF、脉冲压缩（含偏移校正）、
%       MTD和CFAR检测的全流程雷达信号处理仿真与分析脚本。
%       脚本基于对关键偏移量进行了校正，
%       并在每个处理阶段都提供了详细的可视化输出。
%       信号在DBF处理后加噪声
%       加入聚类算法，减少CFAR重复测量点
% 日期: 2025年09月17日

clc; clear; close all;

%% 1. 用户配置区
% =========================================================================
fprintf('--- 1. 开始进行用户参数配置 ---\n');

% --- 1.1 定义仿真的目标 ---
targets = struct(...
    'Range', {10000}, ...         % 目标距离 (m)
    'Velocity', {10}, ...         % 目标速度 (m/s, 正为靠近雷达)
    'RCS', {1}, ...             % 雷达散射截面 (m^2)
    'ElevationAngle', {43.3} ... % 目标俯仰角 (度)
);

SNR_dB = 10; % 信噪比 功率（20dB-100倍，10dB-10倍，3dB-1倍）

prt_selected = 5;   % 选择后续用于画图分析的prt   
beam_selected = 10; % 选择后续用于画图分析的波束

% --- 1.2 文件路径配置 ---
base_path = 'C:\Users\a\Desktop\9.8 问题\Simulation'; 
dbf_coef_path = fullfile(base_path, 'X8数据采集250522_DBFcoef.csv');
% fir_filter_path = fullfile(base_path, 'Simulation', 'FIR.mat'); % FIR系数保存在此

% --- 1.3 CFAR 参数配置 ---
cfar_params.refCells_V = 4;      % 速度维参考单元
cfar_params.guardCells_V = 8;   % 速度维保护单元
cfar_params.refCells_R = 4;      % 距离维参考单元
cfar_params.guardCells_R = 8;    % 距离维保护单元
cfar_params.T_CFAR = 6;          % 检测门限因子
cfar_params.method = 'GOCA';     % 'CA' (平均), 'GOCA' (选大), 'SOCA' (选小)

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
config.Sig_Config.tao = [0.16e-6, 8e-6, 28e-6]; % 脉宽 [窄, 中, 长] (s)，雷达脉冲发射持续时长，此距离为该脉冲盲距
config.Sig_Config.gap_duration = [11.4e-6, 31.8e-6, 153.4e-6];  % 三个脉冲脉宽后的间隔时间（信号为0）
config.Sig_Config.point_prt_segments = [228, 723, 2453]; % 各脉冲段对应的采样波门点数

% --- 2.2 天线阵列参数 ---
config.Array.element_spacing = 0.0138;          % 阵元间距 (m)

% --- 2.3 派生参数 ---
config.Sig_Config.prf = 1/config.Sig_Config.prt;
config.Sig_Config.wavelength = config.Sig_Config.c / config.Sig_Config.fc;
config.Sig_Config.K2 = -config.Sig_Config.B / config.Sig_Config.tao(2); % 中脉冲调频斜率
config.Sig_Config.K3 = config.Sig_Config.B / config.Sig_Config.tao(3);  % 长脉冲调频斜率
ts = 1 / config.Sig_Config.fs;

deltaR = config.Sig_Config.c*ts/2;


%% 3. 生成发射参考波形 (物理真实脉冲)
% =========================================================================
fprintf('--- 3. 生成发射参考波形 ---\n');
% 3.1 各脉冲基本参数
% 窄脉冲参数
tau1 = config.Sig_Config.tao(1);                    % 窄脉冲脉冲宽度
gap_duration1 = config.Sig_Config.gap_duration(1);  % 窄脉冲后脉冲间隔

% 中脉冲参数
tau2 = config.Sig_Config.tao(2);              % 中脉冲脉冲宽度
gap_duration2 = config.Sig_Config.gap_duration(2);  % 中脉冲后脉冲间隔
k2 = -config.Sig_Config.B/tau2;         % 中脉冲调频斜率

% 长脉冲参数
tau3 = config.Sig_Config.tao(3);                % 长脉冲脉冲宽度
gap_duration3 = config.Sig_Config.gap_duration(3);  % 长脉冲后脉冲间隔
k3 = config.Sig_Config.B/tau3;              % 长脉冲调频斜率

% 3.2 生成各脉冲信号
% 计算信号（脉宽内）长度
num_all_prt = round(config.Sig_Config.prt * config.Sig_Config.fs );  % 一个完整PRT信号的总长度
num_samples_1 = round(tau1 * config.Sig_Config.fs);                  % 窄脉冲信号时间长度   
num_samples_2 = round(tau2 * config.Sig_Config.fs);                  % 中脉冲信号时间长度
num_samples_3 = round(tau3 * config.Sig_Config.fs);                  % 长脉冲信号时间长度

num_receiver_1 = round(gap_duration1 * config.Sig_Config.fs);        % 窄脉冲脉冲间隔（窄脉冲发射后能用于接受回波的时间段）
num_receiver_2 = round(gap_duration2 * config.Sig_Config.fs);        % 中脉冲脉冲间隔（中脉冲发射后能用于接受回波的时间段）
num_receiver_3 = round(gap_duration3 * config.Sig_Config.fs);        % 长脉冲脉冲间隔（长脉冲发射后能用于接受回波的时间段）

% 计算信号时间向量
t1 = linspace(-tau1/2, tau1/2, num_samples_1);       
t2 = linspace(-tau2/2, tau2/2, num_samples_2);
t3 = linspace(-tau3/2, tau3/2, num_samples_3);
t_tx = linspace(0, config.Sig_Config.prt, num_all_prt);

% 计算画图用距离轴
R_rx = linspace(24, 24 + 3404*6, 3404);
v_max = config.Sig_Config.wavelength / (2 * config.Sig_Config.prt); % 计算雷达最大不模糊速度
velocity_axis = linspace(-v_max/2, v_max/2, config.Sig_Config.prtNum);
range_axis = (0:config.Sig_Config.point_PRT-1) * (config.Sig_Config.c / (2 * config.Sig_Config.fs));
% range_axis = R_rx-24;


% 生成各脉冲信号
pulse1 = ones(1, num_samples_1);
pulse2 = exp(1j*2*pi*(0.5*k2*(t2.^2)));
pulse3 = exp(1j*2*pi*(0.5*k3*(t3.^2)));

% 3.3 将各子脉冲信号拼接为一个完整的长脉冲信号
tx_pulse = complex(zeros(1, num_all_prt));                              % 构建一个完整脉冲，初始化一个复数矩阵
tx_pulse(1:num_samples_1) = pulse1;                                     % 放置短脉冲

offset1 = round((tau1+gap_duration1) * config.Sig_Config.fs);           % 中脉冲相对于窄脉冲延后的长度
tx_pulse(offset1+1 : offset1+num_samples_2) = pulse2;                   % 放置中脉冲

offset2 = offset1 + round((tau2+gap_duration2) * config.Sig_Config.fs); % 长脉冲相对于中脉冲延后的长度
tx_pulse(offset2+1 : offset2+num_samples_3) = pulse3;                   % 放置长脉冲

% --- 3.4 (可选但建议) 可视化发射信号频谱，检查是否正确 ---
figure('Name', '发射信号频谱检查');
Nfft = 2^nextpow2(num_all_prt);
Tx_fft = fftshift(fft(tx_pulse, Nfft));
freq_axis = linspace(-config.Sig_Config.fs/2, config.Sig_Config.fs/2, Nfft);
plot(freq_axis/1e6, 20*log10(abs(Tx_fft)));
title('发射信号频谱 (应为一个平顶矩形)');
xlabel('频率 (MHz)');
ylabel('幅度 (dB)');
grid on;
ylim([-40, max(20*log10(abs(Tx_fft))) + 5]); % 设置一个好的观察范围

%% 4. 生成16通道目标回波信号
% =========================================================================
fprintf('--- 4. 模拟16通道目标回波 ---\n');
% (此部分代码保持不变，仅修改了噪声和增益)
raw_iq_data = complex(zeros(config.Sig_Config.prtNum, num_all_prt, config.Sig_Config.channel_num));

for m = 1:config.Sig_Config.prtNum % 循环生成每一个prt的信号  
    pulse_echo_all_channels = complex(zeros(num_all_prt, config.Sig_Config.channel_num));
    for k = 1:length(targets)
        range = targets(k).Range;
        % velocity = targets(k).Velocitys;
        % elevation_angle = targets(k).ElevationAngle;
        
        delay = 2 * range / config.Sig_Config.c;
        delay_samples = round(delay / ts); % 计算目标产生的时延点数
        
        doppler_freq = 2 * targets(k).Velocity / config.Sig_Config.wavelength; % 计算目标运动产生的多普勒频率，速度为正时（目标向雷达靠近），fd为正
        doppler_phase_shift = exp(1j * 2 * pi * doppler_freq * (m-1) * config.Sig_Config.prt); % 计算目标运动产生的多普勒相移
        
        % amplitude_gain = 1e11; 
        % lambda_sq = config.Sig_Config.wavelength^2;
        % amplitude = amplitude_gain * sqrt(targets(k).RCS * lambda_sq) / (range^2 * (4*pi)^(3/2));
        
        amplitude = 1;  % 简化幅度
        
        % 1. 先初始化一个全零的回波基带信号
        target_echo_base = complex(zeros(1, num_all_prt));
        % 2. 只有在有效延迟范围内，才将发射波形放入
        if (delay_samples > 0) && (delay_samples < num_all_prt)
            len_echo = min(length(tx_pulse), num_all_prt - delay_samples);
            target_echo_base(delay_samples+1 : delay_samples+len_echo) = tx_pulse(1:len_echo); % 零填充，禁止回卷
            %target_echo_base = circshift(tx_pulse, delay_samples);   % 目标时延，信号循环移位，可能会回卷
        end
        % 3. 对整个基带信号（包含波形和零的部分）乘以幅度和多普勒相移
        target_echo_base = amplitude * target_echo_base * doppler_phase_shift; % 幅度 * 信号 * 多普勒相移
        phase_shifts_rad = deg2rad(calculate_phase_shifts(targets(k).ElevationAngle, config.Array.element_spacing, config.Sig_Config.wavelength)); % 调用子函数计算每个通道的相位偏移
        
        % 计算多通道相位偏移
        channel_phasors = exp(1j * phase_shifts_rad); % 构造阵列流形向量，给每个通道（阵元）叠加一个与入射角相关的相位推进量，来模拟平面波以入射角 θ 打到均匀线阵时，阵元间的相位差。
        target_echo_multichannel = target_echo_base.' * channel_phasors; % 
        pulse_echo_all_channels = pulse_echo_all_channels + target_echo_multichannel;
    end
    raw_iq_data(m, :, :) = reshape(pulse_echo_all_channels, [1, num_all_prt, config.Sig_Config.channel_num]);
end

% 回波信号闭锁期置零
raw_iq_data_reset0 = raw_iq_data; % 闭锁期置零信号

% reset_zero_1_end = num_samples_1;                                  % 窄脉冲信号发射结束点
% reset_zero_2_start = reset_zero_1_end + num_receiver_1 + 1;        % 中脉冲信号发射起始点
% reset_zero_2_end = reset_zero_2_start + num_samples_2 - 1;         % 中脉冲信号发射结束点
% reset_zero_3_start = reset_zero_2_end + num_receiver_2 + 1;        % 长脉冲信号发射起始点
% reset_zero_3_end = reset_zero_3_start + num_samples_3 - 1;         % 长脉冲信号发射结束点
% 
% raw_iq_data_reset0(:, 1:reset_zero_1_end, :) = 0;                  % 窄脉冲发射沿不接收回波信号，置零
% raw_iq_data_reset0(:, reset_zero_2_start:reset_zero_2_end, :) = 0; % 中脉冲发射沿不接收回波信号，置零
% raw_iq_data_reset0(:, reset_zero_3_start:reset_zero_3_end, :) = 0; % 长脉冲发射沿不接收回波信号，置零

% 额外：打印一些自检信息
fprintf('完整PRT点数: %d\n', num_all_prt);
fprintf('预期6 km单程时延≈ %.2f us，对应样点≈ %d 点\n', 2*targets(1).Range/config.Sig_Config.c*1e6, ...
        round((2*targets(1).Range/config.Sig_Config.c)/ts));

% 打印相位步进用于对照
lambda = config.Sig_Config.wavelength;
dphi_ch = 2*pi*config.Array.element_spacing * sind(targets(1).ElevationAngle) / lambda;
fd = -2*targets(1).Velocity / lambda;
dphi_prt = 2*pi*fd*config.Sig_Config.prt;
fprintf('阵元间相位差≈ %.1f°，相邻PRT多普勒步进≈ %.1f°\n', rad2deg(dphi_ch), rad2deg(dphi_prt));

% % DBF前加入高斯白噪声
% I_noise = randn(size(raw_iq_data_reset0));                    % I通道高斯噪声
% Q_noise = randn(size(raw_iq_data_reset0));                    % Q通道高斯噪声
% noise = I_noise + 1j * Q_noise;                               % 合成IQ噪声
% gauss_noise_power = mean(abs(noise).^2,'all');                % 这是多脉冲多波束信号矩阵，求功率应该对每个元素都要处理，结果是IQ信号的总功率
% 
% signal_power = mean(abs(raw_iq_data_reset0(raw_iq_data_reset0 ~= 0)).^2); % 计算信号平均功率
% noise_power = signal_power / (10^(SNR_dB / 10));          % 根据信号平均功率和信噪比计算高斯白噪声功率
% noise_power = noise_power / (gauss_noise_power/2);        % 调整噪声幅度，高斯白噪声功率（方差）大致为1
% noise = noise .* sqrt(noise_power/2);                     % 双路IQ信号功率各占总功率一半，平方根是把功率值转换为幅度值
% raw_iq_data_noise = raw_iq_data_reset0 + noise;           % 在13波束数据上添加既定功率的高斯白噪声


% DBF前对16通道信号数据加入高斯白噪声，循环添加，保证每个波束间的噪声都是相互独立的
raw_iq_data_noise = complex(zeros(332,5819,16));
for c = 1:16
    raw_iq_data_temp = squeeze(raw_iq_data_reset0(:,:,c));
    I_noise = randn(size((raw_iq_data_temp)));                % I通道高斯噪声
    Q_noise = randn(size(raw_iq_data_temp));                  % Q通道高斯噪声
    noise = I_noise + 1j * Q_noise;                           % 合成IQ噪声
    gauss_noise_power = mean(abs(noise).^2,'all');            % 这是多脉冲多波束信号矩阵，求功率应该对每个元素都要处理，结果是IQ信号的总功率

    signal_power = mean(abs(raw_iq_data_temp(raw_iq_data_temp ~= 0)).^2); % 计算信号平均功率
    noise_power = signal_power / (10^(SNR_dB / 10));          % 根据信号平均功率和信噪比计算高斯白噪声功率
    noise_power = noise_power / (gauss_noise_power/2);        % 调整噪声幅度，高斯白噪声功率（方差）大致为1
    noise = noise .* sqrt(noise_power/2);                     % 双路IQ信号功率各占总功率一半，平方根是把功率值转换为幅度值
    raw_iq_data_noise(:,:,c) = raw_iq_data_temp + noise;      % 在13波束数据上添加既定功率的高斯白噪声

end

%% 5. 数字波束形成 (DBF)
% =========================================================================
fprintf('--- 5. 执行数字波束形成 (DBF) ---\n');
% --- 加载DBF系数 ---
DBF_coeffs_data = readmatrix(dbf_coef_path);
DBF_coeffs_data_C = DBF_coeffs_data(:, 1:2:end) + 1j * DBF_coeffs_data(:, 2:2:end);
% --- 逐脉冲DBF ---
iq_data_13beam = complex(zeros(config.Sig_Config.prtNum, num_all_prt, config.Sig_Config.beam_num));
for i = 1:config.Sig_Config.prtNum
    single_pulse_16ch = squeeze(raw_iq_data_noise(i, :, :));
    single_pulse_13beam = single_pulse_16ch * DBF_coeffs_data_C';
    iq_data_13beam(i, :, :) = single_pulse_13beam;
end


% --- DBF结果可视化 ---
figure('Name', '阶段1: DBF后波束功率');
beam_to_check = 10; % 选择一个靠近目标的波束来显示
imagesc(db(abs(squeeze(iq_data_13beam(:,:,beam_to_check)))));
title(sprintf('DBF后波束 #%d 的功率 (dB)', beam_to_check));
xlabel('距离单元 (快时间)'); ylabel('脉冲数 (慢时间)'); colorbar;

% DBF添加噪声，避免噪声幅度相位相关性等被改变
% 在闭锁期置零后的回波信号矩阵上添加符合SNR要求的高斯白噪声

% signal_power_before_DBF = mean(abs(raw_iq_data(raw_iq_data ~= 0)).^2,"all");      % 信号DBF处理前平均功率大小（排除信号为0的点）
% signal_power_after_DBF = mean(abs(iq_data_13beam(iq_data_13beam ~= 0)).^2,"all"); % 信号DBF处理后平均功率大小（排除信号为0的点）
%% 6. 脉冲压缩 (修正版)
% =========================================================================
fprintf('--- 6. 执行分段脉冲压缩 ---\n');

% --- 6.1 定义匹配滤波器 ---
% 窄脉冲: 使用提供的FIR系数
fir_coeffs = [794,1403,2143,2672,2591,1711,-58,-2351,-4592,-5855,-5338,-2389,3005,10341,18410,25779,30907,32768,30907,25779,18410,10341,3005,-2389,-5338,-5855,-4592,-2351,-58,1711,2591,2672,2143,1403,794];
% fir_coeffs = fir1(35, 0.1, 'low'); % 使用低通FIR滤波器
fir_coeffs = 6 * fir_coeffs/max(fir_coeffs); % 归一化

MF_narrow = fir_coeffs; 
fir_delay = round(mean(grpdelay(MF_narrow))); % FIR滤波器的固有延迟

% 中脉冲: LFM信号的时间反转共轭
MF_medium = fliplr(conj(pulse2));

% 长脉冲: LFM信号的时间反转共轭
MF_long = fliplr(conj(pulse3));

% --- 6.2 计算各滤波器的能量，用于增益归一化或提升FIR滤波器增益 ---
% 增益与滤波器系数的能量成正比 (sum of squares of taps)
gain_narrow = sqrt(sum(abs(MF_narrow).^2));
gain_medium = sqrt(sum(abs(MF_medium).^2));
gain_long = sqrt(sum(abs(MF_long).^2));

% --- 6.2 定义拼接参数 ---
N_gate_narrow = config.Sig_Config.point_prt_segments(1); % 228
N_gate_medium = config.Sig_Config.point_prt_segments(2); % 723
N_gate_long = config.Sig_Config.point_prt_segments(3);   % 2453
N_total_gate = sum(config.Sig_Config.point_prt_segments);

% 脉压信号分段切割点数设置
cutFir_start = 4; % 窄脉冲FIR滤波后前端舍弃点数
cutSed_start = 0;  % 中脉冲脉冲压缩后前端舍弃点数
cutLong_start = 0; % 长脉冲脉冲压缩后前端舍弃点数

cutFir_end = 8;    % 窄脉冲FIR滤波后后端舍弃点数
cutSec_end = 502;  % 中脉冲脉冲压缩后后端舍弃点数
cutLong_end = 961; % 长脉冲脉冲压缩后后端舍弃点数

% 对应各脉冲采样波门的起始点
gap_duration1_num = gap_duration1 * config.Sig_Config.fs; % 窄脉冲发射后脉冲间隔点数
gap_duration2_num = gap_duration2 * config.Sig_Config.fs; % 中脉冲发射后脉冲间隔点数
seg_start_narrow = 1;  % 窄脉冲发射起始处 第5点
seg_start_medium = num_samples_1 + gap_duration1_num + 1;  % 中脉冲发射起始处 第290点
seg_start_long = num_samples_1 + gap_duration1_num + num_samples_2 + gap_duration2_num  + 1;  % 长脉冲采样波门起始处 第1285点


% --- 6.3 逐波束进行处理 ---
pc_results_13beam = complex(zeros(config.Sig_Config.prtNum, N_total_gate, config.Sig_Config.beam_num));

for b = 1:config.Sig_Config.beam_num
    % 提取当前波束的所有PRT数据
    beam_data = squeeze(iq_data_13beam(:,:,b));

    % --- 1. 提取用于处理的、连续的回波信号段 ---
    % 窄脉冲处理段: 从窄脉冲发射结束开始
    % seg_start_narrow = num_samples_1 + 1;
    % segment_narrow = beam_data(:, seg_start_narrow:end);
    segment_narrow = beam_data;
    % 中脉冲处理段: 从中脉冲发射结束开始
    % seg_start_medium = offset1 + num_samples_2 + 1;
    % segment_medium = beam_data(:, seg_start_medium:end);
    segment_medium = beam_data;
    
    % 长脉冲处理段: 从长脉冲发射结束开始
    % seg_start_long = offset2 + num_samples_3 + 1;
    % segment_long = beam_data(:, seg_start_long:end);
    segment_long = beam_data;
    
    % --- 2. 对每个连续段进行独立的脉冲压缩 ---
    % 窄脉冲滤波 (时域FIR滤波)
    % filter函数会逐行处理矩阵
    pc_out_narrow_full = filter(MF_narrow, 1, segment_narrow, [], 2);
    % 校正FIR引入的延迟
    pc_out_narrow_full = circshift(pc_out_narrow_full, -fir_delay, 2);

    % 中脉冲脉冲压缩 (时域卷积)
    pc_out_medium_full = complex(zeros(size(segment_medium,1), size(segment_medium,2) + num_samples_2 - 1));
    for i = 1:config.Sig_Config.prtNum
        pc_out_medium_full(i,:) = conv(segment_medium(i,:), MF_medium);
    end
    pc_out_medium_full = pc_out_medium_full(:,1:(size(pc_out_medium_full,2)-cutSec_end));

    
    % 长脉冲脉冲压缩 (时域卷积)
    pc_out_long_full = complex(zeros(size(segment_long,1), size(segment_long,2) + num_samples_3 - 1));
    for i = 1:config.Sig_Config.prtNum
        pc_out_long_full(i,:) = conv(segment_long(i,:), MF_long);
    end
    pc_out_long_full = pc_out_long_full(:,1:(size(pc_out_long_full,2)-cutLong_end));

    % --- 3. 从脉压结果中截取并拼接 ---
    % 窄脉冲结果段
    %piece1 = pc_out_narrow_full(:, 1+cutFir_start : N_gate_narrow+cutFir_start);
    piece1 = pc_out_narrow_full(:, 1 : N_gate_narrow);   % 取窄脉冲处理后信号的5-232点（共计228点）
    % 中脉冲结果段
    % 注意：这里我们截取的是中脉冲处理结果中，紧接着窄脉冲之后的部分
    idx2_start = seg_start_medium + N_gate_narrow;
    %idx2_start = 490;
    idx2_end   = idx2_start + N_gate_medium - 1;
    pc_out_medium_full = pc_out_medium_full(:,490:end);
    piece2 = pc_out_medium_full(:, 1:N_gate_medium);           % 取中脉冲处理后信号的5-232点（共计228点）

    % 长脉冲结果段
    idx3_start = seg_start_long + N_gate_medium;
    %idx3_start = 1985;
    idx3_end   = idx3_start + N_gate_long-1;
    pc_out_long_full = pc_out_long_full(:,1985:end);
    piece3 = pc_out_long_full(:, 1:N_gate_long);

    pc_results_13beam(:,:,b) = [piece1, piece2, piece3];
end

% 脉冲分割拼接时距离补偿
rSysErr_narrow = cutFir_start * deltaR; % 窄脉冲距离补偿
rSysErr_medium = rSysErr_narrow + N_gate_narrow*deltaR; % 中脉冲距离补偿
rSysErr_long = rSysErr_narrow + (N_gate_narrow + N_gate_medium)*deltaR; % 长脉冲距离补偿

rMeasureErr_narrow = 0;
rMeasureErr_medium = 0;
rMeasureErr_long = 0;

% 距离轴拼接
rScale_narrow = (0:N_gate_narrow-1)*deltaR + rSysErr_narrow - rMeasureErr_narrow;      % 短脉冲距离轴
rScale_medium = (0:N_gate_medium-1)*deltaR + rSysErr_medium - rMeasureErr_medium;      % 中脉冲距离轴
rScale_long = (0:N_gate_long-1)*deltaR + rSysErr_long - rMeasureErr_long;              % 长脉冲距离轴

rScale = [rScale_narrow, rScale_medium, rScale_long];

fprintf('脉冲压缩完成。拼接后的总距离门数为: %d\n', size(pc_results_13beam, 2));

%% 7. MTD 处理 (对所有波束)
%=========================================================================
fprintf('--- 7. 执行MTD处理 ---\n');

% 对13个波束的脉压结果，逐一进行MTD
rdm_13beam = complex(zeros(size(pc_results_13beam)));

for b = 1:config.Sig_Config.beam_num
    % 提取单个波束的脉压数据立方体
    pc_data_single_beam = squeeze(pc_results_13beam(:,:,b));

    % 沿慢时间维(第1维)做FFT，从而得到目标速度，并进行fftshift将零频移到中心
    rdm_13beam(:,:,b) = fftshift(fft(pc_data_single_beam, [], 1), 1);
    % rdm_13beam(:,:,b) = fft(pc_data_single_beam, [], 1);
end

fprintf('MTD处理完成，已生成13个波束的距离-多普勒数据。\n');

%% 8. CFAR 检测与目标聚类 (函数调用版)
% =========================================================================
fprintf('--- 8. 执行CFAR检测与目标聚类 ---\n');

% --- 8.1 形成用于检测的和波束 ---
beam_idx_A = 9;
beam_idx_B = 10;
rdm_beamA_power = abs(squeeze(rdm_13beam(:,:,beam_idx_A)));
rdm_beamB_power = abs(squeeze(rdm_13beam(:,:,beam_idx_B)));
rdm_for_cfar = rdm_beamA_power + rdm_beamB_power;

fprintf('已将波束 #%d 和 #%d 合成为用于CFAR检测的和波束。\n', beam_idx_A, beam_idx_B);


% --- 8.2 实现 2D-GOCA-CFAR 检测器 ---
% 获取CFAR参数
P_fa = 1e-6; % 虚警率
alpha = cfar_params.T_CFAR;
guard_R = cfar_params.guardCells_R;
guard_V = cfar_params.guardCells_V;
ref_R = cfar_params.refCells_R;
ref_V = cfar_params.guardCells_V;
[num_V, num_R] = size(rdm_for_cfar);

cfar_detections = zeros(num_V, num_R); % 初始化检测结果矩阵
threshold_map = zeros(num_V, num_R);   % 用于可视化的门限图

% 遍历所有待检测单元 (CUT - Cell Under Test)
for r = (ref_R + guard_R + 1) : (num_R - ref_R - guard_R)
    for v = (ref_V + guard_V + 1) : (num_V - ref_V - guard_V)

        cut_power = rdm_for_cfar(v, r);
 
        % --- 距离维 GOCA ---
        % 前窗
        range_win_leading = rdm_for_cfar(v, r - guard_R - ref_R : r - guard_R - 1);
        % 后窗
        range_win_trailing = rdm_for_cfar(v, r + guard_R + 1 : r + guard_R + ref_R);

        noise_R_leading = mean(range_win_leading);
        noise_R_trailing = mean(range_win_trailing);
        noise_R = max(noise_R_leading, noise_R_trailing);

        % --- 速度维 GOCA ---
        % 前窗 (v-方向)
        doppler_win_leading = rdm_for_cfar(v - guard_V - ref_V : v - guard_V - 1, r);
        % 后窗 (v+方向)
        doppler_win_trailing = rdm_for_cfar(v + guard_V + 1 : v + guard_V + ref_V, r);

        noise_V_leading = mean(doppler_win_leading);
        noise_V_trailing = mean(doppler_win_trailing);
        noise_V = max(noise_V_leading, noise_V_trailing);

        % --- 综合噪声估计 (取两者中的较大值以保证虚警率) ---
        noise_estimate = max(noise_R, noise_V);

        % --- 生成门限并进行比较 ---
        threshold = alpha * noise_estimate;
        threshold_map(v, r) = threshold; % 保存门限

        if cut_power > threshold
            cfar_detections(v, r) = 1;
        end
    end
end

[detected_v_idx, detected_r_idx] = find(cfar_detections);
num_detections = length(detected_v_idx);
fprintf('CFAR检测完成，共发现 %d 个目标点迹。\n', num_detections);

%% 11. 目标聚类与参数精估计
% =========================================================================
fprintf('\n--- 11. 执行目标聚类 ---\n');

if num_detections > 0
    % --- 11.1 聚类参数定义 ---
    max_range_sep = 5;  % 距离维上，点被视为同一簇的最大间隔（单元数）
    max_vel_sep = 3;    % 速度维上，点被视为同一簇的最大间隔（单元数）

    % --- 11.2 简单邻域聚类算法 ---
    detections = [detected_v_idx, detected_r_idx];
    cluster_IDs = zeros(num_detections, 1);
    current_cluster_ID = 0;

    for i = 1:num_detections
        if cluster_IDs(i) == 0 % 如果当前点尚未被分配
            current_cluster_ID = current_cluster_ID + 1; % 创建一个新簇
            
            % 使用一个栈来进行深度优先搜索，找到所有相连的点
            points_to_visit = i; 
            while ~isempty(points_to_visit)
                current_point_idx = points_to_visit(1);
                points_to_visit(1) = []; % 出栈
                
                if cluster_IDs(current_point_idx) == 0
                    cluster_IDs(current_point_idx) = current_cluster_ID; % 分配簇ID
                    
                    % 寻找所有邻近且未被分配的点，并将它们入栈
                    for j = 1:num_detections
                        if cluster_IDs(j) == 0
                            dist_v = abs(detections(current_point_idx, 1) - detections(j, 1));
                            dist_r = abs(detections(current_point_idx, 2) - detections(j, 2));
                            if dist_v <= max_vel_sep && dist_r <= max_range_sep
                                points_to_visit(end+1) = j; % 入栈
                            end
                        end
                    end
                end
            end
        end
    end
    
    num_clusters = current_cluster_ID;
    fprintf('聚类完成，将 %d 个检测点合并为 %d 个目标。\n', num_detections, num_clusters);

    % --- 11.3 计算簇中心 (功率加权平均) 与最终报告 ---
    final_targets = [];
    fprintf('\n--- 最终目标列表 (聚类后) ---\n');
    
    for i = 1:num_clusters
        points_in_cluster_mask = (cluster_IDs == i);
        points_in_cluster_indices = detections(points_in_cluster_mask, :);
        
        % 提取簇内所有点的线形功率
        linear_powers = zeros(size(points_in_cluster_indices, 1), 1);
        for p_idx = 1:size(points_in_cluster_indices, 1)
            v_idx = points_in_cluster_indices(p_idx, 1);
            r_idx = points_in_cluster_indices(p_idx, 2);
            linear_powers(p_idx) = rdm_for_cfar(v_idx, r_idx);
        end
        
        total_power = sum(linear_powers);
        
        % 计算功率加权的平均距离和速度索引
        centroid_r_idx = sum(points_in_cluster_indices(:, 2) .* linear_powers) / total_power;
        centroid_v_idx = sum(points_in_cluster_indices(:, 1) .* linear_powers) / total_power;
        
        % 插值得到更精确的距离和速度
        final_range = interp1(1:num_R, range_axis, centroid_r_idx);
        final_velocity = interp1(1:num_V, velocity_axis, centroid_v_idx);
        
        % (此处可以再次调用单脉冲测角等)
        
        fprintf('目标 %d: 距离 ~%.2f m, 速度 ~%.2f m/s\n', ...
            i, final_range, final_velocity);
    end
    
else
    fprintf('未检测到目标，无需聚类。\n');
end







% %% 9. 单脉冲测角
% % =========================================================================
% fprintf('\n--- 10. 执行单脉冲测角 ---\n');
% 
% if num_detections > 0
%     % 遍历每一个CFAR检测到的目标点
%     for i = 1:num_detections
%         % 获取当前目标在RDM中的索引
%         r_idx = detected_r_idx(i);
%         v_idx = detected_v_idx(i);
% 
%         fprintf('\n--- 正在分析目标 #%d (距离: %.2f m, 速度: %.2f m/s) ---\n', ...
%             i, range_axis(r_idx), velocity_axis(v_idx));
% 
%         % --- 1. 寻找该目标信号最强的"主波束" ---
%         all_beams_power_at_target = zeros(1, config.Sig_Config.beam_num);
%         for b = 1:config.Sig_Config.beam_num
%             rdm = squeeze(rdm_13beam(:,:,b));
%             all_beams_power_at_target(b) = abs(rdm(v_idx, r_idx));
%         end
%         [max_power, main_beam_idx] = max(all_beams_power_at_target);
%         fprintf('目标信号在波束 #%d 中最强。\n', main_beam_idx);
% 
%         % --- 2. 选择相邻波束组成和差通道 ---
%         beam_A_idx = 0;
%         beam_B_idx = 0;
% 
%         % 决定使用左邻波束还是右邻波束
%         % (处理边界情况：如果最强波束是#1或#13)
%         if main_beam_idx == 1
%             beam_A_idx = 1; beam_B_idx = 2;
%         elseif main_beam_idx == config.Sig_Config.beam_num
%             beam_A_idx = config.Sig_Config.beam_num - 1;
%             beam_B_idx = config.Sig_Config.beam_num;
%         else
%             % 比较左右邻波束的功率，选择更强的那个组成最优测角波束对
%             power_left = all_beams_power_at_target(main_beam_idx - 1);
%             power_right = all_beams_power_at_target(main_beam_idx + 1);
%             if power_left > power_right
%                 beam_A_idx = main_beam_idx - 1;
%                 beam_B_idx = main_beam_idx;
%             else
%                 beam_A_idx = main_beam_idx;
%                 beam_B_idx = main_beam_idx + 1;
%             end
%         end
%         fprintf('选择波束对 (#%d, #%d) 进行单脉冲测角。\n', beam_A_idx, beam_B_idx);
% 
%         % --- 3. 提取对应波束在目标位置的复数值 ---
%         rdm_A = squeeze(rdm_13beam(:,:,beam_A_idx));
%         rdm_B = squeeze(rdm_13beam(:,:,beam_B_idx));
% 
%         S_A = rdm_A(v_idx, r_idx);
%         S_B = rdm_B(v_idx, r_idx);
% 
%         % --- 4. 计算和、差信号及单脉冲比 ---
%         S_sum = S_A + S_B;
%         S_diff = S_A - S_B; % 差信号
% 
%         monopulse_ratio = S_diff / S_sum;
% 
%         fprintf('和信号(Σ)幅值: %.2f\n', abs(S_sum));
%         fprintf('差信号(Δ)幅值: %.2f\n', abs(S_diff));
%         fprintf('单脉冲比 (Δ/Σ): %.4f + %.4fi\n', real(monopulse_ratio), imag(monopulse_ratio));
%     end
% else
%     fprintf('未检测到目标，无法进行单脉冲测角。\n');
% end



%% 可视化分析
% 1. 发射信号实部时域图
figure;
plot(t_tx*1e6, real(tx_pulse));
title('发射信号实部时域波形图');
xlabel('时间（us）');
ylabel('信号幅度');

% 2. 正交高斯白噪声时域图
figure;
subplot(3, 1, 1);
plot(real(noise(1,:)));
title('正交高斯白噪声时域信号');
subtitle('I通道（实部）');
xlabel('样本点');
ylabel('幅度');

subplot(3, 1, 2);
plot(imag(noise(1,:)));
subtitle('Q通道（虚部）');
xlabel('样本点');
ylabel('幅度');

subplot(3, 1, 3);
plot(abs(noise(1,:)));
subtitle('幅值');
xlabel('样本点');
ylabel('幅度');

% 3. DBF前回波信号时域图
% 加噪前
figure;
subplot(3,1,1);
plot(t_tx*1e6,real(raw_iq_data(prt_selected,:,beam_selected))); % 画prt 1，通道2的波形图
title(sprintf('加噪前回波信号实部时域波形图 PRT %d BEAM %d',prt_selected,beam_selected));
subtitle('实部');
xlabel('时间（us）');
ylabel('信号幅度');

subplot(3,1,2);
plot(t_tx*1e6,imag(raw_iq_data(prt_selected,:,beam_selected))); % 画prt 1，通道2的波形图
subtitle('虚部');
xlabel('时间（us）');
ylabel('信号幅度');

subplot(3,1,3);
plot(t_tx*1e6,abs(raw_iq_data(prt_selected,:,beam_selected))); % 画prt 1，通道2的波形图
subtitle('幅值');
xlabel('时间（us）');
ylabel('信号幅度');

% 4. DBF前回波信号时域图
% 加噪后
figure;
subplot(3,1,1);
plot(t_tx*1e6,real(raw_iq_data_noise(prt_selected,:,beam_selected))); % 画prt 1，通道2的波形图
title(sprintf('加噪后回波信号实部时域波形图 PRT %d BEAM %d',prt_selected,beam_selected));
subtitle('实部');
xlabel('时间（us）');
ylabel('信号幅度');

subplot(3,1,2);
plot(t_tx*1e6,imag(raw_iq_data_noise(prt_selected,:,beam_selected))); % 画prt 1，通道2的波形图
subtitle('虚部');
xlabel('时间（us）');
ylabel('信号幅度');

subplot(3,1,3);
plot(t_tx*1e6,abs(raw_iq_data_noise(prt_selected,:,beam_selected))); % 画prt 1，通道2的波形图
subtitle('幅值');
xlabel('时间（us）');
ylabel('信号幅度');

% 5. DBF后波束数据时域波形图
figure;
subplot(3,1,1);
plot(t_tx*1e6,real(iq_data_13beam(prt_selected,:,beam_selected))); % 画prt 1，通道2的波形图
title('DBF后波束数据实部时域波形图');
subtitle('实部');
xlabel('时间（us）');
ylabel('信号幅度');

subplot(3,1,2);
plot(t_tx*1e6,imag(iq_data_13beam(prt_selected,:,beam_selected))); % 画prt 1，通道2的波形图
subtitle('虚部');
xlabel('时间（us）');
ylabel('信号幅度');

subplot(3,1,3);
plot(t_tx*1e6,abs(iq_data_13beam(prt_selected,:,beam_selected))); % 画prt 1，通道2的波形图
subtitle('幅值');
xlabel('时间（us）');
ylabel('信号幅度');


% 预期时间标注（μs）
t_delay = 2*targets(1).Range/config.Sig_Config.c;
tau1 = config.Sig_Config.tao(1);  tau2 = config.Sig_Config.tao(2);
gap1 = config.Sig_Config.gap_duration(1);  gap2 = config.Sig_Config.gap_duration(2);
t_mid_start  = t_delay + (tau1 + gap1);
t_mid_end    = t_mid_start + tau2;
t_long_start = t_delay + (tau1 + gap1 + tau2 + gap2);
t_long_end   = t_long_start + config.Sig_Config.tao(3);

figure; 
plot(t_tx*1e6, abs(squeeze(iq_data_13beam(prt_selected,:,beam_selected)))); 
hold on; 
grid on;
xline([t_delay t_mid_start t_mid_end t_long_start t_long_end]*1e6, '--');
title(sprintf('DBF后波束数据实部时域波形图  PRT %d  BEAM %d', prt_selected, beam_selected)), xlabel('时间 (\mus)'), ylabel('Real');


% 6. 画脉冲压缩后距离效果图
figure;
subplot(4,1,1);
plot(R_rx-24, real(pc_results_13beam(prt_selected,:,beam_selected))); 
title('脉压后拼接结果图');
subtitle('实部');
xlabel('距离（m）');
ylabel('信号幅度');
for i=1:length(targets)
    xline(targets(i).Range,'--r','目标','LabelOrientation','horizontal') 
end

subplot(4,1,2);
plot(R_rx-24, imag(pc_results_13beam(prt_selected,:,beam_selected))); 
title('虚部');
xlabel('距离（m）');
ylabel('信号幅度');
for i=1:length(targets)
    xline(targets(i).Range,'--r','目标','LabelOrientation','horizontal') 
end

subplot(4,1,3);
plot(rScale, abs(pc_results_13beam(prt_selected,:,beam_selected))); 
title('幅值');
xlabel('距离（m）');
ylabel('信号幅度');
for i=1:length(targets)
    xline(targets(i).Range,'--r','目标','LabelOrientation','horizontal') 
end

subplot(4,1,4);
plot(rScale, 20*log10(abs(pc_results_13beam(prt_selected,:,beam_selected)))); 
title('幅值（dB）');
xlabel('距离（m）');
ylabel('信号幅度（dB）');
for i=1:length(targets)
    xline(targets(i).Range,'--r','目标','LabelOrientation','horizontal') 
end

% 7. MTD结果可视化
% =========================================================================
fprintf('--- 生成MTD可视化图表 ---\n');

% --- 定义用于绘图的坐标轴 ---
v_max = config.Sig_Config.wavelength / (2 * config.Sig_Config.prt); % 计算雷达最大不模糊速度
velocity_axis = linspace(-v_max/2, v_max/2, config.Sig_Config.prtNum);
% range_axis = (0:config.Sig_Config.point_PRT-1) * (config.Sig_Config.c / (2 * config.Sig_Config.fs));
range_axis = R_rx-24;
% 选择一个波束用于展示，根据仿真，选择其中一个对准目标的波束来展示
beam_to_display = beam_selected;
rdm_single_beam = squeeze(rdm_13beam(:,:,beam_to_display));

% 图 7.1：距离-多普勒二维图 (RDM)
figure('Name', 'MTD结果: 2D距离多普勒图');
imagesc(range_axis, velocity_axis, abs(rdm_single_beam));
colormap('jet');
colorbar;
axis xy; % 将y轴原点置于左下角，符合常规视图
xlabel('距离 (m)');
ylabel('速度 (m/s)');
title(sprintf('MTD结果：距离-多普勒二维图 (RDM) - 波束 #%d', beam_to_display));
grid on;


% 图 7.2：绘制速度单元剖面图
% 选择一个特定的速度单元（例如第 10 个速度单元）
[~, velocity_unit] = min(abs(velocity_axis-targets.Velocity));     % 找到符合设定的目标速度的速度单元
rdm_single_velocity = squeeze(rdm_single_beam(velocity_unit, :));  % 提取该速度单元下所有距离的值

% 绘制速度维剖面图
figure;
plot(range_axis, 20*log10(abs(rdm_single_velocity)));  % 以距离作为X轴，绘制该速度单元下的距离剖面
xlabel('距离 (m)');
ylabel('幅度');
title(sprintf('波束 #%d 速度单元剖面图 (速度: %.2f m/s)', beam_to_display, velocity_axis(velocity_unit)));
grid on;


% 图 7.3：绘制距离单元剖面图
% 选择一个特定的距离单元来绘制剖面图
distance_unit = round(targets.Range/6);  % 选择距离单元，假设为第100个距离点
range_profile = 20*log10(abs(rdm_single_beam(:, distance_unit)));  % 提取对应距离的多普勒信息

% 绘制距离单元剖面图
figure;
plot(velocity_axis, range_profile);  % 将幅度转换为dB，绘制速度和能量的关系
xlabel('速度 (m/s)');
ylabel('幅度 (dB)');
title(sprintf('波束 #%d 距离单元剖面图 (距离: %.2f m)', beam_to_display, targets.Range));
grid on;

% 图 7.4：MTD结果三维图
figure('Name', 'MTD结果: 3D距离多普勒图');
% 为了提高绘图效率和美观度，可以对数据进行适当的降采样
[R_grid, V_grid] = meshgrid(range_axis, velocity_axis);
skip_rate = 4; % 每隔4个点绘制一个，避免图像过于密集
surf(R_grid(1:skip_rate:end, 1:skip_rate:end), ...
     V_grid(1:skip_rate:end, 1:skip_rate:end), ...
     20*log10(abs(rdm_single_beam(1:skip_rate:end, 1:skip_rate:end))), ...
     'EdgeColor', 'none'); % 使用'none'边缘颜色更平滑
shading interp;
colormap('jet');
colorbar;
xlabel('距离 (m)');
ylabel('速度 (m/s)');
zlabel('幅度 (dB)');
title(sprintf('MTD结果三维图 - 波束 #%d', beam_to_display));
view(45, 30); % 设置一个较好的初始观察视角

% 图 7.5：MTD杂波抑制效果对比图
figure('Name', 'MTD结果: 杂波抑制对比');
% MTD处理前 (取第一个PRT的脉压结果)
range_profile_before_mtd = abs(squeeze(pc_results_13beam(1, :, beam_to_display)));

% MTD处理后 (对所有非零多普勒单元的能量进行累加)
rdm_power = abs(rdm_single_beam).^2;
% 将零多普勒通道置零 (或用一个极小值代替以避免log(0)错误)
rdm_power(floor(config.Sig_Config.prtNum/2)+1, :) = 1e-10; 
range_profile_after_mtd = sum(rdm_power, 1);

% 绘图 (转换为dB)
plot(range_axis, 20*log10(range_profile_before_mtd), 'b-', 'LineWidth', 1);
hold on;
plot(range_axis, 10*log10(range_profile_after_mtd), 'r-', 'LineWidth', 2);
grid on;
xlabel('距离 (m)');
ylabel('幅度 (dB)');
title(sprintf('MTD杂波抑制效果对比 - 波束 #%d', beam_to_display));
legend('MTD处理前 (单PRT)', 'MTD处理后 (积累动目标)');
% 动态调整Y轴范围以便观察
max_val = max(20*log10(range_profile_before_mtd));
if isfinite(max_val)
    ylim([max_val - 80, max_val + 10]); 
end




%% 单脉冲全程滤波（脉压）结果
beam_data_temp = squeeze(iq_data_13beam(:,:,beam_selected));

% 窄脉冲 (时域FIR滤波)
% 加窗
window = hanning(length(MF_narrow))';
fir_coeffs_windowed = MF_narrow.*window;
% filter函数会逐行处理矩阵
pc_out_narrow_full_temp = filter(fir_coeffs_windowed, 1, beam_data_temp, [], 2);
% 校正FIR引入的延迟
pc_out_narrow_full_temp = circshift(pc_out_narrow_full_temp, -fir_delay, 2);

% 中脉冲 (时域卷积)
pc_out_medium_full_temp = complex(zeros(size(beam_data_temp,1), size(beam_data_temp,2) + num_samples_2 - 1));
for i = 1:config.Sig_Config.prtNum
    pc_out_medium_full_temp(i,:) = conv(beam_data_temp(i,:), MF_medium);
end

% 长脉冲 (时域卷积)
pc_out_long_full_temp = complex(zeros(size(beam_data_temp,1), size(beam_data_temp,2) + num_samples_3 - 1));
for i = 1:config.Sig_Config.prtNum
    pc_out_long_full_temp(i,:) = conv(beam_data_temp(i,:), MF_long);
end

% 窄脉冲prt上滤波结果
RX1 = linspace(0, length(pc_out_narrow_full_temp)*6, length(pc_out_narrow_full_temp));
figure;
subplot(4,1,1)
plot(RX1, real(pc_out_narrow_full_temp(1,:))); % 输出脉冲压缩结果
title('窄脉冲全程prt上滤波信号幅度');
subtitle('实部');
xlabel('距离 (m)');
ylabel('幅度');

subplot(4,1,2)
plot(RX1, imag(pc_out_narrow_full_temp(1,:))); % 输出脉冲压缩结果
subtitle('虚部');
xlabel('距离 (m)');
ylabel('幅度');

subplot(4,1,3)
plot(RX1, abs(pc_out_narrow_full_temp(1,:))); % 输出脉冲压缩结果
subtitle('幅值');
xlabel('距离 (m)');
ylabel('幅度');

subplot(4,1,4)
plot(RX1, 20*log10(abs(pc_out_narrow_full_temp(1,:)))); % 输出脉冲压缩结果
subtitle('幅值（dB）');
xlabel('距离 (m)');
ylabel('幅度');

% 中脉冲prt上脉压结果
%pc_out_medium_full_temp = cicshift()
RX2 = linspace(0, length(pc_out_medium_full_temp)*6, length(pc_out_medium_full_temp));
figure;
subplot(4,1,1);
plot(RX2, real(pc_out_medium_full_temp(1,:))); % 输出脉冲压缩结果
title('中脉冲全程prt上脉压后信号幅度');
subtitle('实部');
xlabel('距离 (m)');
ylabel('幅度');

subplot(4,1,2);
plot(RX2, imag(pc_out_medium_full_temp(1,:))); % 输出脉冲压缩结果
subtitle('虚部');
xlabel('距离 (m)');
ylabel('幅度');

subplot(4,1,3);
plot(RX2, abs(pc_out_medium_full_temp(1,:))); % 输出脉冲压缩结果
subtitle('幅值');
xlabel('距离 (m)');
ylabel('幅度');

subplot(4,1,4);
plot(RX2, 20*log10(abs(pc_out_medium_full_temp(1,:)))); % 输出脉冲压缩结果
subtitle('幅值（dB）');
xlabel('距离 (m)');
ylabel('幅度');


% 长脉冲prt上脉压结果
RX3 = linspace(0, length(pc_out_long_full_temp)*6, length(pc_out_long_full_temp));
figure;
subplot(4,1,1);
plot(RX3, real(pc_out_long_full_temp(1,:))); % 输出脉冲压缩结果
title('长脉冲全程prt脉压后信号幅度');
subtitle('实部');
xlabel('距离 (m)');
ylabel('幅度');

subplot(4,1,2);
plot(RX3, imag(pc_out_long_full_temp(1,:))); % 输出脉冲压缩结果
subtitle('虚部');
xlabel('距离 (m)');
ylabel('幅度');

subplot(4,1,3);
plot(RX3, abs(pc_out_long_full_temp(1,:))); % 输出脉冲压缩结果
subtitle('幅值');
xlabel('距离 (m)');
ylabel('幅度');


subplot(4,1,4);
plot(RX3, 20*log10(abs(pc_out_long_full_temp(1,:)))); % 输出脉冲压缩结果
subtitle('幅值（dB）');
xlabel('距离 (m)');
ylabel('幅度');

% % 噪声FFT
% NOISE = fft(noise(1,:,1));
% figure;
% subplot(3,1,1)
% plot(real(noise(1,:,1)));
% title('原始噪声信号幅度');
% subplot(3,1,2)
% plot(imag(noise(1,:,1)));
% subplot(3,1,3)
% plot(20*log10(abs(noise(1,:,1))));
% 
% 
% figure;
% subplot(3,1,1)
% plot(real(NOISE));
% title('原始噪声FFT后后信号幅度');
% subplot(3,1,2)
% plot(imag(NOISE));
% subplot(3,1,3)
% plot(20*log10(abs(NOISE)));
% 
% % 窄脉冲 (时域FIR滤波)
% % filter函数会逐行处理矩阵
% noise_filter = filter(fir_coeffs, 1, noise(1,:,1));
% 
% figure;
% subplot(3,1,1)
% plot(real(noise_filter));
% title('原始噪声滤波后后信号幅度');
% subplot(3,1,2)
% plot(imag(noise_filter));
% subplot(3,1,3)
% plot(20*log10(abs(noise_filter)));




%% 9. 最终结果可视化
% =========================================================================
fprintf('--- 9. 生成最终可视化图表 ---\n');

% --- 9.1 和波束的距离-多普勒图 ---
figure('Name', '最终检测结果');
subplot(2,1,1);
v_max = config.Sig_Config.wavelength / (2 * config.Sig_Config.prt);
velocity_axis = linspace(-v_max/2, v_max/2, config.Sig_Config.prtNum);
range_axis_cfar = (0:config.Sig_Config.point_PRT-1) * (config.Sig_Config.c / (2 * config.Sig_Config.fs));
imagesc(range_axis_cfar, velocity_axis, rdm_for_cfar);
colormap('jet'); colorbar; axis xy;
title(sprintf('用于CFAR的和波束 (波束 #%d + #%d) RDM', beam_idx_A, beam_idx_B));
xlabel('距离 (m)');
ylabel('速度 (m/s)');
% 
% --- 9.2 CFAR检测结果叠加显示 ---
subplot(2,1,2);
imagesc(range_axis_cfar, velocity_axis, rdm_for_cfar);
colormap('jet'); colorbar; axis xy;
hold on;

% 在图上标记检测到的目标
if num_detections > 0
    detected_ranges = range_axis_cfar(detected_r_idx);
    detected_velocities = velocity_axis(detected_v_idx);
    plot(detected_ranges, detected_velocities, 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'CFAR 检测点');
end
title('CFAR检测结果叠加图');
xlabel('距离 (m)');
ylabel('速度 (m/s)');
legend;

% --- 额外：打印检测到的目标信息 ---
if num_detections > 0
    fprintf('\n--- 检测到的目标列表 ---\n');
    for i = 1:num_detections
        fprintf('目标 %d: 距离 ~%.2f m, 速度 ~%.2f m/s\n', ...
            i, detected_ranges(i), detected_velocities(i));
    end
end

%% 本地子函数
% =========================================================================
function channel_phase_shifts_deg = calculate_phase_shifts(arrival_angle_deg, element_spacing_m, wavelength_m)
    arrival_angle_rad = deg2rad(arrival_angle_deg);
    delta_phi_rad = (2 * pi * element_spacing_m * sin(arrival_angle_rad)) / wavelength_m;
    channel_indices = 0:15;
    channel_phase_shifts_rad = channel_indices * delta_phi_rad;
    channel_phase_shifts_deg = rad2deg(channel_phase_shifts_rad);
end




