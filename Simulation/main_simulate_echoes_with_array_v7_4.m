% main_simulate_echoes_with_array_v7_4.m
% 版本: v7.4
% 描述: 这是一个集成了高保真回波仿真、DBF、脉冲压缩（含偏移校正）、
%       MTD、CFAR检测、波束归并聚类和单脉冲测角参数测量的全流程雷达信号处理仿真与分析脚本。
%       脚本基于对关键偏移量进行了校正，
%       并在每个处理阶段都提供了详细的可视化输出。
%       信号在DBF处理前加噪声，每个波束的噪声要单独加，这样可以避免波束间产生相关性
%       加入2D-GOCA-CFAR 十字形检测器，进行多波束目标检测，效果较好
%       加入波束归并聚类算法，减少CFAR重复测量点
%       加入单脉冲测角算法，实现测角功能和参数测量
% 日期: 2025年10月20日

%  修改记录
%  date       by      version   modify
%  25/08/24   XZR      v1.0     创建0 
%  25/08/24   XZR      v2.0     创建，增加了对16通道线性天线阵列的模拟
%  25/08/26   XZR      V2.1     在仿真结束后，脚本会自动对生成的数据进行脉冲压缩处理，并绘制距离-幅度图。
%  25/08/27   XZR      V2.2     增添多个绘图指令用于分析模拟生成的数据。
%  25/09/04   XZR      V2.3     统一绘图程序，将脉压
%  25/10/12   XZR      V7.1
%  25/10/17   XZR      V7.2     先对初步检测点参数测量，再聚类归并
%  25/10/20   XZR      V7.3     对脉压和MTD过程加入kaiser窗函数
%  25/10/21   XZR      V7.4     对所有原始目标检测点进行三次样条插值参数精测

clc; clear; close all;

%% 1. 用户配置区
% =========================================================================
fprintf('--- 1. 开始进行用户参数配置 ---\n');

% --- 1.1 定义仿真的目标 ---
targets = struct(...
    'Range', {10000}, ...        % 目标距离 (m)
    'Velocity', {25}, ...       % 目标速度 (m/s, 正为靠近雷达)
    'RCS', {1}, ...             % 雷达散射截面 (m^2)
    'ElevationAngle', {35} ...  % 目标俯仰角 (度)
);

SNR_dB = -10; % 信噪比 功率（20dB-100倍，10dB-10倍，3dB-1倍）

prt_selected = 5;   % 选择后续用于画图分析的prt   
beam_selected = 9; % 选择后续用于画图分析的波束

% --- 1.2 文件路径配置 ---
base_path = 'C:\Users\a\Desktop\9.8 问题\Simulation'; 
dbf_coef_path = fullfile(base_path, 'X8数据采集250522_DBFcoef.csv');
% fir_filter_path = fullfile(base_path, 'Simulation', 'FIR.mat'); % FIR系数保存在此

% --- 1.3 CFAR 参数配置 ---
cfar_params.refCells_V = 5;      % 速度维参考单元
cfar_params.guardCells_V = 10;   % 速度维保护单元
cfar_params.refCells_R = 5;      % 距离维参考单元
cfar_params.guardCells_R = 10;    % 距离维保护单元
cfar_params.T_CFAR = 8;          % 检测门限因子
cfar_params.method = 'GOCA';     % 'CA' (平均), 'GOCA' (选大), 'SOCA' (选小)

% --- 1.4：聚类参数配置 ---
cluster_params.max_range_sep = 50;  % 距离维上，点被视为同一簇的最大间隔 (米)
cluster_params.max_vel_sep = 1.0;   % 速度维上，点被视为同一簇的最大间隔 (米/秒)
cluster_params.max_angle_sep = 5.0;  % 角度维上，点被视为同一簇的最大间隔 (度)

% % --- 1.5：脉压窗函数选择 ---
% PC_win_TYPE = 1; % 快时间维加窗：1-hamming；2-hanning 3- kaiser; 4-blackman
% switch PC_win_TYPE
%     case 1 % hamming
%         pcWh = hamming(prtNum);
%     case 2 % hanning
%         pcWh = hann(prtNum);
%     case 3 % kaiser
%         betaPC = 4.5;
%         pcWh = kaiser(prtNum,betaPC);
%     case 4 % blackman
%         pcWh = blackman(prtNum);
% end
% 
% % --- 1.5：MTD窗函数选择 ---
% MTD_win_TYPE = 1; % 慢时间维加窗：1-hamming；2-hanning 3- kaiser; 4-blackman
% switch MTD_win_TYPE
%     case 1 % hamming
%         mtdWh = hamming(prtNum);
%     case 2 % hanning
%         mtdWh = hann(prtNum);
%     case 3 % kaiser
%         betaMTD = 4.5;
%         mtdWh = kaiser(prtNum,betaMTD);
%     case 4 % blackman
%         mtdWh = blackman(prtNum);
% end

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
config.Sig_Config.K2 = -config.Sig_Config.B / config.Sig_Config.tao(2);    % 中脉冲调频斜率
config.Sig_Config.K3 = config.Sig_Config.B / config.Sig_Config.tao(3);     % 长脉冲调频斜率
ts = 1 / config.Sig_Config.fs;
lamda = config.Sig_Config.c/config.Sig_Config.fc; % 波长 
deltaDoppler = config.Sig_Config.prf/config.Sig_Config.prtNum;
deltaR = config.Sig_Config.c*ts/2; % 距离单元分辨率
deltaV = lamda*deltaDoppler/2; 

%% 3. 生成发射参考波形 (物理真实脉冲)
% =========================================================================
fprintf('--- 3. 生成发射参考波形 ---\n');
% 3.1 各脉冲基本参数
% 窄脉冲参数
tau1 = config.Sig_Config.tao(1);                    % 窄脉冲脉冲宽度
gap_duration1 = config.Sig_Config.gap_duration(1);  % 窄脉冲后脉冲间隔

% 中脉冲参数
tau2 = config.Sig_Config.tao(2);                    % 中脉冲脉冲宽度
gap_duration2 = config.Sig_Config.gap_duration(2);  % 中脉冲后脉冲间隔
k2 = -config.Sig_Config.B/tau2;                     % 中脉冲调频斜率

% 长脉冲参数
tau3 = config.Sig_Config.tao(3);                    % 长脉冲脉冲宽度
gap_duration3 = config.Sig_Config.gap_duration(3);  % 长脉冲后脉冲间隔
k3 = config.Sig_Config.B/tau3;                      % 长脉冲调频斜率

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


% DBF前对16通道信号数据加入高斯白噪声
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
    raw_iq_data_noise(:,:,c) = raw_iq_data_temp + noise;           % 在13波束数据上添加既定功率的高斯白噪声

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

% --- 6.1 定义匹配滤波器并加窗（可选） ---
% 脉压前加窗抑制旁瓣减少伪目标，加窗的本质为对匹配滤波器的脉冲响应函数加权
% 窄脉冲时域匹配滤波器: 使用提供的FIR系数
fir_coeffs = [794,1403,2143,2672,2591,1711,-58,-2351,-4592,-5855,-5338,-2389,3005,10341,18410,25779,30907,32768,30907,25779,18410,10341,3005,-2389,-5338,-5855,-4592,-2351,-58,1711,2591,2672,2143,1403,794];
% fir_coeffs = fir1(35, 0.1, 'low'); % 使用低通FIR滤波器
fir_coeffs = 6 * fir_coeffs/max(fir_coeffs); % 归一化

MF_narrow = fir_coeffs;
%MF_narrow_win = MF_narrow .* hanning(length(fir_coeffs))';
fir_delay = round(mean(grpdelay(MF_narrow))); % FIR滤波器的固有延迟

% 中脉冲时域匹配滤波器: LFM信号的时间反转共轭
MF_medium = fliplr(conj(pulse2));
win_medium = kaiser(length(pulse2), 4.5); % 使用kaiser窗控制旁瓣抑制。beta值越大，旁瓣越低，但主瓣越宽。常用值在4到8之间。
MF_medium_win = fliplr(conj(pulse2 .* win_medium'));

% 长脉冲时域匹配滤波器: LFM信号的时间反转共轭
MF_long = fliplr(conj(pulse3));
win_long = kaiser(length(pulse3), 4.5);
MF_long_win = fliplr(conj(pulse3 .* win_long'));

% --- 6.2 计算各滤波器的能量，用于增益归一化或提升FIR滤波器增益 ---
% 增益与滤波器系数的能量成正比 (sum of squares of taps)
gain_narrow = sqrt(sum(abs(MF_narrow).^2));
gain_medium = sqrt(sum(abs(MF_medium_win).^2));
gain_long = sqrt(sum(abs(MF_long_win).^2));

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
seg_start_narrow = num_samples_1 + 1;  % 窄脉冲采样波门起始处 第5点
seg_start_medium = num_samples_1 + gap_duration1_num + num_samples_2 + 1;  % 中脉冲采样波门起始处 第490点
seg_start_long = num_samples_1 + gap_duration1_num + num_samples_2 + gap_duration2_num + num_samples_3 + 1;  % 长脉冲采样波门起始处 第1985点


% --- 6.3 逐波束进行处理 ---
pc_results_13beam = complex(zeros(config.Sig_Config.prtNum, N_total_gate, config.Sig_Config.beam_num));

for b = 1:config.Sig_Config.beam_num
    % 提取当前波束的所有PRT数据
    beam_data = squeeze(iq_data_13beam(:,:,b));

    % --- 1. 提取用于处理的、连续的回波信号段 ---
    % 窄脉冲处理段: 从窄脉冲发射结束开始
    % seg_start_narrow = num_samples_1 + 1;
    segment_narrow = beam_data(:, seg_start_narrow:end);
    
    % 中脉冲处理段: 从中脉冲发射结束开始
    % seg_start_medium = offset1 + num_samples_2 + 1;
    segment_medium = beam_data(:, seg_start_medium:end);
    
    % 长脉冲处理段: 从长脉冲发射结束开始
    % seg_start_long = offset2 + num_samples_3 + 1;
    segment_long = beam_data(:, seg_start_long:end);
    
    % --- 2. 对每个连续段进行独立的脉冲压缩 ---
    % 窄脉冲滤波 (时域FIR滤波)
    % filter函数会逐行处理矩阵
    pc_out_narrow_full = filter(MF_narrow, 1, segment_narrow, [], 2);
    % 校正FIR引入的延迟
    pc_out_narrow_full = circshift(pc_out_narrow_full, -fir_delay, 2);

    % 中脉冲脉冲压缩 (时域卷积)
    pc_out_medium_full = complex(zeros(size(segment_medium,1), size(segment_medium,2) + num_samples_2 - 1));
    for i = 1:config.Sig_Config.prtNum
        pc_out_medium_full(i,:) = conv(segment_medium(i,:), MF_medium_win);
    end

    % 长脉冲脉冲压缩 (时域卷积)
    pc_out_long_full = complex(zeros(size(segment_long,1), size(segment_long,2) + num_samples_3 - 1));
    for i = 1:config.Sig_Config.prtNum
        pc_out_long_full(i,:) = conv(segment_long(i,:), MF_long_win);
    end

    % --- 3. 从脉压结果中截取并拼接 ---
    % 窄脉冲结果段
    %piece1 = pc_out_narrow_full(:, 1+cutFir_start : N_gate_narrow+cutFir_start);
    piece1 = pc_out_narrow_full(:, 1+cutFir_start : N_gate_narrow+cutFir_start);   % 取窄脉冲处理后信号的5-232点（共计228点）
    % 中脉冲结果段
    % 注意：这里我们截取的是中脉冲处理结果中，紧接着窄脉冲之后的部分
    idx2_start = N_gate_narrow + 1;
    idx2_end   = N_gate_narrow + N_gate_medium;
    piece2 = pc_out_medium_full(:, idx2_start:idx2_end);           % 取中脉冲处理后信号的5-232点（共计228点）

    % 长脉冲结果段
    idx3_start = N_gate_narrow + N_gate_medium + 1;
    idx3_end   = N_gate_narrow + N_gate_medium + N_gate_long;
    piece3 = pc_out_long_full(:, idx3_start:idx3_end);

    pc_results_13beam(:,:,b) = [piece1, piece2, piece3];
end

fprintf('脉冲压缩完成。拼接后的总距离门数为: %d\n', size(pc_results_13beam, 2));

%% 7. MTD 处理 (对所有波束)
%=========================================================================
fprintf('--- 7. 执行MTD处理 ---\n');
% MTD处理前加窗
beta_MTD = 4.5;
MTD_win = kaiser(config.Sig_Config.prtNum, beta_MTD);
% 对13个波束的脉压结果，逐一进行MTD
rdm_13beam = complex(zeros(size(pc_results_13beam)));

for b = 1:config.Sig_Config.beam_num
    % 提取单个波束的脉压数据立方体
    pc_data_single_beam = squeeze(pc_results_13beam(:,:,b));
    % 加窗
    pc_data_single_beam_win = pc_data_single_beam .* MTD_win;
    % 沿慢时间维(第1维)做FFT，从而得到目标速度，并进行fftshift将零频移到中心
    rdm_13beam(:,:,b) = fftshift(fft(pc_data_single_beam_win, [], 1), 1);
    
end

fprintf('MTD处理完成，已生成13个波束的距离-多普勒数据。\n');

%% 8. 多波束CFAR检测 (新)
% =========================================================================
fprintf('--- 8. 执行多波束CFAR检测 ---\n');

% --- 8.1 初始化总的检测列表 ---
% 用于存储在所有和波束中发现的所有原始点迹
% 格式: [速度索引, 距离索引, 波束对索引, 和波束幅值]
all_raw_detections = []; 

% --- 8.2 循环遍历所有12个相邻波束对 ---
for pair_idx = 1:(config.Sig_Config.beam_num - 1)
    
    % --- 动态形成当前用于检测的和波束 ---
    beam_idx_A = pair_idx;
    beam_idx_B = pair_idx + 1;
    rdm_beamA = abs(squeeze(rdm_13beam(:,:,beam_idx_A)));
    rdm_beamB = abs(squeeze(rdm_13beam(:,:,beam_idx_B)));
    rdm_for_cfar = rdm_beamA + rdm_beamB; % 用于当前检测的波束对组成的和波束
    rdm_for_cfar_all(:,:,pair_idx) = rdm_for_cfar; % 保存12组和波束的RD图
    fprintf('正在对波束对 #%d-#%d 进行CFAR检测...\n', beam_idx_A, beam_idx_B);

    % --- 8.3 对当前和波束执行 2D-GOCA-CFAR 检测器 ---
    % 十字形CFAR检测，距离维、速度维各自独立计算噪声估计后选大，顺序无影响
    T_CFAR = cfar_params.T_CFAR;% 门限因子，它决定了门限要比噪声高多少倍
    guard_R = cfar_params.guardCells_R;     % 距离维保护单元数
    guard_V = cfar_params.guardCells_V;     % 速度维保护单元数
    ref_R = cfar_params.refCells_R;         % 距离维参考单元数
    ref_V = cfar_params.refCells_V;         % 速度维参考单元数
    [num_V, num_R] = size(rdm_for_cfar);    % 和波束
    cfar_detections = zeros(num_V, num_R);

    % 遍历所有待检测单元 (CUT - Cell Under Test)
    % 循环的起始和结束点留出了足够的边界，以确保CFAR窗口不会超出矩阵范围
    for r = (ref_R + guard_R + 1) : (num_R - ref_R - guard_R)
        for v = (ref_V + guard_V + 1) : (num_V - ref_V - guard_V)
            % 提取当前待检测单元的幅值
            cut_power = rdm_for_cfar(v, r);
            
            % --- 距离维 GOCA ---
            % 定义前向参考窗（距离更近的单元）
            range_win_leading = rdm_for_cfar(v, r - guard_R - ref_R : r - guard_R - 1);
            % 定义后向参考窗（距离更远的单元）
            range_win_trailing = rdm_for_cfar(v, r + guard_R + 1 : r + guard_R + ref_R);
            % 分别计算两个窗的平均功率
            noise_R_leading = mean(range_win_leading);
            noise_R_trailing = mean(range_win_trailing);
            % GOCA核心：选取两个窗中较大的平均值作为距离维的噪声估计
            noise_R = max(noise_R_leading, noise_R_trailing);

            % --- 速度维 GOCA ---
            % 定义前向参考窗 (速度更小的单元)
            doppler_win_leading = rdm_for_cfar(v - guard_V - ref_V : v - guard_V - 1, r);
            % 定义后向参考窗 (速度更大的单元)
            doppler_win_trailing = rdm_for_cfar(v + guard_V + 1 : v + guard_V + ref_V, r);
            % 分别计算平均幅值
            noise_V_leading = mean(doppler_win_leading);
            noise_V_trailing = mean(doppler_win_trailing);
            % GOCA核心：选取两个窗中较大的平均值作为速度维的噪声估计
            noise_V = max(noise_V_leading, noise_V_trailing);
            
            % --- 综合噪声估计 (取两者中的较大值以保证虚警率) ---
            noise_estimate = max(noise_R, noise_V);
            
            % --- 生成门限并进行比较 ---
            threshold = T_CFAR * noise_estimate;         % 最终的检测门限 = 噪声估计 * 门限因子
            threshold_map(v, r, pair_idx) = threshold;   % 保存当前单元的门限值，方便后续分析
            
            % 如果待检测单元的功率超过了动态计算出的门限，则判定为目标，即在结果矩阵对应位置记为1
            if cut_power > threshold            
                cfar_detections(v, r) = 1;        
            end
        end
    end
    
    % --- 8.4 收集当前和波束的检测结果 ---
    [detected_v_idx, detected_r_idx] = find(cfar_detections);   
    num_detections_current_pair = length(detected_v_idx);
    
    if num_detections_current_pair > 0
        fprintf('  > 在波束对 #%d-#%d 中发现 %d 个点迹。\n', beam_idx_A, beam_idx_B, num_detections_current_pair);
        % 将当前找到的点及其信息（包括所在的波束对索引）添加到总列表中
        for i = 1:num_detections_current_pair
            v_idx = detected_v_idx(i);
            r_idx = detected_r_idx(i);
            rdm_abs = rdm_for_cfar(v_idx, r_idx);
            all_raw_detections(end+1, :) = [v_idx, r_idx, pair_idx, rdm_abs];
        end
    end
end

num_total_raw_detections = size(all_raw_detections, 1);
fprintf('\n多波束CFAR检测完成，在12组和波束中总共发现 %d 个原始点迹。\n', num_total_raw_detections);


%% 9. 逐点参数估计 (含单脉冲测角) 
% =========================================================================
fprintf('--- 9. 对每个原始点迹进行参数估计 (三次样条插值) ---\n');
parameterized_detections = []; % 初始化结构体数组

if num_total_raw_detections > 0
    % --- 定义坐标轴和标定数据 ---
    N_total_gate = size(rdm_13beam, 2);
    range_axis = (0:N_total_gate-1) * (config.Sig_Config.c / (2 * config.Sig_Config.fs));
    v_max = config.Sig_Config.wavelength / (2 * config.Sig_Config.prt);
    % --- 核心修正：确保速度轴范围正确 ---
    velocity_axis = linspace(-v_max/2, v_max/2, config.Sig_Config.prtNum); 
    [num_V, num_R] = size(squeeze(rdm_13beam(:,:,1)));
    beam_angles_deg = [-16, -9.6, -3.2, 3.2, 9.6, 16, 22.6, 29.2, 36.1, 43.3, 51, 59.6, 70.3]; % 13个波束波峰对应的角度值
    k_slopes_LUT = [-4.6391,-4.6888,-4.7578,-4.7891,-4.7214,-4.7513,-5.2343,-5.4529,-5.7323,-6.1685,-7.0256,-8.7612]; % 测角K值

    % --- 插值参数 ---
    extraDots = 2;      % 在峰值左右/上下各取几个点进行插值 (共 2*extraDots + 1 个点)
    rInterpTimes = 8;   % 距离维插值倍数(每个单元格内插多少个点)
    vInterpTimes = 4;   % 速度维插值倍数(每个单元格内插多少个点)
    
    % --- 遍历每一个原始目标检测点迹 ---
    for i = 1:num_total_raw_detections
        v_idx = all_raw_detections(i, 1);   % 峰值速度索引 (整数)
        r_idx = all_raw_detections(i, 2);   % 峰值距离索引 (整数)
        pair_idx = all_raw_detections(i, 3);% 波束对索引
        power = all_raw_detections(i, 4);   % 峰值幅值 (线性) - 注意参考代码是在幅度上插值
        
        % --- 获取用于插值的RDM数据 (和波束幅度) ---
        % !! 注意: rdm_for_cfar_all 必须是在 %% 8. 中保存的包含所有和波束幅度数据的 3D 矩阵
        if ~exist('rdm_for_cfar_all', 'var') || size(rdm_for_cfar_all,3) < pair_idx
            error('变量 rdm_for_cfar_all 未定义或维度不正确，请检查 %% 8. CFAR 检测部分');
        end
        rdm_interp = rdm_for_cfar_all(:,:,pair_idx); % 提取对应的和波束幅度RDM
        
        % --- 距离维峰值样条插值 ---
        rCellsFix = (r_idx - extraDots) : (r_idx + extraDots);      % 选取插值窗口的索引（要检测的距离采样单元）
        % 边界处理：确保窗口索引在有效范围内
        rCellsFix = rCellsFix(rCellsFix >= 1 & rCellsFix <= num_R); 
        if length(rCellsFix) < 3 % 样条插值至少需要3个点
            rCellMax = r_idx; % 点数不足，无法插值，使用原始峰值索引
        else
            rCellsFix = sort(rCellsFix); % 由大到小排列
            % 插值
            mtdDataUsed_r = rdm_interp(v_idx, rCellsFix); % 取出该目标所在位置附近的数据进行插值
            rCellsFixQ = rCellsFix(1): 1/rInterpTimes :rCellsFix(end); % 距离维插值变量
            mtdDataUsedQ_r = interp1(rCellsFix-rCellsFix(1),mtdDataUsed_r,rCellsFixQ-rCellsFix(1),'spline');
            [M1,I1] = max(mtdDataUsedQ_r);
            if ~isempty(I1)
                % rCellMax = rCellsFixQ(I1);     % 幅度值最大的位置
                rCellMax = rCellsFixQ(I1(1));    % 幅度值最大的位置
                est_range = range_axis(r_idx) + (rCellMax-r_idx)*deltaR;   % 目标测量距离
            else
                est_range = [];
            end
        end

       
        % --- 速度维峰值样条插值 ---
        vCellsFix = (v_idx - extraDots) : (v_idx + extraDots); % 选取插值窗口的索引（要检测的速度单元）
        % 边界处理
        vCellsFix = vCellsFix(vCellsFix >= 1 & vCellsFix <= num_V);
        if length(vCellsFix) < 3
            vCellMax = v_idx; % 点数不足，无法插值
        else
            vCellsFix = sort(vCellsFix); % 由大到小排列
            % 插值
            mtdDataUsed_v = rdm_interp(vCellsFix,r_idx);  % 取出目标所在位置附近的数据进行插值
            vCellsFixQ = vCellsFix(1): 1/vInterpTimes :vCellsFix(end);   % 距离维插值变量
            mtdDataUsedQ_v = interp1(vCellsFix-vCellsFix(1),mtdDataUsed_v,vCellsFixQ-vCellsFix(1),'spline');
            [M2,I2] = max(mtdDataUsedQ_v);
            if ~isempty(I2) && ~isempty(est_range)
                % vCellMax = vCellsFixQ(I2);     % 幅度值最大的位置
                vCellMax = vCellsFixQ(I2(1));    % 幅度值最大的位置
                est_velocity = velocity_axis(fix(vCellMax)) + (vCellMax-fix(vCellMax))*deltaV;   % 目标测量速度
            else
                est_velocity = [];
            end
        end

        % --- 角度估计 (单脉冲测角) ---
        % (测角部分不变，仍然使用原始整数索引 v_idx, r_idx 查找复数值)
        beam_A_idx = pair_idx;
        beam_B_idx = pair_idx + 1;
        S_A = rdm_13beam(v_idx, r_idx, beam_A_idx);
        S_B = rdm_13beam(v_idx, r_idx, beam_B_idx);
        monopulse_ratio = (S_A - S_B) / (S_A + S_B + eps);
        k_slope = k_slopes_LUT(beam_A_idx);
        angle_A = beam_angles_deg(beam_A_idx);
        angle_B = beam_angles_deg(beam_B_idx);
        angle_offset = k_slope * real(monopulse_ratio);
        est_angle = (angle_A + angle_B)/2 + angle_offset;
        
        % --- 存储结果 ---
        parameterized_detections(i).Range = est_range;
        parameterized_detections(i).Velocity = est_velocity;
        parameterized_detections(i).Angle = est_angle;
        parameterized_detections(i).Power = power; % 存储的是和波束的幅度值
        parameterized_detections(i).PairIndex = pair_idx;
    end
end
fprintf('已完成对 %d 个原始点迹的参数估计 (含三次样条插值)。\n', num_total_raw_detections);

%% 10. 最终聚类 (基于物理参数)
% =========================================================================
% 基于广度优先搜索（Breadth-First Search, BFS）的聚类归并算法
fprintf('--- 10. 执行最终聚类 (基于插值后的物理参数) ---\n');
% 初始化最终的目标列表，准备存放聚类合并后的结果
final_targets = [];
% CFAR检测到了至少一个原始点迹时才执行聚类
if num_total_raw_detections > 0
    % --- 聚类算法核心 ---
    % 创建一个与原始点迹数量相同的向量，用于存储每个点迹所属的簇ID
    % 0 表示 "尚未分配给任何簇"
    cluster_IDs = zeros(num_total_raw_detections, 1);
    current_cluster_ID = 0;
    % 遍历每一个原始点迹 (i 是点迹的索引)
    for i = 1:num_total_raw_detections
        if cluster_IDs(i) == 0
            current_cluster_ID = current_cluster_ID + 1;
            points_to_visit = i;
            
            % --- 广度优先搜索 (BFS)：只要 "待访问" 列表不为空，就继续扩展当前簇 ---
            while ~isempty(points_to_visit)
                % 1. 出队：从 "待访问" 列表中取出第一个点迹
                current_idx = points_to_visit(1);
                points_to_visit(1) = [];
                % 2. 访问：再次检查这个点是否已被访问 (可能被其他邻居重复加入)，如果未被访问 (ID为0)，则处理它
                if cluster_IDs(current_idx) == 0
                    % 3. 标记：将当前簇的ID分配给这个点
                    cluster_IDs(current_idx) = current_cluster_ID;
                    % 4. 扩展：遍历所有其他点迹 (j)，寻找 current_idx 的 "邻居"
                    for j = 1:num_total_raw_detections
                        if cluster_IDs(j) == 0
                            % --- 核心：定义 "邻居" 关系 ---
                            % 提取两个点 (current_idx 和 j) 的物理参数
                            % 计算两点之间的 "距离" (在三个物理维度上)
                            dist_r = abs(parameterized_detections(current_idx).Range - parameterized_detections(j).Range);       % 维度1：距离 (Range)
                            dist_v = abs(parameterized_detections(current_idx).Velocity - parameterized_detections(j).Velocity); % 维度2：速度 (Velocity)
                            dist_a = abs(parameterized_detections(current_idx).Angle - parameterized_detections(j).Angle);       % 维度3：角度 (Angle)
                            % --- 核心：基于物理参数进行比较 ---
                            % --- 邻居判断 (聚类准则) ：如果两点在所有三个维度上的距离都小于设定的阈值，则认为 j 是 current_idx 的 "邻居"，属于同一个簇 ---
                            if dist_r <= cluster_params.max_range_sep && dist_v <= cluster_params.max_vel_sep && dist_a <= cluster_params.max_angle_sep
                                % 5. 入队：将这个新邻居 j 添加到 "待访问" 列表的末尾，它将在后续的循环中被访问和扩展
                                points_to_visit(end+1) = j;
                            end
                        end
                    end
                end
            end
        end
    end

    % --- 聚类完成，开始合并 ---
    num_clusters = current_cluster_ID;  % 记录找到的总簇数
    fprintf('聚类完成，将 %d 个参数化点迹合并为 %d 个唯一目标。\n', num_total_raw_detections, num_clusters);

    % --- 对每个簇进行合并 (功率（幅值）加权平均) ---
    final_targets = struct('Range', cell(1,num_clusters), 'Velocity', cell(1,num_clusters), 'Angle', cell(1,num_clusters));
    fprintf('\n--- 最终目标列表 (聚类后) ---\n');
    
    % 遍历每一个找到的簇 (从 1 到 num_clusters)
    for i = 1:num_clusters
        cluster_mask = (cluster_IDs == i);
        detections_in_cluster = parameterized_detections(cluster_mask);
        
        powers = [detections_in_cluster.Power]';
        total_power = sum(powers);
        
        % 功率加权平均：(参数1 * 功率1 + 参数2 * 功率2 + ...) / (总功率)
        % 功率越大的点迹，在计算最终平均值时的权重就越大
        final_range = sum([detections_in_cluster.Range]' .* powers) / total_power;
        final_velocity = sum([detections_in_cluster.Velocity]' .* powers) / total_power;
        final_angle = sum([detections_in_cluster.Angle]' .* powers) / total_power;
        
        final_targets(i).Range = final_range;
        final_targets(i).Velocity = final_velocity;
        final_targets(i).Angle = final_angle;
        
        fprintf('目标 %d: 距离 ~%.2f m, 速度 ~%.2f m/s, 估算角度 ~%.2f 度\n', ...
            i, final_range, final_velocity, final_angle);
    end
else
    fprintf('无原始点迹，无需进行聚类。\n');
end

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
plot(R_rx-24, abs(pc_results_13beam(prt_selected,:,beam_selected))); 
title('幅值');
xlabel('距离（m）');
ylabel('信号幅度');
for i=1:length(targets)
    xline(targets(i).Range,'--r','目标','LabelOrientation','horizontal') 
end

subplot(4,1,4);
plot(R_rx-24, 20*log10(abs(pc_results_13beam(prt_selected,:,beam_selected)))); 
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

% 图 7.5：MTD杂波抑制效果对比图（没什么用，可以删了）
% figure('Name', 'MTD结果: 杂波抑制对比');
% % MTD处理前 (取第一个PRT的脉压结果)
% range_profile_before_mtd = abs(squeeze(pc_results_13beam(1, :, beam_to_display)));
% 
% % MTD处理后 (对所有非零多普勒单元的能量进行累加，在一段持续时间（332个PRT）内)
% rdm_power = abs(rdm_single_beam).^2;
% % 将零多普勒通道置零 (或用一个极小值代替以避免log(0)错误)
% rdm_power(floor(config.Sig_Config.prtNum/2)+1, :) = 1e-10; 
% range_profile_after_mtd = sum(rdm_power, 1);
% 
% % 绘图 (转换为dB)
% plot(range_axis, 20*log10(range_profile_before_mtd), 'b-', 'LineWidth', 1);
% hold on;
% plot(range_axis, 10*log10(range_profile_after_mtd), 'r-', 'LineWidth', 2);
% grid on;
% xlabel('距离 (m)');
% ylabel('幅度 (dB)');
% title(sprintf('MTD杂波抑制效果对比 - 波束 #%d', beam_to_display));
% legend('MTD处理前 (单PRT)', 'MTD处理后 (332prt内积累动目标)');
% % 动态调整Y轴范围以便观察
% max_val = max(20*log10(range_profile_before_mtd));
% if isfinite(max_val)
%     ylim([max_val - 80, max_val + 10]); 
% end


%% 8. 最终结果可视化
% =========================================================================
fprintf('--- 8. 生成最终可视化图表 ---\n');
% --- 8.1 和波束的幅值与CFAR动态检测门限图 ---

% figure;
% plot(threshold_map(clusters{1}(1), :,beam_selected), 'r');
% hold on;
% plot(rdm_for_cfar_all(clusters{1}(1), :,beam_selected),'b');
% legend('检测门限','和波束幅度');
% hold off;
% title('和波束幅度与检测门限情况');

% --- 8.2 和波束的距离-多普勒图 ---
figure('Name', '最终检测结果');
subplot(2,1,1);
v_max = config.Sig_Config.wavelength / (2 * config.Sig_Config.prt);
velocity_axis = linspace(-v_max/2, v_max/2, config.Sig_Config.prtNum);
range_axis_cfar = (0:config.Sig_Config.point_PRT-1) * (config.Sig_Config.c / (2 * config.Sig_Config.fs));
imagesc(range_axis_cfar, velocity_axis, rdm_for_cfar_all(:,:,beam_selected));
colormap('jet'); colorbar; axis xy;
title(sprintf('用于CFAR的和波束 (波束 #%d + #%d) RDM', beam_idx_A, beam_idx_B));
xlabel('距离 (m)');
ylabel('速度 (m/s)');
% 
% --- 8.3 CFAR检测结果叠加显示 ---
subplot(2,1,2);
imagesc(range_axis_cfar, velocity_axis, rdm_for_cfar_all(:,:,beam_selected));
colormap('jet'); colorbar; axis xy;
hold on;

% 在图上标记检测到的目标
if num_clusters > 0
    detected_ranges = range_axis_cfar(detected_r_idx); % 这里点的数据应该是聚类后的
    detected_velocities = velocity_axis(detected_v_idx);
    plot(detected_ranges, detected_velocities, 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'CFAR 检测点');
end
title('CFAR检测结果叠加图');
xlabel('距离 (m)');
ylabel('速度 (m/s)');
legend;

%% 打印所有波束的检测点（聚类归并前）

% 绘制所有原始CFAR检测点
figure('Name', '所有原始CFAR检测点');
% 提取点迹的距离速度序号
v_indices = all_raw_detections(:, 1);
r_indices = all_raw_detections(:, 2);
% 把速度距离序号转换为具体速度距离值
detected_ranges_m = range_axis(r_indices);
detected_velocities_ms = velocity_axis(v_indices);
scatter(detected_ranges_m, detected_velocities_ms, 30, 'b', 'filled');
grid on;
xlabel('距离 (m)');
ylabel('速度 (m/s)');
title(sprintf('所有原始CFAR检测点迹图 (%d 个点)', size(all_raw_detections, 1)));
%fprintf('已将 %d 个原始CFAR检测点绘制在图中。\n', size(all_raw_detections, 1));

%% 绘制所有原始CFAR检测点的距离速度波束图
% =========================================================================
figure('Name', '所有原始CFAR检测点 (距离-角度)');

if isempty(parameterized_detections)
    fprintf('没有参数化的点迹，无法绘图。\n');
else
    % --- 1. 正确地从结构体数组中提取数据 ---
    % 使用 [struct_array.FieldName] 语法可以一次性提取所有同名字段的数据
    ranges = [parameterized_detections.Range];
    velocities = [parameterized_detections.Velocity];
    pair_indices = [parameterized_detections.PairIndex];

    % --- 2. 使用scatter函数绘制点迹图 ---
    % X轴: 距离
    % Y轴: 俯仰角度
    % 颜色: 由波束对索引决定
    scatter(ranges, velocities, 50, pair_indices, 'filled');

    % --- 3. 添加标注和颜色条 ---
    grid on;
    xlabel('距离 (m)');
    ylabel('速度 (m/s)');
    title(sprintf('所有原始CFAR检测点迹 (%d 个点)', length(ranges)));
    
    % 添加颜色条，并为其添加标签
    h = colorbar;
    ylabel(h, '波束对 索引');
    
    fprintf('已将 %d 个参数化点迹的距离-速度图绘制完成。\n', length(ranges));
end

%% 绘制所有原始CFAR检测点的距离俯仰角度波束图
% =========================================================================
figure('Name', '所有原始CFAR检测点 (距离-角度)');

if isempty(parameterized_detections)
    fprintf('没有参数化的点迹，无法绘图。\n');
else
    % --- 1. 正确地从结构体数组中提取数据 ---
    % 使用 [struct_array.FieldName] 语法可以一次性提取所有同名字段的数据
    ranges = [parameterized_detections.Range];
    Powers = [parameterized_detections.Angle];
    pair_indices = [parameterized_detections.PairIndex];

    % --- 2. 使用scatter函数绘制点迹图 ---
    % X轴: 距离
    % Y轴: 俯仰角度
    % 颜色: 由波束对索引决定
    scatter(ranges, Powers, 50, pair_indices, 'filled');

    % --- 3. 添加标注和颜色条 ---
    grid on;
    xlabel('距离 (m)');
    ylabel('俯仰角度 (度)');
    title(sprintf('所有原始CFAR检测点迹 (%d 个点)', length(ranges)));
    
    % 添加颜色条，并为其添加标签
    h = colorbar;
    ylabel(h, '波束对 索引');
    
    fprintf('已将 %d 个参数化点迹的距离-俯仰角度图绘制完成。\n', length(ranges));
end

%% 绘制所有原始CFAR检测点的速度俯仰角度波束图
% =========================================================================
figure('Name', '所有原始CFAR检测点 (距离-角度)');

if isempty(parameterized_detections)
    fprintf('没有参数化的点迹，无法绘图。\n');
else
    % --- 1. 正确地从结构体数组中提取数据 ---
    % 使用 [struct_array.FieldName] 语法可以一次性提取所有同名字段的数据
    velocities = [parameterized_detections.Velocity];
    Powers = [parameterized_detections.Angle];
    pair_indices = [parameterized_detections.PairIndex];

    % --- 2. 使用scatter函数绘制点迹图 ---
    % X轴: 距离
    % Y轴: 俯仰角度
    % 颜色: 由波束对索引决定
    scatter(velocities, Powers, 50, pair_indices, 'filled');

    % --- 3. 添加标注和颜色条 ---
    grid on;
    xlabel('速度 (m/s)');
    ylabel('俯仰角度 (度)');
    title(sprintf('所有原始CFAR检测点迹 (%d 个点)', length(ranges)));
    
    % 添加颜色条，并为其添加标签
    h = colorbar;
    ylabel(h, '波束对 索引');
    
    fprintf('已将 %d 个参数化点迹的速度-俯仰角度图绘制完成。\n', length(ranges));
end



%% 绘制所有原始CFAR检测点的距离-信号幅度波束图
% =========================================================================
figure('Name', '所有原始CFAR检测点 (距离-角度)');

if isempty(parameterized_detections)
    fprintf('没有参数化的点迹，无法绘图。\n');
else
    % --- 1. 正确地从结构体数组中提取数据 ---
    % 使用 [struct_array.FieldName] 语法可以一次性提取所有同名字段的数据
    ranges = [parameterized_detections.Range];
    Powers = [parameterized_detections.Power];
    Powers = 20*log10(Powers);
    pair_indices = [parameterized_detections.PairIndex];

    % --- 2. 使用scatter函数绘制点迹图 ---
    % X轴: 距离
    % Y轴: 俯仰角度
    % 颜色: 由波束对索引决定
    scatter(ranges, Powers, 50, pair_indices, 'filled');

    % --- 3. 添加标注和颜色条 ---
    grid on;
    xlabel('距离 (m)');
    ylabel('信号幅值');
    title(sprintf('所有原始CFAR检测点迹 (%d 个点)', length(ranges)));
    
    % 添加颜色条，并为其添加标签
    h = colorbar;
    ylabel(h, '波束对 索引');
    
    fprintf('已将 %d 个参数化点迹的距离-信号幅度图绘制完成。\n', length(ranges));
end



% figure('Name', '所有原始CFAR检测点 (距离-角度)');
% 
% if isempty(parameterized_detections)
%     fprintf('没有参数化的点迹，无法绘图。\n');
% else
%     % --- 1. 正确地从结构体数组中提取数据 ---
%     ranges = [parameterized_detections.Range];
%     angles = [parameterized_detections.Angle];
%     pair_indices = [parameterized_detections.PairIndex];
% 
%     % --- 2. 查找数据中存在哪些唯一的波束对 ---
%     unique_pairs = unique(pair_indices);
%     num_unique_pairs = length(unique_pairs);
% 
%     % --- 3. 创建一个颜色查找表 ---
%     % lines(N) 会生成 N 个视觉上区分度高的颜色
%     colors = lines(num_unique_pairs);
% 
%     % --- 4. 循环绘制每个波束对的点迹 ---
%     hold on; % 打开"保持"开关，准备叠加绘制
% 
%     for i = 1:num_unique_pairs
%         % 获取当前的波束对索引
%         current_pair = unique_pairs(i);
% 
%         % 创建一个逻辑掩码，找到所有属于当前波束对的点
%         mask = (pair_indices == current_pair);
% 
%         % 提取这些点的距离和角度坐标
%         ranges_for_this_pair = ranges(mask);
%         angles_for_this_pair = angles(mask);
% 
%         % 使用scatter函数，为当前波束对的点指定一种特定颜色
%         % 'DisplayName' 用于在图例中显示文字
%         scatter(ranges_for_this_pair, angles_for_this_pair, 70, colors(i, :), 'filled', ...
%             'DisplayName', sprintf('波束对 #%d', current_pair));
%     end
% 
%     hold off; % 关闭"保持"开关
% 
%     % --- 5. 添加标注和普通图例 ---
%     grid on;
%     xlabel('距离 (m)');
%     ylabel('俯仰角度 (度)');
%     title(sprintf('所有原始CFAR检测点迹 (%d 个点)', length(ranges)));
% 
%     % 关键：调用 legend('show') 来显示我们在循环中设置的DisplayName
%     legend('show', 'Location', 'bestoutside');
% 
%     fprintf('已将 %d 个参数化点迹的距离-角度图绘制完成。\n', length(ranges));
% end

%% 单脉冲全程滤波（脉压）结果
% 测试用，不用就注释掉，计算量较大
% beam_data_temp = squeeze(iq_data_13beam(:,:,beam_selected));
% 
% % 窄脉冲 (时域FIR滤波)
% % 非加窗
% fir_coeffs = MF_narrow;
% % filter函数会逐行处理矩阵
% pc_out_narrow_full_temp = filter(fir_coeffs, 1, beam_data_temp, [], 2);
% % 校正FIR引入的延迟
% pc_out_narrow_full_temp = circshift(pc_out_narrow_full_temp, -fir_delay, 2);
% % 加窗
% window = hanning(length(MF_narrow))';
% fir_coeffs_windowed = MF_narrow.*window;
% % filter函数会逐行处理矩阵
% pc_out_narrow_full_temp_win = filter(fir_coeffs_windowed, 1, beam_data_temp, [], 2);
% % 校正FIR引入的延迟
% pc_out_narrow_full_temp_win = circshift(pc_out_narrow_full_temp_win, -fir_delay, 2);
% 
% 
% % 中脉冲 (非加窗时域卷积)
% pc_out_medium_full_temp = complex(zeros(size(beam_data_temp,1), size(beam_data_temp,2) + num_samples_2 - 1));
% for i = 1:config.Sig_Config.prtNum
%     pc_out_medium_full_temp(i,:) = conv(beam_data_temp(i,:), MF_medium);
% end
% % 中脉冲 (加窗时域卷积)
% pc_out_medium_full_temp_win = complex(zeros(size(beam_data_temp,1), size(beam_data_temp,2) + num_samples_2 - 1));
% for i = 1:config.Sig_Config.prtNum
%     pc_out_medium_full_temp_win(i,:) = conv(beam_data_temp(i,:), MF_medium_win);
% end
% 
% 
% % 长脉冲 (非加窗时域卷积)
% pc_out_long_full_temp = complex(zeros(size(beam_data_temp,1), size(beam_data_temp,2) + num_samples_3 - 1));
% for i = 1:config.Sig_Config.prtNum
%     pc_out_long_full_temp(i,:) = conv(beam_data_temp(i,:), MF_long);
% end
% % 长脉冲 (加窗时域卷积)
% pc_out_long_full_temp_win = complex(zeros(size(beam_data_temp,1), size(beam_data_temp,2) + num_samples_3 - 1));
% for i = 1:config.Sig_Config.prtNum
%     pc_out_long_full_temp_win(i,:) = conv(beam_data_temp(i,:), MF_long_win);
% end
% 
% 
% % 非加窗窄脉冲prt上滤波结果
% RX1 = linspace(0, length(pc_out_narrow_full_temp)*6, length(pc_out_narrow_full_temp));
% figure;
% subplot(4,1,1)
% plot(RX1, real(pc_out_narrow_full_temp(1,:))); % 输出脉冲压缩结果
% title('窄脉冲全程prt上滤波信号幅度');
% subtitle('实部');
% xlabel('距离 (m)');
% ylabel('幅度');
% 
% subplot(4,1,2)
% plot(RX1, imag(pc_out_narrow_full_temp(1,:))); % 输出脉冲压缩结果
% subtitle('虚部');
% xlabel('距离 (m)');
% ylabel('幅度');
% 
% subplot(4,1,3)
% plot(RX1, abs(pc_out_narrow_full_temp(1,:))); % 输出脉冲压缩结果
% subtitle('幅值');
% xlabel('距离 (m)');
% ylabel('幅度');
% 
% subplot(4,1,4)
% plot(RX1, 20*log10(abs(pc_out_narrow_full_temp(1,:)))); % 输出脉冲压缩结果
% subtitle('幅值（dB）');
% xlabel('距离 (m)');
% ylabel('幅度');
% % 加窗窄脉冲prt上滤波结果
% RX1 = linspace(0, length(pc_out_narrow_full_temp_win)*6, length(pc_out_narrow_full_temp_win));
% figure;
% subplot(4,1,1)
% plot(RX1, real(pc_out_narrow_full_temp_win(1,:))); % 输出脉冲压缩结果
% title('窄脉冲全程prt上滤波信号幅度（加窗后）');
% subtitle('实部');
% xlabel('距离 (m)');
% ylabel('幅度');
% 
% subplot(4,1,2)
% plot(RX1, imag(pc_out_narrow_full_temp_win(1,:))); % 输出脉冲压缩结果
% subtitle('虚部');
% xlabel('距离 (m)');
% ylabel('幅度');
% 
% subplot(4,1,3)
% plot(RX1, abs(pc_out_narrow_full_temp_win(1,:))); % 输出脉冲压缩结果
% subtitle('幅值');
% xlabel('距离 (m)');
% ylabel('幅度');
% 
% subplot(4,1,4)
% plot(RX1, 20*log10(abs(pc_out_narrow_full_temp_win(1,:)))); % 输出脉冲压缩结果
% subtitle('幅值（dB）');
% xlabel('距离 (m)');
% ylabel('幅度');
% 
% % 非加窗中脉冲prt上脉压结果
% RX2 = linspace(0, length(pc_out_medium_full_temp)*6, length(pc_out_medium_full_temp));
% figure;
% subplot(4,1,1);
% plot(RX2, real(pc_out_medium_full_temp(1,:))); % 输出脉冲压缩结果
% title('中脉冲全程prt上脉压后信号幅度');
% subtitle('实部');
% xlabel('距离 (m)');
% ylabel('幅度');
% 
% subplot(4,1,2);
% plot(RX2, imag(pc_out_medium_full_temp(1,:))); % 输出脉冲压缩结果
% subtitle('虚部');
% xlabel('距离 (m)');
% ylabel('幅度');
% 
% subplot(4,1,3);
% plot(RX2, abs(pc_out_medium_full_temp(1,:))); % 输出脉冲压缩结果
% subtitle('幅值');
% xlabel('距离 (m)');
% ylabel('幅度');
% 
% subplot(4,1,4);
% plot(RX2, 20*log10(abs(pc_out_medium_full_temp(1,:)))); % 输出脉冲压缩结果
% subtitle('幅值（dB）');
% xlabel('距离 (m)');
% ylabel('幅度');
% % 加窗中脉冲prt上脉压结果
% RX2 = linspace(0, length(pc_out_medium_full_temp_win)*6, length(pc_out_medium_full_temp_win));
% figure;
% subplot(4,1,1);
% plot(RX2, real(pc_out_medium_full_temp_win(1,:))); % 输出脉冲压缩结果
% title('中脉冲全程prt上脉压后信号幅度（加窗后）');
% subtitle('实部');
% xlabel('距离 (m)');
% ylabel('幅度');
% 
% subplot(4,1,2);
% plot(RX2, imag(pc_out_medium_full_temp_win(1,:))); % 输出脉冲压缩结果
% subtitle('虚部');
% xlabel('距离 (m)');
% ylabel('幅度');
% 
% subplot(4,1,3);
% plot(RX2, abs(pc_out_medium_full_temp_win(1,:))); % 输出脉冲压缩结果
% subtitle('幅值');
% xlabel('距离 (m)');
% ylabel('幅度');
% 
% subplot(4,1,4);
% plot(RX2, 20*log10(abs(pc_out_medium_full_temp_win(1,:)))); % 输出脉冲压缩结果
% subtitle('幅值（dB）');
% xlabel('距离 (m)');
% ylabel('幅度');
% 
% % 非加窗长脉冲prt上脉压结果
% RX3 = linspace(0, length(pc_out_long_full_temp)*6, length(pc_out_long_full_temp));
% figure;
% subplot(4,1,1);
% plot(RX3, real(pc_out_long_full_temp(1,:))); % 输出脉冲压缩结果
% title('长脉冲全程prt脉压后信号幅度');
% subtitle('实部');
% xlabel('距离 (m)');
% ylabel('幅度');
% 
% subplot(4,1,2);
% plot(RX3, imag(pc_out_long_full_temp(1,:))); % 输出脉冲压缩结果
% subtitle('虚部');
% xlabel('距离 (m)');
% ylabel('幅度');
% 
% subplot(4,1,3);
% plot(RX3, abs(pc_out_long_full_temp(1,:))); % 输出脉冲压缩结果
% subtitle('幅值');
% xlabel('距离 (m)');
% ylabel('幅度');
% 
% subplot(4,1,4);
% plot(RX3, 20*log10(abs(pc_out_long_full_temp(1,:)))); % 输出脉冲压缩结果
% subtitle('幅值（dB）');
% xlabel('距离 (m)');
% ylabel('幅度');
% % 加窗长脉冲prt上脉压结果
% RX3 = linspace(0, length(pc_out_long_full_temp_win)*6, length(pc_out_long_full_temp_win));
% figure;
% subplot(4,1,1);
% plot(RX3, real(pc_out_long_full_temp_win(1,:))); % 输出脉冲压缩结果
% title('长脉冲全程prt脉压后信号幅度（加窗后）');
% subtitle('实部');
% xlabel('距离 (m)');
% ylabel('幅度');
% 
% subplot(4,1,2);
% plot(RX3, imag(pc_out_long_full_temp_win(1,:))); % 输出脉冲压缩结果
% subtitle('虚部');
% xlabel('距离 (m)');
% ylabel('幅度');
% 
% subplot(4,1,3);
% plot(RX3, abs(pc_out_long_full_temp_win(1,:))); % 输出脉冲压缩结果
% subtitle('幅值');
% xlabel('距离 (m)');
% ylabel('幅度');
% 
% 
% subplot(4,1,4);
% plot(RX3, 20*log10(abs(pc_out_long_full_temp_win(1,:)))); % 输出脉冲压缩结果
% subtitle('幅值（dB）');
% xlabel('距离 (m)');
% ylabel('幅度');


%% 噪声FFT结果图
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

%% 本地子函数
% =========================================================================
function channel_phase_shifts_deg = calculate_phase_shifts(arrival_angle_deg, element_spacing_m, wavelength_m)
    arrival_angle_rad = deg2rad(arrival_angle_deg);
    delta_phi_rad = (2 * pi * element_spacing_m * sin(arrival_angle_rad)) / wavelength_m;
    channel_indices = 0:15;
    channel_phase_shifts_rad = channel_indices * delta_phi_rad;
    channel_phase_shifts_deg = rad2deg(channel_phase_shifts_rad);
end




