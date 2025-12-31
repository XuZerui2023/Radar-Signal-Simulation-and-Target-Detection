% main_simulate_echoes_with_array_v2.m
% 雷达目标信号回波模拟仿真程序
% 本程序用于对雷达目标信号回波进行模拟仿真，并产生数据来测试雷达信号处理算法
% --- 工作流程 ---
% 1. 配置雷达系统参数和天线阵列参数。
% 2. 定义一个或多个点目标，包括其距离、速度、RCS和俯仰角。
% 3. 三脉冲发射波形生成和目标信号回波模拟。
% 4. 在生成单个目标回波后，为其计算16个通道的相位偏置。
% 5. 将带有相位偏置的回波信号复制到16个通道中。
% 6. 累加所有目标在所有通道的回波。
% 7. 添加噪声并保存为与项目兼容的 .mat 文件。
% 8. 脉冲压缩和MTD处理
% 9. 对以上结果可视化分析。
%
%  修改记录
%  date       by      version   modify
%  25/08/24   XZR      v1.0     创建0 
%  25/08/24   XZR      v2.0     创建，增加了对16通道线性天线阵列的模拟
%  25/08/26   XZR      V2.1     在仿真结束后，脚本会自动对生成的数据进行脉冲压缩处理，并绘制距离-幅度图。
%  25/08/27   XZR      V2.2     增添多个绘图指令用于分析模拟生成的数据。
%  25/09/04   XZR      V2.3     统一绘图程序，将脉压

%  MTD内容提前，在脉压中加入汉明窗函数，在MTD中加入kaiser窗函数，显著提高算法性能。
clc; clear; close all;

%% 1. 用户配置区
% =========================================================================
% 描述: 这是修改仿真参数的核心区域。您可以在这里定义需要仿真的目标特性。

% --- 1.1 定义仿真的目标 ---
% targets 是一个结构体数组，每一列代表一个独立的目标。
targets = struct(...
    'Range', {1000}, ...            % 目标的初始距离 (单位: 米)
    'Velocity', {20}, ...           % 目标的速度 (单位: 米/秒, 正值表示朝向雷达，负值表示远离)
    'RCS', {2}, ...                 % 目标的雷达散射截面 (单位: 平方米), 这个值决定了目标回波的强度
    'ElevationAngle', {5} ...        % 目标的俯仰角/入射角 (单位: 度)，相对于阵面法线方向
);

% --- 1.2 仿真控制 ---
frame_index_to_save = 1; % 定义要保存的仿真数据帧的编号
output_folder = fullfile(pwd, 'simulated_data_array'); % 定义存放仿真结果的文件夹名称

%% 2. 雷达与天线阵列参数配置
% =========================================================================
% 描述: 此处定义了雷达系统和天线阵列的物理参数。
% 注意: 这些参数必须与项目中 `main_integrated_processing.m` 的配置完全一致，
%       以确保生成的仿真数据能够被后续处理流程正确解析。

fprintf('--- 开始进行雷达与天线阵列参数配置 ---\n');

% --- 雷达系统参数 ---
config.Sig_Config.c = 2.99792458e8;             % 光速 (m/s)
config.Sig_Config.fs = 25e6;                    % 信号采样率 (Hz)
config.Sig_Config.fc = 9450e6;                  % 雷达工作中心频率 (Hz)
config.Sig_Config.f0 = 0;                       % 模拟信号起始频率，给基带LFM信号引入一个频率偏移。
                                                % 但要注意匹配滤波器（即
                                                % tx_pulse）的参考波形也必须包含完全相同的f0，即发射波形和匹配滤波器波形中心都应该在f0
config.Sig_Config.prtNum = 332;                 % 每帧的脉冲数 (慢时间维度的样本数)
config.Sig_Config.point_PRT = 3404;             % 每个脉冲的采样点数 (快时间维度的样本数)
config.Sig_Config.channel_num = 16;             % 接收通道数
config.Sig_Config.prt = 232.76e-6;              % 脉冲重复时间(PRT周期) (单位: s)
config.Sig_Config.B = 20e6;                     % 信号带宽 (Hz)
config.Sig_Config.tao = [0.16e-6, 8e-6, 28e-6]; % 脉宽 [窄, 中, 长] (s)
config.Sig_Config.tao_interval = [11.4e-6, 31.8e-6, 153.4e-6]; % 三脉冲脉冲间隔（脉冲前信号为0时刻）
config.Sig_Config.point_prt_segments = [228, 723, 2453]; % 各脉冲段对应的采样波门采样点数

% --- 天线阵列参数 ---
config.Array.element_spacing = 0.0138;          % 阵元间距 (13.8mm), 单位: 米

% --- 由上述参数计算得出的其他参数 ---
config.Sig_Config.prf = 1/config.Sig_Config.prt; % 脉冲重复频率 (Hz)
config.Sig_Config.wavelength = config.Sig_Config.c / config.Sig_Config.fc; % 信号波长 (m)
config.Sig_Config.K1 = config.Sig_Config.B / config.Sig_Config.tao(1);     % 短脉冲调频斜率
config.Sig_Config.K2 = -config.Sig_Config.B / config.Sig_Config.tao(2);    % 中脉冲调频斜率
config.Sig_Config.K3 = config.Sig_Config.B / config.Sig_Config.tao(3);     % 长脉冲调频斜率
ts = 1 / config.Sig_Config.fs; % 采样时间间隔 (s)

%% 3. 生成发射参考波形 (参考 fun_MTD_produce.m)
% =========================================================================
% 描述: 此部分根据雷达参数，生成一个完整的、由三段不同类型脉冲拼接而成的发射波形。这个波形是后续生成目标回波的基础。
% --- 3.1 根据实际脉冲宽度和采样率计算信号点数 ---
num_samples_1 = round(config.Sig_Config.tao(1) * config.Sig_Config.fs); % 窄脉冲, 0.16us * 25MHz = 4 points
num_samples_2 = round(config.Sig_Config.tao(2) * config.Sig_Config.fs); % 中脉冲, 8us * 25MHz = 200 points
num_samples_3 = round(config.Sig_Config.tao(3) * config.Sig_Config.fs); % 长脉冲, 28us * 25MHz = 700 points

% --- 3.2 为每个脉冲段创建正确长度的时间向量 ---
t1 = linspace(-config.Sig_Config.tao(1)/2, config.Sig_Config.tao(1)/2, num_samples_1);
t2 = linspace(-config.Sig_Config.tao(2)/2, config.Sig_Config.tao(2)/2, num_samples_2);
t3 = linspace(-config.Sig_Config.tao(3)/2, config.Sig_Config.tao(3)/2, num_samples_3);

% --- 3.3 生成各段脉冲的波形 ---
pulse1 = sin(2*pi*t1 + pi/2); % 窄脉冲: 长度为 4
pulse2 = exp(1j*2*pi*(config.Sig_Config.f0*t2 + 0.5*config.Sig_Config.K2*(t2.^2))); % 中脉冲: 长度为 200
pulse3 = exp(1j*2*pi*(config.Sig_Config.f0*t3 + 0.5*config.Sig_Config.K3*(t3.^2))); % 长脉冲: 长度为 700

% --- 3.4 构建一个完整PRT脉冲波形，以窄脉冲上升沿为零点，依次在对应位置填充中脉冲和长脉冲 ---
num_samples_prt = config.Sig_Config.prt * config.Sig_Config.fs; % 一个完整PRT信号长度（点数）
t_prt = linspace(-config.Sig_Config.prt/2, config.Sig_Config.prt/2, num_samples_prt); % 完整PRT的时间向量
tx_pulse = zeros(1, num_samples_prt); % 构建一个完整PRT的发射信号

% 3.4.1 填充窄脉冲段
% 将长度为 4 的 pulse1 放置在整段PRT的开头
tx_pulse(1:num_samples_1) = pulse1;

% 3.4.2 填充中脉冲段
% 计算中脉冲的起始点
offset_m = (config.Sig_Config.tao(1) + config.Sig_Config.tao_interval(1)) * config.Sig_Config.fs; % 中脉冲偏移点数：在窄脉冲基础上右移窄脉冲的脉宽和窄中脉冲间的脉冲间隔
% 将长度为 200 的 pulse2 右移 offset_m 点
tx_pulse(offset_m+1 : offset_m+num_samples_2) = pulse2;

% 3. 填充长脉冲段
% 计算长脉冲的起始点
offset_l = offset_m + (config.Sig_Config.tao(2) + config.Sig_Config.tao_interval(2)) * config.Sig_Config.fs; % 长脉冲偏移点数：在中脉冲基础上右移中脉冲的脉宽和中长脉冲间的脉冲间隔
% 将长度为 700 的 pulse3 右移 offset_l 点
tx_pulse(offset_l+1 : offset_l+num_samples_3) = pulse3;

fprintf('  > 物理真实的发射波形已生成。\n');
figure;
plot(t_prt, abs(tx_pulse));  % 绘制幅值
ylim([-1.5,1.5]);
xlabel('时间');
ylabel('幅值');
title('复数信号幅值时域图');

figure;
plot(t_prt, real(tx_pulse));  % 绘制幅值
ylim([-1.5,1.5]);
xlabel('时间');
ylabel('信号实部');
title('发射信号实部波形图');

figure;
plot(t_prt, unwrap(angle(tx_pulse))); % 绘制相位，使用unwrap避免相位跳变
xlabel('时间');
ylabel('弧度');
title('发射信号相位图');

%% 4. 模拟多通道回波信号
% =========================================================================
% 程序循环遍历每一个脉冲和每一个目标，计算并合成最终包含所有目标回波、噪声和通道间相位差的16通道数据。

fprintf('--- 开始模拟16通道目标回波信号 ---\n');

% 初始化一帧的三维复数数据矩阵 (脉冲数 x 距离点数 x 通道数)
raw_iq_data = complex(zeros(config.Sig_Config.prtNum, num_samples_prt, config.Sig_Config.channel_num));

% --- 慢时间维循环 (速度维，逐个脉冲进行模拟) ---
for m = 1:config.Sig_Config.prtNum
    % m = 1; % 调试用
    % 初始化当前脉冲在所有通道的回波矩阵
    pulse_echo_all_channels = complex(zeros(num_samples_prt, config.Sig_Config.channel_num));
    
    % --- 遍历所有定义的目标 ---
    for k = 1:length(targets)
        % k = 1; % 调试用
        % 1. 计算当前脉冲时刻，目标的瞬时距离
        %range = targets(k).Range + targets(k).Velocity * (m-1) * config.Sig_Config.prt;  % 模拟目标在匀速直线运动
        range = targets(k).Range;  % 模拟目标在某点“静止”，但具有速度。
        % 2. 计算来回时延 (Round-trip delay) 并转换为采样点数
        delay = 2 * range / config.Sig_Config.c;
        delay_samples = round(delay / ts);
        
        % 3. 计算多普勒频移和当前脉冲的多普勒相移
        doppler_freq = 2 * targets(k).Velocity / config.Sig_Config.wavelength;
        doppler_phase_shift = exp(1j * 2 * pi * doppler_freq * (m-1) * config.Sig_Config.prt); % 目标速度引起的多普勒效应，对整个回波信号施加一个相位偏移，相当于相位旋转因子
         
        % 4. 根据简化的雷达方程计算回波幅度
        % 原始的雷达方程会导致信号幅度过小，完全被噪声淹没。为了验证，在此添加一个较大的增益，以确保目标可被检测。
        amplitude_gain = 1e8; % 引入信号增益，提高信噪比 1e2时检测会被干扰，1e3及以上有效
        lambda_sq = config.Sig_Config.wavelength^2;
        amplitude = amplitude_gain * sqrt(targets(k).RCS * lambda_sq) / (range^2 * (4*pi)^(3/2));
        amplitude = 1; % 调试用，信号幅度直接设为 1。
        
        % 5. 生成无方向性的基础回波信号 (所有通道接收到的信号主体是相同的)
        target_echo_base = zeros(1, num_samples_prt);
        if (delay_samples > 0) && (delay_samples < num_samples_prt)
            % 通过截取和移位发射脉冲来模拟延迟后的回波
            len_echo = min(length(tx_pulse), num_samples_prt - delay_samples);
            target_echo_base(delay_samples+1 : delay_samples+len_echo) = tx_pulse(1:len_echo);
        end
        % 在目标回波信号上应用幅度和多普勒相移
        target_echo_base = amplitude * target_echo_base * doppler_phase_shift;
        
        % 6. 计算并施加16个通道的相位偏移
        % 调用辅助函数，根据目标的俯仰角计算16个通道的相位差
        phase_shifts_deg = calculate_phase_shifts(targets(k).ElevationAngle, config.Array.element_spacing, config.Sig_Config.wavelength);
        phase_shifts_rad = deg2rad(phase_shifts_deg);
        channel_phasors = exp(1j * phase_shifts_rad); % 将相位差转换为复数相量 e^(j*phi)
        
        % 7. 使用外积 (outer product) 将相位偏置高效地应用到所有通道
        %    target_echo_base.' 是一个列向量 (point_PRT x 1)
        %    channel_phasors 是一个行向量 (1 x 16)
        %    结果是一个矩阵 (point_PRT x 16)，每一列都乘以了对应的复数相量，从而施加了相位
        target_echo_multichannel = target_echo_base.' * channel_phasors;
        
        % 8. 将当前目标的多通道回波累加到总回波中
        pulse_echo_all_channels = pulse_echo_all_channels + target_echo_multichannel;
    end
    
    % 将当前脉冲的所有通道回波存入最终的数据帧矩阵中
    
    % raw_iq_data(m, :, :) = pulse_echo_all_channels.'; % 需要转置以匹配 (prtNum, point_PRT, channel_num) 的维度格式
    raw_iq_data(m, :, :) = reshape(pulse_echo_all_channels, [1, num_samples_prt, config.Sig_Config.channel_num]);
end

% --- 添加高斯白噪声到所有通道中 ---
noise_power = 1e-3; % 可调噪声功率
noise = sqrt(noise_power/2) * (randn(size(raw_iq_data)) + 1j * randn(size(raw_iq_data)));
% raw_iq_data = raw_iq_data*10;
raw_iq_data_noise = raw_iq_data + noise;

% --- 模拟伺服角度 ---
% 假设雷达在采集一帧数据的过程中匀速扫描
start_angle = 10; % 假设起始方位角为10度
time_per_frame = config.Sig_Config.prt * config.Sig_Config.prtNum;
antangle_per_Time = 72; % 伺服角每秒转72°
scan_rate_per_frame = time_per_frame * antangle_per_Time;
scan_rate_per_pulse = (scan_rate_per_frame / config.Sig_Config.prtNum); % 假设每帧扫描0.5度，计算每个脉冲对应的角度增量
servo_angle = start_angle + (0:config.Sig_Config.prtNum-1) * scan_rate_per_pulse;
% 转换为雷达数据格式中定义的原始整数值 (乘以0.01后为度)
% servo_angle = servo_angle / 0.01; 

% 目标回波信号时域波形
figure;
hold on;
plot((0:length(tx_pulse)-1), real(raw_iq_data_noise(1,:,1)));
% plot((0:length(tx_pulse)-1)*ts, real(raw_iq_data_noise(1,:,2)));
grid on;
title('回波信号时域波形图 (实部)');
xlabel('采样点数');
ylabel('幅度');
hold off;

%% 5. 可视化检查
% =========================================================================
% 描述: 提供两个简单的图表，用于快速检查仿真结果是否符合预期。

figure('Name', '多通道仿真回波信号检查');
% 图1: 显示通道1的回波幅度图 (距离-脉冲图)
subplot(1, 2, 1);
plot(db(abs(raw_iq_data(:,:,1))));
xlabel('距离单元 (快时间)');
ylabel('脉冲数 (慢时间)');
title('通道 1 幅度');
colorbar;
% 图2: 显示第10个脉冲在所有16个通道上的相位分布
% 正常情况下，应该能看到相位随通道号线性变化
subplot(1, 2, 2);
imagesc(angle(squeeze(raw_iq_data(10,:,:)))); % squeeze移除单维度
xlabel('通道号');

ylabel('距离单元');
title('单脉冲相位分布');
colorbar;

%% 6. 分割三脉冲，拼接成3404点的信号数据
% 窄脉冲采样点数 228（采样取83 - 312点）
raw_iq_data_noise_sample_s = raw_iq_data_noise(:, (83:310), :);

% 中脉冲采样点数 723（采样取311 - 1033点）
raw_iq_data_noise_sample_m = raw_iq_data_noise(:, (311:1033), :);

% 长脉冲采样点数 2453（采样取1034 - 3486点）
raw_iq_data_noise_sample_l = raw_iq_data_noise(:, (1034:3486), :);

% 拼接三脉冲 
raw_iq_data_noise_sample = [raw_iq_data_noise_sample_s, raw_iq_data_noise_sample_m, raw_iq_data_noise_sample_l];

% 目标回波信号时域波形
figure;
hold on;
plot((1:config.Sig_Config.point_PRT), real(raw_iq_data_noise_sample(1,:,1)));
% plot((0:length(tx_pulse)-1)*ts, real(raw_iq_data_noise(1,:,2)));
grid on;
title('拼接后回波信号时域波形图 (实部)');
xlabel('采样点数');
ylabel('幅度');
hold off


%% 7. 保存仿真数据
% =========================================================================
% 描述: 将最终生成的仿真数据保存为 .mat 文件，以便后续处理脚本调用。

if ~exist(output_folder, 'dir'), mkdir(output_folder); end
output_filename = fullfile(output_folder, sprintf('frame_sim_array_%d.mat', frame_index_to_save));

% 保存的变量名 'raw_iq_data' 和 'servo_angle' 与真实数据文件中的变量名保持一致
save(output_filename, 'raw_iq_data_noise_sample', 'servo_angle');

fprintf('--- 仿真完成 ---\n');
fprintf('16通道仿真数据已成功保存到: %s\n', output_filename);

%% 8. 脉冲压缩
% =========================================================================
% 描述: 此部分通过对仿真数据进行脉冲压缩，来验证目标是否被正确地放置在了预设的距离上。
% 对所有脉冲进行脉冲压缩（快时间维，脉冲压缩），提取距离信息
fprintf('--- 开始进行仿真验证 ---\n');

% --- 7.1 提取单脉冲回波 ---
% 提取单通道完整脉冲回波用于分析，从16个通道中选择第6个通道。
echo_for_pc = squeeze(raw_iq_data_noise_sample(:, :, 6)); 

% --- 7.2 执行脉冲压缩 (匹配滤波) ---
% 使用与 fun_pulse_compression.m 相同的原理
% 1. 生成匹配滤波器 (发射信号的时间反转共轭)
matched_filter = conj(fliplr(tx_pulse));

% 2. 距离维加窗 (用于抑制距离旁瓣)
window_pc = hamming(length(matched_filter))';
matched_filter_windowed = matched_filter .* window_pc;

% 4. 使用FFT进行快速卷积
pc_result = ifft(fft(echo_for_pc, [], 2) .* repmat(fft(matched_filter(1:3404)), 332, 1), [], 2);                     % 不加窗脉压结果
pc_result_windowed = ifft( fft(echo_for_pc, [], 2) .* repmat(fft(matched_filter_windowed(1:3404)), 332, 1), [], 2 ); % 加窗脉压结果

%% 9. MTD
% --- 速度维多普勒FFT (已加入窗函数) ---
% 速度维加窗 (用于抑制速度旁瓣/频谱泄漏)
window_mtd = kaiser(config.Sig_Config.prtNum, 8);

% 不加窗做MTD（用于对比）
rdm = fftshift(fft(pc_result, [], 1), 1); 

% 对每个距离单元的脉冲串数据加窗
pc_all_pulses_windowed = pc_result_windowed .* repmat(window_mtd, 1, size(pc_result, 2));

% 进行加窗MTD处理
rdm_windowed = fftshift(fft(pc_all_pulses_windowed, [], 1), 1);


%% 可视化

% --- 9.2 接收与处理结果分析 ---
% 子图4: 脉冲压缩结果 (距离剖面)
figure;
subplot(2,1,1);
deltaR = config.Sig_Config.c / (2 * config.Sig_Config.fs);
range_axis_pc = (0:length(pc_result)-1) * deltaR;
plot(range_axis_pc, abs(pc_result));
grid on; 
hold on;
for i = 1:length(targets)
    xline(targets(i).Range, '--r', sprintf('目标 %d', i));
end
title('脉冲压缩后的距离-幅度图（未加窗）');
xlabel('距离 (m)');
ylabel('幅度');
xlim([0, max([targets.Range])*1.5]);
hold off;

subplot(2,1,2);
range_axis_pc_windowed = (0:length(pc_result_windowed)-1) * deltaR;
plot(range_axis_pc_windowed, abs(pc_result_windowed));
grid on; 
hold on;
for i = 1:length(targets)
    xline(targets(i).Range, '--r', sprintf('目标 %d', i));
end
title('脉冲压缩后的距离-幅度图（加窗）');
xlabel('距离 (m)');
ylabel('幅度');
xlim([0, max([targets.Range])*1.5]);
hold off;

% 子图5: 距离-多普勒图 (RDM)
figure;
% 绘制未加窗RDM
subplot(2,1,1);
% 创建速度轴
velocity_axis = linspace(-config.Sig_Config.prf/2, config.Sig_Config.prf/2, config.Sig_Config.prtNum) * config.Sig_Config.wavelength / 2;
imagesc(range_axis_pc, velocity_axis, 20*log10(abs(rdm)));
title('距离-多普勒图 (RDM)（未加窗）');
xlabel('距离 (m)');
ylabel('速度 (m/s)');
colorbar;
axis xy; % 将坐标原点置于左下角

% 绘制加窗后RDM
subplot(2,1,2);
velocity_axis = linspace(-config.Sig_Config.prf/2, config.Sig_Config.prf/2, config.Sig_Config.prtNum) * config.Sig_Config.wavelength / 2;
imagesc(range_axis_pc_windowed, velocity_axis, 20*log10(abs(rdm_windowed)));
title('距离-多普勒图 (RDM)（加窗）');
xlabel('距离 (m)');
ylabel('速度 (m/s)');
colorbar;
axis xy; % 将坐标原点置于左下角







% %% 9. 仿真信号可视化分析
% % =========================================================================
% fprintf('--- 开始进行仿真信号可视化分析 ---\n');
% 
% % --- 9.1 发射信号分析 ---
% % 子图1: 发射信号时域波形 (实部和虚部)
% figure;
% subplot(2,1,1);
% hold on;
% plot((0:length(tx_pulse)-1)*ts, real(raw_iq_data_noise(1,:,1)));
% % plot((0:length(tx_pulse)-1)*ts, real(raw_iq_data(1,:,2)));
% grid on;
% title('发射信号时域波形 (实部)');
% xlabel('时间 (s)');
% ylabel('幅度');
% hold off;
% 
% subplot(2,1,2);
% hold on;
% plot((0:length(tx_pulse)-1)*ts, imag(raw_iq_data_noise(1,:,1)));
% % plot((0:length(tx_pulse)-1)*ts, imag(raw_iq_data(1,:,2)));
% grid on;
% title('发射信号时域波形 (虚部)');
% xlabel('时间 (s)');
% ylabel('幅度');
% hold off;
% % 子图2: 发射信号的频谱
% figure;
% N_fft = 2^nextpow2(length(tx_pulse));
% freq_axis = (-N_fft/2:N_fft/2-1)*(config.Sig_Config.fs/N_fft)/1e6;
% tx_spectrum = fftshift(fft(tx_pulse, N_fft));
% plot(freq_axis, 20*log10(abs(tx_spectrum)/max(abs(tx_spectrum))));
% grid on;
% title('发射信号频谱');
% xlabel('频率 (MHz)');
% ylabel('归一化幅度 (dB)');
% ylim([-100,10]);
% 
% % 子图3: 发射信号时频谱 (Spectrogram)
% figure;
% % 分析中频调脉冲段(pulse2)以获得清晰的图像
% spectrogram(pulse2, kaiser(128, 5), 120, 1024, config.Sig_Config.fs, 'yaxis');
% title('中频调脉冲(pulse2)的时频谱');
% 
% 
% % --- 9.2 接收与处理结果分析 ---
% % 子图4: 脉冲压缩结果 (距离剖面)
% figure;
% subplot(2,1,1);
% deltaR = config.Sig_Config.c / (2 * config.Sig_Config.fs);
% range_axis_pc = (0:length(pc_result)-1) * deltaR;
% plot(range_axis_pc, abs(pc_result));
% grid on; 
% hold on;
% for i = 1:length(targets)
%     xline(targets(i).Range, '--r', sprintf('目标 %d', i));
% end
% title('脉冲压缩后的距离-幅度图（未加窗）');
% xlabel('距离 (m)');
% ylabel('幅度');
% xlim([0, max([targets.Range])*1.5]);
% hold off;
% 
% subplot(2,1,2);
% range_axis_pc_windowed = (0:length(pc_result_windowed)-1) * deltaR;
% plot(range_axis_pc_windowed, abs(pc_result_windowed));
% grid on; 
% hold on;
% for i = 1:length(targets)
%     xline(targets(i).Range, '--r', sprintf('目标 %d', i));
% end
% title('脉冲压缩后的距离-幅度图（加窗）');
% xlabel('距离 (m)');
% ylabel('幅度');
% xlim([0, max([targets.Range])*1.5]);
% hold off;
% 
% % 子图5: 距离-多普勒图 (RDM)
% figure;
% % 绘制未加窗RDM
% subplot(2,1,1);
% % 创建速度轴
% velocity_axis = linspace(-config.Sig_Config.prf/2, config.Sig_Config.prf/2, config.Sig_Config.prtNum) * config.Sig_Config.wavelength / 2;
% imagesc(range_axis_pc, velocity_axis, 20*log10(abs(rdm)));
% title('距离-多普勒图 (RDM)（未加窗）');
% xlabel('距离 (m)');
% ylabel('速度 (m/s)');
% colorbar;
% axis xy; % 将坐标原点置于左下角
% 
% % 绘制加窗后RDM
% subplot(2,1,2);
% velocity_axis = linspace(-config.Sig_Config.prf/2, config.Sig_Config.prf/2, config.Sig_Config.prtNum) * config.Sig_Config.wavelength / 2;
% imagesc(range_axis_pc_windowed, velocity_axis, 20*log10(abs(rdm_windowed)));
% title('距离-多普勒图 (RDM)（加窗）');
% xlabel('距离 (m)');
% ylabel('速度 (m/s)');
% colorbar;
% axis xy; % 将坐标原点置于左下角
% 
% % 速度维


%% 本地子函数
% =========================================================================
% 该函数用于引入16通道相位差
function channel_phase_shifts_deg = calculate_phase_shifts(arrival_angle_deg, element_spacing_m, wavelength_m)
    % 将入射角从度转换为弧度
    arrival_angle_rad = deg2rad(arrival_angle_deg);
    
    % 计算相邻阵元间的相位差 (弧度)
    delta_phi_rad = (2 * pi * element_spacing_m * sin(arrival_angle_rad)) / wavelength_m;
    
    % 计算每个通道相对于第一个通道(索引为0)的总相位差
    channel_indices = 0:15;
    channel_phase_shifts_rad = channel_indices * delta_phi_rad;
    
    % 将结果从弧度转换为度
    channel_phase_shifts_deg = rad2deg(channel_phase_shifts_rad);
end
 