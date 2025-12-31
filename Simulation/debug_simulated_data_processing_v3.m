clc,clear;close all;

%% 1 数据加载与参数定义
% =========================================================================
fprintf('--- STAGE 1: 正在加载数据与系数 ---\n');
% --- 目标设定 --- 
targets = struct(...
    'Range', {1000}, ...            % 目标的初始距离 (单位: 米)
    'Velocity', {20}, ...           % 目标的速度 (单位: 米/秒, 正值表示朝向雷达，负值表示远离)
    'RCS', {2}, ...                 % 目标的雷达散射截面 (单位: 平方米), 这个值决定了目标回波的强度
    'ElevationAngle', {5} ...        % 目标的俯仰角/入射角 (单位: 度)，相对于阵面法线方向
);

% --- 1.1 加载仿真数据 ---
base_path = 'C:\Users\a\Desktop\9.8 问题\Simulation';
sim_data_folder = fullfile(base_path, 'simulated_data_array');
sim_frame_file = 'frame_sim_array_1.mat';
sim_data_path = fullfile(sim_data_folder, sim_frame_file);
if ~exist(sim_data_path, 'file'), error('仿真数据文件不存在!'); end
sim_data = load(sim_data_path);
raw_iq_data = sim_data.raw_iq_data_noise_sample;
angle = sim_data.servo_angle;

% --- 1.2 加载DBF系数 ---
dbf_coef_path = fullfile(base_path, 'X8数据采集250522_DBFcoef.csv')
if ~exist(dbf_coef_path, 'file'), error('DBF系数文件不存在!'); end
DBF_coeffs_data = readmatrix(dbf_coef_path);
DBF_coeffs_data_I = double(DBF_coeffs_data(:, 1:2:end));
DBF_coeffs_data_Q = double(DBF_coeffs_data(:, 2:2:end));
DBF_coeffs_data_C = DBF_coeffs_data_I + 1j * DBF_coeffs_data_Q;

% --- 1.3 加载测角K值 ---
k_value_path = fullfile(base_path, 'R9-DMX3-2024001_Angle_k.csv')
if ~exist(k_value_path, 'file'), error('DBF系数文件不存在!'); end
k_value_data = readmatrix(k_value_path);

% --- 1.4 雷达系统参数 ---
% 伺服角（方位角）修正系数
config.corrected.North = - 242;                % 雷达指北角（见雷达系统设置文件 SysSet.ini.bak）
config.corrected.FixAngle = 35;              % 雷达固定角（见初始化参数文件 InitPara.ini）
% 俯仰角修正系数 
config.corrected.ELeAngleSettingValue = -10; % 雷达俯仰设置值（见雷达系统设置文件 SysSet.ini.bak）

% 雷达系统基础参数
MHz = 1e6; % 定义MHz单位
config.Sig_Config.c = 2.99792458e8;            % 光速
config.Sig_Config.pi = pi;                     % 圆周率
config.Sig_Config.fs = 25e6;                   % 采样率 (Hz)
config.Sig_Config.fc = 9450e6;                 % 中心频率 (Hz)
config.Sig_Config.timer_freq = 200e6;          % 时标计数频率 
config.Sig_Config.prtNum = 332;                % 定义每帧信号的脉冲数，每帧信号包含 332 个脉冲
config.Sig_Config.point_PRT = 3404;            % 定义每个PRT中的采样点数（距离单元）
config.Sig_Config.channel_num = 16;            % 通道数（阵元数目）
config.Sig_Config.beam_num = 13;               % 波束数
config.Sig_Config.prt = 232.76e-6;             % 脉冲重复时间 (s)
config.Sig_Config.prf = 1/config.Sig_Config.prt;
config.Sig_Config.B = 20e6;                    % 带宽 (Hz)
config.Sig_Config.bytesFrameHead = 64;         % 每个PRT帧头字节数
config.Sig_Config.bytesFrameEnd = 64;          % 每个PRT帧尾字节数
config.Sig_Config.bytesFrameRealtime = 128;    % 实时参数的字节数 
config.Sig_Config.tao = [0.16e-6, 8e-6, 28e-6];       % 脉宽 [窄, 中, 长]
config.Sig_Config.point_prt = [3404, 228, 723, 2453]; % 采集点数 [总采集点数，窄脉冲采集点数，中脉冲采集点数，长脉冲采集点数]   
config.Sig_Config.wavelength = config.Sig_Config.c / config.Sig_Config.fc;   % 信号波长
config.Sig_Config.deltaR = config.Sig_Config.c / (2 * config.Sig_Config.fs); % 距离分辨率由采样率决定
config.Sig_Config.tao1 = config.Sig_Config.tao(1);    % 窄脉宽 
config.Sig_Config.tao2 = config.Sig_Config.tao(2);    % 中脉宽
config.Sig_Config.tao3 = config.Sig_Config.tao(3);    % 长脉宽
config.Sig_Config.K1   = config.Sig_Config.B/config.Sig_Config.tao1;   % 短脉冲调频斜率
config.Sig_Config.K2   = -config.Sig_Config.B/config.Sig_Config.tao2;  % 中脉冲调频斜率（负）
config.Sig_Config.K3   = config.Sig_Config.B/config.Sig_Config.tao3;   % 长脉冲调频斜率

config.Sig_Config.beam_angles_deg = [-12.5, -7.5, -2.5, 2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5]; % 13个波束的标称俯仰角 (Nominal Elevation Angles)
config.Sig_Config.beam_angles_deg = config.Sig_Config.beam_angles_deg - config.corrected.ELeAngleSettingValue;

config.Sig_Config.debug.show_PC = 0;              % 脉冲压缩结果显示
config.Sig_Config.debug.show_FFT = 0;             % 速度维显示
config.Sig_Config.debug.graph = 0;                % 是否画图 

% --- 1.5 MTD处理参数 ---
config.mtd.win_size = 4;                          % MTD 窗口切片数  
config.mtd.c = config.Sig_Config.c;               % 光速
config.mtd.pi = config.Sig_Config.pi;             % 圆周率
config.mtd.prtNum = config.Sig_Config.prtNum;     % 每帧信号的脉冲数
config.mtd.fs = config.Sig_Config.fs;             % 采样频率 (Hz)
config.mtd.fc = config.Sig_Config.fc;             % 中心频率 (Hz)
config.mtd.beam_num = config.Sig_Config.beam_num; % 雷达总波束数量
config.mtd.prt = config.Sig_Config.prt;           % 脉冲重复时间 (s)
config.mtd.B = config.Sig_Config.B;               % 带宽 (Hz)
config.mtd.tao = config.Sig_Config.tao;           % 脉宽 [窄, 中, 长]
config.mtd.point_prt = config.Sig_Config.point_prt;   % 采集点数 [总采集点数，窄脉冲采集点数，中脉冲采集点数，长脉冲采集点数]   
config.mtd.prf = 1 / config.mtd.prt;                  
config.mtd.wavelength = config.Sig_Config.wavelength; % 信号波长
config.mtd.deltaR = config.Sig_Config.deltaR;         % 距离分辨率由采样率决定
config.mtd.tao1 = config.Sig_Config.tao1;             % 窄脉宽 
config.mtd.tao2 = config.Sig_Config.tao2;             % 中脉宽
config.mtd.tao3 = config.Sig_Config.tao3;             % 长脉宽
config.mtd.K1   = config.Sig_Config.B/config.Sig_Config.tao1;   % 短脉冲调频斜率
config.mtd.K2   = -config.Sig_Config.B/config.Sig_Config.tao2;  % 中脉冲调频斜率（负）
config.mtd.K3   = config.Sig_Config.B/config.Sig_Config.tao3;   % 长脉冲调频斜率

% --- 1.6 CFAR处理参数 ---
% --- CFAR核心参数 ---
config.cfar.T_CFAR = 3;                    % 恒虚警标称化因子
config.cfar.MTD_V = 3;                     % 杂波区速度范围，速度在 -3 m/s 到 +3 m/s 范围内的区域都当作是地杂波区域，在CFAR检测中忽略掉。

% --- 速度维参数 ---
% config.cfar.refCells_V = 5;                % 速度维 参考单元数
config.cfar.refCells_V = 4;                % 速度维 参考单元数
% config.cfar.saveCells_V = 14;              % 速度维 保护单元数
config.cfar.saveCells_V = 7;              % 速度维 保护单元数
config.cfar.T_CFAR_V = config.cfar.T_CFAR; % 速度维恒虚警标称化因子
config.cfar.CFARmethod_V = 0;              % 0--选大；1--选小

% --- 距离维参数 ---
config.cfar.rCFARDetect_Flag = 1;          % 距离维CFAR检测操作标志。 0-否； 1-是
config.cfar.refCells_R = 5;                % 距离维 参考单元数
config.cfar.saveCells_R = 14;              % 距离维 保护单元数
config.cfar.T_CFAR_R = config.cfar.T_CFAR; % 距离维恒虚警标称化因子7,越低，门限越低，虚警率越高
config.cfar.CFARmethod_R = 0;              % 0--选大；1--选小

% --- 计算杂波区对应的速度单元数 ---
config.cfar.deltaDoppler = config.Sig_Config.prf/config.Sig_Config.prtNum;      % 计算多普勒频率分辨率
config.cfar.deltaV = config.Sig_Config.wavelength*config.cfar.deltaDoppler/2;   % 计算速度分辨率：计算出了一个频率单元（deltaDoppler）等效于多少米/秒（m/s）的速度。这个 deltaV 就是雷达能分辨的最小速度差。
config.cfar.MTD_0v_num = floor(config.cfar.MTD_V/config.cfar.deltaV);           % 计算杂波区的宽度（以单元数计），在进行CFAR检测时，需要以零速为中心，向两侧各跳过 MTD_0v_num 个速度单元，以避开强大的地杂波对噪声估计的干扰。

% --- 画图参数（绘图坐标轴参数）---
config.cfar.graph = 0;
config.cfar.prtNum = config.Sig_Config.prtNum;            % 每帧信号prt数量
config.cfar.point_prt = config.Sig_Config.point_prt(1);   % 3个脉冲的PRT总采集点数
config.cfar.R_point = 6;                                  % 每个距离单元长度（两点间距6m）
config.cfar.r_axis = 0 : config.cfar.R_point: config.cfar.point_prt*config.cfar.R_point-config.cfar.R_point;  % 距离轴
config.cfar.fd = linspace(-config.Sig_Config.prf/2, config.Sig_Config.prf/2, config.Sig_Config.prtNum);           
config.cfar.v_axis = config.cfar.fd*config.Sig_Config.wavelength/2;                                           % 速度轴
% v_axis = v_axis(691:845);
% (其他如refCells, saveCells等参数也应在此定义)

% --- 目标参数测量参数与插值参数 --- 
config.interpolation.extra_dots = 2;            % 插值时在峰值两侧各取几个点
config.interpolation.range_interp_times = 8;    % 距离维插值倍数
config.interpolation.velocity_interp_times = 4; % 速度维插值倍数

%% 2 信号处理流程
% =========================================================================
% --- 2.1 DBF处理 ---
fprintf('--- STAGE 2: 正在执行DBF处理 ---\n');
raw_iq_16ch = raw_iq_data;
iq_data_13beam = complex(zeros(config.Sig_Config.prtNum, config.Sig_Config.point_PRT, config.Sig_Config.beam_num));
for i = 1:config.Sig_Config.prtNum
    single_pulse_16ch = squeeze(raw_iq_16ch(i, :, :));
    single_pulse_13beam = single_pulse_16ch * DBF_coeffs_data_C';
    iq_data_13beam(i, :, :) = single_pulse_13beam; % DBF后13波束三维信号矩阵 332*3404*13
end

% --- 可视化信号数据分析 ---
% 目标回波信号时域波形
figure;
hold on;
plot((1:config.Sig_Config.point_PRT), real(iq_data_13beam(1,:,1)));
% plot((0:length(tx_pulse)-1)*ts, real(raw_iq_data_noise(1,:,2)));
grid on;
title('拼接后回波信号时域波形图 (实部)');
xlabel('采样点数');
ylabel('幅度');
hold off



% --- 可视化DBF结果分析 ---
% 选择一个波束进行显示
% 预设目标在 5度，我们检查第2、3波束 (标称7.5度)
beam_to_check_dbf = 3; 
figure('Name', '调试 - DBF后波束功率');
imagesc(db(abs(squeeze(iq_data_13beam(:,:,beam_to_check_dbf)))));
title(sprintf('DBF后波束 #%d 的功率 (dB)', beam_to_check_dbf));
xlabel('距离单元 (快时间)');
ylabel('脉冲数 (慢时间)');
colorbar;
% fprintf('  > DBF可视化完成。请检查Figure 1中是否有明显的水平能量条带。\n');

% --- 2.2 执行脉冲压缩测距和MTD测速 ---
% % 窄脉冲段（采样点数 228）
% raw_iq_data_noise_sample_s = iq_data_13beam(:, (1:228), :);
% % 中脉冲段（采样点数 723）
% raw_iq_data_noise_sample_m = iq_data_13beam(:, (229:951), :);
% % 长脉冲段（采样点数 2453）
% raw_iq_data_noise_sample_l = iq_data_13beam(:, (952:3404), :);

% 调用MTD执行函数进行脉冲压缩和MTD处理
[MTD_results, PC_results] = process_stage2_mtd(iq_data_13beam, angle, config);

% --- 可视化分析脉冲压缩和MTD过程结果 --- 
% 脉冲压缩
% 选择可视化的波束号
beam_to_plot_pc = 7;
% 从三维数据中提取指定波束的二维数据
pc_data_2d = PC_results(:, :, beam_to_plot_pc);
% 2. 数据处理与绘图
% 通常我们关心的是信号的幅度，所以取绝对值
% 为了能同时看清强信号和弱信号，通常会转换为dB单位
pc_magnitude_db = 20 * log10(abs(pc_data_2d) + eps); % eps防止log10(0)
figure;
imagesc(pc_magnitude_db);
title(['脉冲压缩结果图 (波束 ' num2str(beam_to_plot_pc) ')'], 'FontSize', 14);
xlabel('距离单元 (Range Bin)', 'FontSize', 12);
ylabel('脉冲序号 (Pulse Index)', 'FontSize', 12);
colorbar; % 显示颜色条
axis xy; % 将坐标原点(1,1)设置在左下角，这是常规做法

figure;
plot(abs(PC_results(5,:,beam_to_plot_pc)))

figure;
plot(20*log10(abs(PC_results(5,:,beam_to_plot_pc))));

figure;
deltaR = config.Sig_Config.c / (2 * config.Sig_Config.fs);
range_axis_pc = (0:length(PC_results)-1) * deltaR;
plot(range_axis_pc, abs(PC_results(5,:,beam_to_plot_pc)));
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

% --- 2.2 CFAR处理 ---




































