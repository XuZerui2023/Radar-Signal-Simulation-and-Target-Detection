% debug_simulated_data_processing_v2.m
%
% v2.0 (CFAR门限验证版):
% - 【核心功能】修改了CFAR相关函数，使其能够返回计算出的动态检测门限矩阵。
% - 【新增可视化】创建了一个新的 Figure 4，用于将信号强度与CFAR门限进行直接对比。
%
% --- 工作流程 ---
% 1. 加载仿真数据和DBF系数。
% 2. 执行DBF处理并可视化波束功率 (Figure 1)。
% 3. 执行MTD处理并可视化RDM图 (Figure 2)。
% 4. 执行CFAR检测，并同时获取检测结果和检测门限。
% 5. 可视化CFAR的输入和输出 (Figure 3)。
% 6. 【新增】可视化最强目标的信号剖面及其对应的CFAR门限 (Figure 4)。
% 修改记录
%  date       by      version   modify
%  25/09/01   XZR      v2.0     加入CFAR目标信号剖面图分析并返回其计算的门限
clc; clear; close all;

%% 1. 用户配置区
% =========================================================================
% --- 1.1 数据与算法参数 ---
sim_data_folder = 'C:\Users\a\Desktop\9.8 问题\Simulation\simulated_data_array';
sim_frame_file = 'frame_sim_array_0.mat';

amplitude_gain = 1e8; % 信号增益，如果RDM图中目标不明显，请增大此值

% --- 1.2 DBF系数路径 ---
% 注意: 请确保此路径指向您项目中的DBF系数文件
dbf_coef_path = 'C:\Users\a\Desktop\9.8 问题\Simulation\X8数据采集250522_DBFcoef.csv';

% --- 1.3 雷达与CFAR参数 (与main_integrated_processing.m保持一致) ---
% (此处省略了大部分参数定义，直接使用您项目中的config结构体)
% 为了独立运行，我们在这里重新定义必要的参数
config.base_path = 'C:\Users\a\Desktop\9.8 问题\Simulation';
config.DBF_coef_path = fullfile(config.base_path, 'X8数据采集250522_DBFcoef.csv');          % 用于 DDC 数据波束成形的 DBF 系数，先处理成CSV（逗号分隔值）文件
config.angle_k_path = fullfile(config.base_path,'R9-DMX3-2024001_Angle_k.csv');             % 用于和差比幅的K值矩阵

% --- 0.0 实际测量时一些修正参数 ---
% 伺服角（方位角）修正系数
config.corrected.North = - 242;                % 雷达指北角（见雷达系统设置文件 SysSet.ini.bak）
config.corrected.FixAngle = 35;              % 雷达固定角（见初始化参数文件 InitPara.ini）
% 俯仰角修正系数 
config.corrected.ELeAngleSettingValue = -10; % 雷达俯仰设置值（见雷达系统设置文件 SysSet.ini.bak）

% --- 1.4 雷达系统基础参数 ---
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


% --- 仿真目标定义 (用于绘图标记) ---
targets = struct(...
    'Range', {1500, 2500, 4000}, ...
    'Velocity', {20, -15, 10}, ...
    'RCS', {10, 5, 8}, ...
    'ElevationAngle', {-5, 8.2, 15} ...
);

%% 2. 加载数据与系数
% =========================================================================
fprintf('--- STAGE 1: 正在加载数据与系数 ---\n');
% --- 加载仿真数据 ---
sim_data_path = fullfile(sim_data_folder, sim_frame_file);
if ~exist(sim_data_path, 'file'), error('仿真数据文件不存在!'); end
sim_data = load(sim_data_path);
raw_iq_data = sim_data.raw_iq_data_noise_sample;

% --- 加载DBF系数 ---
if ~exist(dbf_coef_path, 'file'), error('DBF系数文件不存在!'); end
DBF_coeffs_data = readmatrix(dbf_coef_path);
DBF_coeffs_data_I = double(DBF_coeffs_data(:, 1:2:end));
DBF_coeffs_data_Q = double(DBF_coeffs_data(:, 2:2:end));
DBF_coeffs_data_C = DBF_coeffs_data_I + 1j * DBF_coeffs_data_Q;

%% 3. 逐级信号处理与可视化
% =========================================================================
% --- 步骤1: DBF处理 ---
fprintf('--- STAGE 2: 正在执行DBF处理 ---\n');
% 选取一个有代表性的波束进行分析，例如第5个
beam_to_check_dbf = 5; 
dbf_output = zeros(config.Sig_Config.prtNum, config.Sig_Config.point_PRT);
for m = 1:config.Sig_Config.prtNum
    % 提取当前脉冲的16通道数据
    single_pulse_16ch = squeeze(raw_iq_data(m, :, :));
    % 执行DBF
    dbf_output(m, :) = (single_pulse_16ch * DBF_coeffs_data_C(beam_to_check_dbf, :).').';
end

% --- Figure 1: 可视化DBF结果 ---
figure('Name', 'Figure 1: 调试 - DBF后波束功率');
imagesc(1:config.Sig_Config.point_PRT, 1:config.Sig_Config.prtNum, db(abs(dbf_output)));
xlabel('距离单元 (快时间)');
ylabel('脉冲数 (慢时间)');
title(sprintf('DBF后波束 #%d 的功率 (dB)', beam_to_check_dbf));
colorbar;

% --- 步骤2: MTD处理 ---
fprintf('--- STAGE 3: 正在执行MTD处理 ---\n');
% 为了简化，我们直接对DBF处理后的单波束数据进行MTD
% 注意：这与完整流程中的拼接窗口方式略有不同，但足以用于验证
[mtd_rdm, ~] = fun_MTD_produce(dbf_output, config);

% --- Figure 2: 可视化MTD结果 ---
figure('Name', 'Figure 2: 调试 - MTD后RDM');
range_axis = (0:size(mtd_rdm, 2)-1) * (config.Sig_Config.c / (2*config.Sig_Config.fs));
velocity_axis = linspace(-config.Sig_Config.prf/2, config.Sig_Config.prf/2, size(mtd_rdm, 1)) * config.Sig_Config.wavelength / 2;
imagesc(range_axis, velocity_axis, db(abs(mtd_rdm)));
title(sprintf('和波束对 #%d 的距离-多普勒图 (RDM)', beam_to_check_dbf));
xlabel('距离 (m)');
ylabel('速度 (m/s)');
colorbar;
axis xy;

% --- 步骤3: CFAR检测 ---
fprintf('--- STAGE 4: 正在执行CFAR检测 ---\n');
cfar_input = abs(mtd_rdm);
% 【关键】调用修正后的CFAR函数，同时接收检测结果和门限矩阵
[cfar_output, cfar_threshold] = local_execute_cfar(cfar_input, config.cfar, config);  

% --- Figure 3: 可视化CFAR结果 ---
figure('Name', 'Figure 3: 调试 - CFAR检测');
subplot(1, 2, 1);
imagesc(range_axis, velocity_axis, db(cfar_input));
title(sprintf('送入CFAR的信号 (和波束 #%d)', beam_to_check_dbf));
xlabel('距离 (m)'); ylabel('速度 (m/s)'); colorbar; axis xy;
subplot(1, 2, 2);
imagesc(range_axis, velocity_axis, cfar_output);
title('CFAR检测输出 (1=目标)');
xlabel('距离 (m)'); ylabel('速度 (m/s)'); colorbar; axis xy;

% --- 【新增】步骤4: 可视化CFAR门限 ---
fprintf('--- STAGE 5: 正在可视化CFAR门限 ---\n');
% 找到信号最强的点，以该点为中心进行剖面分析
[~, max_idx] = max(cfar_input(:));
[max_v_idx, max_r_idx] = ind2sub(size(cfar_input), max_idx);

% --- Figure 4: 信号与门限对比图 ---
figure('Name', 'Figure 4: 调试 - 信号与CFAR门限对比');
% 提取最强目标所在距离单元的速度剖面
signal_slice = cfar_input(:, max_r_idx);
threshold_slice = cfar_threshold(:, max_r_idx);
plot(velocity_axis, signal_slice, 'b-', 'LineWidth', 1.5, 'DisplayName', '信号强度');
hold on;
plot(velocity_axis, threshold_slice, 'r--', 'LineWidth', 1.5, 'DisplayName', 'CFAR检测门限');
grid on;
legend;
title(sprintf('信号与门限对比 @ 距离 %.2f m', range_axis(max_r_idx)));
xlabel('速度 (m/s)');
ylabel('幅度');
% 标记出预设的目标速度
for i = 1:length(targets)
    xline(targets(i).Velocity, ':k', sprintf('目标%d 速度', i));
end
hold off;

%% 本地子函数区
% =========================================================================
% (此处省略了 fun_MTD_produce, local_execute_cfar, executeCFAR_2D, 
%  以及修正后的 Function_CFAR1D_sub 的代码，因为它们很长。
%  请将您项目中这些函数的最新版本复制到此处，
%  特别是要使用我们之前修正过的、能返回 threshold_matrix 的版本。)

% 示例：您需要将下面这个函数替换为您项目中的版本
function [MTD_Signal, pc_signal] = fun_MTD_produce(echoData, params)
%% 1.参数传递
% 1.1 控制与调试开关
%((abs(V)<1)||(R<5))如果速度太小或距离太近就不仿真了，当作背景杂波
% 只在这里控制是否画图
show_PC = params.Sig_Config.debug.show_PC;   % 脉冲压缩结果显示
show_FFT = params.Sig_Config.debug.show_FFT; % 速度维显示
graph = params.Sig_Config.debug.graph;       % 是否画图 
         
% 2.2 常数与单位定义
cj = sqrt(-1);
c = 2.99792458e8;           % 电磁波传播速度
PI2 = 2*pi;
us  = 1e-6;                 % time unit(us)
ns  = 1e-9;                 % time unit(ns)
MHz = 1e+6;                 % frequency unit(MHz)
KHz = 1e+3;                 % frequency unit(KHz)
GHz = 1e+9;                 % frequency unit(GHz)
fileTotalNums = 380;        % 数据bin文件个数 19*20=380
framesEachFile = 10;        % 每个新文件存储的帧数10个332PRT

% 2.3 雷达系统参数
prt   = params.Sig_Config.prt;             % X波段雷达PRT
prf   = params.Sig_Config.prf;
prtNum = params.Sig_Config.prtNum;         % 每帧信号的脉冲数
fs = params.Sig_Config.fs;                 % 生成原始信号的采样频率，每秒25,000,000次采样
ts = 1/fs;                                 % 采样时间间隔
deltaR = params.Sig_Config.deltaR;         % 距离分辨率单元

f0 = 0*MHz;                                % 起始频率
fc = params.Sig_Config.fc;                 % X波段雷达中心频率
wavelength = params.Sig_Config.wavelength; % 信号波长

B = params.Sig_Config.B;                   % 带宽
tao1  = params.Sig_Config.tao(1);          % 脉冲1脉宽
tao2  = params.Sig_Config.tao(2);          % 脉冲2脉宽
tao3  = params.Sig_Config.tao(3);          % 脉冲3脉宽
K1    = params.Sig_Config.K1;              % 短脉冲调频斜率（简单脉冲不需要）
K2    = params.Sig_Config.K2;              % 中脉冲负线性调频斜率
K3    = params.Sig_Config.K3;              % 长脉冲线性调频斜率

point_prt  = params.Sig_Config.point_prt(1); % 总PRT采集点数
point_prt1 = params.Sig_Config.point_prt(2); % 1脉冲区间点数
point_prt2 = params.Sig_Config.point_prt(3); % 2脉冲区间点数 
point_prt3 = params.Sig_Config.point_prt(4); % 3脉冲区间点数


% 2.生成用于脉冲压缩的发射参考信号

% 画出chirp信号时域波形
t1 = -tao1/2 : ts : tao1/2-ts;  % 时间变量
t2 = -tao2/2 : ts : tao2/2-ts;  
t3 = -tao3/2 : ts : tao3/2-ts; 
t123 = linspace(0,prt,point_prt);

% 理想仿真情况
pulse1 = sin(2*pi*t1+pi/2);                     % 窄脉冲：简单脉冲   4个点
pulse2 = exp(cj*2*pi*(f0*t2+0.5*K2*(t2.^2)));   % 中脉冲：负线性调频 200个点
pulse3 = exp(cj*2*pi*(f0*t3+0.5*K3*(t3.^2)));   % 长脉冲：正线性调频 700个点

% 手动输入实际情况，目前未知
%pulse2_real = [8.73485e-18 ,-0.0124099 ,-0.0104873 ,0.00915782 ,0.0237637 ,0.0113691 ,-0.0197015 ,-0.0368711 ,-0.0181063 ,0.0238755 ,0.0526078 ,0.0408238 ,-0.00647206 ,-0.0568911 ,-0.0755893 ,-0.0487868 ,0.00971216 ,0.0695179 ,0.101737 ,0.0925585 ,0.0466803 ,-0.0181685 ,-0.0805229 ,-0.123282 ,-0.137906 ,-0.124546 ,-0.0895256 ,-0.0419417 ,0.0092644 ,0.0571391 ,0.097408 ,0.128325 ,0.150067 ,0.163996 ,0.172 ,0.175996 ,0.17761 ,0.177992 ,0.177734 ,0.17684 ,0.174725 ,0.17026 ,0.16186 ,0.147675 ,0.125902 ,0.0952795 ,0.0557189 ,0.00900592 ,-0.0406418 ,-0.0864689 ,-0.119894 ,-0.132301 ,-0.117856 ,-0.0767001 ,-0.0172411 ,0.0441251 ,0.0871372 ,0.095372 ,0.0648793 ,0.00902167 ,-0.0450937 ,-0.0694989 ,-0.0520121 ,-0.00588109 ,0.0368519 ,0.0471478 ,0.0212279 ,-0.0159561 ,-0.0321687 ,-0.0169928 ,0.00967583 ,0.0199053 ,0.00752329 ,-0.00840642 ,-0.00962958];
%pulse2_imag = [-0.010182 ,-0.00426338 ,0.0125694 ,0.0176974 ,-0.00109573 ,-0.0255354 ,-0.0257313 ,0.00450108 ,0.0380632 ,0.0409565 ,0.00552931 ,-0.0420387 ,-0.0641617 ,-0.041699 ,0.0131109 ,0.0671493 ,0.0888244 ,0.0658318 ,0.00897237 ,-0.0565867 ,-0.104846 ,-0.119569 ,-0.0981699 ,-0.0494093 ,0.0121622 ,0.0719068 ,0.118977 ,0.147825 ,0.157798 ,0.151695 ,0.134071 ,0.109833 ,0.0833203 ,0.0578807 ,0.0358073 ,0.0184979 ,0.0066989 ,0.000745574 ,0.000744496 ,0.00666986 ,0.0183643 ,0.0354451 ,0.0571271 ,0.0819922 ,0.107759 ,0.131141 ,0.147924 ,0.153396 ,0.143243 ,0.114915 ,0.0692207 ,0.0116679 ,-0.0472347 ,-0.0935092 ,-0.113466 ,-0.0991065 ,-0.0532723 ,0.00841107 ,0.0614392 ,0.0825094 ,0.0620661 ,0.0120545 ,-0.0381229 ,-0.058303 ,-0.0379486 ,0.00495544 ,0.0364146 ,0.0335431 ,0.00392703 ,-0.0221936 ,-0.0217323 ,-0.000917822 ,0.0145387 ,0.0100754 ,-0.00330821];
%pulse2 = pulse2_real+1j*pulse2_imag; % 75个点

%pulse3_real = [0.00347321 ,0.00769471 ,0.00119518 ,-0.00865314 ,-0.00821249 ,0.00383747 ,0.0130384 ,0.00670282 ,-0.00935889 ,-0.0164105 ,-0.00467536 ,0.0143311 ,0.0196654 ,0.00388581 ,-0.0181421 ,-0.0239369 ,-0.00612985 ,0.0197439 ,0.0296144 ,0.0131487 ,-0.0167403 ,-0.0350039 ,-0.0256129 ,0.0053683 ,0.0347717 ,0.0402595 ,0.0167782 ,-0.0202684 ,-0.0461937 ,-0.0436764 ,-0.0135048 ,0.0267132 ,0.0537673 ,0.0523094 ,0.0227437 ,-0.0199996 ,-0.0549552 ,-0.0656411 ,-0.0472785 ,-0.00790095 ,0.0361765 ,0.0676954 ,0.0751575 ,0.0564085 ,0.0183415 ,-0.0264459 ,-0.0643147 ,-0.0847896 ,-0.083097 ,-0.0606393 ,-0.0237243 ,0.0187109 ,0.0575581 ,0.0855287 ,0.0983183 ,0.0948924 ,0.0770564 ,0.0485885 ,0.0142245 ,-0.0212723 ,-0.0538106 ,-0.0803855 ,-0.0992217 ,-0.109702 ,-0.11216 ,-0.107604 ,-0.0974397 ,-0.0832246 ,-0.0664766 ,-0.0485472 ,-0.0305497 ,-0.013335 ,0.00249855 ,0.0165738 ,0.0286919 ,0.0387858 ,0.0468762 ,0.0530328 ,0.0573421 ,0.0598832 ,0.060711 ,0.0598457 ,0.0572702 ,0.0529331 ,0.0467587 ,0.0386642 ,0.0285839 ,0.016501 ,0.002486 ,-0.0132596 ,-0.0303575 ,-0.0482109 ,-0.0659736 ,-0.0825414 ,-0.096577 ,-0.106581 ,-0.111021 ,-0.108516 ,-0.0980836 ,-0.07941 ,-0.0531215 ,-0.0209855 ,0.014023 ,0.0478671 ,0.0758589 ,0.0933513 ,0.096652 ,0.084018 ,0.0564998 ,0.0183531 ,-0.0232531 ,-0.0593895 ,-0.0813211 ,-0.082912 ,-0.0628401 ,-0.0258185 ,0.0178915 ,0.0549782 ,0.0731889 ,0.0658645 ,0.0351666 ,-0.00767338 ,-0.0458739 ,-0.0636303 ,-0.0532198 ,-0.0193486 ,0.0219808 ,0.0505018 ,0.0518531 ,0.0257335 ,-0.0129946 ,-0.0419767 ,-0.0443417 ,-0.0194313 ,0.0160642 ,0.0384939 ,0.0332 ,0.00511814 ,-0.024382 ,-0.0332687 ,-0.0158839 ,0.0124542 ,0.0279986 ,0.0186306 ,-0.00577235 ,-0.022492 ,-0.0170077 ,0.00363389 ,0.018342 ,0.013329 ,-0.00433515 ,-0.015166 ,-0.00861794 ,0.00614772 ,0.0119063 ,0.00348724 ,-0.00742221 ,-0.00777207 ,0.00106588 ,0.00680575];
%pulse3_imag = [-0.00601577 ,0.0015725 ,0.00872502 ,0.00461064 ,-0.00708884 ,-0.0113048 ,-0.000957481 ,0.0125797 ,0.0123299 ,-0.00335368 ,-0.0174487 ,-0.0131091 ,0.00684821 ,0.0219256 ,0.0153313 ,-0.00812548 ,-0.0261348 ,-0.0204812 ,0.00532846 ,0.028919 ,0.0289951 ,0.00395725 ,-0.0267089 ,-0.0384408 ,-0.0210585 ,0.0136663 ,0.0411138 ,0.0416486 ,0.0139467 ,-0.0247613 ,-0.0504005 ,-0.0471194 ,-0.0162333 ,0.0254565 ,0.0557316 ,0.0589169 ,0.0332821 ,-0.00916684 ,-0.0493015 ,-0.0698879 ,-0.0626595 ,-0.0307793 ,0.013523 ,0.054378 ,0.0781994 ,0.0779072 ,0.0543502 ,0.015027 ,-0.0289374 ,-0.0662922 ,-0.0885401 ,-0.0915578 ,-0.07583 ,-0.0455722 ,-0.00722005 ,0.0322116 ,0.0665133 ,0.0911896 ,0.103841 ,0.104092 ,0.0932026 ,0.0735308 ,0.0479656 ,0.0194422 ,-0.00941832 ,-0.0365266 ,-0.0604152 ,-0.0802288 ,-0.0956473 ,-0.106774 ,-0.114013 ,-0.117955 ,-0.119279 ,-0.11868 ,-0.116812 ,-0.114259 ,-0.111514 ,-0.108975 ,-0.106943 ,-0.105628 ,-0.105155 ,-0.105562 ,-0.106809 ,-0.10877 ,-0.111235 ,-0.113901 ,-0.116373 ,-0.118159 ,-0.11868 ,-0.117288 ,-0.113296 ,-0.106034 ,-0.0949235 ,-0.0795702 ,-0.0598803 ,-0.0361795 ,-0.00932269 ,0.019232 ,0.0474155 ,0.0726385 ,0.0920092 ,0.102688 ,0.102371 ,0.0898357 ,0.0654796 ,0.0316885 ,-0.00709768 ,-0.0447673 ,-0.0744358 ,-0.0898073 ,-0.0867818 ,-0.0649259 ,-0.028319 ,0.0146942 ,0.053104 ,0.076059 ,0.076281 ,0.0529992 ,0.0131688 ,-0.0299468 ,-0.0609103 ,-0.067875 ,-0.0478369 ,-0.00888603 ,0.032231 ,0.0569993 ,0.0538624 ,0.0245768 ,-0.0156554 ,-0.0453913 ,-0.0484963 ,-0.0237976 ,0.0133876 ,0.0399284 ,0.0393641 ,0.0130669 ,-0.0201066 ,-0.0366495 ,-0.0254254 ,0.00376108 ,0.0275117 ,0.0273915 ,0.00503775 ,-0.0193262 ,-0.0246105 ,-0.007635 ,0.0143726 ,0.0205041 ,0.00638737 ,-0.0121924 ,-0.016179 ,-0.00309935 ,0.0113537 ,0.0115379 ,-0.000874347 ,-0.0102731 ,-0.00640668 ,0.00414119 ,0.00778108 ,0.00139083];
%pulse3 = pulse3_real+1j*pulse3_imag; % 160个点

% 3.ISTC 调制补偿
% [stc_cure, echo] = fun_iSTC(echo); % 缺失数据


% 4.脉冲压缩

[Echo_0] = fun_lss_pulse_compression(echoData,params,show_PC, pulse1, pulse2, pulse3,point_prt1,point_prt2,point_prt3);   % (1536，1031)

% 删除重复的区间
% Echo_0=fun_lss_range_concate(prtNum,Echo_0);%(1536,868)

% point_prt长度改变了
[~,point_prt] = size(Echo_0);

% 将脉冲压缩的结果赋值给新的输出变量
pc_signal = Echo_0;

% 5.MTD (多普勒处理)
[m,n] = size(Echo_0);
MTD_Signal = fun_Process_MTD(Echo_0, n, m); % 做MTD，多普勒处理


% 6.杂波抑制
MTD_Signal = fun_0v_pressing(MTD_Signal, params);   % 0速抑制，压制0速附近一定区域的峰值


% 7. 调试与可视化
if(show_FFT==1)
    for i=1:point_prt
        figure(1);
        plot(20*log10(abs(MTD_Signal(:,i)))),title('FFT');

        pause(0.05);
    end
end

% 8.画出MTD三维图
if(graph==1)
    
    R_range = [500,4000];
    V_range = [-20,20];
    R_range_point = round(R_range/6);
    V_range_point = [691,845];
    
    figure(8);
    R_point = 6;%两点间距6m
    r0 = 0:R_point:point_prt*R_point-R_point;%距离轴
    fd = linspace(-prf/2,prf/2,prtNum);
    v0 = fd*wavelength/2;%速度轴
    MTD_Signal_simu_log = 20*log10(abs(MTD_Signal)/max(max(abs(MTD_Signal))));
    mesh(r0,v0,MTD_Signal_simu_log);
    xlabel('距离');
    ylabel('速度m/s');
    zlabel('幅度dB');
    title('0波束仿真信号MTD');
    xlim(R_range);
    ylim(V_range);
    
    MTD_local = MTD_Signal_simu_log(V_range_point(1) : V_range_point(2),R_range_point(1) : R_range_point(2));
    r0_local = r0(R_range_point(1) : R_range_point(2));
    v0_local = v0(V_range_point(1) : V_range_point(2));
    

%     % 画出速度维
%     MTD_max=max(max(MTD_local));
%     [vindex,rindex]=find(MTD_local(:,:)==MTD_max);
% 
%     figure(9);
%     plot(v0_local,(MTD_local(:,rindex))),xlabel('速度m/s'),ylabel('幅度dB');title('速度维');
%     
% 
%     
%     % 画出距离维
%     figure(11);
%     plot(r0_local,(MTD_local(vindex,:))),xlabel('距离'),ylabel('幅度dB');title('距离维');  
end

end


function ref_pulse = tx_pulse_placeholder(config)
    % 简化的发射脉冲生成函数 (占位符)
    ts = 1/config.Sig_Config.fs;
    t2 = linspace(-config.Sig_Config.tao(2)/2, config.Sig_Config.tao(2)/2, config.Sig_Config.point_prt_segments(2));
    pulse2 = exp(1j*pi*(-config.Sig_Config.B / config.Sig_Config.tao(2))*(t2.^2));
    ref_pulse = zeros(1, config.Sig_Config.point_PRT);
    ref_pulse(229:229+723-1) = pulse2;
end


% --- 【关键】需要修改的CFAR函数 ---
function [cfar_flag, threshold_matrix] = local_execute_cfar(mtd_amplitude_map, cfar_params, config)
    point_prt_narrow = config.Sig_Config.point_prt(2);
    point_prt_medium = config.Sig_Config.point_prt(3);
    
    MTD_p0 = mtd_amplitude_map(:, 1:point_prt_narrow);
    MTD_p1 = mtd_amplitude_map(:, point_prt_narrow+1 : point_prt_narrow+point_prt_medium);
    MTD_p2 = mtd_amplitude_map(:, point_prt_narrow+point_prt_medium+1 : end);

    [cfar_0, th_0] = executeCFAR_2D(MTD_p0, config.cfar);
    [cfar_1, th_1] = executeCFAR_2D(MTD_p1, config.cfar);
    [cfar_2, th_2] = executeCFAR_2D(MTD_p2, config.cfar);

    cfar_flag = zeros(size(mtd_amplitude_map));
    cfar_flag(:, 1:point_prt_narrow) = cfar_0;
    cfar_flag(:, point_prt_narrow+1 : point_prt_narrow+point_prt_medium) = cfar_1;
    cfar_flag(:, point_prt_narrow+point_prt_medium+1 : end) = cfar_2;
    
    threshold_matrix = zeros(size(mtd_amplitude_map));
    threshold_matrix(:, 1:point_prt_narrow) = th_0;
    threshold_matrix(:, point_prt_narrow+1 : point_prt_narrow+point_prt_medium) = th_1;
    threshold_matrix(:, point_prt_narrow+point_prt_medium+1 : end) = th_2;
end

function [cfarResultFlag, threshold_matrix] = executeCFAR_2D(echo_mtd_amp, cfar_params)
    [vCellNum_org, rCellNum_org] = size(echo_mtd_amp);
    MTD_0v_num = cfar_params.MTD_0v_num;
    
    valid_vel_mask = true(vCellNum_org, 1);
    zero_v_center = round(vCellNum_org / 2) + 1;
    suppress_start = max(1, zero_v_center - MTD_0v_num);
    suppress_end = min(vCellNum_org, zero_v_center + MTD_0v_num);
    valid_vel_mask(suppress_start:suppress_end) = false;
    
    echo_mtd_abs_used = echo_mtd_amp(valid_vel_mask, :);
    
    % 调用能返回门限的CFAR函数
    [cfarResult_R, threshold_R] = Function_CFAR1D_sub(echo_mtd_abs_used, cfar_params.refCells_R, cfar_params.saveCells_R, cfar_params.T_CFAR, cfar_params.CFARmethod_R);
    
    cfarResultFlag = zeros(vCellNum_org, rCellNum_org);
    cfarResultFlag(valid_vel_mask, :) = cfarResult_R;
    
    threshold_matrix = zeros(vCellNum_org, rCellNum_org);
    threshold_matrix(valid_vel_mask, :) = threshold_R;
end



% --- 【核心修改】Function_CFAR1D_sub ---
function [PeakDetectionoutput, threshold_matrix]=Function_CFAR1D_sub(datamatrix,refCellNum,saveCellNum,T_CFAR,CFARmethod)
    [RowNum_ASR,ColNum_ASR]=size(datamatrix);
    PeakDetectionoutput=zeros(RowNum_ASR,ColNum_ASR);
    threshold_matrix = zeros(RowNum_ASR,ColNum_ASR); % << 新增: 初始化门限矩阵
    
    for y=1:ColNum_ASR
        % ... (此处省略了原函数中计算 ref_average_used 的代码)
        if CFARmethod==0
            refL_average = zeros(RowNum_ASR,1);  % 初始化
            refR_average = zeros(RowNum_ASR,1);  % 初始化
        else
            refL_average = Inf*ones(RowNum_ASR,1);  % 初始化
            refR_average = Inf*ones(RowNum_ASR,1);  % 初始化
        end
        refL1=y-(saveCellNum+refCellNum);% 左参考单元的左边界
        refL2=y-saveCellNum-1;           % 左参考单元的右边界
        refR1=y+saveCellNum+1;           % 右参考单元的左边界
        refR2=y+saveCellNum+refCellNum;  % 右参考单元的右边界

        if refL1>=1
            refL_average=mean(datamatrix(:,refL1:refL2),2);  % 左参考单元的平均
        else % 左参考单元点数不够,则用检测点右边的数据估计杂波水平
            refL_average = mean(datamatrix(:,refR1:refR2),2);
        end

        if refR2<=ColNum_ASR
            refR_average=mean(datamatrix(:,refR1:refR2),2);    % 右参考单元的平均
        else  % 右参考单元点数不够,则用检测点左边的数据估计杂波水平
            refR_average=mean(datamatrix(:,refL1:refL2),2);
        end

        if CFARmethod==0
            ref_average_used = max(refL_average,refR_average);     % 选大
        else
            ref_average_used = min(refL_average,refR_average);     % 选小
        end
        % ref_average_used = max(refL_average,refR_average);
        
        threshold_CFAR = ref_average_used.*T_CFAR;
        threshold_matrix(:,y) = threshold_CFAR; % << 新增: 保存当前列的门限
        
        flag=datamatrix(:,y) >= threshold_CFAR;
        PeakDetectionoutput(:,y) = flag; % 简化输出
    end
end

function flag = check_single_point_cfar(data_row, cut_idx, ref_cells, guard_cells, T_cfar, method)
    % 简化版1D-CFAR，只检测单个点
    flag = 0;
    len = length(data_row);
    
    % 定义参考窗和保护窗的边界
    left_guard_end = cut_idx - guard_cells - 1;
    left_ref_start = left_guard_end - ref_cells + 1;
    
    right_guard_start = cut_idx + guard_cells + 1;
    right_ref_end = right_guard_start + ref_cells - 1;
    
    noise_L = -1; % 初始值，表示左侧噪声未计算
    noise_R = -1; % 初始值，表示右侧噪声未计算
    
    % 计算左侧参考窗的噪声平均值
    if left_ref_start >= 1
        noise_L = mean(data_row(left_ref_start:left_guard_end));
    end
    
    % 计算右侧参考窗的噪声平均值
    if right_ref_end <= len
        noise_R = mean(data_row(right_guard_start:right_ref_end));
    end
    
    % 噪声估计策略
    if noise_L == -1 && noise_R == -1
        return; % 边缘情况：左右都无法估计噪声，直接返回
    elseif noise_L == -1
        noise_est = noise_R; % 左侧无参考窗，使用右侧噪声
    elseif noise_R == -1
        noise_est = noise_L; % 右侧无参考窗，使用左侧噪声
    else
        % 正常情况：两侧都有参考窗
        if method == 0 % 选大 (CA-GO)
            noise_est = max(noise_L, noise_R);
        else % 选小 (CA-SO)
            noise_est = min(noise_L, noise_R);
        end
    end
    
    % 比较门限
    threshold = noise_est * T_cfar;
    if data_row(cut_idx) > threshold
        flag = 1; % 超过门限，判定为目标
    end
end