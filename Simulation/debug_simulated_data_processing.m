% debug_simulated_data_processing.m
%
% 这是一个专门用于调试和验证信号处理流程的脚本。
% 它会加载仿真生成的数据，并逐级执行处理算法，
% 在每个关键步骤后都生成可视化图表，以便于分析信号的变化。
%
% --- 调试步骤 ---
% 1. 运行此脚本。
% 2. 观察Figure 1 (DBF后波束功率): 确认在目标所在距离门附近有能量条带。
% 3. 观察Figure 2 (MTD后RDM图): 确认在预设目标的距离和速度处有清晰的亮点。
%    这是最关键的一步，如果这里没有亮点，说明MTD处理有问题或信噪比太低。
% 4. 观察Figure 3 (CFAR检测结果): 
%    - 左图是送入CFAR检测器的信号，应该有亮点。
%    - 右图是CFAR的输出，如果左图有亮点而右图没有，则说明CFAR门限(T_CFAR)太高。
% 修改记录
%  date       by      version   modify
%  25/08/27   XZR      v1.0     创建

clc; clear; close all;

%% 1. 加载配置、仿真数据和DBF系数
% =========================================================================
% --- 0.0 实际测量时一些修正参数 ---
% 伺服角（方位角）修正系数
config.corrected.North = - 242;                % 雷达指北角（见雷达系统设置文件 SysSet.ini.bak）
config.corrected.FixAngle = 35;              % 雷达固定角（见初始化参数文件 InitPara.ini）
% 俯仰角修正系数 
config.corrected.ELeAngleSettingValue = -10; % 雷达俯仰设置值（见雷达系统设置文件 SysSet.ini.bak）

% --- 1.1 流程控制 --- 
config.frame_range = 0 : 600;                         % 指定要处理的帧范围
config.save_options.save_frameheads_mat = true;       % 开关：是否保存读取文件时检测的文件帧头信息
config.save_options.save_iq_mat_before_DBF = true;    % 开关：是否保存原始的16通道I/Q数据
config.save_options.save_iq_mat_after_DBF = true;     % 开关：是否保存第一阶段的.mat格式I/Q数据
config.save_options.save_pc_mat = true;               % 开关：是否保存脉冲压缩(PC)结果
config.save_options.save_mtd_mat = true;              % 开关：是否保存第二阶段的MTD结果
config.save_options.save_cfar_mat = true;             % 开关：是否保存第三阶段的CFAR结果矩阵 (CFARflag)
config.save_options.save_beam_sum_cfar_mat = true;    % 开关：是否保存第三阶段和波束CFAR目标检测阶段结果
config.save_options.save_final_log = true;            % 开关：是否保存第四阶段差波束参数测量的结果（按帧保存）   
config.save_options.save_cumulative_log = true;       % 开关：是否保存所有帧累积的测量结果
config.save_options.save_to_bin = true;               % 开关：是否保存目标检测点信息为二进制.bin文件

% --- 1.2 路径配置 ---

% --- 1.3 实验数据路径 ---
config.base_path = 'C:\Users\a\Desktop\MatlabProcess_xuzerui_latest\XZR\new_simulation\Simulation';
config.DBF_coef_path = fullfile(config.base_path, 'X8数据采集250522_DBFcoef.csv');          % 用于 DDC 数据波束成形的 DBF 系数，先处理成CSV（逗号分隔值）文件
config.angle_k_path = fullfile(config.base_path,'R9-DMX3-2024001_Angle_k.csv');             % 用于和差比幅的K值矩阵

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
config.cfar.T_CFAR = 5;                    % 恒虚警标称化因子
config.cfar.MTD_V = 3;                     % 杂波区速度范围，速度在 -3 m/s 到 +3 m/s 范围内的区域都当作是地杂波区域，在CFAR检测中忽略掉。

% --- 速度维参数 ---
config.cfar.refCells_V = 5;                % 速度维 参考单元数
config.cfar.saveCells_V = 14;              % 速度维 保护单元数
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


% --- 加载DBF系数 ---
try
    DBF_coeffs_data = readmatrix(config.DBF_coef_path);
    DBF_coeffs_data_I = double(DBF_coeffs_data(:, 1:2:end));
    DBF_coeffs_data_Q = double(DBF_coeffs_data(:, 2:2:end));
    DBF_coeffs_data_C = DBF_coeffs_data_I + 1j * DBF_coeffs_data_Q;
catch ME
    error('读入DBF系数文件失败: %s', ME.message);
end

% --- 加载仿真数据 ---
sim_data_folder = 'simulated_data_array'; 
sim_frame_file = 'frame_sim_array_0.mat';
sim_data_path = fullfile('C:\Users\a\Desktop\MatlabProcess_xuzerui_latest\XZR\new_simulation\Simulation', sim_data_folder, sim_frame_file);
if ~exist(sim_data_path, 'file'), error('仿真数据文件不存在！'); end
sim_data = load(sim_data_path);

%% 2. 执行并可视化 Stage 1: 数字波束形成 (DBF)
% =========================================================================
fprintf('--- 步骤2: 执行并可视化DBF ---\n');
raw_iq_16ch = sim_data.raw_iq_data;
iq_data_13beam = complex(zeros(config.Sig_Config.prtNum, config.Sig_Config.point_PRT, config.Sig_Config.beam_num));
for i = 1:config.Sig_Config.prtNum
    single_pulse_16ch = squeeze(raw_iq_16ch(i, :, :));
    single_pulse_13beam = single_pulse_16ch * DBF_coeffs_data_C';
    iq_data_13beam(i, :, :) = single_pulse_13beam;
end

% --- 可视化DBF结果 ---
% 我们选择一个应该有强信号的波束进行显示
% 预设目标在 -5度, 8.2度, 15度。我们检查第5个波束 (标称7.5度)
beam_to_check_dbf = 5; 
figure('Name', '调试 - DBF后波束功率');
imagesc(db(abs(squeeze(iq_data_13beam(:,:,beam_to_check_dbf)))));
title(sprintf('DBF后波束 #%d 的功率 (dB)', beam_to_check_dbf));
xlabel('距离单元 (快时间)');
ylabel('脉冲数 (慢时间)');
colorbar;
fprintf('  > DBF可视化完成。请检查Figure 1中是否有明显的水平能量条带。\n');

%% 3. 执行并可视化 Stage 2: MTD处理
% =========================================================================
fprintf('--- 步骤3: 执行并可视化MTD ---\n');
% 准备MTD输入
iq_data1 = iq_data_13beam;
angle1 = sim_data.servo_angle;
iq_data2 = iq_data1; % 用同一帧数据
angle2 = angle1;
% 执行MTD
[mtd_results, ~, ~] = process_stage2_mtd(iq_data1, iq_data2, angle1, angle2, config);

% --- 可视化RDM图 ---
% 检查一个包含目标的和差波束对，例如第4对 (波束4和5)
beam_pair_to_check_mtd = 4;
slice_to_check_mtd = 4; % 检查第一个切片
RDM_beam1 = squeeze(mtd_results{beam_pair_to_check_mtd}(slice_to_check_mtd, :, :));
RDM_beam2 = squeeze(mtd_results{beam_pair_to_check_mtd+1}(slice_to_check_mtd, :, :));
RDM_sum_amp = abs(RDM_beam1) + abs(RDM_beam2);  % 和波束

deltaR = config.Sig_Config.c / (2 * config.Sig_Config.fs);
range_axis = (0:config.Sig_Config.point_PRT-1) * deltaR;
velocity_axis = linspace(-config.Sig_Config.prf/2, config.Sig_Config.prf/2, config.Sig_Config.prtNum) * config.Sig_Config.wavelength / 2;

figure('Name', '调试 - MTD后RDM');
imagesc(range_axis, velocity_axis, db(RDM_sum_amp));
title(sprintf('和波束对 #%d 的距离-多普勒图 (RDM)', beam_pair_to_check_mtd));
xlabel('距离 (m)');
ylabel('速度 (m/s)');
colorbar;
axis xy;
fprintf('  > MTD可视化完成。请检查Figure 2中是否有对应目标的亮点。\n');

%% 4. 执行并可视化 Stage 3: CFAR检测
% =========================================================================
fprintf('--- 步骤4: 执行并可视化CFAR ---\n');
current_frame_num = 0;
[prelim_log, cfar_flags] = process_stage3_detection(mtd_results, angle1, config, current_frame_num);

% --- 可视化CFAR的输入和输出 ---
% 我们仍然使用之前选择的那个和差波束对
cfar_input = RDM_sum_amp;
% cfar_flags 是一个cell数组，每个cell对应一个波束对
% 每个cell内部是 (vel, range, slice) 的三维矩阵
cfar_output = squeeze(cfar_flags{beam_pair_to_check_mtd}(:,:,slice_to_check_mtd)); 

figure('Name', '调试 - CFAR检测');
subplot(1, 2, 1);
imagesc(range_axis, velocity_axis, db(cfar_input));
title(sprintf('送入CFAR的信号 (和波束对 #%d)', beam_pair_to_check_mtd));
xlabel('距离 (m)'); ylabel('速度 (m/s)'); colorbar; axis xy;

subplot(1, 2, 2);
imagesc(range_axis, velocity_axis, cfar_output);
title(sprintf('CFAR检测输出 (1=目标)'));
xlabel('距离 (m)'); ylabel('速度 (m/s)'); colorbar; axis xy;
fprintf('  > CFAR可视化完成。请对比Figure 3的左右两图。\n');

%% 5. 执行 Stage 4 并显示最终结果
% =========================================================================
% fprintf('--- 步骤5: 执行参数测量并显示最终结果 ---\n');
% frame_headers1(1).freq_no = 7;
% final_log = process_stage4_measurement(prelim_log, mtd_results, config, frame_headers1);
% 
% fprintf('\n--- 最终检测结果 ---\n');
% if isempty(final_log)
%     fprintf('未检测到任何目标。\n');
% else
%     fprintf('共检测到 %d 个目标。\n', length(final_log));
%     disp(final_log);
% end
