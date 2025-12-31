% main_test_with_simulated_data.m
%
%
%
% 修改记录
% date       by      version   modify
% 25/06/25   XZR      v1.0     创建
% 25/07/15   XZR      v2.0     增添波束成形部分，利用和波束做CFAR目标检测，差波束用来单脉冲测角
% 25/07/27   XZR      v3.0     更新了对第四阶段返回值的处理，以适应新的两层数据结构
% 25/08/28   XZR      v4.0     修改为处理模拟仿真的信号数据

clc; clear; close all;

%% 1. 全局参数配置区 (所有参数在此统一设置 Config)
fprintf('--- 开始进行全局参数配置 ---\n');

% --- 0.0 实际测量时一些修正参数 ---
% 伺服角（方位角）修正系数
config.corrected.North = - 242;                % 雷达指北角（见雷达系统设置文件 SysSet.ini.bak）
config.corrected.FixAngle = 35;              % 雷达固定角（见初始化参数文件 InitPara.ini）
% 俯仰角修正系数 
config.corrected.ELeAngleSettingValue = -10.3; % 雷达俯仰设置值（见雷达系统设置文件 SysSet.ini.bak）

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
% config.DBF_coef_path = fullfile(config.base_path, 'DBF_data', 'X8数据采集250522_DBFcoef.csv');          % 用于 DDC 数据波束成形的 DBF 系数，先处理成CSV（逗号分隔值）文件
% config.angle_k_path = fullfile(config.base_path, 'K_value', 'R9-DMX3-2024001_Angle_k.csv');             % 用于和差比幅的K值矩阵

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
config.cfar.T_CFAR = 10;                   % 恒虚警标称化因子
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

% --- 1.7 创建所有输出目录 ---
% config.output_paths.framehead = fullfile(config.base_path, num2str(config.n_exp), 'Framehead_information');
% config.output_paths.iq_before_DBF = fullfile(config.base_path, num2str(config.n_exp), 'iq_data_before_DBF');
% config.output_paths.iq_after_DBF = fullfile(config.base_path, num2str(config.n_exp), 'BasebandRawData_mat');
% config.output_paths.pc = fullfile(config.base_path, num2str(config.n_exp), 'pulse_compressed_data');
% config.output_paths.mtd = fullfile(config.base_path, num2str(config.n_exp), ['MTD_data_win', num2str(config.mtd.win_size)]);
% config.output_paths.cfar = fullfile(config.base_path, num2str(config.n_exp), ['cfarFlag4_T', num2str(config.cfar.T_CFAR)]);
% config.output_paths.beam_sum_cfar = fullfile(config.base_path, num2str(config.n_exp), ['beam_sum_cfarFlag4_T', num2str(config.cfar.T_CFAR)]);
% config.output_paths.beam_diff_estimation = fullfile(config.base_path, num2str(config.n_exp), ['beam_diff_estimation_cfarFlag4_T', num2str(config.cfar.T_CFAR)]);
% config.output_paths.beam_diff_estimation_cumulative = fullfile(config.base_path, num2str(config.n_exp), 'beam_diff_estimation_cumulative');
% config.output_paths.bin_output = fullfile(config.base_path, num2str(config.n_exp), 'Save_bin');
% 
% if config.save_options.save_frameheads_mat,    mkdir(config.output_paths.framehead);     end   % 开关：是否保存读取文件时检测的文件帧头信息
% if config.save_options.save_iq_mat_before_DBF, mkdir(config.output_paths.iq_before_DBF); end
% if config.save_options.save_iq_mat_after_DBF,  mkdir(config.output_paths.iq_after_DBF);  end
% if config.save_options.save_pc_mat, mkdir(config.output_paths.pc);                       end
% if config.save_options.save_mtd_mat,           mkdir(config.output_paths.mtd);           end
% if config.save_options.save_cfar_mat,          mkdir(config.output_paths.cfar);          end
% if config.save_options.save_beam_sum_cfar_mat, mkdir(config.output_paths.beam_sum_cfar); end
% if config.save_options.save_final_log,         mkdir(config.output_paths.beam_diff_estimation);            end
% if config.save_options.save_cumulative_log,    mkdir(config.output_paths.beam_diff_estimation_cumulative); end
% if config.save_options.save_to_bin,            mkdir(config.output_paths.bin_output);    end

%% 2. 初始化与预加载

DBF_coef_path = 'C:\Users\a\Desktop\MatlabProcess_xuzerui_latest\XZR\new_simulation\Simulation\X8数据采集250522_DBFcoef.csv'; % 用于 DDC 数据波束成形的 DBF 系数，先处理成CSV（逗号分隔值）文件
   
% --- 加载DBF系数 ---
try
    DBF_coeffs_data = readmatrix(DBF_coef_path);
    DBF_coeffs_data_I = double(DBF_coeffs_data(:, 1:2:end));
    DBF_coeffs_data_Q = double(DBF_coeffs_data(:, 2:2:end));
    DBF_coeffs_data_C = DBF_coeffs_data_I + 1j * DBF_coeffs_data_Q;
catch ME
    error('读入DBF系数文件失败: %s', ME.message);
end

%% 3. 加载仿真数据
% =========================================================================
tic;
fprintf('--- 开始处理流程 ---\n');
fprintf('--- 正在加载仿真生成的数据 ---\n');

% --- 2.1 定义仿真数据的文件路径 ---
sim_data_folder = 'simulated_data_array'; 
sim_frame_file = 'frame_sim_array_0.mat';
sim_data_path = fullfile('C:\Users\a\Desktop\MatlabProcess_xuzerui_latest\XZR\new_simulation\Simulation', sim_data_folder, sim_frame_file);

if ~exist(sim_data_path, 'file')
    error('仿真数据文件不存在！请先运行 main_simulate_echoes_with_array.m 生成数据。');
end

% --- 2.2 加载 .mat 文件 ---
sim_data = load(sim_data_path);
fprintf('  > 仿真数据 %s 加载成功。\n', sim_frame_file);

% --- 2.3 准备送入处理流程的变量 ---
% MTD 处理 (Stage 2) 需要两个连续的帧 (iq_data1, iq_data2) 来进行拼接。
% 为了测试，将同一帧的仿真数据使用两次。

% 准备第一帧数据
% 注意: 仿真脚本生成的 raw_iq_data 是16通道的原始数据。
% 在真实流程中，iq_data_13beam 是经过DBF处理的，但在这里我们可以
% 暂时用16通道数据的前13个通道来代替，或者直接用原始数据进行测试。

% 逐脉冲DBF处理
iq_data_before_DBF = sim_data.raw_iq_data;    % DBF处理前的16通道原始I/Q数据
current_prtnum = 1;
while current_prtnum <= config.Sig_Config.prtNum
    iq_data_temp = squeeze(iq_data_before_DBF(current_prtnum, :, :));      % 提取当前脉冲的16通道数据 (3404 x 16)
    iq_data1(current_prtnum,:,:) = iq_data_temp * DBF_coeffs_data_C.';     % 与DBF系数矩阵相乘
    current_prtnum = current_prtnum + 1;
end
angle1 = sim_data.servo_angle;
fprintf('  > DBF处理完成，已生成13波束数据。\n');


% Stage 4 需要帧头信息来获取频点号，我们在这里模拟一个简单的帧头
frame_headers1(1).freq_no = 6; % 使用您在仿真中设置的频点号

% 准备第二帧数据 (使用相同的数据)
iq_data2 = iq_data1;
angle2 = angle1;
frame_headers2 = frame_headers1;

%% 3. 执行后续处理流程 (从 Stage 2 开始)
% =========================================================================
% 从这里开始，后续的代码几乎不需要修改，因为它接收的数据格式
% (iq_data1, iq_data2, angle1, angle2) 与真实数据完全一致。

current_frame_num = 0; % 我们正在处理第0帧仿真数据

fprintf('\n================== 正在处理仿真帧: %d ==================\n', current_frame_num);

% --- STAGE 2: MTD处理 ---
fprintf('STAGE 2: 正在对仿真数据进行MTD处理...\n');
[mtd_results, ~, ~] = process_stage2_mtd(iq_data1, iq_data2, angle1, angle2, config);

% --- STAGE 3: 波束成形及和波束CFAR检测 ---
fprintf('STAGE 3: 正在进行CFAR检测...\n');
[prelim_log, ~] = process_stage3_detection(mtd_results, angle1, config, current_frame_num);

% --- STAGE 4: 差波束单脉冲测角（参数测量） ---
fprintf('STAGE 4: 正在进行参数测量...\n');
final_log = process_stage4_measurement(prelim_log, mtd_results, config, frame_headers1);

%% 4. 分析最终结果
% =========================================================================
fprintf('\n--- 信号处理流程完成 ---\n');
fprintf('共检测到 %d 个目标。\n', length(final_log));
disp('检测到的目标参数为:');
disp(final_log);

% --- 结果验证 ---
if ~isempty(final_log)
    fprintf('\n--- 正在将检测结果与预设目标进行对比 ---\n');
    preset_targets = targets; % 从您的仿真脚本中复制 targets 结构体到这里
    for i = 1:length(final_log)
        fprintf('--- 检测到的目标 #%d ---\n', i);
        fprintf('  - 距离: %.2f m\n', final_log(i).range_m);
        fprintf('  - 速度: %.2f m/s\n', final_log(i).velocity_ms);
        fprintf('  - 俯仰角: %.2f 度\n', final_log(i).elevation_deg);
    end
end

    




