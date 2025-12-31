function [MTD_results, PC_results] = process_stage2_mtd(iq_data, angle, config)
% PROCESS_STAGE2_MTD - 对输入的I/Q数据执行MTD处理
%
% 本函数是原 main_produce_dataset_win_xzr_v2.m 脚本的功能化改造版本。
% 它接收一个数据帧并对其进行完整的脉压 MTD处理流程。
%
% 输入参数:
%   iq_data  - (complex 3D matrix) 当前帧的I/Q数据。
%   config   - (struct) 包含所有MTD处理参数的配置结构体。
%
% 输出参数:
%   mtd_results - (cell) MTD处理结果。每个单元格包含一个波束的所有切片结果。

%% 1. 从config结构体中获取参数
beam_num = config.mtd.beam_num;
win_size = config.mtd.win_size;
% 将核心的雷达系统参数打包，传递给下一层函数
params_for_produce = config; 
% params_for_produce.debug = config.debug; % 传递调试开关

%% 2. 提取13波束信号数据
% 提取信号数据
echo_win_beams = iq_data;

% 提取角度数据
servo_angle_win = angle;

%% 3. 对所有波束和切片进行MTD处理
MTD_results = zeros(332, 3404, 13); % 初始化用于保存最终结果的Cell数组
PC_results = zeros(332, 3404, 13);
% --- 波束循环 ---
MTD_data_for_one_beam = []; % 初始化临时变量
for b = 1:beam_num
    current_beam_echo_win = echo_win_beams(:,:,b);
    [total_prts, ~] = size(current_beam_echo_win);
    prts_per_slice = total_prts;
   
    start_row = 1;
    end_row = total_prts;
    echo_segment = current_beam_echo_win(start_row:end_row, :);

    % --- 调用核心MTD处理链函数 ---
    % 注意：这里 fun_MTD_produce 也需要被相应修改，以接收params结构体
    [MTD_signal, PC_signal] = fun_MTD_produce(echo_segment, params_for_produce);
    MTD_data_for_one_beam = MTD_signal;
    PC_signal_for_one_beam = PC_signal;
    
    MTD_results(:,:,b) = MTD_data_for_one_beam;
    PC_results(:,:,b) = PC_signal_for_one_beam;
end

end
