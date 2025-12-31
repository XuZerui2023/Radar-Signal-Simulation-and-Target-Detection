clc; clear; close all;
%% 1. 参数设置
%%% 工作频率
c = 3e8;
freq = 10e9;
lambda = c/freq;    % 波长
k = 2*pi/lambda;    % 波数

%%% 阵列参数 (URA)
Nx = 8;                 % x 轴阵元数量
Ny = 8;                 % y 轴阵元数量
N = Nx * Ny;            % 总阵元数量
dx = 0.5*lambda;        % x 轴阵元间隔 
dy = 0.5*lambda;        % y 轴阵元间隔

% 生成 URA 阵元坐标 (阵列位于 x-y 平面)
[x_idx, y_idx] = meshgrid(0:Nx-1, 0:Ny-1);
x_coord = x_idx(:) * dx;   % 所有阵元的 x 坐标 (N x 1)
y_coord = y_idx(:) * dy;   % 所有阵元的 y 坐标 (N x 1)

%%% 信号源参数 (目标角度变为 2D)
% 定义：方位角 phi (x-y平面内，从x轴正向逆时针)
%      俯仰角 theta (从 x-y 平面向上，0 表示在 x-y 平面，90 表示在 z 轴正向)
phi = [30, -60]' * pi/180;      % 信号1: 30°, 信号2: -60°
theta = [20, 45]' * pi/180;     % 信号1: 20°, 信号2: 45°
M = length(phi);               % 信号源数目

%%% 仿真参数
SNR = 10;             % 信噪比(dB)
K = 1000;             % 采样点数

%% 2. 阵列接收信号仿真模拟
%%% 关键：构建 2D 流型矩阵 S (N x M)
% 2D URA 的导向矢量 (Steering Vector)
% a(phi, theta) = exp(1j * k * (x * cos(theta) * cos(phi) + y * cos(theta) * sin(phi)))
S = zeros(N, M);
for m = 1:M
    ph = phi(m);
    th = theta(m);
    % 对应 N 个阵元的相位延迟
    phase_delay = k * (x_coord * cos(th) * cos(ph) + y_coord * cos(th) * sin(ph));
    S(:, m) = exp(1j * phase_delay);
end

% 模拟信号源 (使用复高斯信号)
Alpha = (randn(M, K) + 1j*randn(M, K)) / sqrt(2);
X = S*Alpha;                        % 阵列接收信号
X1 = awgn(X, SNR, 'measured');      % 加载高斯白噪声

%% 3. MUSIC 算法 (核心算法)
%%% 阵列接收信号的协方差矩阵的特征分解
R = X1*X1'/K;    % 阵列接收信号的协方差矩阵 (N x N)
[EV, D] = eig(R);       % 特征值分解
EVA = diag(D);          % 提取特征值
[EVA, I] = sort(EVA, 'descend');   % 降序排序
Q = EV(:, I);           % 特征向量构成的矩阵
Q_n = Q(:, M+1:N);      % 噪声子空间 (N x (N-M))

%% 4. 计算 2D MUSIC 谱估计 (关键：1D 搜索变为 2D 搜索)
%%% 建立 2D 搜索网格
search_phi = linspace(-90, 90, 181);      % 方位角搜索范围
search_theta = linspace(0, 90, 91);       % 俯仰角搜索范围
P_MUSIC = zeros(length(search_theta), length(search_phi));

%%% 遍历 2D 网格计算空间谱
% % 循环方式 (易于理解，但较慢)
% for i = 1:length(search_theta)
%     th = search_theta(i) * pi/180;
%     for j = 1:length(search_phi)
%         ph = search_phi(j) * pi/180;
%         
%         % 构建当前搜索角度的导向矢量 s_search
%         phase_delay = k * (x_coord * cos(th) * cos(ph) + y_coord * cos(th) * sin(ph));
%         s_search = exp(1j * phase_delay);
%         
%         % 计算 MUSIC 谱
%         P_MUSIC(i, j) = 1 / (s_search' * Q_n * Q_n' * s_search);
%     end
% end

%%% 矢量化方式 (速度快)
[Phi_grid, Theta_grid] = meshgrid(search_phi * pi/180, search_theta * pi/180);
Phi_vec = Phi_grid(:);
Theta_vec = Theta_grid(:);

% 构建 N x (NumPhi * NumTheta) 的搜索流型矩阵 S1
A = cos(Theta_vec) .* cos(Phi_vec);
B = cos(Theta_vec) .* sin(Phi_vec);
S1 = exp(1j * k * (x_coord * A' + y_coord * B'));

% 一次性计算所有谱值
P_MUSIC_vec = 1./sum(abs(Q_n'*S1).^2, 1);
P_MUSIC = reshape(P_MUSIC_vec, length(search_theta), length(search_phi));

%%% 转换为 dB
P_MUSIC = abs(P_MUSIC);
P_MUSIC_max = max(P_MUSIC(:));
P_MUSIC_dB = 10*log10(P_MUSIC/P_MUSIC_max);

%% 5. 结果提取与绘图 (关键改动：绘制 2D 谱)
%%% 绘图
figure;
imagesc(search_phi, search_theta, P_MUSIC_dB);
colorbar;
axis xy; % 确保 y 轴 (俯仰角) 方向正确
xlabel('方位角 \phi (deg)');
ylabel('俯仰角 \theta (deg)');
title('2D MUSIC 空间谱 (dB)');
grid on;
hold on;

% 标注真实信号位置
plot(phi*180/pi, theta*180/pi, 'r+', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', '真实方向');
legend;

%%% 提取峰值 (使用 2D 峰值查找)
% 注意：2D 峰值查找比 1D 复杂，可能需要图像处理工具箱
% 这里我们使用 imregionalmax 来查找局部最大值
peaks_mask = imregionalmax(P_MUSIC_dB);
[peak_rows, peak_cols] = find(peaks_mask);

% 提取峰值并排序
peak_vals = P_MUSIC_dB(sub2ind(size(P_MUSIC_dB), peak_rows, peak_cols));
[sorted_peaks, sort_idx] = sort(peak_vals, 'descend');

% 提取前 M 个峰值
num_peaks_to_find = min(M, length(sorted_peaks));
phi_e = zeros(num_peaks_to_find, 1);
theta_e = zeros(num_peaks_to_find, 1);

disp('信号源估计方向 (方位角, 俯仰角):');
for m = 1:num_peaks_to_find
    idx = sort_idx(m);
    r = peak_rows(idx);
    c = peak_cols(idx);
    
    phi_e(m) = search_phi(c);
    theta_e(m) = search_theta(r);
    
    fprintf('  信号 %d: (%0.1f°, %0.1f°)\n', m, phi_e(m), theta_e(m));
    
    % 在图上标注估计的峰值
    plot(phi_e(m), theta_e(m), 'wx', 'MarkerSize', 10, 'LineWidth', 2);
end