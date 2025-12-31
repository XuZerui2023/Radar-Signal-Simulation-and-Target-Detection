% written by Patrick_xie
% https://zhuanlan.zhihu.com/p/457093041
% 本程序用于学习雷达目标信号回波生成仿真
% 本程序用于单脉冲单通道多目标雷达回波信号仿真处理
% 可以学习LFM信号如何通过脉冲压缩和MTD处理来同时测量目标的距离和速度

clc;
clear;
close all;
%% 1. 信号与系统参数设置
f = 10e6;         % 调频信号基础频率 10MHz
fs = 100e6;       % 采样频率
ts = 1 / fs;      % 采样周期
B = 10^7;         % 调频范围（带宽）
T = 10^(-5);      % 脉冲宽度
u = B / T;        % 调频率，斜率
c = 3e8;          % 光速
PRT = 500e-6;     % 脉冲重复周期

%% 2. 生成LFM信号与回波
r1 = 10000;             % 目标1距离
tao1 = 2*r1/c;          % 计算目标1时延（电磁波往返时间）
n1 = round(tao1 / ts);  % 把目标1时延转换成采样点数
r2 = 30000;             % 计算目标2时延（电磁波往返时间）
tao2 = 2*r2/c;          % 把目标2时延转换成采样点数
n2 = round(tao2 / ts);  % 目标2时延转换为采样点数
r3 = 35000;             % 目标3距离
tao3 = 2*r3/c;          % 计算目标3时延（电磁波往返时间）
n3 = round(tao3 / ts);  % 把目标3时延转换成采样点数
r4 = 50000;
tao4 = 2*r4/c;
n4 = round(tao4 / ts);
NN = PRT / ts;          % 一个完整PRT内的采样点数

t1 = 0 : ts : T-ts;
y = sin(2*pi*(f*t1+0.5*u*t1.^2));
N = length(y);

y1 = [zeros(1, n1), y, zeros(1, NN-n1-N)];  % 通过前后补零的方式，模拟信号的时延，将三个目标回波信号长度补到NN
y2 = [zeros(1, n2), y, zeros(1, NN-n2-N)];
y3 = [zeros(1, n3), y, zeros(1, NN-n3-N)];
y4 = [zeros(1, n4), y, zeros(1, NN-n4-N)];
xt = y1+y2+y3+y4;                              % 将三个回波信号叠加得到总的回波信号

% --- 在信号上添加高斯白噪声 ---
noise_power = 1e-2; % 可调噪声功率
noise = sqrt(noise_power/2) * randn(size(xt));
xt = xt + noise;




%% 3. 信号时域和频域分析
figure;
subplot(2, 1, 1);
plot((0:length(xt)-1)*ts*c/2/1000, xt);  % 注意横轴坐标将回波信号的采样点索引转换为对应的实际距离，单位：千米
xlabel('距离/km');
subtitle('雷达回波信号的时域图');
xtfft = abs(fft(xt, 50000));
subplot(2, 1, 2);
fx = (0:length(xt)/2-1)*fs/length(xt);   % 计算公式 `k * (fs / N)` 是**傅里叶变换频率分辨率**的通用公式，其中 `k` 是频率索引，`fs` 是采样频率，`N` 是信号长度。
plot(fx/1e6, xtfft(1:length(xt)/2));
xlabel('频率/MHz');
subtitle('雷达回波信号的频域图');

figure;
xtfft_shift = fftshift(xtfft);
fx_shift = (-length(xt)/2 : length(xt)/2-1) * fs/length(xt);
plot(fx_shift/1e6,xtfft_shift);
xlabel('频率/MHz');
title('频谱搬移后雷达回波信号的频域图');


%% 4. 混频解调与滤波
% 对回波信号采样
fs1 = 100e6;
ts1 = 1/fs1;
t2 = 0 : ts1 : (length(xt)-1)*ts1;
xrt = xt .* sin(2*pi*f*t2);             % 混频（下变频）。混频后，高频信号被转移到基带，其频率信息可以被有效提取。

fx1 = (0:length(xt)/2-1)*fs1/length(xt);
figure;
subplot(2, 1, 1);
plot((0:length(xt)-1)*ts*c/2/1000, xrt);
xlabel('距离/km');
title('解调后雷达回波信号的时域图');
xrtfft = abs(fft(xrt, 50000));
subplot(2, 1, 2);
fx1=(0:length(xrt)/2-1)*fs1/length(xrt);
plot(fx1/1e6, xrtfft(1:length(xrt)/2));
xlabel('频率/MHz');
title('解调后雷达回波信号的频域图');

load FIR.mat                      % 滤波，除去混频所产生的高频分量，只保留基带有用的信号
firxrt = filter(Num, 1,  xrt);    
firfft = abs(fft(firxrt, 50000));

figure;
plot(fx1/1e6, firfft(1:length(firxrt)/2));
xlabel('频率/MHz');
title('解调信号滤波后的信号频域图')
xrtdown = downsample(firxrt, 4);    % 下采样，下采样可以减小数据量，降低后续处理的计算负担。
value0 = abs(fft(xrtdown));           % FFT分析

figure;
fx1 = (0:length(xrtdown)/2-1)*(fs1/4)/length(xrtdown);
%plot((fx1(1:length(value0)))/1e6, value0);
plot(fx1/1e6, value0(1:length(fx1)));
xlabel('频率/MHz');
title('解调滤波后对信号抽取（下采样）的频域图')

%% 5. 脉冲压缩（匹配滤波）
T = 10^(-5);
u = B / T;   % 调频斜率
c = 3e8;
fs = 25e6;
ts = 1 / fs;
t1 = 0 : ts : T-ts;
hdt = sin(2*pi*(0*t1 + 0.5*u*t1.^2));  % 匹配滤波器

figure;
plot(hdt);
replica = xrtdown;
% y = fliplr(hdt);        % 时域卷积的匹配滤波器 
% out = conv(replica, y);
replica1= [replica,zeros(1,16384-length(replica))];
y1 = [hdt,zeros(1,16384-length(hdt))];
rfft = fft(replica1);
yfft = fft(y1);
out = abs(ifft((rfft.*conj(yfft))));   % 频域相乘
title('脉冲压缩后信号频域图');

figure;
t = (0:length(xrtdown)-1)*ts*c/2/1000;
plot(t, out(1:length(t)));
xlabel('距离/km');
title('脉冲压缩后信号时域图');