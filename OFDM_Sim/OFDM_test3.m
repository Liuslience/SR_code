clc;
clear all;

%% 参数设置
N_frame = 1024;         %OFDM符号数
N_FFT = 64;             %每个符号FFT长度
N_cp = 16;              %循环前缀长度
N_symbol = 48;          %每个OFDM符号的长度
M = 16;K = 4;           %M：调制阶数
sr = 250000;            %符号速率
br = sr.*K;             %每个载波的比特速率

EbN0 = 0:1:30;                      % 设出比特信噪比(dB)
SNR = EbN0 + 10 * log10(K);         % 由公式推出snr(dB)表达式

BER_CSI = zeros(1, length(SNR));        % 初始化误码率
BER_lslinear = zeros(1, length(SNR));
BER_lsspline = zeros(1, length(SNR));
BER_MMSE = zeros(1, length(SNR));


%% 衰落参数初始化
PowerdB=[0 -8 -17 -21 -25]; % 信道抽头功率特性
% PowerdB=[0 -2 -3 -4 -5];    % 信道抽头功率特性
Delay=[0 3 5 6 8];          % 信道时延,示例
% Delay=[0 3 5 56 78];        % 信道时延
Power=10.^(PowerdB/10);     % 信道抽头功率特性 '线性'
Ntap=length(PowerdB);       % 信道抽头数
Lch=Delay(end)+1;           % 信道长度

%% 生成需要传输的信息，并进行QAM调制
xn = randi([0,15],1,N_symbol * N_frame);                       %生成随机信息
xn_modulated = qammod(xn, M,'gray','UnitAveragePower',true);%进行QAM调制
% xn_modulated1 = xn_modulated;
avgPower = mean(abs(xn_modulated).^2);
xn_modulated = reshape(xn_modulated, N_symbol, N_frame);

%% 插入导频信号
N_pilot_interval = 4;   %导频信号间隔
start_pilot = 1;        %导频信号起始载波位置
N_pilot_number = fix(N_FFT / N_pilot_interval);%计算导频信号的个数
xn_add_pilot = add_pilot(xn_modulated, N_pilot_interval, N_pilot_number, start_pilot, N_symbol, N_frame);
xn_pilot_signal = ones(1, N_pilot_number);
pilot_loc = zeros(1, N_pilot_number);%生成导频信号的位置信息
for i = 1 : N_pilot_number
    pilot_loc(1, i) = 1 + (i-1) * N_pilot_interval;
end

%% 进行ifft,并添加cp
xn_ifft = ifft(xn_add_pilot);
xn_add_cp = add_cp(xn_ifft, N_frame, N_FFT, N_cp);

%% 瑞利衰落信道
channel = (randn(1,Ntap) + 1j * randn(1,Ntap)).*sqrt(Power/2);
h = zeros(1,Lch);
h(Delay+1) = channel;
xn_fading = conv(xn_add_cp,h);
xn_fading = xn_fading(:, 1:length(xn_add_cp));
% figure(1);
% stem(0:length(h_frenqency)-1, abs(h_frenqency), '.');
% title('多径时延信道的频域图像');
% grid on;

%% 使信号通过AWGN
for i = 1:length(SNR)
    
    y_received = awgn(xn_fading,SNR(i),'measured');

    y_received = reshape(y_received, N_FFT + N_cp, N_frame);
    
    %%%%%%去掉cp,并进行fft%%%%%%
    y_remove_cp = fft(remove_cp(y_received, N_frame, N_FFT, N_cp));

    %%%%%%信道估计%%%%%%
    for j = 1:3
        if j ==1
            H_est_ls_linear = LS_CE(y_remove_cp,xn_pilot_signal.',pilot_loc,N_FFT,N_pilot_interval,'linear');
        elseif j ==2
            H_est_ls_spline = LS_CE(y_remove_cp,xn_pilot_signal.',pilot_loc,N_FFT,N_pilot_interval,'spline');
        else
            H_est_MMSE = MMSE_CE(y_remove_cp,xn_pilot_signal.',pilot_loc,N_FFT,N_pilot_interval,h,SNR(i));
        end
    end

    %%%%%%已知CSI下的信道均衡%%%%%%
    H = fft([h,zeros(1,N_FFT-Lch)].');
    y_equalization_CSI = y_remove_cp./H;

    %%%%%%LS_linear下的信道均衡%%%%%%
    y_equalization_LS_linear = y_remove_cp./H_est_ls_linear.';

    %%%%%%LS_spline下的信道均衡%%%%%%
    y_equalization_LS_spline = y_remove_cp./H_est_ls_spline.';

    %%%%%%MMSE下的信道均衡%%%%%%
    y_equalization_MMSE = y_remove_cp./H_est_MMSE.';
    
    %%%%%%去除导频信号%%%%%%
    r_remove_pilot_CSI = remove_pilot(y_equalization_CSI, N_pilot_interval, N_pilot_number, start_pilot, N_FFT, N_frame);
    r_remove_pilot_LS_linear = remove_pilot(y_equalization_LS_linear, N_pilot_interval, N_pilot_number, start_pilot, N_FFT, N_frame);
    r_remove_pilot_LS_spline = remove_pilot(y_equalization_LS_spline, N_pilot_interval, N_pilot_number, start_pilot, N_FFT, N_frame);
    r_remove_pilot_MMSE = remove_pilot(y_equalization_MMSE, N_pilot_interval, N_pilot_number, start_pilot, N_FFT, N_frame);

    %%%%%%QAM解调%%%%%%
    r_remove_pilot_CSI = reshape(r_remove_pilot_CSI, 1, N_frame * N_symbol);
    y_demodulated_CSI = qamdemod(r_remove_pilot_CSI, M,'gray','UnitAveragePower',true);

    r_remove_pilot_LS_linear = reshape(r_remove_pilot_LS_linear, 1, N_frame * N_symbol);
    y_demodulated_LS_linear = qamdemod(r_remove_pilot_LS_linear, M,'gray','UnitAveragePower',true);

    r_remove_pilot_LS_spline = reshape(r_remove_pilot_LS_spline, 1, N_frame * N_symbol);
    y_demodulated_LS_spline = qamdemod(r_remove_pilot_LS_spline, M,'gray','UnitAveragePower',true);

    r_remove_pilot_MMSE = reshape(r_remove_pilot_MMSE, 1, N_frame * N_symbol);
    y_demodulated_MMSE = qamdemod(r_remove_pilot_MMSE, M,'gray','UnitAveragePower',true);

    %%%%%%误比特率计算%%%%%%
    Ntb = N_symbol * N_frame * K;                              %仿真的总比特数

    Neb_CSI = sum(sum(int2bit(y_demodulated_CSI,K) ~= int2bit(xn,K))); %转为2进制，计算具体有几bit错误
    BER_CSI(i) = Neb_CSI / Ntb;

    Neb_lslinear = sum(sum(int2bit(y_demodulated_LS_linear,K) ~= int2bit(xn,K))); %转为2进制，计算具体有几bit错误
    BER_lslinear(i) = Neb_lslinear / Ntb;

    Neb_lsspline = sum(sum(int2bit(y_demodulated_LS_spline,K) ~= int2bit(xn,K))); %转为2进制，计算具体有几bit错误
    BER_lsspline(i) = Neb_lsspline / Ntb;

    Neb_MMSE = sum(sum(int2bit(y_demodulated_MMSE,K) ~= int2bit(xn,K))); %转为2进制，计算具体有几bit错误
    BER_MMSE(i) = Neb_MMSE / Ntb;
end

%% 画图
figure(2);
semilogy(BER_CSI, 'r-o');
hold on;
semilogy(BER_lslinear, 'b-*');
hold on;
semilogy(BER_lsspline, 'g-+');
hold on;
semilogy(BER_MMSE, 'c-o');
xlabel('EbN0(比特信噪比)');ylabel('BER(误比特率)');
title('多径衰落信道下误码率仿真曲线');
legend('已知CSI的曲线', 'lslinear曲线', 'lsspline曲线', 'MMSE曲线');
grid on;







