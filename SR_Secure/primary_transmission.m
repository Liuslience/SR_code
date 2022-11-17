clc;
clear all;

%% 参数设置
N_frame = 64;         %OFDM符号数
N_FFT = 64;             %每个符号FFT长度
N_cp = 16;              %循环前缀长度
N_symbol = 60;          %每个OFDM符号的长度
N_pilot_number = 4;     %导频信号的个数
M = 4;K = 2;           %M：调制阶数
sr = 250000;            %符号速率
br = sr.*K;             %每个载波的比特速率

EbN0 = 1:1:15;                      % 设出比特信噪比(dB)
% SNR = EbN0 + 10 * log10(K);         % 由公式推出snr(dB)表达式
SNR = 0:1:30;

BER_CSI_Sim_p = zeros(1e6, length(SNR));        % 初始化误码率
BER_CSI_theo_p = zeros(1e6, length(SNR));

BER_lslinear_Sim_p = zeros(1e6, length(SNR));
BER_lslinear_theo_p = zeros(1e6, length(SNR));

BER_lsspline_Sim_p = zeros(1e6, length(SNR));

BER_CSI_Sim_s = zeros(1e6, length(SNR));
BER_CSI_theo_s = zeros(1e6, length(SNR));

BER_MMSE = zeros(1, length(SNR));


%% 衰落参数初始化
PowerdB_d=[0 -8 -17 -21]; % 信道抽头功率特性
% PowerdB_b=[0 -3 -18 -23];
PowerdB_b=[0 -8 -17 -21];
% PowerdB=[0 -2 -3 -4];    % 信道抽头功率特性
Delay=[0 1 2 3];          % 信道时延,示例
% Delay=[0 3 5 56];        % 信道时延
Power_d=10.^(PowerdB_d/10);     % 信道抽头功率特性 '线性'
Power_b=10.^(PowerdB_b/10);
Ntap=length(PowerdB_d);       % 信道抽头数
Lch=Delay(end)+1;           % 信道长度

number_realization = 1e2;

for jj = 1 : number_realization
%% 生成需要传输的信息，并进行QPSK调制
xn = randi([0,1], N_symbol * N_frame, K);                       %生成随机信息
xn_modulated = qpskmod(xn, N_symbol, N_frame);
% avgPower = mean(abs(xn_modulated).^2);
% xn_modulated1 = xn_modulated;
xn_modulated = reshape(xn_modulated, N_symbol, N_frame);

%% 插入导频信号
xn_pilot_signal = ones(1, N_pilot_number);
pilot_loc = [11, 25, 39, 53];
xn_add_pilot = add_pilot(xn_modulated, xn_pilot_signal, pilot_loc, N_pilot_number, N_symbol, N_frame);


%% 进行ifft,并添加cp
xn_ifft = ifft(xn_add_pilot) .* 8;
% avgPower_1 = mean(abs(xn_ifft).^2);
xn_add_cp = add_cp(xn_ifft, N_frame, N_FFT, N_cp);

P_t = mean(abs(xn_add_cp).^2);%计算发射功率



%% 瑞利衰落信道(直射信道)
channel_d = (randn(1,Ntap) + 1j * randn(1,Ntap)).*sqrt(Power_d/2);
h_d = zeros(1,Lch);
h_d(Delay+1) = channel_d;
h_d = h_d ./ sqrt(mean(abs(h_d) .^ 2)) ./ 2;
H_d = fft([h_d,zeros(1,N_FFT-Lch)].');

xn_fading_d = conv(xn_add_cp,h_d);%通过H_d信道
xn_fading_d = xn_fading_d(:, 1:length(xn_add_cp));

%% 生成c(n)信号
c_n = randi([0, 1], 1, N_frame - 2);
c_n = [1, 0, c_n] .* 2 - 1;

%% 用c（n）调制反射信道的信号
xn_add_cp_1 = reshape(xn_add_cp, N_cp + N_FFT, N_frame);
xn_add_cp_2 = xn_add_cp_1 .* c_n;
xn_add_cp_3 = reshape(xn_add_cp_2, 1, (N_cp + N_FFT) * N_frame);

%% 瑞利衰落信道(反射信道)
channel_b = (randn(1,Ntap) + 1j * randn(1,Ntap)).*sqrt(Power_b/2);
h_b = zeros(1,Lch);
h_b(Delay+1) = channel_b;
h_b = h_b ./ sqrt(mean(abs(h_b) .^ 2)) ./ 2;
H_b = fft([h_b,zeros(1,N_FFT-Lch)].');

xn_fading_b = conv(xn_add_cp_3,h_b);%通过H_b信道
xn_fading_b = xn_fading_b(:, 1:length(xn_add_cp_3));

xn_fading = xn_fading_d + xn_fading_b;%合并信号

%% 生成组合信道
H_b = H_b .* c_n;
H = repmat(H_d, 1, N_frame) + H_b; %组合信道的CFR

%% 生成W矩阵
W = zeros(N_FFT, N_FFT);
for p = 1 : N_FFT
    for q = 1 : N_FFT
        W(p, q) = exp(-1i * 2 * pi / N_FFT * (p-1) * (q-1));
    end
end

F_L = W(:, 1 : Lch);
F_p = F_L(pilot_loc, :);

% a = F_L(1, :) * inv(F_p' * F_p) * F_L(1, :)';

% c_estimated = arg_min(H, H_d, H_b, c_n, N_frame);
H_est_ls_linear = zeros(N_FFT, N_FFT);

for i = 1:length(SNR)
%% 使信号通过AWGN
    y_received = add_awgn(xn_fading, SNR(i), P_t, H_d);


%% BER Performance with perfect CSI(Primary transmission)
%     N = sqrt(10 .^ (-SNR(i) / 10) * avgPower);
    N = 10 .^ (-SNR(i) / 10) * P_t;     %计算噪声功率
    BER_CSI_theo_p(jj, i) = 1 / (N_FFT) * sum(qfunc(sqrt(P_t .* abs(H_d) .^ 2 ./ (P_t .* abs(H_b(:, 1)) .^ 2 + N))));

%% BER Performance with Estimated CSI(Second transmission)
%     BER_CSI_theo_s(jj, i) = qfunc(sqrt(2 * P_t .* sum(abs(H_b(:, 1).^ 2))  / N));

%% BER Performance with perfect CSI(Primary transmission)
    
    result = 0;
    for ii = 1 : N_FFT
        a = F_L(ii, :) * inv(F_p' * F_p) * F_L(ii, :)' + 1;
        result = result + 1 / (2 * N_FFT)...
            * (qfunc(sqrt(P_t .*abs(H_d(ii) + H_b(ii, 1)) .^ 2 ./ N / abs(a)))...
            + qfunc(sqrt(P_t .*abs(H_d(ii) - H_b(ii, 1)) .^ 2 ./ N / abs(a))));
    end
    BER_lslinear_theo_p(jj, i) = result;


    y_received = reshape(y_received, N_FFT + N_cp, N_frame);
    
    %%%%%%去掉cp,并进行fft%%%%%%
    y_remove_cp = fft(remove_cp(y_received, N_frame, N_FFT, N_cp)) ./ 8;

    %%%%%%信道估计%%%%%%
    for j = 1
        if j ==1
            H_est_ls_linear = LS_CE(y_remove_cp,xn_pilot_signal.',pilot_loc,N_FFT,N_pilot_number,'linear');
        elseif j ==2
            H_est_ls_spline = LS_CE(y_remove_cp,xn_pilot_signal.',pilot_loc,N_FFT,N_pilot_number,'spline');
        else
            H_est_MMSE = MMSE_CE(y_remove_cp,xn_pilot_signal.',pilot_loc,N_FFT,N_pilot_interval,h,SNR(i));
        end
    end

    for iii = 1 : N_FFT
        h_est_ls(1:4, iii) = inv(P_t * F_p' * F_p) * sqrt(P_t) * F_p' * y_remove_cp(pilot_loc, iii);
    end
    h_est_ls(5:64, :) = 0;
    H_est_ls_linear2 = fft(h_est_ls);

%     H_est_ls_linear2 = estimate_1(xn_add_pilot, y_remove_cp, P_t, F_p, N_FFT, Lch);


%     H_est_ls_linear1 = ifft(H_est_ls_linear.');
%     H_est_ls_linear1(5:64, :) = 0;
%     H_est_ls_linear2 = fft(H_est_ls_linear1);

    %%%%%%已知CSI下的信道均衡%%%%%%
    y_equalization_CSI = y_remove_cp ./ repmat(H_d, 1, N_frame);

    %%%%%%LS_linear下的信道均衡%%%%%%
    y_equalization_LS_linear = y_remove_cp./H_est_ls_linear2;

    %%%%%%LS_spline下的信道均衡%%%%%%
%     y_equalization_LS_spline = y_remove_cp./H_est_ls_spline.';

%     %%%%%%MMSE下的信道均衡%%%%%%
%     y_equalization_MMSE = y_remove_cp./H_est_MMSE.';
    
    %%%%%%去除导频信号%%%%%%
    r_remove_pilot_CSI = remove_pilot(y_equalization_CSI, pilot_loc, N_frame, N_symbol);
    r_remove_pilot_LS_linear = remove_pilot(y_equalization_LS_linear, pilot_loc, N_frame, N_symbol);
%     r_remove_pilot_LS_spline = remove_pilot(y_equalization_LS_spline, pilot_loc, N_frame, N_symbol);
%     r_remove_pilot_MMSE = remove_pilot(y_equalization_MMSE, N_pilot_interval, N_pilot_number, start_pilot, N_FFT, N_frame);

    %%%%%%QPSK解调%%%%%%
    r_remove_pilot_CSI = reshape(r_remove_pilot_CSI, 1, N_frame * N_symbol);
    y_demodulated_CSI = qpskdemod(r_remove_pilot_CSI, N_symbol, N_frame);

    r_remove_pilot_LS_linear = reshape(r_remove_pilot_LS_linear, 1, N_frame * N_symbol);
    y_demodulated_LS_linear = qpskdemod(r_remove_pilot_LS_linear, N_symbol, N_frame);

%     r_remove_pilot_LS_spline = reshape(r_remove_pilot_LS_spline, 1, N_frame * N_symbol);
%     y_demodulated_LS_spline = qpskdemod(r_remove_pilot_LS_spline, N_symbol, N_frame);

%     r_remove_pilot_MMSE = reshape(r_remove_pilot_MMSE, 1, N_frame * N_symbol);
%     y_demodulated_MMSE = qamdemod(r_remove_pilot_MMSE, M,'gray','UnitAveragePower',true);

    %%%%%%误比特率计算%%%%%%
    Ntb = N_symbol * N_frame * K;                              %仿真的总比特数

%     Neb_CSI = sum(sum(int2bit(y_demodulated_CSI,K) ~= int2bit(xn,K))); %转为2进制，计算具体有几bit错误
    Neb_CSI_p = sum(sum(y_demodulated_CSI ~= xn));
    BER_CSI_Sim_p(jj, i) = Neb_CSI_p / Ntb;

    Neb_lslinear = sum(sum(y_demodulated_LS_linear ~= xn));       %转为2进制，计算具体有几bit错误
    BER_lslinear_Sim_p(jj, i) = Neb_lslinear / Ntb;

%     Neb_lsspline = sum(sum(y_demodulated_LS_spline ~= xn)); %转为2进制，计算具体有几bit错误
%     BER_lsspline_Sim_p(i) = BER_lsspline_Sim_p(i) + Neb_lsspline / Ntb;

%     c_estimated = arg_min(H, H_d, H_b, c_n, N_frame);
%     Neb_CSI_s = sum(sum(c_estimated ~= c_n));
%     BER_CSI_Sim_s(i) = Neb_CSI_s / Ntb;

%     Neb_MMSE = sum(sum(int2bit(y_demodulated_MMSE,K) ~= int2bit(xn,K))); %转为2进制，计算具体有几bit错误
%     BER_MMSE(i) = Neb_MMSE / Ntb;
end
jj
end

BER_CSI_Sim_p = sum(BER_CSI_Sim_p, 1) ./ number_realization;
BER_CSI_theo_p = sum(BER_CSI_theo_p, 1) ./ number_realization;
BER_lslinear_Sim_p = sum(BER_lslinear_Sim_p, 1) ./ number_realization;
BER_lslinear_theo_p = sum(BER_lslinear_theo_p, 1) ./ number_realization;

%% 画图
figure(2);
semilogy(SNR, BER_CSI_Sim_p, 'r-o');
hold on;
semilogy(SNR, BER_CSI_theo_p, 'k-^');
hold on;
% semilogy(SNR, BER_lslinear_Sim_p, 'b-*');
% hold on;
% semilogy(SNR, BER_lslinear_theo_p, 'c-v');
% hold on;
% semilogy(SNR, BER_lsspline_Sim_p, 'g-+');
% hold on;
% semilogy(BER_MMSE, 'c-o');
xlabel('SNR(比特信噪比)');ylabel('BER(误比特率)');
title('多径衰落信道下误码率仿真曲线');
% lgd = legend('sim, with perfect CSI', 'theo, with perfect CSI', 'sim, with lslinear CSI', 'theo, with lsspline CSI');
lgd = legend('sim, with perfect CSI', 'theo, with perfect CSI');
position = lgd.Location;
lgd.Location = 'southwest';
grid on;






