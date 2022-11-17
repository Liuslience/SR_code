clc;
clear all;

%% 参数设置
N_frame = 64;            % OFDM符号数
N_FFT = 64;              % 每个符号FFT长度
N_cp = 16;               % 循环前缀长度
N_symbol = 60;           % 每个OFDM符号的长度
N_pilot_number = 4;      % 导频信号的个数
M = 4;K = 2;             % M：调制阶数
sr = 250000;             % 符号速率
br = sr.*K;              % 每个载波的比特速率
N_realization = 1e2;     % 程序实现次数
pilot_loc = [11, 25, 39, 53];% 导频信号的位置
N = 1;                   % 噪声功率

EbN0 = 1:1:15;                      % 设出比特信噪比(dB)
% SNR = EbN0 + 10 * log10(K);         % 由公式推出snr(dB)表达式
SNR = 0:1:30;

BER_CSI_Sim_p = zeros(N_realization, length(SNR));        % 初始化误码率
BER_CSI_theo_p = zeros(N_realization, length(SNR));

BER_LSE_Sim_p = zeros(N_realization, length(SNR));
BER_LSE_theo_p = zeros(N_realization, length(SNR));

BER_CSI_Sim_s = zeros(N_realization, length(SNR));
BER_CSI_theo_s = zeros(N_realization, length(SNR));


%% 衰落参数初始化
PowerdB_d=[0 -8 -17 -21]; % 信道抽头功率特性
PowerdB_b=[0 -8 -17 -21];
Delay=[0 1 2 3];          % 信道时延,示例
Power_d=10.^(PowerdB_d/10);     % 信道抽头功率特性 '线性'
Power_b=10.^(PowerdB_b/10);
Ntap=length(PowerdB_d);       % 信道抽头数
Lch=Delay(end)+1;           % 信道长度

%% generate W matrix
W = zeros(N_FFT, N_FFT);
for p = 1 : N_FFT
    for q = 1 : N_FFT
        W(p, q) = exp(-1i * 2 * pi / N_FFT * (p-1) * (q-1));
    end
end

F_L = W(:, 1 : Lch);
F_p = F_L(pilot_loc, :);

H_est_ls_linear = zeros(N_FFT, N_FFT);

for jj = 1 : N_realization
%% generate x(n), thus modulate them with QPSK
xn = randi([0,1], N_symbol * N_frame, K);
xn_modulated = qpskmod(xn, N_symbol, N_frame);
xn_modulated = reshape(xn_modulated, N_symbol, N_frame);

%% add pilot signal
xn_pilot_signal = ones(1, N_pilot_number);
xn_add_pilot = add_pilot(xn_modulated, xn_pilot_signal, pilot_loc, N_pilot_number, N_symbol, N_frame);

%% generate Rayleigh Channel(direct link)
channel_d = (randn(1,Ntap) + 1j * randn(1,Ntap)).*sqrt(Power_d/2);
hd = zeros(1,Lch);
hd(Delay+1) = channel_d;
hd = hd ./ sqrt(mean(abs(hd) .^ 2)) ./ 2;
Hd = fft([hd,zeros(1,N_FFT-Lch)].');

%% generate Rayleigh Channel(backscatter link)
channel_b = (randn(1,Ntap) + 1j * randn(1,Ntap)).*sqrt(Power_b/2);
hb = zeros(1,Lch);
hb(Delay+1) = channel_b;
hb = hb ./ sqrt(mean(abs(hb) .^ 2)) ./ 2;
Hb = fft([hb,zeros(1,N_FFT-Lch)].');

%% generate c(n)
cn = randi([0, 1], 1, N_frame - 2);
cn = [1, 0, cn] .* 2 - 1;

%% generate the composite Channel
Hb = Hb .* cn;
H = repmat(Hd, 1, N_frame) + Hb; %组合信道的CFR

for i = 1:length(SNR)
%% ifft, then add cp
xn_ifft = ifft(xn_add_pilot) .* sqrt(N_FFT);
xn_add_cp = add_cp(xn_ifft, N_frame, N_FFT, N_cp);

% P = mean(abs(xn_add_cp).^2);

P_t = 10 ^ (SNR(i) / 10) * N / mean(abs(Hd(:, 1) .^ 2));    %calculate the translate power
% P_t = 10 ^ (SNR(i) / 10) * N / sum(abs(Hd(:, 1) .^ 2));

xn_add_cp = xn_add_cp .* sqrt(P_t);
% P1 = mean(abs(xn_add_cp).^2);

%% modulate signal with c(n)
xn_add_cp_1 = reshape(xn_add_cp, N_cp + N_FFT, N_frame);
xn_add_cp_2 = xn_add_cp_1 .* cn;
xn_add_cp_3 = reshape(xn_add_cp_2, 1, (N_cp + N_FFT) * N_frame);

xn_fading_d = conv(xn_add_cp,hd);% pass through the hd
xn_fading_d = xn_fading_d(:, 1:length(xn_add_cp));

xn_fading_b = conv(xn_add_cp_3,hb);% pass through the hb
xn_fading_b = xn_fading_b(:, 1:length(xn_add_cp_3));

xn_fading = xn_fading_d + xn_fading_b;%combine the signal


%% add AWGN to the received signal
    y_received = add_awgn(xn_fading, SNR(i), P_t, Hd);


%% BER Performance with perfect CSI(Primary transmission)
%     N = sqrt(10 .^ (-SNR(i) / 10) * avgPower);
    N = 10 .^ (-SNR(i) / 10) * P_t * mean(abs(Hd(:, 1) .^ 2));     % calculate the N
    BER_CSI_theo_p(jj, i) = 1 / (2 * N_FFT) * sum(qfunc(sqrt(P_t .*abs(Hd + Hb(:, 1)) .^ 2 ./ N)) + qfunc(sqrt(P_t .*abs(Hd - Hb(:, 1)) .^ 2 ./ N)));

%% BER Performance with Estimated CSI(Second transmission)
%     BER_CSI_theo_s(jj, i) = qfunc(sqrt(2 * P_t .* sum(abs(H_b(:, 1).^ 2))  / N));

%% BER Performance with perfect CSI(Primary transmission)
    
    result = 0;
    for ii = 1 : N_FFT
%         a = 2;
        a = F_L(ii, :) * inv(F_p' * F_p) * F_L(ii, :)' + 1;
        result = result + 1 / (2 * N_FFT)...
            * (qfunc(sqrt(P_t .*abs(Hd(ii) + Hb(ii, 1)) .^ 2 ./ N / abs(a)))...
            + qfunc(sqrt(P_t .*abs(Hd(ii) - Hb(ii, 1)) .^ 2 ./ N / abs(a))));
    end
    BER_LSE_theo_p(jj, i) = result;


    y_received = reshape(y_received, N_FFT + N_cp, N_frame);
    
    %%%%%%去掉cp,并进行fft%%%%%%
    y_remove_cp = fft(remove_cp(y_received, N_frame, N_FFT, N_cp)) ./ 8;

    %%%%%%信道估计%%%%%%
%     for j = 1
%         if j ==1
%             H_est_ls_linear = LS_CE(y_remove_cp,xn_pilot_signal.',pilot_loc,N_FFT,N_pilot_number,'linear');
%         else 
%             H_est_ls_spline = LS_CE(y_remove_cp,xn_pilot_signal.',pilot_loc,N_FFT,N_pilot_number,'spline');
%         end
%     end

    H_est_LSE = primary_H_est(y_remove_cp, pilot_loc, P_t, F_p, N_FFT, N_frame, F_L);

%     H_est_ls_linear2 = estimate_1(xn_add_pilot, y_remove_cp, P_t, F_p, N_FFT, Lch);

    %%%%%%已知CSI下的信道均衡%%%%%%
    y_equalization_CSI = y_remove_cp ./ H;

    %%%%%%LS_linear下的信道均衡%%%%%%
    y_equalization_LS_linear = y_remove_cp./H_est_LSE;

    %%%%%%去除导频信号%%%%%%
    r_remove_pilot_CSI = remove_pilot(y_equalization_CSI, pilot_loc, N_frame, N_symbol);
    r_remove_pilot_LS_linear = remove_pilot(y_equalization_LS_linear, pilot_loc, N_frame, N_symbol);

    %%%%%%QPSK解调%%%%%%
    r_remove_pilot_CSI = reshape(r_remove_pilot_CSI, 1, N_frame * N_symbol);
    y_demodulated_CSI = qpskdemod(r_remove_pilot_CSI, N_symbol, N_frame);

    r_remove_pilot_LS_linear = reshape(r_remove_pilot_LS_linear, 1, N_frame * N_symbol);
    y_demodulated_LS_linear = qpskdemod(r_remove_pilot_LS_linear, N_symbol, N_frame);

    %%%%%%误比特率计算%%%%%%
    Ntb = N_symbol * N_frame * K;                              %仿真的总比特数

%     Neb_CSI = sum(sum(int2bit(y_demodulated_CSI,K) ~= int2bit(xn,K))); %转为2进制，计算具体有几bit错误
    Neb_CSI_p = sum(sum(y_demodulated_CSI ~= xn));
    BER_CSI_Sim_p(jj, i) = Neb_CSI_p / Ntb;

    Neb_lslinear = sum(sum(y_demodulated_LS_linear ~= xn));       %转为2进制，计算具体有几bit错误
    BER_LSE_Sim_p(jj, i) = Neb_lslinear / Ntb;

%     c_estimated = arg_min(H, H_d, H_b, c_n, N_frame);
%     Neb_CSI_s = sum(sum(c_estimated ~= c_n));
%     BER_CSI_Sim_s(i) = Neb_CSI_s / Ntb;
end
jj
end

BER_CSI_Sim_p = sum(BER_CSI_Sim_p, 1) ./ N_realization;
BER_CSI_theo_p = sum(BER_CSI_theo_p, 1) ./ N_realization;

BER_LSE_Sim_p = sum(BER_LSE_Sim_p, 1) ./ N_realization;
BER_LSE_theo_p = sum(BER_LSE_theo_p, 1) ./ N_realization;

BER_CSI_Sim_s = sum(BER_CSI_Sim_s, 1) ./ N_realization;
BER_CSI_theo_s = sum(BER_CSI_theo_s, 1) ./ N_realization;

%% 画图
figure(1);
semilogy(SNR, BER_CSI_Sim_p, 'r-o');
hold on;
semilogy(SNR, BER_CSI_theo_p, 'k-^');
hold on;
semilogy(SNR, BER_LSE_Sim_p, 'b-*');
hold on;
semilogy(SNR, BER_LSE_theo_p, 'c-v');
hold on;
xlabel('SNR(比特信噪比)');ylabel('BER(误比特率)');
title('多径衰落信道下误码率仿真曲线');
lgd = legend('sim, with perfect CSI', 'theo, with perfect CSI', 'sim, with lslinear CSI', 'theo, with lsspline CSI');
position = lgd.Location;
lgd.Location = 'southwest';
grid on;






