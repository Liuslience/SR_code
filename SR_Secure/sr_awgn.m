%% parameter settings
SNR = 0:1:30;
N_signal = 64;%number of signal
N_block = 64;
K = 2;
number_realization = 1;

BER_CSI_Sim_p = zeros(number_realization, length(SNR));        % 初始化误码率
BER_CSI_theo_p = zeros(number_realization, length(SNR));

%% 衰落参数初始化
PowerdB_d=-8; % 信道抽头功率特性
PowerdB_b=-8;
Power_d=10.^(PowerdB_d/10);     % 信道抽头功率特性 '线性'
Power_b=10.^(PowerdB_b/10);
Ntap=length(PowerdB_d);       % 信道抽头数
channel_d = (randn(1,Ntap) + 1j * randn(1,Ntap)).*sqrt(Power_d/2);%直接链路的瑞利衰落信道
channel_b = (randn(1,Ntap) + 1j * randn(1,Ntap)).*sqrt(Power_b/2);%反射链路的瑞利衰落信道

for jj = 1 : number_realization
    
xn = randi([0,1], N_signal * N_block, K);                       %generate random signal
xn_modulated = qpskmod(xn, N_signal * N_block);
P_t = mean(abs(xn_modulated).^2);

%% 生成c(n)信号
c_n = randi([0, 1], 1, N_block);
c_n = c_n .* 2 - 1;

%% 用c（n）调制反射信道的信号
xn_modulated = reshape(xn_modulated, N_signal, N_block);
xn_modulated1 = xn_modulated .* c_n;

%     xn_modulated = xn_modulated .* 
    
for i = 1:length(SNR)
    y_received = awgn(xn_modulated, SNR(i), 'measured');
    N = 10 .^ (-SNR(i) / 10) * P_t;
%         BER_CSI_theo_p(jj, i) = 1 / (N_FFT) * sum(qfunc(sqrt(P_t .* abs(H_d) .^ 2 ./ (P_t .* abs(H_b(:, 1)) .^ 2 + N))));
%         y_equalization_CSI = y_received ./ 
    
    
    end
end




















