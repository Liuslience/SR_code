clc;

%%%%%%参数设置%%%%%%
N_frame = 1024;         %OFDM符号数
N_FFT = 64;             %每个符号FFT长度
N_cp = 16;              %循环前缀长度
N_symbol = N_FFT + N_cp;%每个OFDM符号的长度
M = 16;K = 4;           %M：调制阶数
sr = 250000;            %符号速率
br = sr.*K;             %每个载波的比特速率

EbN0 = 0:1:50;                      % 设出比特信噪比(dB)
snr = EbN0 + 10 * log10(K);         % 由公式推出snr(dB)表达式
BER = zeros(1, length(snr));          % 初始化误码率

%%%%%%Fading initialization%%%%%%
PowerdB=[0 -8 -17 -21 -25]; % 信道抽头功率特性
% PowerdB=[0 -2 -3 -4 -5]; % 信道抽头功率特性
Delay=[0 3 5 6 8];          % 信道时延,示例
% Delay=[0 3 5 56 78];        % 信道时延
Power=10.^(PowerdB/10);     % 信道抽头功率特性 '线性'
Ntap=length(PowerdB);       % 信道抽头数
Lch=Delay(end)+1;           % 信道长度

%%%%%%生成需要传输的信息，并进行QAM调制%%%%%%
xn = randi([0,15],1,N_FFT * N_frame);                       %生成随机信息
xn_modulated = qammod(xn, M,'gray','UnitAveragePower',true);%进行QAM调制
avgPower = mean(abs(xn_modulated).^2);

%%%%%%进行ifft,并添加cp%%%%%%
xn_modulated1 = reshape(xn_modulated, N_FFT, N_frame);
x_ifft = ifft(xn_modulated1);
x_add_cp = add_cp(x_ifft, N_frame, N_FFT, N_cp);

%%%%%%瑞利衰落信道%%%%%%
channel = (randn(1,Ntap) + 1j * randn(1,Ntap)).*sqrt(Power/2);
h = zeros(1,Lch);
h(Delay+1) = channel;
x_fading = conv(x_add_cp,h);
x_fading1 = x_fading(:, 1:length(x_add_cp));
h_frenqency = fft(h);
% figure(1);
% stem(0:length(h_frenqency)-1, abs(h_frenqency), '.');
% title('多径时延信道的频域图像');
% grid on;

% H = fft([h,zeros(1,N_FFT-Lch)].');
% y_equalization = zeros(N_FFT, N_frame);

%%%%%%使信号通过AWGN%%%%%%
for i = 1:length(snr)
    
    y_received = awgn(x_fading1,snr(i),'measured');

    y_received1 = reshape(y_received, N_FFT + N_cp, N_frame);
    
    %%%%%%去掉cp,并进行fft%%%%%%
    y_remove_cp = fft(remove_cp(y_received1, N_frame, N_FFT, N_cp));

    %%%%%%信道均衡%%%%%%
    H = fft([h,zeros(1,N_FFT-Lch)].');
    y_equalization = zeros(N_FFT, N_frame);
    for number = 1:N_frame
        y_equalization(:, number) = y_remove_cp(:, number)./H;
    end
    
    %%%%%%QAM解调%%%%%%
    y_equal_1 = reshape(y_equalization, 1, N_frame * N_FFT);
    y_demodulated = qamdemod(y_equal_1, M,'gray','UnitAveragePower',true);

    %%%%%%误比特率计算%%%%%%
    Neb = sum(sum(int2bit(y_demodulated,K) ~= int2bit(xn,K))); %转为2进制，计算具体有几bit错误
    Ntb = N_FFT * N_frame * K;                             %仿真的总比特数
    BER(i) = Neb / Ntb;
end

figure(2);
semilogy(BER, 'r-o');
xlabel('EbN0(比特信噪比)');ylabel('BER(误比特率)');
title('多径衰落信道下误码率仿真曲线');
grid on;










