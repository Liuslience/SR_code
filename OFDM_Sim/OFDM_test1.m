clc;

%%%%%%参数设置%%%%%%
N_frame = 512;          %OFDM符号数
N_FFT = 64;             %每个符号FFT长度
N_cp = 16;              %循环前缀长度
N_symbol = N_FFT + N_cp;%每个OFDM符号的长度
M = 16;K = 4;           %M：调制阶数
sr = 250000;  %Symbol rate 符号速率
br = sr.*K;  %Bit rate per carrier

EbN0 = 0:1:20;                      % 设出比特信噪比(dB)
snr = EbN0 + 10 * log10(K);         % 由公式推出snr(dB)表达式，公式见信噪比补充部分
BER(1 : length(EbN0)) = 0;          % 初始化误码率

%%%%%%生成需要传输的信息，并进行QAM调制%%%%%%
xn = randi([0,15],1,N_FFT * N_frame);                       %生成随机信息
xn_modulated = qammod(xn, M,'gray','UnitAveragePower',true);%进行QAM调制
avgPower = mean(abs(xn_modulated).^2);

%%%%%%添加cp,并进行ifft%%%%%%
xn_modulated1 = reshape(xn_modulated, N_FFT, N_frame);
x_ifft = ifft(xn_modulated1);
x_add_cp = add_cp(x_ifft, N_frame, N_FFT, N_cp);


%%%%%%使信号通过AWGN%%%%%%
parfor i = 1:length(snr)
    y_received = awgn(x_add_cp,snr(i),'measured');

    y_received1 = reshape(y_received, N_FFT + N_cp, N_frame);

    %%%%%%去掉cp,并进行fft%%%%%%
    y_remove_cp = fft(remove_cp(y_received1, N_frame, N_FFT, N_cp));

    %%%%%%QAM解调%%%%%%
    y_demodulated = qamdemod(y_remove_cp, M,'gray','UnitAveragePower',true);

    %%%%%%误比特率计算%%%%%%
    Neb = sum(sum(de2bi(y_demodulated,K) ~= de2bi(xn,K))); % 转为2进制，计算具体有几bit错误
    Ntb = N_FFT * N_frame * K;  % 仿真的总比特数
    BER(i) = Neb / Ntb;
end

semilogy(BER, 'r-o');
xlabel('EbN0(比特信噪比)');ylabel('BER(误比特率)');
grid on;














