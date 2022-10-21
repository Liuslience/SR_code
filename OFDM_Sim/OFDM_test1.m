clc;

%%%%%%��������%%%%%%
N_frame = 512;          %OFDM������
N_FFT = 64;             %ÿ������FFT����
N_cp = 16;              %ѭ��ǰ׺����
N_symbol = N_FFT + N_cp;%ÿ��OFDM���ŵĳ���
M = 16;K = 4;           %M�����ƽ���
sr = 250000;  %Symbol rate ��������
br = sr.*K;  %Bit rate per carrier

EbN0 = 0:1:20;                      % ������������(dB)
snr = EbN0 + 10 * log10(K);         % �ɹ�ʽ�Ƴ�snr(dB)���ʽ����ʽ������Ȳ��䲿��
BER(1 : length(EbN0)) = 0;          % ��ʼ��������

%%%%%%������Ҫ�������Ϣ��������QAM����%%%%%%
xn = randi([0,15],1,N_FFT * N_frame);                       %���������Ϣ
xn_modulated = qammod(xn, M,'gray','UnitAveragePower',true);%����QAM����
avgPower = mean(abs(xn_modulated).^2);

%%%%%%���cp,������ifft%%%%%%
xn_modulated1 = reshape(xn_modulated, N_FFT, N_frame);
x_ifft = ifft(xn_modulated1);
x_add_cp = add_cp(x_ifft, N_frame, N_FFT, N_cp);


%%%%%%ʹ�ź�ͨ��AWGN%%%%%%
parfor i = 1:length(snr)
    y_received = awgn(x_add_cp,snr(i),'measured');

    y_received1 = reshape(y_received, N_FFT + N_cp, N_frame);

    %%%%%%ȥ��cp,������fft%%%%%%
    y_remove_cp = fft(remove_cp(y_received1, N_frame, N_FFT, N_cp));

    %%%%%%QAM���%%%%%%
    y_demodulated = qamdemod(y_remove_cp, M,'gray','UnitAveragePower',true);

    %%%%%%������ʼ���%%%%%%
    Neb = sum(sum(de2bi(y_demodulated,K) ~= de2bi(xn,K))); % תΪ2���ƣ���������м�bit����
    Ntb = N_FFT * N_frame * K;  % ������ܱ�����
    BER(i) = Neb / Ntb;
end

semilogy(BER, 'r-o');
xlabel('EbN0(���������)');ylabel('BER(�������)');
grid on;














