clc;

%%%%%%��������%%%%%%
N_frame = 1024;         %OFDM������
N_FFT = 64;             %ÿ������FFT����
N_cp = 16;              %ѭ��ǰ׺����
N_symbol = N_FFT + N_cp;%ÿ��OFDM���ŵĳ���
M = 16;K = 4;           %M�����ƽ���
sr = 250000;            %��������
br = sr.*K;             %ÿ���ز��ı�������

EbN0 = 0:1:50;                      % ������������(dB)
snr = EbN0 + 10 * log10(K);         % �ɹ�ʽ�Ƴ�snr(dB)���ʽ
BER = zeros(1, length(snr));          % ��ʼ��������

%%%%%%Fading initialization%%%%%%
PowerdB=[0 -8 -17 -21 -25]; % �ŵ���ͷ��������
% PowerdB=[0 -2 -3 -4 -5]; % �ŵ���ͷ��������
Delay=[0 3 5 6 8];          % �ŵ�ʱ��,ʾ��
% Delay=[0 3 5 56 78];        % �ŵ�ʱ��
Power=10.^(PowerdB/10);     % �ŵ���ͷ�������� '����'
Ntap=length(PowerdB);       % �ŵ���ͷ��
Lch=Delay(end)+1;           % �ŵ�����

%%%%%%������Ҫ�������Ϣ��������QAM����%%%%%%
xn = randi([0,15],1,N_FFT * N_frame);                       %���������Ϣ
xn_modulated = qammod(xn, M,'gray','UnitAveragePower',true);%����QAM����
avgPower = mean(abs(xn_modulated).^2);

%%%%%%����ifft,�����cp%%%%%%
xn_modulated1 = reshape(xn_modulated, N_FFT, N_frame);
x_ifft = ifft(xn_modulated1);
x_add_cp = add_cp(x_ifft, N_frame, N_FFT, N_cp);

%%%%%%����˥���ŵ�%%%%%%
channel = (randn(1,Ntap) + 1j * randn(1,Ntap)).*sqrt(Power/2);
h = zeros(1,Lch);
h(Delay+1) = channel;
x_fading = conv(x_add_cp,h);
x_fading1 = x_fading(:, 1:length(x_add_cp));
h_frenqency = fft(h);
% figure(1);
% stem(0:length(h_frenqency)-1, abs(h_frenqency), '.');
% title('�ྶʱ���ŵ���Ƶ��ͼ��');
% grid on;

% H = fft([h,zeros(1,N_FFT-Lch)].');
% y_equalization = zeros(N_FFT, N_frame);

%%%%%%ʹ�ź�ͨ��AWGN%%%%%%
for i = 1:length(snr)
    
    y_received = awgn(x_fading1,snr(i),'measured');

    y_received1 = reshape(y_received, N_FFT + N_cp, N_frame);
    
    %%%%%%ȥ��cp,������fft%%%%%%
    y_remove_cp = fft(remove_cp(y_received1, N_frame, N_FFT, N_cp));

    %%%%%%�ŵ�����%%%%%%
    H = fft([h,zeros(1,N_FFT-Lch)].');
    y_equalization = zeros(N_FFT, N_frame);
    for number = 1:N_frame
        y_equalization(:, number) = y_remove_cp(:, number)./H;
    end
    
    %%%%%%QAM���%%%%%%
    y_equal_1 = reshape(y_equalization, 1, N_frame * N_FFT);
    y_demodulated = qamdemod(y_equal_1, M,'gray','UnitAveragePower',true);

    %%%%%%������ʼ���%%%%%%
    Neb = sum(sum(int2bit(y_demodulated,K) ~= int2bit(xn,K))); %תΪ2���ƣ���������м�bit����
    Ntb = N_FFT * N_frame * K;                             %������ܱ�����
    BER(i) = Neb / Ntb;
end

figure(2);
semilogy(BER, 'r-o');
xlabel('EbN0(���������)');ylabel('BER(�������)');
title('�ྶ˥���ŵ��������ʷ�������');
grid on;










