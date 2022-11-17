function H_MMSE = MMSE_CE(Y,Xp,pilot_loc,N_FFT,N_pilot_interval,h,SNR)
%function H_MMSE = MMSE_CE(Y,Xp,pilot_loc,Nfft,Nps,h,ts,SNR)
% MMSE channel estimation function
% Inputs:
%       Y         = Frequency-domain received signal
%       Xp        = Pilot signal
%       pilot_loc = Pilot location
%       Nfft      = FFT size
%       Nps       = Pilot spacing
%       h         = Channel impulse response
%       ts        = Sampling time
%       SNR       = Signal-to-Noise Ratio[dB]
% output:
%      H_MMSE     = MMSE channel estimate

%MIMO-OFDM Wireless Communications with MATLAB¢ç   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
%?2010 John Wiley & Sons (Asia) Pte Ltd

%H = fft(h,N);
%%
snr = 10^(SNR*0.1);
Np=N_FFT/N_pilot_interval; 
k=1:Np; 
H_tilde = (Y(pilot_loc(k), 1)./Xp(k)).';    % LS estimate
k=0:length(h)-1;                       %多径信道时延矩阵
hh = h*h';                             %多径信道的平均接收功率的和
tmp = h.*conj(h).*k;                   %多径信道接收功率的时延加权
r = sum(tmp)/hh;                       %平均附加时延（多径信道接收功率的时延加权/多径信道的平均接收功率的和）
r2 = tmp*k.'/hh;                       %二阶中心矩的开根号/多径信道的平均接收功率的和
tau_rms = sqrt(r2-r^2);                %rms时延
df = 1/N_FFT;                          %子载波间隔，即FFT长度的倒数
j2pi_tau_df = 1i*2*pi*tau_rms*df;

%%
K1 = repmat((0:N_FFT-1).',1,Np); 
K2 = repmat(0:Np-1,N_FFT,1);
rf = 1./(1+j2pi_tau_df*(K1-K2*N_pilot_interval));

K3 = repmat((0:Np-1).',1,Np); 
K4 = repmat(0:Np-1,Np,1);
rf2 = 1./(1+j2pi_tau_df*N_pilot_interval*(K3-K4));

Rhp = rf;
Rpp = rf2 + eye(length(H_tilde),length(H_tilde))/snr;
% H_MMSE = transpose(Rhp*Rpp\H_tilde.');  % MMSE channel estimate
H_MMSE = transpose(Rhp*inv(Rpp)*H_tilde.');  % MMSE channel estimate