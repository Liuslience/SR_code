% N = sqrt(10 .^ (-SNR / 10) * P_t / 2);
N_awgn = randn(1, 64) + 1i * randn(1, 64);
p1 = mean(abs(N_awgn) .^ 2)
% N_awgn = N_awgn ./ sqrt(mean(abs(N_awgn) .^ 2));
% p2 = mean(abs(N_awgn) .^ 2);
% x_add_awgn = xn + N_awgn;

















