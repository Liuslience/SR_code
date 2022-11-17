function  x_add_awgn = add_awgn(xn, SNR, P_t, H_d)

N = sqrt(10 .^ (-SNR / 10) * P_t * mean(abs(H_d(:, 1) .^ 2)) / 2);
N_awgn = N .* (randn(size(xn)) + 1i * randn(size(xn)));
x_add_awgn = xn + N_awgn;

end



























