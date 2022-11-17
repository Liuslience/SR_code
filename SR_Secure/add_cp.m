function x_add_cp = add_cp(xn_modulated, N_frame, N_FFT, N_cp)

% xn_modulated = reshape(xn_modulated, N_frame, N_FFT);
% x_add_cp = reshape(xn_modulated, N_frame, N_FFT);
x_add_cp = [xn_modulated(N_FFT - N_cp +1: N_FFT, :); xn_modulated];
x_add_cp = reshape(x_add_cp, 1, N_frame * (N_FFT + N_cp));

end