function H_est_LS = estimate_2(xn_add_pilot, y_remove_cp, P_t, F_L, N_FFT, Lch, N_frame)

S_n = zeros(N_FFT);
H_est_LS = zeros(N_FFT, N_frame);

for i = 1 : N_frame
    for j = 1 : N_FFT
        S_n(j, j) = xn_add_pilot(j, i);
    end
    H_est_LS(:, i) = F_L * inv(P_t .* F_L' * S_n' * S_n * F_L) * sqrt(P_t) * F_L' * S_n' * y_remove_cp(:, i);
end






























end