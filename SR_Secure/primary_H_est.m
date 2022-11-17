function H_est = primary_H_est(y_remove_cp, pilot_loc, P_t, F_p, N_FFT, N_frame, F_L)

H_est = zeros(N_FFT, N_frame);

for i = 1 : N_frame
    H_est(:, i) = F_L * inv(P_t * F_p' * F_p) * sqrt(P_t) * F_p' * y_remove_cp(pilot_loc, i);
end









end










