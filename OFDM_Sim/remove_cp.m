function y_remove_cp = remove_cp(y_received, N_frame, N_FFT, N_cp)

% y_received = reshape(y_received, N_frame, N_FFT + N_cp);
y_remove_cp = y_received(N_cp + 1 : N_FFT + N_cp, :);
% y_remove_cp = reshape(y_remove_cp, 1, N_frame * N_FFT);

end