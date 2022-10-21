function r_remove_pilot = remove_pilot(y_remove_cp, N_pilot_interval, N_pilot_number, start_pilot, N_FFT, N_frame)

for j = 1 : N_pilot_number
        r_remove_pilot(1 + (j-1) * (N_pilot_interval - 1) : j * (N_pilot_interval - 1), :) = y_remove_cp(start_pilot + 1 + (j-1) * N_pilot_interval : j * N_pilot_interval, :);
end

end