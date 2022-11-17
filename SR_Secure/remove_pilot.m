function r_remove_pilot = remove_pilot(y_remove_cp, pilot_loc, N_frame, N_FFT)

r_remove_pilot = zeros(N_FFT, N_frame);
pilot_loc = [0, pilot_loc, 65];
pilot_interval_1 = [1, 11, 24, 37, 50];
pilot_interval_2 = [10, 23, 36, 49, 60];
for i = 1 : length(pilot_loc) - 1
    r_remove_pilot(pilot_interval_1(i) : pilot_interval_2(i), :) = y_remove_cp(1 + pilot_loc(i) : pilot_loc(i + 1) - 1, :);
end

end