function xn_add_pilot = add_pilot(xn_modulated, N_pilot_interval, N_pilot_number, start_pilot, N_symbol, N_frame)

xn_pilot_signal = ones(1, N_pilot_number);
xn_add_pilot = zeros(N_symbol + N_pilot_number, N_frame);
for i = 1 : N_pilot_number
    xn_add_pilot(start_pilot + (i-1) * N_pilot_interval, :) = xn_pilot_signal(1, i);
    xn_add_pilot(start_pilot + 1 + (i-1) * N_pilot_interval : i * N_pilot_interval, :) = xn_modulated(1 + (i-1) * (N_pilot_interval - 1) : i * (N_pilot_interval - 1), :);
end

end






















