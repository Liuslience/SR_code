function xn_add_pilot = add_pilot(xn_modulated, xn_pilot_signal, pilot_loc, N_pilot_number, N_symbol, N_frame)

xn_add_pilot = zeros(N_symbol + N_pilot_number, N_frame);
pilot_interval_1 = [1, 11, 24, 37, 50];
pilot_interval_2 = [10, 23, 36, 49, 60];

%% 添加导频信号
for i = 1 : N_pilot_number
    xn_add_pilot(pilot_loc(i), :) = xn_pilot_signal(1, i);
end

pilot_loc = [0, pilot_loc, 65];

%% 搬移有用信号
for i = 1 : length(pilot_loc) - 1
    xn_add_pilot(1 + pilot_loc(i) : pilot_loc(i + 1) - 1, :) = xn_modulated(pilot_interval_1(i) : pilot_interval_2(i), :);
end

end






















