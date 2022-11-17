function c = arg_min(H, Hd, Hb, N_frame)

c = zeros(1, N_frame);

for i = 1 : N_frame

    c_0 = abs(1 / sum(abs(Hb(:, 1)) .^ 2) * Hb(:, 1)' * (H(:, i) - Hd) + 1);
    c_1 = abs(1 / sum(abs(Hb(:, 1)) .^ 2) * Hb(:, 1)' * (H(:, i) - Hd) - 1);
    
    if c_0 <= c_1
        c(i) = -1;
    else
        c(i) = 1;
    end

end

end































