function xn_modulated = qpskmod(xn, N_signal)

xn_modulated = zeros(N_signal, 1);
for i = 1 : N_signal

    if xn(i, 1) == 1 && xn(i, 2) == 1
        xn_modulated(i, 1) = (1 + 1i)/sqrt(2);
    elseif xn(i, 1) == 0 && xn(i, 2) == 1
        xn_modulated(i, 1) = (-1 + 1i)/sqrt(2);
    elseif xn(i, 1) == 0 && xn(i, 2) == 0
        xn_modulated(i, 1) = (-1 - 1i)/sqrt(2);
    else
        xn_modulated(i, 1) = (1 - 1i)/sqrt(2);
    end
end

end








