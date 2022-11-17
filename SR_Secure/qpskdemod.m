function xn_demodulated = qpskdemod(xn, N_symbol, N_frame)

xn_demodulated = zeros(N_symbol * N_frame, 2);
for i = 1 : N_symbol * N_frame
    if real(xn(i)) >= 0 && imag(xn(i)) > 0
        xn_demodulated(i, :) = [1, 1];
    elseif real(xn(i)) < 0 && imag(xn(i)) >= 0
        xn_demodulated(i, :) = [0, 1];
    elseif real(xn(i)) <= 0 && imag(xn(i)) < 0
        xn_demodulated(i, :) = [0, 0];
    else
        xn_demodulated(i, :) = [1, 0];
    end
end
% xn_demodulated = xn_demodulated .* sqrt(2);

end