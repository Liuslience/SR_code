function H_LS = LS_CE(Y, Xp, pilot_loc, N_FFT, Np, int_opt)
 
k=1:Np; 
if  lower(int_opt(1))=='l'
    method='linear'; 
else
    method='spline';
end

for i = 1 : N_FFT

    LS_est(k) = Y(pilot_loc(k), i)./Xp(k);  % LS channel estimation
    H_LS(i, :) = interpolate(LS_est,pilot_loc,N_FFT,method); % Linear/Spline interpolation

end

% LS_est(k) = Y(pilot_loc(k), 1)./Xp(k);  % LS channel estimation
% H_LS = interpolate(LS_est,pilot_loc,N_FFT,method); % Linear/Spline interpolation

end










