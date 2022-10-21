function H_LS = LS_CE(Y,Xp,pilot_loc,N_FFT,N_pilot_interval,int_opt)

Np=N_FFT/N_pilot_interval; 
k=1:Np; 
LS_est(k) = Y(pilot_loc(k))./Xp(k);  % LS channel estimation
if  lower(int_opt(1))=='l'
    method='linear'; 
else
    method='spline';
end
H_LS = interpolate(LS_est,pilot_loc,N_FFT,method); % Linear/Spline interpolation