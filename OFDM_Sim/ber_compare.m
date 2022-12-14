%%%%%%%%%%%%%%%%%%%%% 比较不同场景下的误码率 %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% ber_comppare.m  %%%%%%%%%
%%%%%%%%%  data:2020年11月7日  author:飞蓬大将军 %%%%%%%%%%

ebn0 = 3:1:10;
awgn_theory = [0.0228784075610853,0.0125008180407376,0.00595386714777866,0.00238829078093281,0.000772674815378444,0.000190907774075993,3.36272284196176e-05,3.87210821552205e-06];
awgn_no_compensation = [3.698496e-02,2.254329e-02,1.226654e-02,5.823633e-03,2.305339e-03,7.492187e-04,1.757812e-04,3.170573e-05];
rayleign_one_path_theory = [0.125000000000000,0.100000000000000,0.0833333333333333,0.0714285714285715,0.0625000000000000,0.0555555555555556,0.0500000000000000,0.0454545454545455]; 

semilogy(ebn0,awgn_theory,'-*',ebn0,awgn_no_compensation,'-+');
xlabel('比特信噪比');
ylabel('误码率');
title('OFDM在AWGN中采用QPSK调制不同信噪比下误码率仿真曲线');
legend('理论曲线','实验曲线');
grid on;