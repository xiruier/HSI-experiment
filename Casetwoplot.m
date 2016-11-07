%在Casetwo中，我们需要对Pavia city center和Washington DC这两个数据集都做实验
load Var_res_RLRTA_Pavia; % 数据格式是，[snr,alpha,rho_s, var_res],其中[snr,alpha,rho_s]是casetwo的试验参数，分别是snr=20dB, alpha=1, rho_s=0.1
                          % var_res是一个长度为n3 (n3 is the number of spectral
                          % bands).是对应于每个spectral band的residual noise
                          % variance
load Var_res_Hewei_Pavia;  
load Var_res_parafac_Pavia % Liu Xuefeng的文章

figure(1)
plot(Var_res_RLRTA_Pavia(3:end), 'b--o');
hold on
plot(Var_res_Hewei_Pavia(3:end), 'c-*');
hold on
plot(Var_res_parafac_Pavia(3:end), 'g-'); % Liu的实验还没做，先放那里

hold off

load Var_res_RLRTA_DC; % 数据格式是，[snr,alpha,rho_s, var_res],其中[snr,alpha,rho_s]是casetwo的试验参数，分别是snr=20dB, alpha=1, rho_s=0.1
                          % var_res是一个长度为n3 (n3 is the number of spectral
                          % bands).是对应于每个spectral band的residual noise
                          % variance
load Var_res_Hewei_DC;  
load Var_res_parafac_DC % Liu Xuefeng的算法



figure (2)
plot(Var_res_RLRTA_DC(3:end), 'b--o');
hold on
plot(Var_res_Hewei_DC(3:end), 'c-*');
hold on
plot(Var_res_parafac_DC(3:end), 'g-');

hold off


