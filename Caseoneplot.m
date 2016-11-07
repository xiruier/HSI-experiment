%% Caseone,可以参见文章"Reduction of Signal-Dependent Noise From Hyperspectral Images for Target Detection"
% written by Liu Xuefeng.  在Caseone实验结果的展示主要借鉴的就是她的文章
% 如 文章里的图5， 我们需要显示不同算法在不同(snr, alpha)组合中，得到的residual noise, 即 F-Fhat.
% 在她的文章里选择的是， snr=20dB,alpha=0.1; snr=30dB, alpha=1; snr=40, alpha=10; 这三组值
% 我们也先暂定这三组，可以根据具体实验结果来选择更适合显示的组合

load Var_res_RLRTA_Caseone; % Var_res_RLRTA 是一个size为(snr_len*alpha_len, 2+n3)的矩阵
% Var_res_RLRTA_Caseone 是用Caseone_simu_RLRTA.m 得到，并保存下来的
load Var_res_Hewei_Caseone; %  Var_res_Hewei
% Var_res_Hewei_Caseone 是用 Casetwo_Hewei.m 得到， 并保存下来的
load Var_res_parafac_Caseone; %Liu Xuefeng的parafac 算法，
table=[20,0.1;30,1;40,10]; %选了三组不同的snr+alpha组合，也可以选其他的
num1=size(table,1);
index=zeros(num1,1);
[m1,m2]=size(Var_res_RLRTA);
for Outloop=1:numl
    for j=1:m1
    if Var_res_RLRTA(j,1)==table(Outloop,1)&&Var_res_RLRTA(j,2)==table(Outloop,2)
       index(Outloop)=j;
       break;
    end
    end
end

%% 显示结果

for Outloop=1:numl
    
    figure(Outloop)
    plot(Var_res_RLRTA(index(Outloop), 3:end), 'b--o');
    hold on
    plot(Var_res_Hewei(index(Outloop), 3:end), 'c-*');
    hold on
    plot(Var_res_parafac(index(Outloop), 3:end), 'g-'); 
    
    hold off
    
    
end