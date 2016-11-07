addpath(genpath(pwd)); % 是为了把同一目录下的tensor_toolbox添加进来
% A mixture of random noise N(F) and sparse noise S

% 在Casetwo中，需要同时用到Pavia 和 DC Mall这两个数据集
%
k=1;
%k=2;

load Pavia_80;
F=OriData3;

% load Dc_Mall_data_info;
% F=D; 

F_max=max(reshape(F,numel(F),1));
F=F./F_max;
F=max(F,0);
maxP=max(F(:));
[n1,n2,n3]=size(F);

snr=20; alpha=1;
rho_s=0.1; % the sparse parameter
ratio = rho_s*ones(1,n3);

%%-----Noisy Image generation--------------------------
[N_sd, N_si]=Noise_simu (F, snr, alpha) ;
      Y=F+N_sd+N_si;  %Y is the noisy Hyper-sepctral data by adding elementwise noise N(F) to clean data F
for i =1:n3
     Y(:,:,i)=imnoise(Y(:,:,i),'salt & pepper',ratio(i));
%        mixed(vv,:)=addnoise(mixed(vv,:),'mg',ratio);
end
e1=1e-6;
e2=1e-5;
e3=1e-7;
[Fhat,Ksi,Ksd] = parafac_si_sd(Y,e1,e2,e3);

 Fhat = max(Fhat,0);
 Fhat = min(Fhat,maxP);
 psnr= PSNR(F,Fhat,maxP);
 PSNRvec=[snr,alpha,rho_s,psnr];
 
 Nresidual=F-Fhat;
 var_res=zeros(1,n3);
 for k=1:n3
     temp=Nresidual(:,:,k).*Nresidual(:,:,k);
     temp=sum(temp(:));
     var_res(k)=temp/(n1*n2);
 end %this loop computes the residual noise variance in each spectral band
 Var_res=[snr,alpha,rho_s, var_res]; %保存residual noise variance in each spectral band
 
 % 如果是用Pavia做实验
 PSNRvec_parafac_Pavia=PSNRvec;
 Var_res_parafac_Pavia=Var_res;
 Fhat_parafac_Pavia=Fhat;
 

 save PSNRvec_parafac_Pavia  PSNRvec_parafac_Pavia;
 save Var_res_parafac_Pavia    Var_res_parafac_Pavia;
 save Fhat_parafac_Pavia      Fhat_parafac_Pavia;
 
%---------------------------------------------------------------------------
 
% %  % 如果是用DC做数据集
% PSNRvec_parafac_DC=PSNRvec;
% Var_res_parafac_DC=Var_res;
% Fhat_parafac_DC=Fhat;


%  save PSNRvec_parafac_DC  PSNRvec_parafac_DC;
%  save Var_res_parafac_DC    Var_res_parafac_DC;
%  save Fhat_parafac_DC     Fhat_parafac_DC;
