% Generate the SD and SI noise according to the given SNR and $\alpha$
% assuming that sigma_w and sigma_u are contants for all bands
function [N_sd, N_si]=Noise_simu (F, snr, alpha)
% F is the input clean Hyper-spectral data
% snr is the SNR=10log_{10}\frac{\|\mathcal{F} \|_F^2}{\|\mathcal{N}(\mathcal{X}) \|_F^2} (dB)
% alpha=sigma2_sd/sigma2_w
% N_sd is the SD noise
% N_si is the SI noise
[n1,n2,n3]=size(F);
v = reshape(F, numel(F), 1);
Power_F= norm(v)^2/(n1*n2*n3); % the mean power of Hyper-spectral data
sigma2_total=Power_F*10^(-snr/10);
sigma2_sd=sigma2_total*alpha/(1+alpha);
sigma2_w=sigma2_total/(1+alpha);
sigma_w=sqrt(sigma2_w);
mu_n3=mean(mean(F,1),2);
mu=mean(mu_n3);
sigma2_u=sigma2_sd/mu;
sigma_u=sqrt(sigma2_u);
% N_sd=zeros(n1,n2,n3);
% N_si=zeros(n1,n2,n3);
% W=zeros(n1,n2,n3);
% U=zeros(n1,n2,n3);
W= sigma_w*randn(n1,n2,n3);
U=sigma_u*randn(n1,n2,n3);

% for i=1:n3
%   W(:,:,i)= sigma_w*randn(n1,n2);
%   U(:,:,i)=sigma_u*randn(n1,n2);
%    
% end
N_sd=sqrt(F).*U;
N_si=W;

