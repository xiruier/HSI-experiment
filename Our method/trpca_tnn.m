function [F,S,outputs] = trpca_tnn(Y, mu1, mu2, mu3,opts)
% Solve the HSIs-Rec problem based on Low-rank Tensor decomposition via variant of ADMM
% min_{F,U,N,S} ||F||_* + mu1*||sqrt(F).*U||_2^2 + mu2*||N||_2^2 + mu3*||S||_1 ,
% s.t. Y = F + sqrt(F).*U + N + S
% ---------------------------------------------
% reformulated as:
% min_{F,A,U,N,S} ||F||_* + mu1*||A.*U||_2^2 + mu2*||N||_2^2 + mu3*||S||_1 ,
% s.t. Y = F + A.*U + N + S, A.*A = F
% ---------------------------------------------
% Input:
%       Y             -    d1*d2*d3 tensor
%       mu1, mu2, m3  -    >0, parameter
%       opts          -    Structure value in Matlab. The fields are
%         opts.tol        -   termination tolerance
%         opts.max_iter   -   maximum number of iterations
%         opts.beta       -   stepsize for dual variable updating in ADMM
%         opts.max_beta   -   maximum dual variable stepsize
%         opts.rho        -   rho>=1, ratio used to increase beta
%         opts.DEBUG      -   0 or 1
% ---------------------------------------------
% Output:
%       F       -    d1*d2*d3 tensor
%       A       -    d1*d2*d3 tensor
%       U       -    d1*d2*d3 tensor
%       N       -    d1*d2*d3 tensor
%       S       -    d1*d2*d3 tensor
%       obj     -    objective function value
%       err     -    residual 
%       iter    -    number of iterations
outputs = [];
objs = [];
errs = [];
iters= [];

rho      = 1.1;
tol      = 1e-8; 
beta     = 1e-4;
max_beta = 1e10;
max_iter = 500;
DEBUG    = 0;

if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'beta');        beta = opts.beta;            end
if isfield(opts, 'max_beta');    max_beta = opts.max_beta;    end
if isfield(opts, 'DEBUG');       DEBUG = opts.DEBUG;          end

dim     = size(Y);
F       = zeros(dim);
A       = zeros(dim);
U       = zeros(dim);
N       = zeros(dim);
S       = zeros(dim);
Lambda1 = zeros(dim);
Lambda2 = zeros(dim);

iter = 0;
chg_pre = inf;
for iter = 1 : max_iter
    Fk = F;
    Sk = S;
    
    % update F
    t_var = 0.5*(A.*A + Lambda1/beta + Y - (A.*U + N + S) + Lambda2/beta);
    %t_var(isnan(t_var) | isinf(t_var)) = 0;
    %sum(sum(sum(isnan(t_var) + isinf(t_var))))    
    [F,tnnL] = prox_tnn(t_var, 1/(2*beta));   
    % update S
    S = prox_l1(Y - (F + A.*U + N) + Lambda2/beta, mu3/beta);   
    
    % update N
    N = (beta/(2*mu2+beta))*(Y + Lambda2/beta - (F + A.*U + S));      
    
%     % update A
%     A = (beta*U.*(Y + Lambda2/beta - (F + N + S)))...
%         ./(2*(beta*(A.*A - F) + Lambda1) + (2*mu1+beta)*(U.*U));     
%     A(isnan(A) | isinf(A)) = 0;

    % update A
    A = sqrt(F);
    
    % update U
    U = (beta/(2*mu1+beta))*(Y + Lambda2/beta - (F + N + S))./A;
    %U(isnan(U) | isinf(U)) = 0;
    
    % update Lambda1
    Lambda1 = Lambda1 + beta*(A.*A - F);  
  
    dY = Y - (F + A.*U + N + S);
    Lambda2 = Lambda2 + beta*dY;
    
    chgL = max(abs(Fk(:)-F(:)));
    chgS = max(abs(Sk(:)-S(:)));
    chg = max([ chgL chgS max(abs(dY(:))) ]);
    if DEBUG
        if iter == 1 || mod(iter, 10) == 0
            obj = tnnL + mu1*norm(sqrt(F(:)).*U(:)) + mu2*norm(N(:)) + mu3*norm(S(:),1);
            err = norm(dY(:));
            objs = [objs obj];
            errs = [errs err];
            iters= [iters iter];
            disp(['iter ' num2str(iter) ', beta=' num2str(beta) ...
                    ', obj=' num2str(obj) ', err=' num2str(err)...
                    ', chgL=' num2str(chgL) ', chgS=' num2str(chgS)...
                    ', chg=' num2str(chg)]...
                    ); 
        end
    end
    
    if abs(chg - chg_pre) < tol
        break;
    end 
    chg_pre = chg;
    beta = min(rho*beta,max_beta);
    
    %imshow(F/max(F(:)));
end

outputs.objs = objs;
outputs.errs = errs;
outputs.iters = iters;

