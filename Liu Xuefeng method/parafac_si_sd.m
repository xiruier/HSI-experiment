% [Fhat] = parafac_si_sd(Y,e1,e2,e3)

function [Fhat,Ksi,Ksd] = parafac_si_sd(Y,e1,e2,e3)

[n1,n2,n3]=size(Y);
n=min(min(n1,n2),n3);
n0=n/2;
muY=mean(mean(Y,1),2);
muY=muY(:);
invmuY=1./muY;
diaginvmu=diag(invmuY);

opts.maxit = 200; % max number of iterations
opts.tol = 1e-4; % stopping tolerance

esty = zeros(n1*n2,n3);

%step one, remove the SI noise
for R=n0:n
    [A,Out] = ncp(Y,R,opts);
    % removed noise
    Nre=Y-A;
    
    for i= 1:n3
        mid = Nre(:,:,i);
        esty(:,i) = mid(:);
    end
    esty=max(esty,0);
    RR=cov(esty);
    
    rr=diag(RR);
    rr2=sum(rr.*rr);
    Q=RR.*RR;
    Q0=sum(Q(:));
    if abs(Q0-rr2)<e2
        break
    end
    
    
end

Ksi=R;
Yhat=A;

%step two, SD noise reduction

for R=n0:n
    [A,Out] = ncp(Yhat,R,opts);
    % removed noise
    Nre=Yhat-A;
    
    for i= 1:n3
        mid = Nre(:,:,i);
        esty(:,i) = mid(:);
    end
    esty=max(esty,0);
    RR=cov(esty);
    RR=RR*diaginvmu;
    rr=diag(RR);
    rr2=sum(rr.*rr);
    Q=RR.*RR;
    Q0=sum(Q(:));
    if abs(Q0-rr2)<e3
        break
    end
    
    
end

Ksd=R;
Fhat=A;



    
    
    
    
    
    
    
    
    
    
    
    %step two, remove the SD noise