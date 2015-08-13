% This is the implementation of nonconvex RPCA method in ICDM 2015 paper
% 'Robust PCA via Nonconvex Rank Approximation'
% Zhao Kang, August 2015. Questions? Zhao.Kang@siu.edu;

clear all
close all

name='subject5';
load([name,'.mat']);
[m,n]=size(X);


muzero=.5;    % the only tuning parameter

lambda=1e-3;  % model parameter
type=21;   %different modeling of Sparsity.
rate=1.1;   %update rate of \mu
gamma=1e-2;     %gamma parameter in the rank approximation
tol=1e-3;  % stopping criterion

%initializations
S=zeros(m,n);
Y=zeros(m,n);
L=X;
sig=zeros(min(m,n),1); % for DC 
mu=muzero;

tic;
for ii=1:500
  
    D=X-S-Y/mu;
    [ L,sig] = DC(D,mu/2,sig,gamma);
    [S]=errorsol(Y,X,L,lambda,mu,type);
    Y=Y+mu*(L-X+S);
    mu=mu*rate;
    
    sigma=norm(X-S-L,'fro');
    RRE=sigma/norm(X,'fro');
    
    if RRE<tol
        break
    end
    
end
time_cost = toc;
rk=rank(L);
disp(['  relative err ',num2str(RRE),' rank ',num2str(rk),' time ',num2str(time_cost)])






