function [a_x,b_x,k_t,fval,exitflag]=FitLC_simplex(lnM, a_start, b_start, k_start,D)
%FITLC_GDO fit the extended GDP Lee Carter model: ln m(x,t)=a_x+b_x*k_t+y_x*g_t+e
    %D optinal weights    
    %see Fit_wLC.m for more details.
    %here i use the simplex method instead

n=size(lnM);
x=n(1);
t=n(2);
%if no weights given set to 1 except for where 42 (M<=0)
if nargin<=4
    D=ones(n);
    B=lnM==log(42);
    D(B)=0;
end

%options = optimset('MaxFunEvals',4000000,'MaxIter',8000000, 'TolFun',1e-10,'TolX',1e-10);
options = optimset('MaxFunEvals',40000,'MaxIter',5, 'TolFun',10e-10,'TolX',10e-10);
FUN = @(z) sum(sum(D.*(lnM-repmat(z(1:x),1,t)-z(x+1:2*x)*z(2*x+1:2*x+t)').^2   ))   ; 
x0 = [a_start;b_start;k_start];

[z,fval,exitflag] = fminsearch(FUN,x0,options);
a_x=z(1:x);
b_x=z(x+1:2*x);
k_t=z(2*x+1:2*x+t);

        %normalize such that sum b =1
c=sum(b_x);
b_x=b_x/c;
k_t=k_t*c;

%sum(k_t)=0
c=-sum(k_t)/t;
a_x=a_x-b_x*c;
k_t=k_t+c;

  
end

