function [a_x,b_x,k_t,fval,exitflag]=PoissonFit_newton(D, E, a_start,b_start,k_start)
%FITLC_GDO fit the extended GDP Lee Carter model: ln m(x,t)=a_x+b_x*k_t+y_x*g_t+e
    %D optinal weights    
    %see Fit_wLC.m for more details.
    %here i use the quasi newton method instead

    
%starting values
n=size(D);
x=n(1);
t=n(2);

if nargin==2 
   a_start=zeros(x,1);
   b_start=ones(x,1);
   k_start=zeros(t,1);
end    
    
options = optimoptions(@fminunc,'Algorithm','quasi-newton','MaxFunEvals',800000,'MaxIter',2000, 'TolFun',1e-6,'TolX', 1e-6);

L = @(z) sum(sum( E.*exp(repmat(z(1:x),1,t)).*exp(z(x+1:2*x)*z(2*x+1:2*x+t)') - D.*(repmat(z(1:x),1,t)+z(x+1:2*x)*z(2*x+1:2*x+t)')  ))   ; 
x0 = [a_start;b_start;k_start];

[z,fval,exitflag] = fminunc(L,x0,options);
fval=-fval;
a_x=z(1:x);
b_x=z(x+1:2*x);
k_t=z(2*x+1:2*x+t);

%ensuring model indetification sum(x) b_x = 1 
%note for any c is a,b*c,k/c also a solution for ln m_xt=a+b*k
c=sum(b_x);
b_x=b_x/c;
k_t=k_t*c;
%ensuring model indetification sum k=0
%note a -bc, b , k+c is also a solution
c=sum(k_t)/t;
k_t=k_t-c;
a_x=a_x+b_x*c;

        
        %hhh=norm(FUN([a_x;b_x;k_t;y_x])-fval);
        %display(hhh); 3.6e-12 so no difference from transforming
end

