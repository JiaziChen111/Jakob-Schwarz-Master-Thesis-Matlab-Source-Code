function [a_x,b_x,k_t,y_x,fval,exitflag]=FitLC_Explanation_newton(lnM,g_t, a_start, b_start, k_start,y_start,W,norming)
%FITLC_GDO fit the extended GDP Lee Carter model: ln m(x,t)=a_x+b_x*k_t+y_x*g_t+e
    %W  weights    
    %see Fit_wLC.m for more details.
    %here i use the quasi newton method instead
    %norming=1 cov(k_t-k_t-1,g_t-g_t-1)=0  else cov(k_t,g_t)=0
    %g_t is a matrix= [g_1(1:n)|...|g_t(1:n)]
    %y_x is a matrix= [y_1:x(1)|...|y_1:x(n)]

n=size(lnM);
x=n(1);
t=n(2);
n=size(y_start);
n=n(2);

y0=reshape(y_start,x*n,1);

options = optimoptions(@fminunc,'Algorithm','quasi-newton','MaxFunEvals',800000,'MaxIter',2000, 'TolFun',1e-10,'TolX', 10e-6);
FUN = @(z) sum(sum(W.*(lnM-repmat(z(1:x),1,t)-z(x+1:2*x)*z(2*x+1:2*x+t)'-reshape(z(2*x+t+1:end),x,n)*g_t).^2   ))   ; 
x0 = [a_start;b_start;k_start;y0];
[z,fval,exitflag] = fminunc(FUN,x0,options);

a_x=z(1:x);
b_x=z(x+1:2*x);
k_t=z(2*x+1:2*x+t);
y_x=z(2*x+t+1:end);
y_x=reshape(y_x,x,n);

%norming
if (norming ==1) %differene independece neu
     dg=diff(g_t');
     CovMatrixdgti= cov(dg);
     covVectorktdgti=zeros(1,n);
     for i=1:n
         temp=cov(diff(k_t),diff(g_t(i,:)));
         covVectorktdgti(i)=temp(2,1);      
     end
         
     e=CovMatrixdgti\covVectorktdgti;      
     
else %normal alt
     CovMatrixgti= cov(g_t');
     covVectorktgti=zeros(n,1);
     for i=1:n
         temp=cov(k_t,g_t(i,:));
         covVectorktgti(i)=temp(2,1);      
     end
     e=CovMatrixgti\covVectorktgti;     
    
end
        
        k_t=k_t-g_t'*e;
        y_x=y_x+b_x*e';     
        h=sum(b_x);
        b_x=b_x/h;
        k_t=k_t*h;
        c=-sum(k_t)/t;
        a_x=a_x-b_x*c;
        k_t=k_t+c;        
    
end

