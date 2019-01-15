function [A_1, A_2, mu, C]=q_RegressionFit(q,x_min_C,x_max_C)
%q_REGRESSIONFIT Estimates the Processes A_1 and A_2 in Cairns Model
%Also estiamtes mu and C to fit A(t+1)-A(t)=N(mu,V)=mu+ N(0,C*C')=mu+C*N(0,I)
%The Model is q(t,x)=exp(A_1(t)+A_2(t)*(x+t))/(1+exp(A_1(t)+A_2(t)*(x+t)))
%Using the funktion f(x)=log(x/(1-x) on both sides we get the linear model
%log(q(x,t)/(1-q(x,t)))=A_1(t)+A_2(t)*(x+t)
%for fixed t this gives a model Y=X*b with b=[A_1, A,2]
%simple linear regression gives the least sqarres estimate

n_C=size(q);
x_C=n_C(1);
t_C=n_C(2);

A_1=zeros(t_C,1);
A_2=zeros(t_C,1);

y=linspace(x_min_C,x_max_C,x_C)';
for i=1:t_C
    Y=log(q(:,i)./(1-q(:,i)));
    X=[ones(x_C,1),y+i-1];
    b=X\Y;
    A_1(i)=b(1);
    A_2(i)=b(2);
end

%A(t+1)-A(t)=N(mu,V)=mu+ N(0,C*C')=mu+C*N(0,I)
h1=A_1(2:end)-A_1(1:(end-1));
h2=A_2(2:end)-A_2(1:(end-1));
h=[h1,h2];
mu=zeros(2,1);
mu(1)=mean(h1);
mu(2)=mean(h2);
V=cov(h)
%C*C'=V lower triangular matrix %why did cairns say upper??
C= chol(V,'lower'); 
end

