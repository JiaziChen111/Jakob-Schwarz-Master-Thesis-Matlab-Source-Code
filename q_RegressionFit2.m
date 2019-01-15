function [A_1, A_2]=q_RegressionFit2(q,x_min_C,x_max_C)
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
    %if q=1 Y becomes Inf, ignore those rows
    H=(isinf(Y)-1)*(-1);
    H=logical(H);
    
    b=X(H,:)\Y(H);
    A_1(i)=b(1);
    A_2(i)=b(2);
end

end

