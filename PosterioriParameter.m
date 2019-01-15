function [mu_hat,V_hat,df]=PosterioriParameter(A_1,A_2)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

h1=A_1(2:end)-A_1(1:(end-1));
h2=A_2(2:end)-A_2(1:(end-1));
H=[h1,h2];
df=length(h1);
mu_hat=mean(H);
V_hat=zeros(2,2);
for i=1:df
    V_hat=V_hat+[H(i,1)-mu_hat(1); H(i,2)-mu_hat(2)]*[H(i,1)-mu_hat(1), H(i,2)-mu_hat(2)];
end
V_hat=V_hat/df;
end

