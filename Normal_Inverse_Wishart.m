function [mu,C]=Normal_Inverse_Wishart(A_1,A_2,k)
%Normal_Inverse_Wishar Puls k realizations from the Jeffreys prior posterior distribution
%  mu, C are cell arrays where each entry mu[i] C[i] is a realization

h1=A_1(2:end)-A_1(1:(end-1));
h2=A_2(2:end)-A_2(1:(end-1));
D=[h1,h2];
n=length(h1);
mu_hat=mean(D);
V_hat=zeros(2,2);
for i=1:n
    V_hat=V_hat+[D(i,1)-mu_hat(1); D(i,2)-mu_hat(2)]*[D(i,1)-mu_hat(1), D(i,2)-mu_hat(2)];
end
V_hat=V_hat/n;

%cell array of posteriori realizations
mu = cell(k, 1) ;
C = cell(k, 1) ;



for i=1:k
    %V_inv=wishrnd(inv(V_hat)/n,n-1);
    %V=inv(V_inv);
    %V_inv has wishard disttribution with parameter n and inv(T)
    %then V has iwishard distributinon with parameter n and T

    V=iwishrnd(V_hat*n,n-1);
    C{i}= chol(V,'lower');
    mu_temp=mvnrnd(mu_hat,V/n,1);
    mu{i}=mu_temp(1,:);
end
end

