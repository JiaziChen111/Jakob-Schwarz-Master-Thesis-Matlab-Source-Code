function [A_future]=A_UncertainSamplePaths(A1_start,A2_start,mu_hat,V_hat,df,t,n)
%A_SAMPLE_PATHS  creates a cell array A_future size n of simulated paths legnth t
%   For each path in the cell array a different mu and C is used given by
%   a realization of the posterior distribution with parameters mu_hat,
%   V_hat and df

A_future = cell(n, 1) ;
for i=1:n
    
    %create random mu and C from posterior distributon
    V=iwishrnd(V_hat*df,df-1);
    C= chol(V,'lower');
    mu=mvnrnd(mu_hat,V/df,1);
    
    %create the random path
    A_future{i}=zeros(t,2);
    Z= mvnrnd([0 0],[1 0;0 1]);
    A_future{i}(1,:)=[A1_start,A2_start]+mu+(C*Z')';
    
    for j=2:t   
        Z= mvnrnd([0 0],[1 0;0 1]);
        A_future{i}(j,:)= A_future{i}(j-1,:)+mu+(C*Z')';
    end
    
end



end