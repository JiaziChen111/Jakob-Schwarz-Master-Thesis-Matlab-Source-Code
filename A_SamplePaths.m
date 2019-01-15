function [A_future]=A_SamplePaths(A1_start,A2_start,mu,C,t,n);
%A_SAMPLE_PATHS  creates a cell array A_future size n of simulated paths legnth t
%   Detailed explanation goes here

A_future = cell(n, 1) ;
for i=1:n
    A_future{i}=zeros(t,2);
    Z= mvnrnd([0 0],[1 0;0 1]);
    A_future{i}(1,:)=[A1_start,A2_start]+mu'+(C*Z')';
    
    for j=2:t   
        Z= mvnrnd([0 0],[1 0;0 1]);
        A_future{i}(j,:)= A_future{i}(j-1,:)+mu'+(C*Z')';
    end
    
end



end

