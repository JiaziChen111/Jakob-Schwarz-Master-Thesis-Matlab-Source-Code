function [A_1_lower,A_1_upper,A_2_lower,A_2_upper]=QuantilPath(A_future,alpha)
%calc alpha/2 and 1-alpha/2 quantil path for A_1 and A_2 sample path
%   Detailed explanation goes here

n=length(A_future);
h=size(A_future{1});
t=h(1);
X=zeros(t,n); %data matrix A1(i,j)
    for j=1:n %which path
        for i=1:t %which period
            X(i,j)=A_future{j}(i,1);
        end
    end
A_1_lower = prctile(X,100*(alpha/2),2);
A_1_upper = prctile(X,100*(1-alpha/2),2);
    
    for j=1:n %which path
        for i=1:t %which period
            X(i,j)=A_future{j}(i,2);
        end
    end

A_2_lower = prctile(X,100*alpha/2,2);
A_2_upper = prctile(X,100*(1-alpha/2),2);
    
end

