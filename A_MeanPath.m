function [A_future_mean] =A_MeanPath(A1_start,A2_start,mu,t)
%A_MEANPATH  Calculates expeted future path of A(t)=[A1(t), A2(t)]
%   given he mean mu and the starting values A1 and A2 for the last period
%   before forecasting shall start. t tells the amount of forecastes
%   years

A_future_mean=zeros(t,2);
A_future_mean(1,:)=[A1_start,A2_start]+mu';

for i=2:t   
    A_future_mean(i,:)= A_future_mean(i-1,:)+mu';
end


end

