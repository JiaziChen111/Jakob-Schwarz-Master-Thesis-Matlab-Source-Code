function Msim=getMsim(M,a_x,b_x,k_t)
%GETMSIM Creates a cell array of M matrixes Fit+Forecast
%   M is the fitted matrix to historical data
%   a_x, b_x the fitted model parameter
%   k_t the forecasted k_t parameter given as a txn matrix
%   t forecast year and n simualted path
%   Msim is a n cell array, each entry [M,M_forecast(n)]
n=size(k_t);
Msim=cell(n(2),1);
x=length(a_x);
M_Forecast=zeros(x,n(1));
%outer path loop
for k=1:n(2)
    for i=1:x
        %forecasted year loop
        for j=1:n(1)
            M_Forecast(i,j)=exp(a_x(i)+b_x(i)*k_t(j,k));
        end
    end
    Msim{k}=[M,M_Forecast]; 
end

end

