function [ Est,EstCov,res, yF,UB,LB, Ysim ,simU, simL, YsimParUnc, U_ParUnc, L_ParUnc,aic,bic] = ArimaForecast(t_future,n_paths,Model,kt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[Est,EstCov,logL,~]=estimate(Model,kt);
res = infer(Est,kt);

parameter=2+length(Model.MA)+length(Model.AR);
[aic,bic] = aicbic(logL, parameter, length(kt));

[yF,yMSE] = forecast(Est,t_future,'Y0',kt);
UB = yF + 1.96*sqrt(yMSE);
LB = yF - 1.96*sqrt(yMSE);

Ysim = simulate(Est,t_future,'NumPaths',n_paths,'Y0',kt,'E0',res);
simU = prctile(Ysim,97.5,2);
simL = prctile(Ysim,2.5,2);

%%%including parameter uncertainty
%model mu=N(mu_hat,var(mu_hat)), mu_hat= mean(dk)=(k_t-k_1)/(t-1)
%var(mu_hat)=sigma_hat/(t-1, sigma_hat=sum(dk-mu_hat)^2/(t-2)
%sigma=var(diff(kt))/(t-1);

%matlab use maximul likelihood to find mu and above mu_hat as startwert

%MC with simulated mu for each path to include parameter uncertainty
YsimParUnc=zeros(t_future,n_paths);
temp=Est;
for i = 1:n_paths
    %get new random mu as N(mu,sigma) from the arima estimates
    mu_temp=normrnd(Est.Constant,EstCov(1,1));
        %create path 
        temp.Constant=mu_temp;
        temp_Ysim = simulate(temp,t_future,'NumPaths',1,'Y0',kt,'E0',res);
    YsimParUnc(:,i)=temp_Ysim;
end
U_ParUnc = prctile(YsimParUnc,97.5,2);
L_ParUnc = prctile(YsimParUnc,2.5,2);
end

