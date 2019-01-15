function [ MaxErr,MinErr,MeanErr,Err2,FVU,FVUx ] = errorfkt_q( q,q_hat)
%ERRORFKT Calculates statistics to compare different fittings
%   Defintions are given in the Masterarbeit chapter 2
n=size(q);
x_C=n(1);
t_C=n(2);

Err=q-q_hat;


MaxErr=max(max(Err));
MinErr=min(min(Err));

Err2=sum(sum(Err.^2))/(x_C*t_C);
MeanErr=sum(sum(Err))/(x_C*t_C);

Rx=zeros(x_C,1);
totalVar=0;
totalVar2=0;

%X=zeros(x,2);
for i=1:x_C
  
    k2=var(Err(i,:));
    k=var(q(i,:));
    Rx(i)= k2/k;
    totalVar=totalVar+k;
    totalVar2=totalVar2+k2;
    %X(i,:)=[k2,k]
end 
FVUx=Rx*100;

FVU=totalVar2/totalVar*100;


end

