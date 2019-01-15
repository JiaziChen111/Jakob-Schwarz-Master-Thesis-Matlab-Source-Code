function [ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,FVUx,maxFVUx,x_maxFVUx,FVU ] = errorfkt( M,D,E,a_x,b_x,k_t )
%ERRORFKT Calculates statistics to compare different fittings
%   Defintions are given in the Masterarbeit chapter 2
n=size(M);
x=n(1);
t=n(2);
M_hat=zeros(x,t);
    for i=1:x
        for j=1:t
            M_hat(i,j)=exp(a_x(i)+b_x(i)*k_t(j));
        end
    end
D_hat=zeros(x,t);
    for i=1:x
        for j=1:t
            D_hat(i,j)=E(i,j)*M_hat(i,j);
        end
    end

Err=M-M_hat;

index42=M(:,:)==42;
Err(index42)=0;
n_42=sum(sum(index42));

MaxErr=max(max(Err));
MinErr=min(min(Err));
[x_MaxErr,t_MaxErr]=find(Err==MaxErr);
[x_MinErr,t_MinErr]=find(Err==MinErr);
Err2=sum(sum(Err.^2))/(x*t-n_42);
MeanErr=sum(sum(Err))/(x*t-n_42);

H=D.*log(D./D_hat);
indexNaN=isnan(H);
H(indexNaN)=0;
L2=2* sum(sum(H ));

index42=(index42-1)*(-1);
Rx=zeros(x,1);
totalVar=0;
totalVar2=0;

%X=zeros(x,2);
for i=1:x
    h=logical(index42(i,:));
    k2=var(Err(i,h));
    k=var(M(i,h));
    Rx(i)= k2/k;
    totalVar=totalVar+k;
    totalVar2=totalVar2+k2;
    %X(i,:)=[k2,k]
end 
FVUx=Rx*100;
maxFVUx=max(FVUx);
[x_maxFVUx]=find(FVUx==maxFVUx);

FVU=totalVar2/totalVar*100;


end

