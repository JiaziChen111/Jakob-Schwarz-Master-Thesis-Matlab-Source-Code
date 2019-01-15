function [ Err2 ] = errorfkt3( M,M_hat )
%ERRORFKT Calculates statistics to compare different fittings
%   Defintions are given in the Masterarbeit chapter 2
n=size(M);
x=n(1);
t=n(2);


Err=M-M_hat;

index42=M(:,:)==42;
Err(index42)=0;
n_42=sum(sum(index42));


Err2=sum(sum(Err.^2))/(x*t-n_42);

end

