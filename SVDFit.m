function [ a_x b_x k_t ] = SVDFit( lnM )
%LeeCarterFit Calculates the least squarres estiamtes of the Lee Carter Fit
%   Model ln m(x,t)= a_x + b_x * k_t
%   lnM is the matrix ln m(x,t) of the strictly positive central death rates
%   a_x as average over all periods t for age x a= 1/T ?(1,t) ln m(x,t)
%   b_x and k_t as least squarres approximation
%   least squarres solution ||A-X||_2 for all matrixes X of rank p is
%   given by ?(1,p) sigma_i * u_i * v_i' of SVD U*?*V', vgl. numerik 3
%   k_t*b_x = X has rank 1 
%   therefor best L2 approxiamtion is given by u_1, sigma_1* v_1'

n=size(lnM);
x=n(1);
t=n(2);

a_x=zeros(x,1);
for i=1:x
   for j=1:t
          a_x(i)=a_x(i)+lnM(i,j);
   end
   a_x(i)=a_x(i)/t;
end

%Z_xt= ln m(x,t)-a_x=b_x*k_t
Z=lnM-repmat(a_x,1,t);
[U,S,V] = svd(Z);
k_t=S(1,1)*V(:,1);
b_x=U(:,1);


%sum(b_x.^2)=1 is true; 
%lee carter had it without .^2
 %re normalize such that sum b =1
c=sum(b_x);
b_x=b_x/c;
k_t=k_t*c;

%sum(k_t)=0 is true somehow anyway else do
c=-sum(k_t)/t;
a_x=a_x-b_x*c;
k_t=k_t+c;

end

