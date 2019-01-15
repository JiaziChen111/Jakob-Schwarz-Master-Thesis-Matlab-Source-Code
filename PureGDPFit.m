function [y_0x,y_1x]=PureGDPFit(lnM,g_t)
%PUREDGPFIT fits the model: ln m(x,t)=y_0x+y_1x*g_t+e_x,t

%this model is for fixed x a simple regresion model: Y= X * b
%with Y=ln m(x,1:t), X=[ones, g_t], b=[y_0x,y_1x]

n=size(lnM);
x=n(1);
t=n(2);
y_0x=zeros(x,1);
y_1x=zeros(x,1);
X=[ones(t,1),g_t];
for i= 1:x
    Y=lnM(i,:)'; %all years data for age x   
    b=X\Y;
    y_0x(i)=b(1);
    y_1x(i)=b(2);  
end

end

