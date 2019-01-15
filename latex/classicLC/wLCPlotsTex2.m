index42=M(:,:)==42;

hhh=figure;
subplot(2,2,1);
plot(linspace(1,x,x),a_wLC(1:x),'r');
hold on;
plot(linspace(1,x,x),a_DwLC(1:x),'k');
hold on;
plot(linspace(1,x,x),a_W(1:x),'c');
title('a_x');

subplot(2,2,2);
plot(linspace(1,x,x),b_wLC,'r');
hold on;
plot(linspace(1,x,x),b_DwLC,'k');
hold on;
plot(linspace(1,x,x),b_W,'c');
title('b_x');

subplot(2,2,3:4);
plot(linspace(1,t,t),k_wLC,'r');
hold on;
plot(linspace(1,t,t),k_DwLC,'k');
hold on;
plot(linspace(1,t,t),k_W,'c');

title('k_t');
legend('1-0-weights','Death-weights','w_{x,t}=t+x_{max}-x');
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhh,'axbxkt_w','-dpdf','-r0')

%%%%%%%
M_wLC=zeros(x,t);
    for i=1:x
        for j=1:t
            M_wLC(i,j)=exp(a_wLC(i)+b_wLC(i)*k_wLC(j));
        end
    end
    
M_D=zeros(x,t);
    for i=1:x
        for j=1:t
            M_D(i,j)=exp(a_DwLC(i)+b_DwLC(i)*k_DwLC(j));
        end
    end    
  M_wLC_W=zeros(x,t);
    for i=1:x
        for j=1:t
            M_wLC_W(i,j)=exp(a_W(i)+b_W(i)*k_W(j));
        end
    end
M_wLC2_W=zeros(x,t);
    for i=1:x
        for j=1:t
            M_wLC2_W(i,j)=exp(a_W2(i)+b_W2(i)*k_W2(j));
        end
    end
%%%%%%%%%%%
errhwLC=figure;
subplot(2,2,1)
Merr= M- M_wLC ;
Merr(index42)=0;
xscale = [1956 2015];
yscale = [0 104];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.5 1]);
xlabel('year t')
ylabel('age x')
title(' ')

subplot(2,2,2)
xscale = [1956 2015];
yscale = [0 104];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.05 0.075]);
xlabel('year t')
ylabel('age x')
title(' ')

subplot(2,2,3)
xscale = [1956 2015];
yscale = [0 104];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.005 0.0075]);
xlabel('year t')
ylabel('age x')

subplot(2,2,4)

xscale = [1956 2015];
yscale = [0 104];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.0005 0.00075]);
xlabel('year t')
ylabel('age x')

suptitle('Fitting error colormap with 1-0-weights')
set(errhwLC,'Units','Inches');
pos = get(errhwLC,'Position');
set(errhwLC,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(errhwLC,'errhwLC','-dpdf','-r0')

%%%%%%%%%%
errhwLC_W=figure;
subplot(2,2,1)
Merr= M- M_wLC_W ;
Merr(index42)=0;
xscale = [1956 2015];
yscale = [0 104];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.5 1]);
xlabel('year t')
ylabel('age x')
title('')

subplot(2,2,2)
xscale = [1956 2015];
yscale = [0 104];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.05 0.075]);
xlabel('year t')
ylabel('age x')
title('')

subplot(2,2,3)
xscale = [1956 2015];
yscale = [0 104];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.005 0.0075]);
xlabel('year t')
ylabel('age x')

subplot(2,2,4)
xscale = [1956 2015];
yscale = [0 104];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.0005 0.00075]);
xlabel('year t')
ylabel('age x')
suptitle('Fitting error colormap with w_{x,t}=t+x_{max}-x')

set(errhwLC_W,'Units','Inches');
pos = get(errhwLC_W,'Position');
set(errhwLC_W,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(errhwLC_W,'errhwLC_W','-dpdf','-r0')

%%%%%%%%%%
errhwLC_D=figure;
subplot(2,2,1)
Merr= M- M_D ;
Merr(index42)=0;
xscale = [1956 2015];
yscale = [0 104];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.5 1]);
xlabel('year t')
ylabel('age x')
title('')

subplot(2,2,2)
xscale = [1956 2015];
yscale = [0 104];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.05 0.075]);
xlabel('year t')
ylabel('age x')
title('')

subplot(2,2,3)
xscale = [1956 2015];
yscale = [0 104];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.005 0.0075]);
xlabel('year t')
ylabel('age x')

subplot(2,2,4)
xscale = [1956 2015];
yscale = [0 104];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.0005 0.00075]);
xlabel('year t')
ylabel('age x')

suptitle('Fitting error colormap with death-weights')
set(errhwLC_D,'Units','Inches');
pos = get(errhwLC_D,'Position');
set(errhwLC_D,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(errhwLC_D,'errhwLC_D','-dpdf','-r0')

%%%%%%%%%%
errhwLC_W=figure;
subplot(2,2,1)
Merr= M- M_wLC2_W ;
Merr(index42)=0;
xscale = [1956 2015];
yscale = [0 104];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.5 1]);
xlabel('year t')
ylabel('age x')
title('')

subplot(2,2,2)
xscale = [1956 2015];
yscale = [0 104];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.05 0.075]);
xlabel('year t')
ylabel('age x')
title('')

subplot(2,2,3)
xscale = [1956 2015];
yscale = [0 104];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.005 0.0075]);
xlabel('year t')
ylabel('age x')

subplot(2,2,4)
xscale = [1956 2015];
yscale = [0 104];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.0005 0.00075]);
xlabel('year t')
ylabel('age x')
suptitle('Fitting error colormap with w_{x,t}=t+x_{max}-x after second step')

set(errhwLC_W,'Units','Inches');
pos = get(errhwLC_W,'Position');
set(errhwLC_W,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(errhwLC_W,'errhwLC2_W','-dpdf','-r0')

%%%%%%%%%%%%
Merr=abs( M- M_wLC_W );
Merr(index42)=0;
max(max(Merr))
sum(sum(Merr))/(x*t)
sqrt(sum(sum(Merr.^2)))/(x*t)

Merr= abs( M- M_D );
Merr(index42)=0;
max(max(Merr))
sum(sum(Merr))/(x*t)
sqrt(sum(sum(Merr.^2)))/(x*t)

Merr=abs( M- M_wLC2_W );
Merr(index42)=0;
max(max(Merr))
sum(sum(Merr))/(x*t)
sqrt(sum(sum(Merr.^2)))/(x*t)

HHH=M_wLC_W-M_wLC2_W;