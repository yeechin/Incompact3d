yp = load('yp.dat');
nu = 1/4200;
%%
addpath legendflex-pkg-master/legendflex/
addpath legendflex-pkg-master/setgetpos_V1.2/
%%
nx = 64;ny= 33;nz = 32;
cnt = 0;
for i = 121:130
    fn = sprintf('nut_liutex-%03d.bin',i);
    fid = fopen(fn,'r');
    nut = fread(fid,nx*ny*nz,'real*8');
    fclose(fid);
    nut = reshape(nut,nx,ny,nz);
    if(i==121)
        meannut = squeeze(mean(mean(nut,3),1));
    else
        meannut = meannut + squeeze(mean(mean(nut,3),1));
    end
    cnt = cnt + 1;
end
meannut = meannut /cnt;
meannut = 0.5*(meannut + meannut(end:-1:1));

%%
figure;plot(meannut,yp)
%%
save('meannut-liutex2','yp','meannut');

%%
fig = figure;
load('meannut-smagorinsky.mat')
plot(meannut/nu,yp,'b-.','LineWidth',2)

%
% '--','LineWidth',2,'color',[0.5 0.5 0.5])
hold on
grid on
grid minor
%
load('meannut-liutex2.mat')
plot(meannut/nu,yp,'r-','LineWidth',2)

%
% load('liutex_les2_results.mat')
% semilogx(ypp,umeanp,'r-','LineWidth',2)
%
load('meannut-wale.mat')
plot(meannut/nu,yp,'k--','LineWidth',2)
%
% load('liutex_les2_results.mat')
% semilogx(ypp,umeanp,'r-','LineWidth',2)
%%
set(gca,'FontSize',20);
t=xlabel('$$\left<\nu_t\right>/\nu$$','Interpreter','Latex','FontSize',30);
t=ylabel('$${y}/{h}$$','Interpreter','Latex','FontSize',30);
% set(get(gca,'YLabel'),'Rotation',0);
% pos=get(t,'position');
% pos(1)=pos(1)-0.1;
% set(t,'position',pos);

axis([0 1.8 0 2])
% yticks([0 5 10 15 20 25])
%%
legendflex(gca,{...
    'Smagorinsky','Liutex-based model','WALE'},...
    'interpreter','Latex',...
    'xscale',2,'box','off','FontSize',18,'anchor',{'nw','nw'},...
    'buffer',[-5 -30])


%%
print(gcf,'nut-comparison-Liutex-name','-dpng')


%%
% a sub figure to y = 0
fig = figure;
% load('meannut-smagorinsky.mat')
% plot(meannut(1:10)/nu,yp(1:10),'b-.','LineWidth',2)

%
% '--','LineWidth',2,'color',[0.5 0.5 0.5])
hold on
grid on
grid minor
%
load('meannut-liutex2.mat')
plot(meannut(1:10)/nu,yp(1:10),'r-','LineWidth',2)

%
% load('liutex_les2_results.mat')
% semilogx(ypp,umeanp,'r-','LineWidth',2)
%
load('meannut-wale.mat')
plot(meannut(1:10)/nu,yp(1:10),'k--','LineWidth',2)
%
% load('liutex_les2_results.mat')
% semilogx(ypp,umeanp,'r-','LineWidth',2)

%%
set(gca,'FontSize',20);
t=xlabel('$$\left<\nu_t\right>/\nu$$','Interpreter','Latex','FontSize',30);
t=ylabel('$${y}/{h}$$','Interpreter','Latex','FontSize',30);
% set(get(gca,'YLabel'),'Rotation',0);
% pos=get(t,'position');
% pos(1)=pos(1)-0.1;
% set(t,'position',pos);

axis([0 1.2 0 0.25])
 yticks([0 0.1 0.2])
 xticks([0 0.4 0.8 1.2])
%%
legendflex(gca,{...
    'Liutex-based model','WALE'},...
    'interpreter','Latex',...
    'xscale',2,'box','off','FontSize',18,'anchor',{'nw','nw'},...
    'buffer',[-5 -30])
%%
print(gcf,'nut-comparison-part','-dpng')