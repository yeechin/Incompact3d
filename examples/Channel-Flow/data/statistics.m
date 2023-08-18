clear; clc;
nx = 256; ny = 129; nz = 128;nvar = 8;
Lx = 8; Ly = 2; Lz = 4; 
%%

count = 0;
% Change the data format to plot3d solution files
for i=50010:10:60000
    i
    fin = sprintf('%s%05d%s','channel-',i,'.q');
    fid = fopen(fin,'r');
    dum = fread(fid,7,'integer*4');
    f = fread(fid,nx*ny*nz*nvar,'real*8');
    fclose(fid);
    
    f(isnan(f))=0;
    
%     if(sum(isnan(f))>0)
%         disp('nan numbers found');
%     end 
    f = reshape(f,nx*ny*nz,nvar);
    
    f = f(:,[1:3,6:8]);
    ff = cat(2,f(:,1:3).^2,f(:,1).*f(:,2),...
            f(:,2).*f(:,3),f(:,3).*f(:,1));
    RR = cat(2,f(:,4:6).^2,f(:,4).*f(:,5),...
            f(:,5).*f(:,6),f(:,6).*f(:,4));
        
    f = reshape(f,nx,ny,nz,6);
    ff = reshape(ff,nx,ny,nz,6);
    RR = reshape(RR,nx,ny,nz,6);
    
    f = squeeze(mean(mean(f,3),1));
    ff = squeeze(mean(mean(ff,3),1));
    RR = squeeze(mean(mean(RR,3),1));
    
    f = 0.5*(f+f(end:-1:1,:));
    f = f(1:(end+1)/2,:);
    ff = 0.5*(ff+ff(end:-1:1,:));
    ff = ff(1:(end+1)/2,:);
    RR = 0.5*(RR+RR(end:-1:1,:));
    RR = RR(1:(end+1)/2,:);

    if i==50010
        fmean = f;
        ffmean = ff;
        RRmean = RR;
        count = count + 1;
    else
        fmean = fmean + f;
        ffmean = ffmean + ff;
        RRmean = RRmean + RR;
        count = count + 1;
    end
end

fmean = fmean/count;
ffmean = ffmean/count;
RRmean = RRmean/count;
%%
load('statistics.mat');
yp = load('yp.dat');
yp = yp(1:(end+1)/2);
Re = 4200;
nu = 1/Re;
%%

ustar = sqrt(nu*(fmean(2,1)-fmean(1,1))/(yp(2)-yp(1)));
ystar = nu/ustar;

figure;semilogx(yp/ystar,fmean(:,1)/ustar)
hold on