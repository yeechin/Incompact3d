function data2plot3d(num, nx, ny, nz)
% filename = num2str(num,'%05d');
fns = {'ux-','uy-','uz-','pp-','critq-','Rx-','Ry-','Rz-'};
nvar = length(fns);
fnout = sprintf('%s%05d%s','channel-',num,'.q');
fidout = fopen(fnout,'w');
fwrite(fidout,[16,nx,ny,nz,nvar,16],'integer*4');
fwrite(fidout,nx*ny*nz*nvar*8,'integer*4');
for i=1:length(fns)
    filename = sprintf('%s%05d%s',fns{i},num,'.bin');
    fid = fopen(filename,'r');
    f = fread(fid,nx*ny*nz,'real*8');
    fclose(fid);
    
    fwrite(fidout,f,'real*8');
end
fwrite(fidout,nx*ny*nz*nvar*8,'integer*4');
fclose(fidout);
end