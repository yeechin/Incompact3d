function [ux,uy,uz,pp] = read_restart_file(filename,nx,ny,nz)
fid = fopen(filename,'r');
ux = fread(fid,nx*ny*nz,'real*8');
uy = fread(fid,nx*ny*nz,'real*8');
uz = fread(fid,nx*ny*nz,'real*8');
pp = fread(fid,nx*(ny-1)*nz,'real*8');
fclose(fid);

ux = reshape(ux,nx,ny,nz);
uy = reshape(uy,nx,ny,nz);
uz = reshape(uz,nx,ny,nz);
pp = reshape(pp,nx,ny-1,nz);
end
