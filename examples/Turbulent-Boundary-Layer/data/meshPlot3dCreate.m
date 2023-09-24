% Periodic in x and z, mesh points in y should be read in
function meshPlot3dCreate(Lx,Lz,nx,nz,meshname) 
% Lx = 8; Lz = 4; nx = 256; nz = 128;
dx = Lx/(nx-1); dz = Lz/nz;
x1d = (0:nx-1)*dx;
z1d = (0:nz-1)*dz;
y1d = load('yp.dat'); y1d = y1d';
ny = length(y1d);
%ny=129;Ly=20;dy=Ly/(ny-1);
%y1d=(0:ny-1)*dy;

[y3d, x3d, z3d] = meshgrid(y1d,x1d,z1d);

fid = fopen(meshname,'w');
fwrite(fid,[12,nx,ny,nz,12],'integer*4');
fwrite(fid,nx*ny*nz*3*8,'integer*4');
fwrite(fid,[x3d(:);y3d(:);z3d(:)],'real*8');
fwrite(fid,nx*ny*nz*3*8,'integer*4');
fclose(fid);

end