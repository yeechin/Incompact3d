function [xxi, yyi, zzi] = create_mesh_inter(Lx,Lz,nx,nz,ypi)

dx = Lx/nx; dz = Lz/nz;
x1d = ((0:nx-1)+0.5)*dx;
z1d = ((0:nz-1)+0.5)*dz;

[yyi, xxi, zzi] = meshgrid(ypi,x1d,z1d);
end