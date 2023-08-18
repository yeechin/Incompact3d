function [xx, yy, zz] = create_mesh(Lx,Lz,nx,nz,yp)

dx = Lx/nx; dz = Lz/nz;
x1d = (0:nx-1)*dx;
z1d = (0:nz-1)*dz;

[yy, xx, zz] = meshgrid(yp,x1d,z1d);
end