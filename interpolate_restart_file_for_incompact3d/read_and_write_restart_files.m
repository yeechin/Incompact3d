Lx = 8; Lz = 4;
nx = 256; nz = 128;

yp = load('yp1.dat'); yp = yp';
ny = length(yp);

[xx1,yy1,zz1] = create_mesh(Lx,Lz,nx,nz,yp);

ypi = load('ypi1.dat');ypi = ypi';
[xxi1,yyi1,zzi1] = create_mesh_inter(Lx,Lz,nx,nz,ypi);

%%
[ux1,uy1,uz1,pp1] = read_restart_file('restart0080000',nx,ny,nz);

%%
F1 = griddedInterpolant(xx1,yy1,zz1,ux1);
F2 = griddedInterpolant(xx1,yy1,zz1,uy1);
F3 = griddedInterpolant(xx1,yy1,zz1,uy1);
F4 = griddedInterpolant(xxi1,yyi1,zzi1,pp1);
%%
Lx = 4;  Lz = 2;
nx = 256; nz = 720;

yp = load('yp2.dat'); yp = yp';
ny = length(yp);

[xx2,yy2,zz2] = create_mesh(Lx,Lz,nx,nz,yp);

ypi = load('ypi2.dat');ypi = ypi';
[xxi2,yyi2,zzi2] = create_mesh_inter(Lx,Lz,nx,nz,ypi);

%%
% interpolatation
ux2 = F1(xx2,yy2,zz2);
uy2 = F2(xx2,yy2,zz2);
uz2 = F3(xx2,yy2,zz2);
pp2 = F4(xxi2,yyi2,zzi2);

%%
write_restart_file('restart0000000',ux2,uy2,uz2,pp2);