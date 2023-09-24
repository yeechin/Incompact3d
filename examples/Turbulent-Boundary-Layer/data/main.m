clear; clc;
nx = 257; ny = 129; nz = 128;
Lx = 10; Ly = 10; Lz = 5; 
%%
% create plot3d mesh file
meshPlot3dCreate(Lx,Lz,nx,nz,'tbl-grid.x');

%%
% Change the data format to plot3d solution files
for i=100:25:200
    data2plot3d(i,nx,ny,nz);
end