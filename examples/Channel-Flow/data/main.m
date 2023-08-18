clear; clc;
nx = 256; ny = 129; nz = 128;
Lx = 8; Ly = 2; Lz = 4; 
%%
% create plot3d mesh file
meshPlot3dCreate(Lx,Lz,nx,nz,'channel-grid.x');

%%
% Change the data format to plot3d solution files
for i=50010:10:60000
    data2plot3d(i,nx,ny,nz);
end