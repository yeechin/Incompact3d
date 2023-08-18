function  write_restart_file(filename,ux,uy,uz,pp)
fid = fopen(filename,'w');
fwrite(fid,ux(:),'real*8');
fwrite(fid,uy(:),'real*8');
fwrite(fid,uz(:),'real*8');
fwrite(fid,pp(:),'real*8');
fclose(fid);
end
