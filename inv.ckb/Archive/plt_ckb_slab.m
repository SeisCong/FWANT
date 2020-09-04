
close all; clear all;
MFILE_ROOT='../../../mfiles';
path([MFILE_ROOT '/fun-spool'],path);
path([MFILE_ROOT '/saclab'],path);

fnm_ckb='block_ckb_slabonly.dat'
mx=60; my=63; mz=29;
%mx=272; my=312; mz=72;

v=load(fnm_ckb);
v=v*100;

m=reshape(v,[mx,my,mz]);

cmin=min(v); cmax=max(v);

figure('Position',[400 400 900 600]); 
subplot(2,3,1), hold on, box on, axis on
set(gcf,'renderer','zbuffer');
imagesc(squeeze(m(:,:,19)));
%shading flat
shading interp;
colormap('jetwr');
caxis([cmin cmax])
colorbar
axis equal; axis tight

subplot(2,3,2), hold on, box on, axis on 
set(gcf,'renderer','zbuffer');
imagesc(squeeze(m(:,:,20)));
%shading flat
shading interp;
colormap('jetwr');
caxis([cmin cmax])
colorbar
axis equal; axis tight

subplot(2,3,3), hold on, box on, axis on
set(gcf,'renderer','zbuffer');
imagesc(squeeze(m(:,:,21)));
%shading flat
shading interp;
colormap('jetwr');
caxis([cmin cmax])
colorbar;
axis equal; axis tight

subplot(2,3,4), hold on, box on, axis on
set(gcf,'renderer','zbuffer');
imagesc(squeeze(m(:,:,22)));
%shading flat
shading interp;
colormap('jetwr');
caxis([cmin cmax])
colorbar;
axis equal; axis tight

subplot(2,3,5), hold on, box on, axis on
set(gcf,'renderer','zbuffer');
imagesc(squeeze(m(:,:,23)));
%shading flat
shading interp;
colormap('jetwr');
caxis([cmin cmax])
colorbar;
axis equal; axis tight

subplot(2,3,6), hold on, box on, axis on
set(gcf,'renderer','zbuffer');
imagesc(squeeze(m(:,:,24)));
%shading flat
shading interp;
colormap('jetwr');
caxis([cmin cmax])
colorbar;
axis equal; axis tight

