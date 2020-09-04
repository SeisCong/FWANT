
close all; clear all;

MFILE_ROOT='../../../mfiles';
path([MFILE_ROOT '/fun-spool'],path);
path([MFILE_ROOT '/saclab'],path);

%fnm_ckb='block_ckb_Moho_checkerboard_10x10_44_45_46.dat'
fnm_ckb='block_ckb_6x6x24.dat'
mx=160; my=106; mz=25;

v=load(fnm_ckb);
v=v*100;
m=reshape(v,[mx,my,mz]);
cmin=-10; cmax=10;
%cmin=min(v); cmax=max(v);

figure('Position',[400 400 900 900]); 
subplot(2,3,1), hold on, box on, axis on
set(gcf,'renderer','zbuffer');
imagesc(squeeze(m(:,:,25)));
axis([0 my 0 mx]);
%shading flat
shading interp;
colormap('jetwr');
caxis([cmin cmax]);
colorbar
title('(a) 7km');

subplot(2,3,2), hold on, box on, axis on 
set(gcf,'renderer','zbuffer');
imagesc(squeeze(m(:,:,23)));
axis([0 my 0 mx]);
%shading flat
shading interp;
colormap('jetwr');
caxis([cmin cmax])
colorbar
title('(a) 15km');

subplot(2,3,3), hold on, box on, axis on
set(gcf,'renderer','zbuffer');
imagesc(squeeze(m(:,:,21)));
axis([0 my 0 mx]);
%shading flat
shading interp;
colormap('jetwr');
caxis([cmin cmax])
colorbar;
title('(a) 25km');

subplot(2,3,4), hold on, box on, axis on
set(gcf,'renderer','zbuffer');
imagesc(squeeze(m(:,:,19)));
axis([0 my 0 mx]);
%shading flat
shading interp;
colormap('jetwr');
caxis([cmin cmax])
colorbar;
title('(a) 36km');

subplot(2,3,5), hold on, box on, axis on
set(gcf,'renderer','zbuffer');
imagesc(squeeze(m(:,:,16)));
axis([0 my 0 mx]);
%shading flat
shading interp;
colormap('jetwr');
caxis([cmin cmax])
colorbar;
title('(a) 52km');

subplot(2,3,6), hold on, box on, axis on
set(gcf,'renderer','zbuffer');
imagesc(squeeze(m(:,:,12)));
axis([0 my 0 mx]);
%shading flat
shading interp;
colormap('jetwr');
caxis([cmin cmax])
colorbar;
title('(a) 76km');