
clear all;
close all;

run ~/FWT/ANT/Proj/cascadia/model_updates/set_netcdf.m

% MFILE_ROOT='../../../mfiles';
% path([MFILE_ROOT '/fun-spool'],path);
% path([MFILE_ROOT '/saclab'],path);

% -- parameter --
mx=76;my=50;mz=40;
fnm_blk='../block.76x50x40.1x1x1.1x1x1.nc';
smot_list=[2];
damp_list=[2];
zindx_list=[37 35 33 31 29 27 25];
res_list={'./result.1th.4x4x20'};
cmp_list={'Vp','Vs'};

% ----------------------------------
num_z=numel(zindx_list);
num_smot=numel(smot_list);
num_damp=numel(damp_list);
num_res=numel(res_list);

X=nc_varget(fnm_blk,'x');
Y=nc_varget(fnm_blk,'y');
Z=nc_varget(fnm_blk,'z');

x=90-reshape(X,[mx,my,mz])/pi*180;
y=reshape(Y,[mx,my,mz])/pi*180;
R=reshape(Z,[mx,my,mz]);
z=(6371-reshape(Z,[mx,my,mz])/1e3);

X=reshape(X,[mx,my,mz]);
Y=reshape(Y,[mx,my,mz]);
% convert to Cartesian coord
[yc,xc,zc]=sph2cart(Y,pi/2-X,R);
xc=xc/1e3; yc=yc/1e3; zc=zc/1e3; %in km

minlon=284; maxlon=292; minlat=39;maxlat=47;

hid=figure('Position',[100 100 1200 500]); hold on, box on, axis on
set(hid,'renderer','zbuffer');
set(gca,'TickDir','out');

figlabel={'(a) ','(b) ','(c) ','(d) ','(e) ','(f) ','(g) ','(h) '};
figwid=0.23; figheight=0.35;
% ----------------------------------
for nres=1
    pnm_result=res_list{nres};
    pnm_fig=[pnm_result '.fig'];
    if ~ isdir(pnm_fig)
       mkdir(pnm_fig);
    end
for nsmot=1:num_smot
    KSMOTNM=num2str(smot_list(nsmot));
for ndamp=1:num_damp
    KDAMPNM=num2str(damp_list(ndamp));

    fnm_try=[pnm_result '/try.damp' KDAMPNM '.smot' KSMOTNM '.st0.ev0.lo0.dat']

    M=load(fnm_try);
    W=reshape(M,mx,my,mz,2);
    W=W*100; %percent

% z slice
for ncmp=2
    KCMPNM=cmp_list{ncmp};
    Vel=squeeze(W(:,:,:,ncmp));

for nk=1:num_z
k=zindx_list(nk);
V=squeeze(Vel(:,:,k));
		
XS=squeeze(x(:,:,k));
YS=squeeze(y(:,:,k));
zs=squeeze(zc(:,:,k));
Vmax=max(max(abs(V)));

if(ncmp == 1)
fnm_fig=['figure_Vp_' int2str(nk) '.eps']
fnm_fig_png=['figure_Vp_' int2str(nk) '.png']
fnm_outfile=['Vp_' int2str(nk) '.txt']
Vmax=5;
VelP=Vel;
end;
if(ncmp == 2)
fnm_fig=['figure_Vs_' int2str(nk) '.eps']
fnm_fig_png=['figure_Vs_' int2str(nk) '.png']
fnm_outfile=['Vs_' int2str(nk) '.txt']
Vmax=5;
VelS=Vel;
end;

fnm_ckb='block_ckb_4x4x20.dat'

v=load(fnm_ckb);
v=v*100;
m=reshape(v,[mx,my,mz]);

subplot(2,4,1)
hold on, box on, axis on,
set(gca,'TickDir','out');
set(gcf,'renderer','zbuffer');
pcolor(YS,XS,squeeze(m(:,:,k)));

state=[];
load ~/FWT/ANT/Proj/cascadia/misc/us_states.mat
for sb=1:length(state)
    plot(state(sb).polygon(:,1)+360, state(sb).polygon(:,2),'color','k','LineWidth',1);
end
daspect([1,cosd((minlat+maxlat)/2),1]);
axis([minlon maxlon minlat maxlat]);

shading interp;
colormap('jetwr');
% ch=colorbar('location','eastoutside','plotboxaspectratio',[0.2 4.0 1]);
caxis([-10 10])
daspect([1,cosd((minlat+maxlat)/2),1]);
title([char(figlabel{1}) 'Input Vs model'],'FontSize',18)

subplot(2,4,nk+1)

hold on, box on, axis on
set(gca,'TickDir','out');

pcolor(YS,XS,V);

state=[];
load ~/FWT/ANT/Proj/cascadia/misc/us_states.mat
for sb=1:length(state)
    plot(state(sb).polygon(:,1)+360, state(sb).polygon(:,2),'color','k','LineWidth',1);
end

shading interp
colormap('jetwr');

daspect([1,cosd((minlat+maxlat)/2),1]);
axis([minlon maxlon minlat maxlat]);

title([char(figlabel{nk+1}) 'Vs at ' num2str(round(z(round(mx/2),round(my/2),k))),' km'],'FontSize',18);
caxis([-8 8]);
% if nk<=3
% caxis([-8 8]);
% else
% caxis([-5 5]);
% end

set(gca,'XTick',[round(minlon):3:round(maxlon)],'XTickLabel',[round(minlon):3:round(maxlon)])
set(gca,'YTick',[round(minlat):2:round(maxlat)],'YTickLabel',[round(minlat):2:round(maxlat)])
set(gca,'TickDir','out');

end
end

end
end
end

ch=colorbar('Position',[0.92 0.3 0.01 0.45]);
caxis([-10 10])
set(ch,'FontSize',14)

%%%% save figure
set(gcf,'PaperPositionMode','auto');   
figname = ['ResTest_checkerboard_0.025deg_5x5x20.eps'];

eval(['print -depsc ' figname])


