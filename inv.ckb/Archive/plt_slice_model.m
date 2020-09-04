
clear all;
close all;

MFILE_ROOT='../../../mfiles';
path([MFILE_ROOT '/fun-spool'],path);
path([MFILE_ROOT '/saclab'],path);

% -- parameter --
mx=124; my=104; mz=48;
fnm_blk='../block.124x104x48.1x1x1.1x1x1.nc'

smot_list=[2];
damp_list=[4];
%zval_list=[ 1 ];
zindx_list=[22 21 20];
res_list={'./result.1th.model'};
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
xc=xc/1e3; yc=yc/1e3;; zc=zc/1e3; %in km

%minlon=235; maxlon=250; minlat=38;maxlat=49;
minlon=236; maxlon=245;minlat=39; maxlat=48.8;

hid=figure('Position',[200 400 600 500]); hold on, box on, axis on
set(hid,'renderer','zbuffer');
set(gca,'TickDir','out');

% figwid=0.4; figheight=0.28;

Vmean=zeros(mz,1);
Vmean(20:28)=[4.3642 4.3220 4.3550 4.4505 4.3827 3.9734 3.7648 3.6179 3.5672];
% ----------------------------------
%for nres=1:num_res
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
% z slice
for ncmp=2
    KCMPNM=cmp_list{ncmp};
	fnm_try=[pnm_result '/try.damp' KDAMPNM '.smot' KSMOTNM '.st0.ev0.lo0.dat']
	M=load(fnm_try);
	W=reshape(M,mx,my,mz,2);
    dVel=squeeze(W(:,:,:,ncmp))*100;
%     for xx=1:mx
%     for yy=1:my
%     Vel(xx,yy,:)=(1+squeeze(dVel(xx,yy,:))).*Vmean;
%     end
%     end
    
for nk=1:num_z
%    for nk=8
k=zindx_list(nk);
%d=zval_list(nk);
V=squeeze(dVel(:,:,k));
% remove the mean, 04/26/2011, Y.Shen
%V=V-mean(mean(V));

XS=squeeze(x(:,:,k));
YS=squeeze(y(:,:,k));
zs=squeeze(zc(:,:,k));
%KDNM=num2str(zval_list(nk));
%KINM=num2str(zindx_list(nk));
Vmax=max(max(abs(V)));

if(ncmp == 1)
fnm_fig=['figure_Vp_' int2str(nk) '.eps']
fnm_fig_png=['figure_Vp_' int2str(nk) '.png']
fnm_outfile=['Vp_' int2str(nk) '.txt']
Vmax=15;
end;
if(ncmp == 2)
fnm_fig=['figure_Vs_' int2str(nk) '.eps']
fnm_fig_png=['figure_Vs_' int2str(nk) '.png']
fnm_outfile=['Vs_' int2str(nk) '.txt']
Vmax=15;
end;

fnm_ckb='block_ckb_LVV_ite04.dat'
v=load(fnm_ckb);

m=reshape(v,[mx,my,mz])*100;

subplot(2,num_z,nk);
%subplot('Position',[0.05 0.05+(nk-1)*(figheight+0.03) figwid figheight])
hold on, box on, axis on,
set(gcf,'renderer','zbuffer');
set(gca,'TickDir','out');
pcolor(YS,XS,squeeze(m(:,:,k)));

state=[];
load ../../../misc/us_states;
for sb=1:length(state)
    plot(state(sb).polygon(:,1)+360, state(sb).polygon(:,2),'color','k','LineWidth',1);
end

[clon,clat]=textread('../../../misc/PNWcoast_revised.dat','%f %f');
plot(clon+360,clat,'k.','MarkerSize',3);

daspect([1,cosd((minlat+maxlat)/2),1]);
%xlabel('Longitude, deg.'); ylabel('Latitude, deg.');
axis([minlon maxlon minlat maxlat]);

shading interp;
colormap('jetwr');
%ch=colorbar('location','eastoutside','plotboxaspectratio',[0.3 4.0 1]);
caxis([-15 15])
daspect([1,cosd((minlat+maxlat)/2),1]);
title(['Input at ' num2str(round(z(round(mx/2),round(my/2),k))),' km '],'FontSize',12);
set(gca,'XTick',[minlon:4:maxlon],'XTickLabel',[minlon:4:maxlon]);
set(gca,'YTick',[39:4:47],'YTickLabel',[39:4:47]);

subplot(2,num_z,num_z+nk)
%subplot('Position',[0.55 0.05+(nk-1)*(figheight+0.03) figwid figheight])
hold on, box on, axis on
set(gca,'TickDir','out');
set(gca,'TickDir','out');

pcolor(YS,XS,V);

state=[];
load ../../../misc/us_states;
for sb=1:length(state)
    plot(state(sb).polygon(:,1)+360, state(sb).polygon(:,2),'color','k','LineWidth',1);
end

[clon,clat]=textread('../../../misc/PNWcoast_revised.dat','%f %f');
plot(clon+360,clat,'k.','MarkerSize',3);

shading interp
colormap('jetwr');
caxis([-10 10])

daspect([1,cosd((minlat+maxlat)/2),1]);
%xlabel('Longitude, deg.'); ylabel('Latitude, deg.');
axis([minlon maxlon minlat maxlat]);
%colorbar

%title([cmp_list{ncmp} ' at depth = ',num2str(round(z(round(mx/2),round(my/2),k))),' km'])
%title(['Recovered at ',num2str(round(z(round(mx/2),round(my/2),k))),' km '],'FontSize',12)
title('Recovered ','FontSize',12)
set(gca,'XTick',[minlon:4:maxlon],'XTickLabel',[minlon:4:maxlon]);
set(gca,'YTick',[39:4:47],'YTickLabel',[39:4:47]);
%ch=colorbar('location','eastoutside','plotboxaspectratio',[0.3 4.0 1]);
%Vmin=3.7; Vmax=4.9;
%caxis([Vmin,Vmax]);
end
end

end
end
end
ch=colorbar('Position',[0.92 0.3 0.02 0.45]);
caxis([-10 10])
set(ch,'FontSize',12)

%%%% save figure
set(gcf,'PaperPositionMode','auto');   
figname = ['ResTest_model_vpvs.eps'];

eval(['print -depsc ' figname])

