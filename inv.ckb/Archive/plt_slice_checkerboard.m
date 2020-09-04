
clear all;
close all;

MFILE_ROOT='../../../mfiles';
path([MFILE_ROOT '/fun-spool'],path);
path([MFILE_ROOT '/saclab'],path);

% -- parameter --
mx=50;my=58;mz=20;
fnm_blk='../block.50x58x20.1x1x1.1x1x1.nc'

smot_list=[4];
damp_list=[4];
%zval_list=[ 1 ];
zindx_list=[11 12 13 14 15];
%zindx_list=[11 13 14 16 17];
res_list={'./result.1th.4x4x10/nosynerr'};
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

minlon=236; maxlon=246; minlat=39;maxlat=48;
%minlon=234; maxlon=251;minlat=39; maxlat=50;

hid=figure('Position',[400 400 800 550]); hold on, box on, axis on
set(hid,'renderer','zbuffer');
set(gca,'TickDir','out');

figwid=0.23; figheight=0.35;
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

    fnm_try=[pnm_result '/try.damp' KDAMPNM '.smot' KSMOTNM '.st0.ev0.lo0.dat']

    M=load(fnm_try);
    W=reshape(M,mx,my,mz,2);
    W=W*100; %percent

% z slice
for ncmp=2
    KCMPNM=cmp_list{ncmp};
    Vel=squeeze(W(:,:,:,ncmp));

for nk=1:num_z
%    for nk=8
k=zindx_list(nk);
%d=zval_list(nk);
V=squeeze(Vel(:,:,k));
% remove the mean, 04/26/2011, Y.Shen
V=V-mean(mean(V));
		
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

%subplot(1,3,ncmp+1); 
% if nk==1&ncmp==1, continue, end
% if ncmp==1, subplot(2,3,nk); 
% else, subplot(2,3,nk+3), end

subplot(2,3,nk+1)
% if nk<=2
% subplot('Position',[0.05+nk*0.3 0.55 figwid figheight])
% elseif nk>=3
% subplot('Position',[0.05+(nk-3)*0.3 0.07 figwid figheight])
% end
hold on, box on, axis on
set(gca,'TickDir','out');

pcolor(YS,XS,V);

[clon,clat]=textread('/home/haiying/Proj/cascadia/misc/PNWcoast_revised.dat','%f %f');
plot(clon+360,clat,'k.','MarkerSize',2);

state=[];
load /home/haiying/Proj/cascadia/misc/us_states.mat
for sb=1:length(state)
    plot(state(sb).polygon(:,1)+360, state(sb).polygon(:,2),'color','k','LineWidth',1);
end

shading interp
colormap('jetwr');

daspect([1,cosd((minlat+maxlat)/2),1]);
%xlabel('Longitude, deg.'); ylabel('Latitude, deg.');
axis([minlon maxlon minlat maxlat]);
%colorbar

title(['Resolved ' cmp_list{ncmp} ' at ',num2str(round(z(round(mx/2),round(my/2),k))),' km'])
caxis([-5,5]);

end
end

end
end
end

% plot input model
%fnm_ckb='./block_ckb_ite_01_dVelS.dat'
%fnm_ckb='./result.1th.4x4x10/block_ckb.dat'
fnm_ckb='block_ckb_4x4x10.dat'

v=load(fnm_ckb);
v=v*100;
m=reshape(v,[mx,my,mz]);

subplot('Position',[0.14 0.52 0.27 0.46]); 
%subplot(2,3,1)
hold on, box on, axis on,
set(gcf,'renderer','zbuffer');
pcolor(YS,XS,squeeze(m(:,:,1)));
[clon,clat]=textread('/home/haiying/Proj/cascadia/misc/PNWcoast_revised.dat','%f %f');
plot(clon+360,clat,'k.','MarkerSize',2);
state=[];
load /home/haiying/Proj/cascadia/misc/us_states.mat
for sb=1:length(state)
    plot(state(sb).polygon(:,1)+360, state(sb).polygon(:,2),'color','k','LineWidth',1);
end
daspect([1,cosd((minlat+maxlat)/2),1]);
%xlabel('Longitude, deg.'); ylabel('Latitude, deg.');
axis([minlon maxlon minlat maxlat]);

shading interp;
colormap('jetwr');
ch=colorbar('location','eastoutside','plotboxaspectratio',[0.2 4.0 1]);
caxis([-5 5])
daspect([1,cosd((minlat+maxlat)/2),1]);
title('input velocity anomaly model');

%%%% save figure
set(gcf,'PaperPositionMode','auto');   
figname = ['ResTest_checkerboard_0.05deg_100km.eps'];

eval(['print -depsc ' figname])


