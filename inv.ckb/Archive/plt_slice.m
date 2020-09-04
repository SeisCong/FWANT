clear all
%close all

MFILE_ROOT='../../../mfiles';
path([MFILE_ROOT '/fun-spool'],path);
path([MFILE_ROOT '/saclab'],path);

% -- parameter --
%mx=49; my=45; mz=26;
%mx=47; my=90; mz=17;
%mx=53;my=92;mz=31;
%mx=106;my=183;mz=21;
%mx=69;my=74;mz=21;
%mx=114; my=122; mz=24;
%mx=69;my=74;mz=23;
mx=69;my=74;mz=24;

%fnm_blk='../block.49x45x26.1x1x1.2x2x1.nc'
%fnm_blk='../block.47x90x17.1x1x1.1x1x1.nc'
%fnm_blk='../block.53x92x31.1x1x1.1x1x1.nc'
%fnm_blk='../block.106x183x21.1x1x1.1x1x1.nc'
%fnm_blk='../block.69x74x21.1x1x1.1x1x1.nc'
fnm_blk=['../block.' num2str(mx) 'x' num2str(my) 'x' num2str(mz) '.1x1x1.1x1x1.nc'];

%smot_list=[64];
smot_list=[8];
damp_list=[8];
%zval_list=[4, 20, 52, 94.7];
%zval_list=[ 1 2 3 4 5 6 7 8];
zval_list=[ 1 ];
%zindx_list=[17,15,13,11,9,7,5,3];
%zindx_list=[27];
zindx_list=[16]; % 117 km
%zindx_list=[23];
%zindx_list=[21]; % 256 km
%zindx_list=[17]; % 443 km 
%zindx_list=[15]; % 548 km
%zindx_list=[11]; % 776 km
%res_list={'./result.1th','./result.2th'};
res_list={'./result.1th'};
%res_list={'./result.2th'};
%res_list={'../Freq123.ENZ.dt.lt.10.model.wg1/result.1th','./result.2th'};
cmp_list={'Vp','Vs'};

% ----------------------------------
num_z=numel(zindx_list);
num_smot=numel(smot_list);
num_damp=numel(damp_list);
num_res=numel(res_list);

X=nc_varget(fnm_blk,'x');
Y=nc_varget(fnm_blk,'y');
Z=nc_varget(fnm_blk,'z');
%x=reshape(X,[mx,my,mz])/1000;
%y=reshape(Y,[mx,my,mz])/1000;
%D=(-reshape(Z,[mx,my,mz])/1e3);

%x=permute(x,[2,1,3]);
%y=permute(y,[2,1,3]);
%D=permute(D,[2,1,3]);

x=90-reshape(X,[mx,my,mz])/pi*180;
y=reshape(Y,[mx,my,mz])/pi*180;
R=reshape(Z,[mx,my,mz]);
z=(6371-reshape(Z,[mx,my,mz])/1e3);

X=reshape(X,[mx,my,mz]);
Y=reshape(Y,[mx,my,mz]);
% convert to Cartesian coord
[yc,xc,zc]=sph2cart(Y,pi/2-X,R);
xc=xc/1e3; yc=yc/1e3;; zc=zc/1e3; %in km

%load /home/yang/Proj/easthemi/misc/topogrid.mat
load /home/yang/FD3Dspher/FD3Dspher.090305.s2fr.PML/run/w.china/topogrid.mat
topo=qt;
Zin=topo;
%lattopo=55-0.2*(13-1)-[1:524]*0.2+0.2; % id==2 start at 13th grid
%lontopo=-30+0.2*(13-1)+[1:912]*0.2-0.2;
ti=66+12*0.12+[1:366]*0.12-0.12; tj=51-12*0.1+[1:344]*(-0.1)-(-0.1); % see config and SeisFD3D.conf
lattopo=tj; lontopo=ti;
flatmin=min(lattopo)*pi/180;flatmax=max(lattopo)*pi/180;
flonmin=min(lontopo)*pi/180;flonmax=max(lontopo)*pi/180;
phi_min=flatmin; phi_max=flatmax;
theta_min=flonmin; theta_max=flonmax;
Zin=double(Zin);theta_min=double(theta_min);theta_max=double(theta_max);
phi_min=double(phi_min);phi_max=double(phi_max);
cR=6371*1.e3;

% plate boundaries
!/bin/cp /home/yang/Proj/Shared/plate_boundaries/Bird_plate_boundaries.xy plate.dat
load plate.dat
plon=plate(:,1)*pi/180; plat=plate(:,2)*pi/180;
pr(1:length(plon))=6371e3; pr=pr';
nindplon=find(plon>flonmax | plon < flonmin);
plon(nindplon)=[]; plat(nindplon)=[]; pr(nindplon)=[];
nindplat=find(plat>flatmax | plat < flatmin);
plon(nindplat)=[]; plat(nindplat)=[];pr(nindplat)=[];
plon=plon*180/pi; plat=plat*180/pi;

% hotspots (in lat, lon)
%!/bin/cp /home/yang/Proj/Shared/hotspot_locations/hotspots.xy hotspot.dat
%!cat /home/yang/Proj/Shared/hotspot_locations/hotspot.lst | awk '{print $3, $2}' > hotspot.dat
%load hotspot.dat

% coastline
!cp /home/yang/Proj/Shared/coastline/coastline_1to5m.dat coast.dat
load coast.dat
coastlon=coast(:,1);coastlat=coast(:,2);
%indx=find(coastlon < 137 | coastlon > 146);
indx=find(coastlon < 66 | coastlon > 113);
coastlon(indx)=[];coastlat(indx)=[];clear indx;
indx=find(coastlat < 14 | coastlat > 51);
coastlon(indx)=[];coastlat(indx)=[];clear indx
clen=length(coastlat); cR(1:clen,1)=6381e3;
[cy,cx,cz]=sph2cart(coastlon*pi/180,coastlat*pi/180,cR);

%scl_xlim=[-250,250];
%scl_ylim=[-200,200];

%x2d=squeeze(x(:,:,end));
%y2d=squeeze(y(:,:,end));
%x1d=squeeze(x(:,1,end));
%y1d=squeeze(y(1,:,end));
%x1d=squeeze(x(1,:,end));
%y1d=squeeze(y(:,1,end));

%scl_xlim=[min(x1d)-1 max(x1d)+1];
%scl_ylim=[min(y1d)-1 max(y1d)+1];
%%%%%%%%%%%%%%%%%%%%%%%%%
%% faults
%%%%%%%%%%%%%%%%%%%%%%%%%
%% stations
%!cat /net/fs01/data/yang/n.cascadia/STinfo/sel_stns_cartcoord.lst | awk '{print $7,$8}' > st.lst
%st=load('st.lst');
%sxlon=[st(:,1)]/1000;
%sylat=[st(:,2)]/1000;
%!/bin/rm st.lst

%!cp /net/fs01/data/yang/n.cascadia/STinfo/sel_stns_cartcoord.lst tmp.dat
%[ntwk,stnm,slat,slon,sele,stny,stnx]= textread('tmp.dat','%s%s%f%f%f%f%f');
%nst=length(stny);
%stny=stny/1000;stnx=stnx/1000;
%!/bin/rm tmp.dat

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
%    W=permute(W,[2,1,3,4]);

    W=W*100; %percent

% z slice
for ncmp=1:2
    KCMPNM=cmp_list{ncmp};
    Vel=squeeze(W(:,:,:,ncmp));

    for nk=1:num_z
%    for nk=8
        k=zindx_list(nk);
        d=zval_list(nk);
        V=squeeze(Vel(:,:,k));
		% remove the mean, 04/26/2011, Y.Shen
		V=V-mean(mean(V));
		%
	XS=squeeze(x(:,:,k));
	YS=squeeze(y(:,:,k));
	zs=squeeze(zc(:,:,k));
        KDNM=num2str(zval_list(nk));
        KINM=num2str(zindx_list(nk));
        Vmax=max(max(abs(V)));

%        % shading interp
%        hid=figure;set(gcf,'renderer','zbuffer');
%        pcolor(x2d,y2d,V);
%        shading interp
%        colormap('jetwr');
%        xlim(scl_xlim); ylim(scl_ylim);
%        daspect([1,1,1]);
%        hold on
%        xlabel('x'); ylabel('y');
%        colorbar
%        title([' Depth ',num2str(D(1,1,k)),'km'])
%        text(300, 250, cmp_list(ncmp))
%     plot(fxlon,fylat, '-k')
%     hold on
%     plot(sxlon,sylat, '*k')

if(ncmp == 1)
fnm_fig=['figure_Vp_' int2str(nk) '.eps']
fnm_fig_png=['figure_Vp_' int2str(nk) '.png']
fnm_outfile=['Vp_' int2str(nk) '.txt']
Vmax=5.0;
end;
if(ncmp == 2)
fnm_fig=['figure_Vs_' int2str(nk) '.eps']
fnm_fig_png=['figure_Vs_' int2str(nk) '.png']
fnm_outfile=['Vs_' int2str(nk) '.txt']
Vmax=5.0;
end;

%        caxis([-Vmax,Vmax]);

%        Vmax=max(max(abs(V)));
        hid=figure;set(hid,'renderer','zbuffer');
%       imagesc(x1d,y1d,V);hold on
%	imagesc(y1d,x1d,V);hold on
	pcolor(YS,XS,V);hold on
%        surf(squeeze(permute(ys,[2 1 3])), ...
%            squeeze(permute(xs,[2 1 3])), ...
%            squeeze(permute(zs,[2 1 3])), ...
%            squeeze(permute(V,[2 1 3])));hold on
%sphere3d1c(Zin,theta_min,theta_max,phi_min,phi_max,cR,1,'contour','nearest',0);
%view(-10,5);
%axis off; axis tight;
hold on

        contour(lontopo,lattopo,topo,[10 10],'Color',[0.2 0.8 1.0]);hold on
        contour(lontopo,lattopo,topo,[2000 2000],'Color',[0.8 0.4 0.2]);hold on
	plot(plon,plat,'+','MarkerSize',1,'Color',[0.2 0.2 0.2]);hold on
%	plot(hotspot(:,2),hotspot(:,1),'ko','MarkerSize',5,'LineWidth',1); hold on
%	contour(x1d,y1d,D(:,:,end),[-1.0 -0.5 0],'Color',[0 1 0]);hold on
%	contour(x1d,y1d,D(:,:,end),[-1.5 -1.5],'Color',[1 0 0]); hold on
%	for ii=1:nst
%	plot(stnx(ii),stny(ii),'k^'); hold on
%	end
        shading interp
        colormap('jetwr');
       % xlim(scl_xlim); ylim(scl_ylim);
        daspect([1,1,1]);
        hold on
        xlabel('Longitude, deg.'); ylabel('Latitude, deg.');
        colorbar
%        shading interp
%        colormap('jetwr');
       % xlim(scl_xlim); ylim(scl_ylim);
%        daspect([1,1,1]);
%        hold on
%        xlabel('x'); ylabel('y');
%        colorbar
%        set(gca,'ydir','normal');
        title([cmp_list{ncmp} ' at depth = ',num2str(z(round(mx/2),round(my/2),k),4),' km'])

     caxis([-Vmax,Vmax]);
%      caxis([-Vmax,Vmax]*100);

%cR=6371*1.e3;
%sphere3d(Zin,theta_min,theta_max,phi_min,phi_max,cR,4.0*pi/180,'contour','spline',0);
%Zin=double(Zin);theta_min=double(theta_min);theta_max=double(theta_max);
%phi_min=double(phi_min);phi_max=double(phi_max);
%sphere3d1c(Zin,theta_min,theta_max,phi_min,phi_max,cR,1,'contour','nearest',0);
%view(-10,5);
%axis off; axis tight;

%     plot(fxlon,fylat, '-k')
%     hold on
%     plot(sxlon,sylat, '*k')

%     print('-depsc',fnm_fig);
%     print('-dpng',fnm_fig_png);

%minlon=-124.85; maxlon=-118.75; minlat=45.5;  maxlat=49.3;
% center of the xy coordinates in latitude and longitude
%lat0=(minlat+maxlat)/2; lon0=(minlon+maxlon)/2; 
%x0=0;y0=0;alpha=90;
%fid=fopen(char(fnm_outfile),'w+');
%kk=1;
%for ii=1:length(x1d)
%for jj=1:length(y1d)
%x=x1d(ii)*1000;
%y=y1d(jj)*1000;
%[xlat,ylon]=cart2geo(x,y,x0,y0,lat0,lon0,alpha);
%fprintf(fid,'%f %f %f\n',xlat,ylon,V(kk));
%kk=kk+1;
%end
%end
%fclose(fid);

    end
end
%close all

end
end
end

%plot([55 55],[23,43],'k')
%plot([50.2 50.2],[23,43],'k')
%plot([40, 65],[37.1,37.1],'k')
%plot([40, 65],[29.9,29.9],'k')


