
% generate model of velocity anomaly for resolution test

clear all; close all;

run ~/FWT/ANT/Proj/cascadia/model_updates/set_netcdf.m

% load Vs from the updated model
set_mfiles_path

ite_nm = ['OBS_0.025deg_11'];
fnm_conf=['../../../model_updates/SeisFD3D.conf_' ite_nm];

idx = 0; % 0, absolute velocity; 1, velocity perturbation

id=[];subs=[];subc=[];subt=[];indxem=[];indxkp=[];

%read updated model (same conf and coord)
dir_media=['../../../model_updates/updated_input_' ite_nm];
dir_coord=['../../../model_updates/updated_input_' ite_nm];

disp(['Read updated model... ']);

id{end+1} = 0; subs{end+1}=[1,1,1];subc{end+1}=[-1,-1,-1];subt{end+1}=[1,1,1];
               indxem{end+1}=[];
               indxkp{end+1}=[];
n=1;
[snapinfo{n}]=locate_snap(fnm_conf,id{n},'start',subs{n},'count',subc{n},'stride',subt{n});
[XSIM{n},YSIM{n},ZSIM{n}]=gather_coord(snapinfo{n},'coorddir',dir_coord);
% convert from radian to degrees
XSIM{n}=90-XSIM{n}*180/pi; %latitude
YSIM{n}=YSIM{n}*180/pi;
%define the area of plot (exclude pmls)
npml=12; %number of pml layers
minlat=XSIM{1}(end-npml,1,end);maxlat=XSIM{1}(1+npml,1,end);
minlon=YSIM{1}(1,1+npml,end);maxlon=YSIM{1}(1,end-npml,end);

mrh{n}=gather_media(snapinfo{n},'rho','mediadir',dir_media);
mmu{n}=gather_media(snapinfo{n},'mu','mediadir',dir_media);
mla{n}=gather_media(snapinfo{n},'lambda','mediadir',dir_media);
mvp{n}=((mla{n}+2*mmu{n})./mrh{n}).^0.5;
mvs{n}=(mmu{n}./mrh{n}).^0.5;
mvs{n}=smooth3(mvs{n},'box',[19 19 1]);

% define vertical grid in the inversion model
nzbgrid=[46 45 44 43 42 41 40 39 38 37];
% define the corresponding simulation grid
nzsimgrid=15+(nzbgrid-1)*3; %starting vertical grid + saved every other grid

% dimensions of the inversion blocks
mx=124; my=104; mz=48;

% index of first ckb block
nx1=1; ny1=1; nz1=1;
% index of last ckb block
nx2=mx; ny2=my; nz2=mz;

% Checkerboard with 5x5 inversion cells in lat and lon 
dmx=5; dmy=5; dmz=3;

vmodel = zeros(mx,my,mz);
for zz=1:length(nzbgrid)

v=squeeze(mvs{1}(npml:dmx:npml+mx*dmx-1,npml:dmy:npml+my*dmy-1,nzsimgrid(zz)));
v=double(v)/1000;

% vmeans=mean(mean(v))
% dv = (v-vmeans)./vmeans; % vs perturbation

vmodel(:,:,nzbgrid(zz))=v;
end

m=vmodel(:);

fid = fopen('block_ckb_model.dat','w');
for ii=1:mx*my*mz
    fprintf(fid,'%3.3f \n',m(ii));
end
fclose(fid);
