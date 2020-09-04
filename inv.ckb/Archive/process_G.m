clear all; close all;

[Kap, Kbp, dt, weight, stno,evno,xlocn,xlocv,ylocn,ylocv,zlocn,zlocv,RL]=textread('inv_Gd_list_observed.ori','%s %s %f %f %d %d %d %d %d %d %d %d %s');

% standard deviation of travel times in second
ntstd=textread('dataerror_list.ori','%f');

minlon = 236; maxlon = 246; minlat = 40; maxlat = 48;

[nw, sta, lon, lat] = textread('../../../STinfo/cascadia.stns.Z.0.05deg.lst','%s %s %f %f');
idx = find(lon<minlon | lon>maxlon | lat<minlat | lat>maxlat);
nwtmp = nw(idx);
statmp = sta(idx);
for cc=1:length(idx)
stalist{cc}=[char(nwtmp(cc)) '.' char(statmp(cc))];
end

fid1=fopen('inv_Gd_list_observed','w');
fid2=fopen('dataerror_list','w');
nt = length(Kap);

for ii=1:nt
tt=strfind(Kap{ii},'_');
sta1=Kap{ii}(12:tt(2)-1)
sta2=Kap{ii}(tt(2)+1:tt(3)-1)

pp1=find(strcmp(stalist, sta1))
pp2=find(strcmp(stalist, sta2))
if ~isempty(pp1), continue, end
if ~isempty(pp2), continue, end

fprintf(fid1,'%s %s %f %f %d %d %d %d %d %d %d %d %s\n',Kap{ii}, Kbp{ii}, dt(ii), weight(ii), stno(ii),evno(ii),xlocn(ii),xlocv(ii),ylocn(ii),ylocv(ii),zlocn(ii),zlocv(ii),RL{ii});
fprintf(fid2,'%f\n',ntstd(ii));
end
fclose(fid1);
fclose(fid2);


