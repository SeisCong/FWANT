
% generate model of velocity anomaly for resolution test

clear all; close all;

% load Vs from the updated model
set_mfiles_path

% dimensions of the inversion blocks
mx=51; my=59; mz=20;

% maximum perturbation
dvmax=5

% Checkerboard with 5x5x3 inversion cells in lat and lon 
dmx=4; dmy=4; dmz=1;

vmodel = zeros(mx,my,mz);

intX = floor(mx/dmx);
intY = floor(my/dmy);
intZ = floor(mz/dmz);

% generate a sine function for the velocity perturbation along Z
wavelength=10*dmz;
sinvz1 = dvmax*sin(2*pi*(0:intZ-1)./wavelength);
sinvz2 = dvmax*sin(2*pi*(0:intZ-1)./wavelength+pi/2);

for ii = 1:2:intX
for jj = 1:intY
xx =    
if mod(jj,2)
vmodel(ii:(ii+dmx-1),jj:(jj+dmy-1),:)=sinvz1;
else
vmodel(ii:(ii+dmx-1),jj:(jj+dmy-1),:)=sinvz2;
end

end
end

m=vmodel(:);

fid = fopen(['block_ckb_' num2str(dmx) 'x' num2str(dmy) '.dat'],'w');
for ii=1:mx*my*mz
    fprintf(fid,'%3.3f \n',m(ii));
end
fclose(fid);
