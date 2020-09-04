
% generate std list for dataerror_list

[nstd]=textread('../../measure/tmp_ANT.dat','%*s %*s %*f %*f %*f %*f %*f %f %*s');

fid = fopen('dataerror_list','w');
for ii=1:length(nstd)
    fprintf(fid,'%f \n',nstd(ii));
end
fclose(fid);
