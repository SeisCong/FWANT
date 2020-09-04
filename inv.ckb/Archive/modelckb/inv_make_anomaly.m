function vout=inv_make_anomaly(outfile,nx,ny,nz,varargin)
%Generate velocity anomalies for resolution tests.
%USAGE: vout=inv_make_anomaly(outfile,nx,ny,nz,anomaly1,anomaly2,...,anomalyN)
%Arguments in:
%   outfile: file name storing the computed anomalies;
%   nx,ny,nz: number of blocks in the model domain. nx: number of N-S grid;
%       ny: number of E-W grid;
%   anomaly1,anomaly2,...,anomalyN: a series of vector with 7 elements specifying
%       anomaly bodies. Format: [xstart xend ystart yend zstart zend dv],where dv 
%       is the relative perturbation (0~1) of the anomaly body, other values are the
%       start and end block number. Here x,y,z MUST be the same direction
%       (meaning) as nx,ny,nz.
%
%Arguments out:
%   vout: column vector with length of nx*ny*nz. This is the same as in the
%       outfile.

nanomaly=length(varargin);
vout=zeros(nx*ny*nz,1);
for n=1:nanomaly
    disp(['anomaly: ' num2str(n)]);
    clear anom;
    anom=varargin{n};
    for k=anom(5):anom(6)
        for j=anom(3):anom(4)
            for i=anom(1):anom(2)
                idx0=nx*(j-1) + i + nx*ny*(k-1);
                vout(idx0)=anom(7);
            end
        end
    end
end
fidout=fopen(outfile,'w');
fprintf(fidout,'%f\n',vout);
fclose(fidout);
end