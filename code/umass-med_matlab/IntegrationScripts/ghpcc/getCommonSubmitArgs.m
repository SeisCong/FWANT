function commonSubmitArgs = getCommonSubmitArgs(cluster, numWorkers)
% Get any additional submit arguments for the LSF bsub command
% that are common to both independent and communicating jobs.

% Copyright 2016-2017 The MathWorks, Inc.

% http://wiki.umassrc.org/wiki/index.php/GHPCC

commonSubmitArgs = '';

% Number of nodes/cores
ppn = validatedPropValue(cluster, 'ProcsPerNode', 'double');
if ppn>0
    % Don't request more cores/node than workers
    ppn = min(numWorkers,ppn);
    commonSubmitArgs = sprintf('%s -R "span[ptile=%d]"',commonSubmitArgs,ppn);
end


%% REQUIRED

% Physical Memory used by a single core
mu = validatedPropValue(cluster, 'MemUsage', 'char');
if isempty(mu)
    emsg = sprintf(['\n\t\t>> %% Must set MemUsage.  E.g. 2GB per core\n\n', ...
                       '\t\t>> c = parcluster;\n', ...
                       '\t\t>> c.AdditionalProperties.MemUsage = ''2048'';\n', ...
                       '\t\t>> c.saveProfile\n\n']);
    error(emsg) %#ok<SPERR>
else
    % 1 => 1 MB
    commonSubmitArgs = sprintf('%s -R "rusage[mem=%s]"',commonSubmitArgs,mu);
end

% Walltime
wt = validatedPropValue(cluster, 'WallTime', 'char');
if isempty(wt)
    emsg = sprintf(['\n\t\t>> %% Must set WallTime. E.g. 1 hour\n\n', ...
                       '\t\t>> c = parcluster;\n', ...
                       '\t\t>> c.AdditionalProperties.WallTime = ''01:00'';\n', ...
                       '\t\t>> c.saveProfile\n\n']);
    error(emsg) %#ok<SPERR>
else
    % -W [hour:]minute
    commonSubmitArgs = [commonSubmitArgs ' -W ' wt];
end


%% OPTIONAL

% Project name
pn = validatedPropValue(cluster, 'ProjectName', 'char');
if ~isempty(pn)
    commonSubmitArgs = [commonSubmitArgs ' -P ' pn];
end

% Queue name
ngpus = validatedPropValue(cluster, 'GpusPerNode', 'double');
if ngpus>0
    if ngpus==1
        % If ngpus = 1, need GPU card
        gcard = validatedPropValue(cluster, 'GpuCard', 'string');
        VALID_CARDS = {'tesla','geforce'};
        if ~any(strcmpi(gcard,VALID_CARDS))
            fmt = repmat('%s, ',1,numel(VALID_CARDS));
            fmt(end-1:end) = [];
            error(['When requesting 1 GPU/node, must provide a valid GPU card: ' fmt ' \n\n\t>> c.AdditionalProperties.GpuCard = ''a-gpu-card'';'],VALID_CARDS{:})
        end
        switch lower(gcard)
          case 'tesla'
            qn = 'gpu';
          case 'geforce'
            qn = 'gpu_geforceGTX';
        end
    else
        % Only support GPUs/node > 1 on "gpu" queue
        qn = 'gpu';
    end
    commonSubmitArgs = sprintf('%s -R "rusage[ngpus_excl_p=%d]"', commonSubmitArgs, numWorkers);
else
    qn = validatedPropValue(cluster, 'QueueName', 'char');
end
if ~isempty(qn)
    commonSubmitArgs = [commonSubmitArgs ' -q ' qn];
end

% Email notification
ea = validatedPropValue(cluster, 'EmailAddress', 'char');
if ~isempty(ea)
    commonSubmitArgs = [commonSubmitArgs ' -B -N -u ' ea];
end

% Need to request a certain number of MDCS licenses.
%
% See Tech Note to setup MDCS with elim
%
% Every job is going to require a certain number of MDCS licenses.
% Specification of MDCS licenses which must be allocated to this
% job.  See Tech Note to setup MDCS with elim
%
%    http://www.mathworks.com/matlabcentral/answers/94481
%
commonSubmitArgs = sprintf('%s -R "rusage[mdcs=%d:duration=1]"', commonSubmitArgs, numWorkers);

% Catch-all
asa = validatedPropValue(cluster, 'AdditionalSubmitArgs', 'char');
if ~isempty(asa)
    commonSubmitArgs = [commonSubmitArgs ' ' asa];
end

commonSubmitArgs = strtrim(commonSubmitArgs);
