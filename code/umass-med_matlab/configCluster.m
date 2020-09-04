function configCluster
% Configure MATLAB to submit to the cluster.

% Copyright 2013-2017 The MathWorks, Inc.

% The version of MATLAB being supported
release = ['R' version('-release')];

% Import cluster definitions
def = clusterDefinition();

% Cluster list
cluster_dir = fullfile(fileparts(mfilename('fullpath')),'IntegrationScripts');
% Listing of setting file(s).  Derive the specific one to use.
cluster_list = dir(cluster_dir);
% Ignore . and .. directories
cluster_list = cluster_list(arrayfun(@(x) x.name(1), cluster_list) ~='.');
len = length(cluster_list);
if len==0
    error('No cluster directory exists.')
elseif len==1
    cluster = cluster_list.name;
else
    cluster = lExtractPfile(cluster_list);
end

% Determine the name of the cluster profile
% If using R2015b or older, spaces in the profile name is not supported
if verLessThan('matlab','8.6')
    profile = [cluster '_' def.Type '_' release];
else
    profile = [cluster ' ' def.Type ' ' release];
end        

% Delete the old profile (if it exists)
profiles = parallel.clusterProfiles();
idx = strcmp(profiles,profile);
ps = parallel.Settings;
ws = warning;
warning off
ps.Profiles(idx).delete
warning(ws)

% User's local machine's hostname
if strcmp(def.Type, 'local')
    hostname = '';
else
    [~, hostname] = system('hostname');
    hostname = strtrim(hostname);
end

% Skip this for local
if ~strcmp(def.Type,'local')
    % If multiple releases were specified in the mdcs.rc
    % select the correct one to use.
    releaseBreakDown = strsplit(def.ClusterMatlabRoot,',');
    matchingRelease = (~cellfun(@isempty,regexp(releaseBreakDown,release,'once')));
    if ~matchingRelease
        emsg = sprintf(['\n\t\t The version of MATLAB you are running is not installed on the cluster.\n', ...
            '\t\t Contact your cluster administrator for further assistance. \n']);
        error(emsg) %#ok<SPERR>
    end
    releaseToUse = releaseBreakDown{matchingRelease};
    releaseToUse = strsplit(releaseToUse,':');
    def.ClusterMatlabRoot = releaseToUse{2};
end

% Create the user's local Job Storage Location folder
if strcmp(def.Type, 'remote')
    rootd = lGetLocalRoot();
elseif strcmp(def.Type, 'local')
    if isempty(def.LocalJobStorageLocation)
        rootd = lGetLocalRoot();
    else
        user = lGetRootUsername();
        rootd = [def.LocalJobStorageLocation user];
    end
else
    user = lGetRootUsername();
    rootd = [def.RemoteSubmissionClientDir user];
end

jfolder = fullfile(rootd,'MdcsDataLocation',cluster,hostname,release);

if exist(jfolder,'dir')==false
    [status,err,eid] = mkdir(jfolder);
    if status==false
        error(eid,err)
    end
end

% Configure the user's remote storage location and assemble
% the cluster profile.
if strcmp(def.Type, 'remote')
    user = lGetRemoteUserName(cluster);
    % rootd = [def.RemoteJobStorageLocation user];
    rootd = lGetRemoteJobStorageLocationRoot(cluster,user);
    rjsl = [rootd '/MdcsDataLocation/' cluster '/' hostname '/' release '/' def.Type];
elseif strcmp(def.Type, 'remotesubmission')
    if ispc
        rootd = [def.RemoteJobStorageLocation user];
        rjsl = [rootd '/MdcsDataLocation/' cluster '/' hostname '/' release '/' def.Type];
    else
        rjsl = '';
    end
else
    rjsl = '';
    user = '';
    def.ClusterHost = '';
    def.ClusterMatlabRoot = '';
end
assembleClusterProfile(jfolder, rjsl, cluster, user, profile, def);

lNotifyUserOfCluster(upper(cluster))

% % Validate if you want to
% profiles = parallel.clusterProfiles();
% idx = strcmp(profiles,profile);
% ps = parallel.Settings;
% ps.Profiles(idx).validate

end


function cluster_name = lExtractPfile(cl)
% Display profile listing to user to select from
len = length(cl);
for pidx = 1:len
    name = cl(pidx).name;
    names{pidx,1} = name; %#ok<AGROW>
end

selected = false;
while selected==false
    for pidx = 1:len
        fprintf('\t[%d] %s\n',pidx,names{pidx});
    end
    idx = input(sprintf('Select a cluster [1-%d]: ',len));
    selected = idx>=1 && idx<=len;
end
cluster_name = cl(idx).name;

end


function r = lGetLocalRoot()

if isunix
    % Some Mac user's have noticed that the [/private]/tmp
    % directory gets cleared when the system is reboot, so for UNIX
    % in general, let's just use the user's local home directory.
    uh = java.lang.System.getProperty('user.home');
    
    r = char(uh);
else
    % If this returns an empty string (some how user name is not defined),
    % it's a no-op for FULLFILE, so there's no strong need to error out.
    un = java.lang.System.getProperty('user.name');
    un = char(un);
    
    r = fullfile(tempdir,un);
end

end


function un = lGetRemoteUserName(cluster)
un = input(['Username on ' upper(cluster) ' (e.g. joe): '],'s');
if isempty(un)
    error(['Failed to configure cluster: ' cluster])
end

end


function rootd = lGetRemoteJobStorageLocationRoot(cluster,user)
default_rootd = ['/home/' user];
rootd = input(['JobStorageLocation on ' upper(cluster) ' [' default_rootd ']: '],'s');
if isempty(rootd)
    rootd = default_rootd;
end

end


function user = lGetRootUsername()

user = char(java.lang.System.getProperty('user.name'));

end


function assembleClusterProfile(jfolder, rjsl, cluster, user, profile, def)

% Create generic cluster profile
c = parallel.cluster.Generic;

% Required mutual fields
% Location of the Integration Scripts
c.IntegrationScriptsLocation = fullfile(fileparts(mfilename('fullpath')),'IntegrationScripts', cluster);
c.NumWorkers = str2num(def.NumWorkers); %#ok<ST2NM>
c.OperatingSystem = 'unix';

% Depending on the submission type, populate cluster profile fields
if strcmp(def.Type, 'local')
    c.HasSharedFilesystem = true;
else
    % Set common properties for remote and remoteSubmission
    c.AdditionalProperties.UserNameOnCluster = user;
    c.AdditionalProperties.ClusterHost = def.ClusterHost;
    c.ClusterMatlabRoot = def.ClusterMatlabRoot;
    if strcmp(def.Type, 'remote')
        c.AdditionalProperties.RemoteJobStorageLocation = rjsl;
        c.HasSharedFilesystem = false;
    else
        if ispc
            jfolder = struct('windows',jfolder,'unix',rjsl);
        end
        c.HasSharedFilesystem = true;
    end
end
c.JobStorageLocation = jfolder;

% AdditionalProperties for the cluster:
% username, queue, walltime, e-mail, etc.
c.AdditionalProperties.AdditionalSubmitArgs = '';
% Set the debug log to false by default. Enable if troubleshooting.
c.AdditionalProperties.DebugMessagesTurnedOn = false;
c.AdditionalProperties.EmailAddress = '';
c.AdditionalProperties.GpuCard = '';
c.AdditionalProperties.GpusPerNode = 0;
c.AdditionalProperties.IdentityFile = '';
c.AdditionalProperties.MemUsage = '';
c.AdditionalProperties.PrivateKeyFile = '';
c.AdditionalProperties.PrivateKeyFileHasPassPhrase = false;
c.AdditionalProperties.ProcsPerNode = 0;
c.AdditionalProperties.ProjectName = '';
c.AdditionalProperties.QueueName = '';
c.AdditionalProperties.UseIdentityFile = true;
c.AdditionalProperties.WallTime = '';

% Save Profile
c.saveAsProfile(profile);
if verLessThan('matlab','9.1')
    c.saveProfile()
else
    c.saveProfile('Description', profile)
end

% Set as default profile
parallel.defaultClusterProfile(profile);

end


function lNotifyUserOfCluster(cluster)

fprintf(['\n\t\t>> %% Must set MemUsage and WallTime before submitting jobs to %s.  E.g.\n\n', ...
    '\t\t>> c = parcluster;\n', ...
    '\t\t>> %% Request 2GB per core\n', ...
    '\t\t>> c.AdditionalProperties.MemUsage = ''2048'';\n', ...
    '\t\t>> %% 1 hour\n', ...
    '\t\t>> c.AdditionalProperties.WallTime = ''01:00'';\n', ...
    '\t\t>> c.saveProfile\n\n'], cluster)

end
