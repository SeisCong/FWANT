function id = schedID(job)

% Copyright 2014-2015 The MathWorks, Inc.
    
narginchk(1,1)
if numel(job)>1
    error('Must only supply one job.')
end

if ~isa(job,'parallel.job.CJSIndependentJob') ...
        && ~isa(job,'parallel.job.CJSCommunicatingJob')
    error('Must provide Independent or Communicating Job')
end

jcd = job.Parent.getJobClusterData(job);
if isempty(jcd)
    error('Job has not been submitted.')
end

id = jcd.ClusterJobIDs;
if length(id)==1
    id = id{:};
end

end
