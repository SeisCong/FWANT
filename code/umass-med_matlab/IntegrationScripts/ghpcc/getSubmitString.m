function submitString = getSubmitString(jobName, quotedLogFile, quotedCommand, ...
    additionalSubmitArgs)
%GETSUBMITSTRING Gets the correct bsub command for an LSF cluster

% Copyright 2010-2011 The MathWorks, Inc.

% Submit to LSF using bsub.  Note the following:
% "-J " - specifies the job name
% "-o" - specifies where standard output goes to (and standard error, when -e is not specified)
% Note that extra spaces in the bsub command are permitted
submitString = sprintf('bsub -J %s -o %s %s %s', jobName, quotedLogFile, additionalSubmitArgs, quotedCommand);
