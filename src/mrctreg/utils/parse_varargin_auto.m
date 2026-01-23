function parse_varargin_auto(vargs,defaults)

% parse_varargin_auto.m
% 
% Purpose: parse a varargin cell array and apply defaults where values are
%           not specified by varargin
%
% Usage: process_varargin(vargs,defaults)
%
% Created: 03-09-22 by Tyler Seibert
% See Change Log at end of file for modifications
% Based on parse_varargin.m but automatically creates variable output
%   names. Instead of varargout, use assignin.m to set the variables in
%   the calling function. 
% 
% Required Arguments:
%   vargs           =   [cell 1xn] varargin cell array from function calling
%                       process_varargin.m
%   defaults    =   [cell mx2] cell array with variable/value pairs
% 

% check defaults and get varnames:
if size(defaults,2) ~= 2
    error('%s: ''defaults'' should be an nx2 cell array with rows corresponding to variable/default pairs\n',mfilename);
else
    varnames = {defaults{:,1}}; % variable names
    defaultvals = {defaults{:,2}}; % default values
end

% check for obvious errors:
if numel(vargs) == 0
    % nothing. This means no optional arguments were provided
elseif rem(numel(vargs),2)
    error('%s: optional arguments should come in pairs (key,value)\n',mfilename);
elseif numel(vargs)/2 > numel(varnames)
    fprintf('%s: more optional arguments than possibilities in ''varnames''\n',mfilename);
    error('%s: %d pairs in varargin, but only %d varnames given\n',mfilename,numel(vargs)/2,numel(varnames));
% elseif numel(varnames) ~= nargout
%     fprintf('%s: number of outputs must match number of variable names in ''varnames''\n',mfilename);
%     error('%s: %d elements in nargout, but %d varnames\n',mfilename,nargout,numel(varnames));
end

% assign defaults for all (in workspace for function calling parse_varargin_auto:
for k=1:numel(varnames)
    assignin('caller',varnames{k},defaultvals{k}); % assign variable value in function calling parse_varargin_auto
end

% now overwrite all key/value pairs as variables:
for i=1:2:numel(vargs)
    if ~ismember(vargs{i},varnames)
        error('%s: varargin key ''%s'' not found among varnames of ''defaults''\n',mfilename,vargs{i});
    else
        % eval([vargs{i} '= vargs{i+1};']);
        assignin('caller',vargs{i},vargs{i+1}); % assign variable value in function calling parse_varargin_auto
    end
end


% CHANGE LOG:
% 03-09-22 --> Created by Tyler Seibert
% 
