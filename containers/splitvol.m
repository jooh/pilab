% obscure helper function for splitting a set of volume instances into cell
% arrays as defined by the first input (split). You can enter any number of
% volumes and you will get the same number of cell arrays back.
%
% This is used in e.g. vol2glm_batch and roiglm2rdmvol.
%
% varargout = splitvol(split,varargin)
function varargout = splitvol(split,varargin)

assert(isa(varargin{1},'BaseVolume'),...
    'all inputs must be BaseVolume or subclasses thereof');
nchunk = varargin{1}.desc.samples.nunique.chunks;
% configure split
if isempty(split)
    % just one global split
    split = ones(1,nchunk);
elseif ischar(split)
    split = eval(split);
end
usplit = unique(split);
% support skipping some chunks completely
usplit(isnan(usplit)) = [];
nsplit = length(usplit);

% process each volume input
for n = 1:nargin-1
    thisvol = varargin{n};
    assert(isa(thisvol,'BaseVolume'),...
        'all inputs must be BaseVolume or subclasses thereof');
    assert(thisvol.desc.samples.nunique.chunks==nchunk,...
        'mismatched chunks across vols');
    % apply the split
    for s = 1:nsplit
        splitind = split==usplit(s);
        varargout{n}{s} = thisvol.selectbymeta('chunks',...
            thisvol.desc.samples.unique.chunks(splitind));
    end
end
