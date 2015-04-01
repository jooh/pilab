% convolve a row vector of onset times (in volume units) to generate a
% design matrix.
%
% INPUTS:
% onsets: vector with onset times in units TR. Note that for consistency
%   with SPM we use 0-based indexing, so the start of the scan is time 0,
%   not time 1.
% conind: vector with regressor indices for each onset time
% tr: acquisition time in s
% n: run length in volumes
%
% NAMED, OPTIONAL VARARGIN:
% hrffun: (spm_hrf) function handle that gets fevaled with a tr as input
% nbin: (32) number of bins to subdivide TR into (limits precision of onset
%   times). If all onsets are on the volume trigger (ie integers) this can
%   safely be set to 1 for speed.
% p: ([]) matrix defining parametric modulators on each onset (default means
%   no modulation)
% dur: (0) duration of each response (default means a single bin duration)
%
% [X,onserr] = convolveonsets(onsets,conind,tr,n,varargin)
function [X,onserr] = convolveonsets(onsets,conind,tr,n,varargin)

ntrials = numel(onsets);
getArgs(varargin,{'hrffun','spm_hrf','nbin',32,'p',[],...
    'dur',zeros(ntrials,1)});

% check parametric modulator
assert(isempty(p) || size(p,1)==ntrials,'p does not match ntrials');
npara = size(p,2);

% overall logic: first super-sample to achieve sufficient precision to
% capture onsets with reasonable precision. Then downsample to correct
% length.

% set up super-sampled design matrix
% (we add an extra bin to accommodate very late onsets)
nsuper = n*nbin+1;
ureg = unique(conind);
nreg = length(ureg);
superx = zeros(nsuper,nreg+npara,class(onsets));

% check duration
if isscalar(dur)
    dur = repmat(dur,[ntrials 1]);
end
superdur = floor(dur*nbin);
assert(numel(superdur)==ntrials,'dur does not match ntrials');


% generate super-sampled HRF
supertr = tr/nbin;
hrf = feval(hrffun,supertr);
hrf = hrf / max(hrf);

superons = round(onsets * nbin);
% track errors in onset times introduced by rounding to nearest bin
onserr = onsets-(superons/nbin);

assert(numel(onsets)==numel(conind),'onsets and conind must be equal size');

% plug in each trial's onset in the right column
for t = 1:ntrials
    % add an extra bin here to convert 0-based to 1-based indexing
    ind = 1+[superons(t):superons(t)+superdur(t)];
    superx(ind,conind(t)) = 1;
    if npara
        superx(ind,nreg+1:end) = p(...
            repmat(t,[superdur(t)+1,1]),:);
    end
end

% convolve with HRF
convx = conv2(superx,hrf);

% downsample to TR
X = downsample(convx(1:nsuper-1,:),nbin);
