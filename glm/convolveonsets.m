% convolve a row vector of onset times (in volume units) to generate a
% design matrix.
%
% INPUTS:
% onsets: vector with onset times in units TR
% conind: vector with regressor indices for each onset time
% tr: acquisition time in s
% n: run length in volumes
%
% NAMED, OPTIONAL VARARGIN:
% hrffun: (spm_hrf) function handle that gets fevaled with a tr as input
% nbin: (32) number of bins to subdivide TR into (limits precision of onset
% times)
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
nsuper = n*nbin;
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

% scale up onsets to bin time (and scale up by 1 to since we are using
% 1-indexing in Matlab (while time is 0-indexed, sort of)
superons = round(onsets * nbin)+1;
% track errors in onset times introduced by rounding to nearest bin
onserr = onsets-(superons/nbin-1);

% plug in each trial's onset in the right column
for t = 1:ntrials
    superx(superons(t):superons(t)+superdur(t),conind(t)) = 1;
    if npara
        superx(superons(t):superons(t)+superdur(t),nreg+1:end) = p(...
            repmat(t,[superdur(t)+1,1]),:);
    end
end

% convolve with HRF
convx = conv2(superx,hrf);

% downsample to TR (and chop off over-shoot for free)
inds = floor(linspace(1,nsuper,n));
X = convx(inds,:);

% 
