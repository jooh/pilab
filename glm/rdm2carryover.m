% onsets: vector with onset times in units TR
% conind: vector with regressor indices for each onset time
% rdm: something compatible with asrdmmat. Conditions should be in same
%   order as conind
% nstep: (default 1) carryover degree
%
% NB, we simply plug in the values from the RDM and mean correct. If you
% want to encode specific predictions (e.g. linear/log) you must transform
% the RDM before calling this.
%
% p = rdm2carryover(onsets,conind,rdm,nstep)
function p = rdm2carryover(onsets,conind,rdm,nstep)

assert(issorted(onsets),'onsets must be in ascending time order');

if ieNotDefined('nstep')
    nstep = 1;
end

rdm = asrdmmat(rdm);

ntrials = numel(onsets);
p = zeros([ntrials nstep],class(onsets));
stepsize = 1:nstep;

% nstep rather than 1 here to avoid indexing <1 in rdm
for t = (nstep+1):ntrials
    p(t,:) = rdm(conind(t),conind(t-stepsize));
end

% de-mean (we usually want to use this as a modulator on a primary
% response)
p = bsxfun(@minus,p,mean(p,1));

% Aguirre does something like this. This won't work unless the carry-over
% takes multiple levels, however. I guess the fundamental question is: do
% we want to model 0? If so, we should mean-correct the entire carry-over
% predictor. This is usually how we think of things with RSA, where 0 is
% the more similar and maximally adapting case. If 0 instead reflects 'no
% modulation', we may prefer something like the below, which would insure
% that these 0 modulators remain 0 even after mean correction.
%zeroind = p==0;
% mean center
%for r = 1:nstep
    %valid = p(p(:,r)~=0,r);
    %p(p(:,r)~=0,r) = valid - mean(valid);
%end
