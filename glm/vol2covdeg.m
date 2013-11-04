function covariatedeg = vol2covdeg(epivol)
nchunks = epivol.desc.samples.nunique.chunks;
covdegs = NaN([1 nchunks]);
for c = 1:nchunks
    nvol = sum(epivol.meta.samples.chunks == ...
        epivol.desc.samples.unique.chunks(c));
    % Kay's rule for deciding on n covariates based on run
    % duration.
    covdegs(c) = round(nvol * epivol.frameperiod / 60 / 2);
end
assert(~any(isnan(covdegs)),...
    'failed to identify covariatedeg');
% if we found a single deg that works for all runs, life is
% easy
if all(covdegs(1)==covdegs)
    covariatedeg = covdegs(1);
    %fprintf('adaptively selected covariate degree: %d\n',...
        %covariatedeg);
else
    % otherwise we need to add new functionality to CovGLM to
    % support different covariatedeg for different runs.
    error(['Adaptive covariate deg selection failed. ' ...
        'Different covariates selected per run: ' ...
        mat2str(covdegs)]);
end
