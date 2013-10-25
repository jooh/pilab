% Convert glmdenoise results to SPM format by writing out the cleaned epis
% to 3D nifti and updating the file references in the entered SPM.mat
% SPM = glmdenoise2spm(SPM,epis,[models],[mask],[rootdir])
function SPM = glmdenoise2spm(SPM,epis,models,mask,rootdir)

if ~exist('models','var') || isempty(models)
    domodels = false;
else
    domodels = true;
end

if ~exist('mask','var') || isempty(mask)
    assert(ndims(epis)==4,'must enter mask unless epis are 4D')
    es = size(epis{1});
    mask = true(es(1:3));
    indexepi = @(ru,vol) epis{ru}(:,:,:,vol);
else
    indexepi = @(ru,vol) epis{ru}(:,vol);
end


if ~exist('rootdir','var') || isempty(rootdir)
    rootdir = pwd;
end

mask = mask>0;
nruns = length(epis);

assert(nruns == length(SPM.Sess),'mismatched nruns')
c = 0;
for r = 1:nruns
    % optionally update the design matrix (useful if you do not
    % assumeconvolved)
    if domodels
        SPM.xX.X(SPM.Sess(r).row,SPM.Sess(r).col) = models{r};
    end
    fprintf('writing volumes for run %d of %d...\n',r,nruns);
    for v = 1:size(epis{r},2)
        c = c+1;
        d = zeros(size(mask));
        % put the voxels in their 3D location
        d(mask) = indexepi(r,v);
        [path,fn,ext ] = fileparts(SPM.xY.VY(c).fname);
        [rootdir,dirname] = fileparts(path);
        fileoutdir = fullfile(rootdir,[dirname '_denoised']);
        mkdirifneeded(fileoutdir);
        SPM.xY.VY(c).fname = fullfile(fileoutdir,[fn ext]);
        spm_write_vol(SPM.xY.VY(c),d);
    end
end
% make sure SPM.xY.P matches {SPM.xY.VY.fname}
SPM.xY.P = char({SPM.xY.VY.fname});
