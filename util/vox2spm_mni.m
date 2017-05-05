% Find the MNI coordinate for a given native-space voxel index. Write out a
% nifti that stores the transformation for future use (this gets stored in the
% same dir as seg_sn).
%
% Implements SPM list email from John Ashburner:
% https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=spm;e71a6654.1204
%
% INPUTS
% vox   xyz voxel indices (1-based, so add 1 if obtaining from fslview...)
% nativenifti   path to a native-space nifti file for subject
% seg_sn        path to _seg_sn.mat file for subject
%
% 2017 05 05 J Carlin
%
% mni = vox2spm_mni(vox,nativenifti,seg_sn)
function mni = vox2spm_mni(vox,nativenifti,seg_sn)

outdir = fileparts(seg_sn);
% NB, NOT the specified ofname because SPM is insane
outfn = 'y_vox2spm_mni.nii';
outpath = fullfile(outdir,outfn);

if exist(outpath,'file')
    logstr('using existing mapping from %s...\n',outpath);
else
    %-----------------------------------------------------------------------
    % Job configuration created by cfg_util (rev $Rev: 4252 $)
    %-----------------------------------------------------------------------
    matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.sn2def.matname = {seg_sn};
    matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.sn2def.vox = [NaN NaN NaN];
    matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.sn2def.bb = [NaN NaN NaN
                                                                  NaN NaN NaN];
    matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {nativenifti};
    matlabbatch{1}.spm.util.defs.ofname = 'vox2spm_mni';
    matlabbatch{1}.spm.util.defs.fnames = '';
    matlabbatch{1}.spm.util.defs.savedir.saveusr = {outdir};
    matlabbatch{1}.spm.util.defs.interp = 0;

    spm_jobman('run',matlabbatch);
end

nii = nifti(outpath);
mni = reshape(nii.dat(vox(1),vox(2),vox(3),1,:),[3,1]);
