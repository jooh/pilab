% A convenience wrapper around spm_write_vol - insures that data is written
% in correct precision to prevent quantizing.
% datamat2nifti(data,outpath,V)
function datamat2nifti(data,outpath,V)

assert(isequal(size(data),V.dim),'data does not match mask dims');

V.dt = [matlabclass2spmnifti(data) spm_platform('bigend')];
V.fname = outpath;
spm_write_vol(V,data);
