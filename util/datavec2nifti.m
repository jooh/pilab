% convert vector to matrix form and pass on to datamat2nifti.
%
% data: 1D row or column vector 
% mask: 3D logical array form apping data vector to 3D
% outpath: char array where nifti gets written
% V: volume header (see spm_vol)
%
% datavec2nifti(data,mask,outpath,V)
function datavec2nifti(data,mask,outpath,V)

datamat = datavec2mat(data,mask);
datamat2nifti(datamat,outpath,V);
