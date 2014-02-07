% convert vector to matrix form and pass on to datamat2nift.
% NB, mask can be a logical matrix or linear indices
%
% datavec2nifti(data,mask,outpath,V)
function datavec2nifti(data,mask,outpath,V)

if isa(data,'logical')
    false(V.dim);
else
    datamat = zeros(V.dim,class(data));
end
datamat(mask) = data;

datamat2nifti(datamat,outpath,V);
