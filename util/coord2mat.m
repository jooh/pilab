% convert a set of coordinates (dimensions in rows, features in columns) to
% a logical matrix.
%
% mat = coord2mat(coords,outsize)
function mat = coord2mat(coords,outsize)

assert(ndims(coords)==2,'coords must be 2d dimensions by features');
coordcell = mat2cell(coords,ones(size(coords,1),1),size(coords,2));
% We could check that the mat has the correct dimensionality (corresponding
% to number of rows in coords), but this would disable some of matlab's
% 'flexibility' in terms of squeezing off singleton end dimensions (the
% final row is often all 1 to enable coordinate transforms).
mat = false(outsize);
mat(sub2ind(outsize,coordcell{:})) = true;
