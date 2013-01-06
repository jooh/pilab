% return mm size of voxels for a given SPM volume handle.
% mm = vox2mm(V)
function mm = vox2mm(V)

mm = sqrt(sum(V(1).mat(1:3,1:3).^2));
