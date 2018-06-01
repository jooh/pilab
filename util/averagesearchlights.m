% combine dissimilarity estimates for a feature by averaging RDMs for
% all searchlights that include it.
%
% inputs:
%
% rdms: ndissimilarity by nfeatures matrix
% spheres: logical nfeatures by nfeatures matrix where the rows index each
% searchlight centre and the columns index the features that are included
% in the searchlight (this will be symmetrical if using radius-based
% mapping but not if using nvox-based mapping).
%
% meanrdms = averagesearchlights(rdms,spheres)
function meanrdms = averagesearchlights(rdms,spheres)

meanrdms = NaN(size(rdms));

parfor v = 1:size(spheres,2);
    meanrdms(:,v) = mean(rdms(:,full(spheres(:,v)~=0)),2);
end
