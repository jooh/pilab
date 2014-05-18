% Convert a 2D vector to matrix form by means of the logical indices in
% mask.
%
% outsidemask: (default 'zeros') 'NaN' is also possible
%
% datamat = datavec2mat(data,mask,[outsidemask])
function datamat = datavec2mat(data,mask,outsidemask)

if ieNotDefined('outsidemask')
    outsidemask = 'zeros';
end

dims = size(mask);
if strcmp(outsidemask,'zeros')
    if isa(data,'logical')
        datamat = false(dims);
    else
        datamat = zeros(dims,class(data));
    end
else
    % here you may crash if you e.g. enter logical data
    datamat = feval(outsidemask,dims,class(data));
end

datamat(mask) = data;
