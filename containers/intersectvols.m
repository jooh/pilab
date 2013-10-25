function varargout = intersectvols(varargin)

if ~nargin
    return;
end

varargout = varargin;
maskmat = cell(1,nargin);
for n = 1:nargin
    % first remove any nans
    nanmask = any(isnan(varargout{n}.data),1);
    if any(nanmask)
        varargout{n} = varargout{n}(:,~nanmask);
    end
    % then place the de-naned mask in a cell
    maskmat{n} = varargout{n}.mask;
end
% setup a combined mask
fullmat = cat(4,maskmat{:});
standardmask = all(fullmat,4);
allok = find(standardmask);

% and pass over again, dropping any inconsistent features
if ~isequal(maskmat{:})
    for n = 1:nargin
        inds = varargout{n}.linind2featind(allok);
        varargout{n} = varargout{n}(:,varargout{n}.linind2featind(allok));
    end
end
