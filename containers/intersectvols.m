% insure that a set of volume instances contain the same features. Also
% removes nans. Supports MriVolumes and other BaseVolume sub-classes with the mask
% property. If the mask field is not available (e.g., if it's a BaseVolume)
% we attempt to intersect by using meta.features.names instead.
%
% varargout = intersectvols(varargin)
function varargout = intersectvols(varargin)

if ~nargin
    return;
end

varargout = varargin;
if isprop(varargin{1},'mask')
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
elseif isfield(varargin{1}.meta.features,'names') && ~isempty(varargin{1}.meta.features.names)
    % intersect by name. This is a little hairy.
    % first remove any NaN ROIs
    for n = 1:nargin
        nanmask = any(isnan(varargout{n}.data),1);
        if any(nanmask)
            varargout{n} = varargout{n}(:,~nanmask);
        end
        allrois{n} = varargout{n}.meta.features.names;
    end
    % then get the unique surviving ROIs
    [urois,nvalid] = count_unique(horzcat(allrois{:}));
    % drop the ones that have less than full sample size
    urois(nvalid~=nargin) = [];
    % reconstruct each individual disvol with only the right ROIs in a
    % consistent order
    for n = 1:nargin
        % find the indices for this subject
        [~,inds] = intersect(allrois{n},urois);
        varargout{n} = varargout{n}(:,inds);
    end
else
    error(['intersectvols requires an input vol with a mask property ' ...
        'or entries in meta.features.names']);
end
