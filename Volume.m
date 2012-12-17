classdef Volume < handle
    % vol = Volume([data],[mask],[varargin])
    % Main data container class for pattern analysis. Very flexible input
    % handling: data and mask can be chars with paths to volumes, 3D
    % matrices or another Volume instance (for data, 2D n by nvox matrices
    % and 4D volume matrices are also permissible). The design of this
    % class owes no small debt to PyMVPA's Dataset.
    %
    % This class overloads Matlab's standard concatenation
    % functionality so it is possible to do e.g. volboth = [vola; volb] to
    % obtain a sensibly combined Volume instance.
    %
    % Similarly, direct indexing is overloaded so that e.g.
    % vol(1,:) returns a volume consisting only of the first data point.
    %
    % Optional named input arguments:
    % labels: n by 1 cell array with condition names
    % chunks: n by 1 vector with data split rules (e.g. run indices)
    % order: n by 1 vector providing the correct sequence for datapoints
    % names: n by 1 cell array of names for each datapoint
    % featuregroups: 1 by nfeatures matrix
    % V: struct containing header info from spm_vol 
    properties
        mask % binary mask for analysis
        lininds % linear indices into mask
        V = struct('dim',[]) % mapped volume header for mask (optional)
        voxsize = [1 1 1] % size of voxels in mm
        nfeatures % number of in-mask voxels (columns in data)
        featuregroups
        ndata % number of data points (rows in data)
        data % ndata by nvox matrix
        names = {} % ndata by 1 cell array
        labels = {}; % ndata by 1 cell array
        labelinds = []; % ndata by 1 matrix with indices into uniquelabels
        uniquelabels = {}
        nlabels = {}
        chunks = [] % ndata by 1 vector (default ones(n,1))
        uniquechunks = []
        nchunks = []
        order = [] % ndata by 1 vector (default 1:n)
    end

    methods
        function vol = Volume(data,mask,varargin)
            % empty return for potential sub-classing
            if nargin==0
                return
            end
            if ieNotDefined('data')
                data = {[]};
            end
            if ieNotDefined('mask')
                mask = [];
            end
            % assign named inputs (as column vectors)
            getArgs(varargin,{'labels',{},'chunks',[],'order',[],...
                'names',{},'V',struct('dim',[]),'featuregroups',[]});
            [vol.labels,vol.chunks,vol.order,vol.names,vol.V,...
                vol.featuregroups] = deal(...
                labels(:),chunks(:),order(:),names(:),V(:),...
                featuregroups(:)');
            if isfield(vol.V,'mat')
                vol.voxsize = vox2mm(vol.V);
            end
            if isempty(vol.names)
                donames = 1;
            else
                donames = 0;
            end
            % parse mask input
            if isa(mask,'char')
                % assume it's a SPM-compatible volume
                vol.V = spm_vol(mask);
                vol.voxsize = vox2mm(vol.V);
                vol.mask = spm_read_vols(vol.V) > 0;
                vol.nfeatures = sum(vol.mask(:)>0);
            elseif isa(mask,'Volume')
                % interesting!
                % but here we need to be careful not to include data
                [vol.mask,vol.V,vol.nfeatures,vol.voxsize] = deal(...
                    mask.mask,mask.V,mask.nfeatures,mask.voxsize);
            elseif ~isempty(mask)
                % assume header-less mask
                vol.V = struct('dim',size(mask));
                vol.mask = mask > 0;
                vol.nfeatures = sum(vol.mask(:)>0);
            end
            % transpose because find returns a column vector and we want
            % our data in rows
            vol.lininds = find(vol.mask)';
            % wrap in cell to allow iteration
            if ~iscell(data)
                data = {data};
            end
            for d = 1:length(data)
                thisdata = data{d};
                volnames = {};
                if isa(thisdata,'char')
                    % assume it's a SPM-compatible volume
                    dV = spm_vol(thisdata);
                    % check that you aren't trying to combine data with
                    % a weird mask
                    if isfield(vol.V,'mat')
                        assert(spm_check_orientations([vol.V dV]),...
                            'data header does not match mask');
                    else
                        % set header by data instead of mask
                        vol.V = dV;
                        vol.voxsize = vox2mm(vol.V);
                    end
                    % if this work we can safely forget about dV
                    datamat = spm_read_vols(dV);
                    for v = 1:length(dV)
                        [p,volnames{v,1},ext] = fileparts(dV(v).fname);
                    end
                elseif isa(thisdata,'Volume')
                    % interesting! 
                    if isfield(vol.V,'mat')
                        assert(spm_check_orientations(...
                            [thisdata vol.V]),...
                            'data header does not match mask');
                    end
                    datamat = thisdata.data;
                    volnames = thisdata.names;
                    vol.labels = [vol.labels; thisdata.labels];
                    vol.chunks = [vol.chunks; thisdata.chunks];
                    vol.order = [vol.order; thisdata.order];
                    % also bring along the mask if possible and needed
                    if isempty(vol.mask) && ~isempty(thisdata.mask)
                        vol.mask = thisdata.mask;
                        vol.lininds = thisdata.lininds;
                        vol.V = thisdata.V;
                        if isfield(vol.V,'mat')
                            vol.voxsize = vox2mm(vol.V);
                        end
                    end
                else
                    % assume headerless data
                    datamat = thisdata;
                end
                % short-circuit if there's no data
                if isempty(datamat)
                    continue
                end
                % parse the data
                tdsize = size(datamat);
                switch ndims(datamat);
                    case {3,4}
                        % it's a 3D/4D volume
                        if ~isempty(vol.V.dim)
                            assert(all(tdsize(1:3)==vol.V.dim),...
                                'data dims do not match mask');
                        else
                            % read dims from data
                            vol.V.dim = size(datamat);
                            % avoid 4D case
                            vol.V.dim = vol.V.dim(1:3);
                        end
                        % extra call to size necessary because matlab
                        % squeezes singletons
                        for n = 1:size(datamat,4)
                            thisvol = datamat(:,:,:,n);
                            % if you have no mask
                            if isempty(vol.lininds)
                                % let's assume you want implicit masking
                                vol.lininds = find(thisvol)';
                                vol.mask = logical(zeros(tdsize(1:3)));
                                vol.mask(vol.lininds) = 1;
                            end
                            vol.data = [vol.data; ...
                                thisvol(vol.lininds)];
                        end
                    case 2
                        % it had best be a n by nfeatures matrix (if we
                        % know nfeatures)
                        if isempty(vol.nfeatures) || vol.nfeatures==0
                            % otherwise, set nfeatures by data
                            vol.nfeatures = tdsize(2);
                        else
                            assert(tdsize(2)==vol.nfeatures,...
                            'received 2D data with bad shape');
                        end
                        vol.data = [vol.data; datamat];
                    otherwise
                        error('could not parse input data')
                end
                % add any names
                if donames
                    vol.names = [vol.names; volnames];
                end
            end
            % analyse the meta data
            vol.ndata = size(vol.data,1);
            assert(isempty(vol.labels) || length(vol.labels)==vol.ndata,...
                'number of labels does not match ndata');
            [vol.uniquelabels,junk,vol.labelinds] = unique(vol.labels);
            vol.nlabels = length(vol.uniquelabels);
            if isempty(vol.chunks)
                % assume everything is one big chunk
                vol.chunks = ones(vol.ndata,1);
            else
                assert(length(vol.chunks) == vol.ndata,...
                    'number of chunks does not match ndata');
            end
            vol.uniquechunks = unique(vol.chunks);
            vol.nchunks = length(vol.uniquechunks);
            if isempty(vol.order)
                vol.order = (1:vol.ndata)';
            else
                assert(length(vol.order) == vol.ndata,...
                    'order does not match ndata');
            end
            % Finally, resort according to order
            if ~issorted(vol.order)
                [junk,inds] = sort(vol.order);
                vol = vol(inds,:);
            end
        end

        function varargout = subsref(a,s)
            switch s(1).type
                case '()'
                    assert(length(s)==1,'cannot subindex Volume instance');
                    % pull out data field 
                    dat = builtin('subsref',a.data,s);
                    % handle special cases
                    datind = s.subs{1};
                    if strcmp(datind,':')
                        datind = 1:a.ndata;
                    end
                    labels = {};
                    if ~isempty(a.labels)
                        labels = a.labels(datind);
                    end
                    chunks = [];
                    if ~isempty(a.chunks)
                        chunks = a.chunks(datind);
                    end
                    order = [];
                    if ~isempty(a.order)
                        order = a.order(datind);
                    end
                    names = {};
                    if ~isempty(a.names)
                        names = a.names(datind);
                    end
                    featuregroups = a.featuregroups;
                    mask = a.mask;
                    if length(s.subs) > 1 
                        % update mask if there is one and if we are
                        % selecting a subset of features
                        if ~isempty(a.mask) && strcmp(s.subs{2},':')
                            mask = logical(zeros(a.V.dim));
                            % map featinds to lininds
                            mask(a.lininds(s.subs{2})) = 1;
                        end
                        % update featuregroups
                        if ~isempty(a.featuregroups)
                            featuregroups = a.featuregroups(s.subs{2});
                        end
                    end
                    % make a new Volume instance
                    varargout{1} = Volume(dat,mask,'labels',labels,...
                        'chunks',chunks,'order',order,'names',names,...
                        'V',a.V,'featuregroups',featuregroups);
                otherwise
                    % revert to builtin behaviour
                    try
                        varargout{1} = builtin('subsref',a,s);
                    catch
                        builtin('subsref',a,s);
                    end
            end
        end

        function o = horzcat(varargin)
            error('concatenation in feature dimension is not supported')
            o = Volume(varargin);
        end

        function o = vertcat(varargin)
            % construct a volume where the mask is read from the first
            % entry and all instances are passed as is. Should concatenate
            % nicely
            o = Volume(varargin);
        end

        function o = cat(dim,varargin)
            % Identical to horzcat - effectively Volume is a 1D class
            assert(dim==1,'concatenation is only supported in data dim')
            o = Volume(varargin);
        end

        function xyz = linind2coord(self,linind)
        % xyz = linind2coord(self,linind)
            [x,y,z] = ind2sub(self.V.dim,linind);
            xyz = [x,y,z];
        end

        function linind = coord2linind(self,coords)
        % linind = coord2linind(self,coords)
            linind = sub2ind(self.V.dim,coords(:,1),coords(:,2),...
                coords(:,3));
        end

        function mat = data2mat(self,datavec)
        % mat = data2mat(self,datavec)
            assert(length(datavec)==self.nfeatures,...
                'datavec length does not match nfeatures')
            mat = zeros(self.V.dim);
            mat(self.lininds) = datavec;
        end

        function data2file(self,datavec,outpath)
        % data2file(self,datavec,outpath)
            outV = self.V;
            outV.fname = outpath;
            % now write out in appropriate numeric class (avoid quantizing)
            if isnumeric(datavec)
                outV.dt = [spm_type('float64') spm_platform('bigend')];
            else
                outV.dt = [spm_type('int32') spm_platform('bigend')];
            end
            mat = self.data2mat(datavec);
            spm_write_vol(outV,mat);
        end

        function featind = linind2featind(self,linind)
        % featind = linind2featind(self,linind)
            [junk,featind] = intersect(self.lininds,linind);
        end

    end
end
