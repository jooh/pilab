classdef Volume < handle
    % vol = Volume([data],[mask],[varargin])
    % Main data container class for pattern analysis. Very flexible input
    % handling: data and mask can be chars with paths to volumes, 3D
    % matrices or another Volume instance (for data, 2D n by nvox matrices
    % and 4D volume matrices are also permissible). The design of this
    % class owes no small debt to PyMVPA's Dataset.
    %
    % This class overloads Matlab's standard concatenation
    % functionality so it is possible to do e.g. volboth = [vola volb] to
    % obtain a sensibly combined Volume instance.
    %
    % Similarly, direct indexing and assignment are overloaded so that e.g.
    % vol(1,:) returns the first data point, or vol(:,1) = 0 sets the first
    % feature to 0 across all data points.
    %
    % Optional named arguments:
    % labels: n by 1 cell array with condition names
    % chunks: n by 1 vector with data split rules (e.g. run indices)
    % order: n by 1 vector providing the correct sequence for datapoints
    % names n by 1 cell array of names for each feature
    properties
        mask % binary mask for analysis
        lininds % linear indices into mask
        V = struct('dim',[]) % mapped volume header for mask (optional)
        voxsize = [1 1 1] % size of voxels in mm
        nfeatures % number of in-mask voxels (columns in data)
        ndata % number of data points (rows in data)
        data % ndata by nvox matrix
        names = {} % ndata by 1 cell array of raw file names for data
        labels = {}; % ndata by 1 cell array
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
                'names',{},'V',struct('dim',[])});
            [vol.labels,vol.chunks,vol.order,vol.names,vol.V] = deal(...
                labels(:),chunks(:),order(:),names(:),V(:));
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
                [vol.mask,vol.V] = deal(mask.mask,mask.V);
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
                vol.names = [vol.names; volnames];
            end
            % analyse the meta data
            vol.ndata = size(vol.data,1);
            assert(isempty(vol.labels) || length(vol.labels)==vol.ndata,...
                'number of labels does not match ndata');
            vol.uniquelabels = unique(vol.labels);
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
        end

        function sref = subsref(a,s)
            switch s(1).type
                case '()'
                    assert(length(s)==1,'cannot subindex Volume instance');
                    % pull out data field 
                    dat = builtin('subsref',a.data,s);
                    % handle special cases
                    ind = s.subs{1};
                    if strcmp(ind,':')
                        ind = 1:a.ndata;
                    end
                    labels = {};
                    if ~isempty(a.labels)
                        labels = a.labels(ind);
                    end
                    chunks = [];
                    if ~isempty(a.chunks)
                        chunks = a.chunks(ind);
                    end
                    order = [];
                    if ~isempty(a.order)
                        order = a.order(ind);
                    end
                    names = {};
                    if ~isempty(a.names)
                        names = a.names(ind);
                    end
                    % update mask if there is one and if we are selecting a
                    % subset of features
                    if length(s.subs) > 1 && ~isempty(a.mask) && ...
                            ~strcmp(s.subs{2},':')
                        mask = logical(zeros(a.V.dim));
                        mask(s.subs{2}) = 1;
                    else
                        % just use the same mask again
                        mask = a.mask;
                    end
                    % make a new Volume instance
                    sref = Volume(dat,mask,'labels',labels,'chunks',...
                        chunks,'order',order,'names',names);
                    % and make sure volume header follows
                    sref.V = a.V;
                otherwise
                    % revert to builtin behaviour
                    sref = builtin('subsref',a,s);
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

        function xyz = ascoord(self,linind)
        % xyz = ascoord(self,linind)
            [x,y,z] = ind2sub(self.V.dim,linind);
            xyz = [x,y,z];
        end

        function linind = aslinind(self,coords)
        % linind = aslinind(self,coords)
            linind = sub2ind(self.V.dim,coords(:,1),coords(:,2),...
                coords(:,3));
        end

        function mat = asmat(self,datavec)
        % mat = asmat(self,datavec)
            assert(length(datavec)==self.nfeatures,...
                'datavec length does not match nfeatures')
            mat = zeros(self.vol.V.dim);
            mat(self.lininds) = datavec;
        end

        function asvol(self,datavec,outpath)
        % asvol(self,datavec,outpath)
            outV = self.V;
            outV.fname = outpath;
            mat = self.asmat(datavec);
            spm_write_vol(outV,mat);
        end

    end
end
