classdef MriVolume < BaseVolume
    % vol = MriVolume([data],[mask],[varargin])
    % Main data container class for pattern analysis. Very flexible input
    % handling: data and mask can be chars with paths to volumes, 3D
    % matrices or another MriVolume instance (for data, 2D n by nfeatures
    % matrices and 4D volume matrices are also permissible). The design of
    % this class owes no small debt to PyMVPA's Dataset.
    %
    % This class overloads Matlab's standard concatenation
    % functionality so it is possible to do e.g. volboth = [vola; volb] to
    % obtain a sensibly combined MriVolume instance.
    %
    % Similarly, direct indexing is overloaded so that e.g.
    % vol(1,:) returns a volume consisting only of the first data point.
    %
    % Optional named input arguments:
    % header: struct containing header info from spm_vol 
    % metasamples: struct of numeric or cell arrays with nsamples by 1 shape
    % metafeatures: struct numeric or cell arrays with 1 by nfeatures shape
    properties
        mask % binary mask for analysis
        linind % 1 by nfeatures linear indices into mask
        xyz % 3 by nfeatures coordinates into mask
        header = struct('dim',[]) % mapped volume header for mask (optional)
        voxsize = [1 1 1] % size of voxels in mm
        limitnvol = Inf; % how many volumes to attempt to load at once
    end

    methods
        function vol = MriVolume(data,mask,varargin)
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
            % handle base meta inputs
            vol = vol.initialisemeta(varargin{:});
            % handle non-base inputs
            getArgs(varargin,{'header',struct('dim',[])},'verbose=0',...
                'suppressUnknownArgMessage=1');
            vol.header = header;
            if isfield(vol.header,'mat')
                vol.voxsize = vox2mm(vol.header);
            end
            % Can also be unset at a later stage if meta.name already
            % exists
            if isempty(vol.meta.samples.names)
                donames = 1;
            else
                donames = 0;
            end
            % parse mask input
            if isa(mask,'char')
                % assume it's a SPM-compatible volume
                vol.header = spm_vol(mask);
                vol.voxsize = vox2mm(vol.header);
                vol.mask = spm_read_vols(vol.header) > 0;
                vol.nfeatures = sum(vol.mask(:)>0);
            elseif isa(mask,'MriVolume')
                % interesting!
                % but here we need to be careful not to include data
                [vol.mask,vol.header,vol.nfeatures,vol.voxsize] = deal(...
                    mask.mask,mask.header,mask.nfeatures,mask.voxsize);
            elseif ~isempty(mask)
                % assume header-less mask
                % (but only if we haven't passed a header input argument)
                if ~isfield(vol.header,'mat')
                    vol.header = struct('dim',size(mask));
                end
                vol.mask = mask > 0;
                vol.nfeatures = sum(vol.mask(:)>0);
            end
            % transpose because find returns a column vector and we want
            % our samples in rows
            vol.linind = find(vol.mask)';
            vol.xyz = vol.linind2coord(vol.linind);
            % wrap in cell to allow iteration
            if ~iscell(data)
                if ischar(data) && ~isinf(vol.limitnvol)
                    % split into manageable segments of n
                    % volumes
                    n = size(data,1);
                    indstart = 1:vol.limitnvol:n;
                    indend = indstart+vol.limitnvol-1;
                    indend(end) = n;
                    nsplit = length(indstart);
                    celldata = cell(1,nsplit);
                    for sp = 1:nsplit
                        celldata{sp} = data(indstart(sp):indend(sp),:);
                    end
                    data = celldata;
                    clear celldata;
                else
                    data = {data};
                end
            end
            for d = 1:length(data)
                thisdata = data{d};
                volnames = {};
                if isa(thisdata,'char')
                    % assume it's a SPM-compatible volume
                    dV = spm_vol(thisdata);
                    % check that you aren't trying to combine data with
                    % a weird mask
                    if isfield(vol.header,'mat')
                        % make sure the headers match
                        assert(spm_check_orientations(...
                            [vol.header; dV(:)]),'header mismatch');
                    else
                        % set header by data instead of mask
                        vol.header = dV;
                        vol.voxsize = vox2mm(vol.header);
                    end
                    datamat = spm_get_data(dV,vol.xyz);
                    volnames = {dV.fname}';
                elseif isa(thisdata,'MriVolume')
                    % interesting! 
                    if isfield(vol.header,'mat')
                        assert(spm_check_orientations(...
                            [thisdata.header; vol.header]),...
                            'data header does not match mask');
                    end
                    datamat = thisdata.data;
                    % update meta
                    vol.appendmetasamples(thisdata.meta.samples);
                    vol.appendmetafeatures(thisdata.meta.features);
                    % also bring along the mask if possible and needed
                    if isempty(vol.mask) && ~isempty(thisdata.mask)
                        vol.mask = thisdata.mask;
                        vol.linind = thisdata.linind;
                        vol.xyz = vol.linind2coord(vol.linind);
                        vol.header = thisdata.header;
                        if isfield(vol.header,'mat')
                            vol.voxsize = vox2mm(vol.header);
                        end
                    end
                    if isempty(vol.frameperiod)
                        vol.frameperiod = thisdata.frameperiod;
                    else
                        assert(vol.frameperiod == thisdata.frameperiod,...
                            'data with different frameperiod');
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
                dtype = class(datamat);
                % initialise with class of first data. Helps conserve
                % memory if you are using a lower precision datatype
                % (otherwise the concatenation operation upcasts all data
                % to double)
                if isempty(vol.data)
                    vol.data = feval(dtype,[]);
                end
                switch ndims(datamat);
                    case {3,4}
                        % it's a 3D/4D volume
                        if ~isempty(vol.header.dim)
                            assert(all(tdsize(1:3)==vol.header.dim),...
                                'data dims do not match mask');
                        else
                            % read dims from data
                            vol.header.dim = size(datamat);
                            % avoid 4D case
                            vol.header.dim = vol.header.dim(1:3);
                        end
                        % deal with no mask case
                        if isempty(vol.linind)
                            % let's assume you want implicit masking
                            vol.linind = find(thisvol(:,:,:,1))';
                            vol.mask = true(tdsize(1:3));
                        end
                        % preallocate to speed things up a bit
                        nvol = size(datamat,4);
                        nsamples = size(vol.data,1);
                        vol.data(nsamples+1:nvol,1:length(vol.linind)) = NaN;
                        % in principle we could do reshape here but in my
                        % experience it is too memory-heavy
                        for v = 1:nvol
                            vtemp = datamat(:,:,:,v);
                            vol.data(nsamples+v,:) = vtemp(vol.mask);
                        end
                        clear datamat
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
                    vol.meta.samples.names = [vol.meta.samples.names; volnames];
                end
            end
            vol.nsamples = size(vol.data,1);
            % analyse the standard descriptives
            vol.checkmeta;
        end

        function xyz = linind2coord(self,linind)
        % convert linear indices (1 by n shape) to xyz voxel coordinates (3
        % by n shape).
        % xyz = linind2coord(linind)
            [x,y,z] = ind2sub(self.header.dim,linind);
            xyz = [x; y; z];
        end

        function linind = coord2linind(self,coords)
        % convert a set of 3 by n voxel coordinates to linear indices
        % linind = coord2linind(coords)
            linind = sub2ind(self.header.dim,coords(1,:),coords(2,:),...
                coords(3,:));
        end

        function mat = data2mat(self,datavec)
        % convert a 1 by n data vector to 3D voxel space.
        % mat = data2mat(datavec)
            assert(length(datavec)==self.nfeatures,...
                'datavec length does not match nfeatures')
            % deal with Matlab's failure to recognise logical as a possible
            % class for zeros...
            if isa(datavec,'logical')
                mat = false(self.header.dim);
            else
                mat = zeros(self.header.dim,class(datavec));
            end
            mat(self.linind) = datavec;
        end

        function data2file(self,datavec,outpath)
        % export a 1 by n data vector to an spm_write_vol-compatible file
        % format (typically nii,hdr/img).
        % data2file(datavec,outpath)
            outV = self.header;
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
        % return only those of the input linear indices that are inside the
        % mask.
        % featind = linind2featind(linind)
            [junk,featind] = intersect(self.linind,linind);
        end

        function smooth(self,fwhm)
        % smooth(fwhm)
        % smooth the data in place along the feature dimension (space)
        % using a gaussian filter with fwhm mm. Smoothing is done in the
        % original 3D shape of the volume and the mask is re-applied after
        % smoothing.
            sigma = fwhm / sqrt(8*log(2));
            % configure filter
            % convert fwhm to standard deviation of gaussian
            % get sigma in voxel size units
            sigmavox = sigma ./ self.voxsize;
            % make a filter
            % filter size should be big enough to cover gaussian
            fsize = ceil(3*sigmavox)*2+1;
            % 1D filter placeholder
            dsize = ones(1,3);
            for dat = 1:self.nsamples
                datmat = self.data2mat(self.data(dat,:));
                % take advantage of the fact that smoothing thrice with 3 1D
                % gaussians is equivalent to but faster than smoothing once
                % with a 3D gaussian (also gets around the issue that
                % smooth3 requires a symmetric 3D kernel, with a scalar
                % sigma)
                for dim = 1:3
                    dimsize = dsize;
                    dimsize(dim) = fsize(dim);
                    datmat = smooth3(datmat,'gaussian',dimsize,...
                        sigmavox(dim));
                end
                % return to data
                self.data(dat,:) = datmat(self.mask);
            end
        end

        function vol = copy(self,dat,meta);
            vol = MriVolume(dat,self.mask,'metasamples',meta.samples,...
                'metafeatures',meta.features,'header',self.header,...
                'frameperiod',self.frameperiod);
        end

        function varargout = subsref(a,s)
        % overloading of round bracket operator.
        % varargout = subsref(a,s)
            switch s(1).type
                case '()'
                    % basic parsing
                    [dat,meta] = basesubsref(a,s);
                    % MriVolume specific parsing
                    mask = a.mask;
                    % if we are also indexing in feature dimension
                    if (length(s.subs) > 1) && (~isempty(a.mask))
                        mask = false(a.header.dim);
                        % map featinds to linind
                        mask(a.linind(s.subs{2})) = 1;
                    end
                    % make a new instance
                    varargout{1} = MriVolume(dat,mask,'metasamples',meta.samples,...
                        'metafeatures',meta.features,'header',a.header,...
                        'frameperiod',a.frameperiod);
                otherwise
                    % revert to builtin behaviour
                    try
                        varargout{1} = builtin('subsref',a,s);
                    catch
                        % maybe no output?
                        builtin('subsref',a,s);
                    end
            end
        end
    end
end
