% vol = MriVolume([data],[mask],[varargin])
% Main data container class. Very flexible input
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
% named input arguments (used here):
% header: struct containing header info from spm_vol 
% voxsize: [1 by 3] size of voxels in mm
%
% named input arguments (used in super-class):
% metasamples: struct of numeric or cell arrays with nsamples by 1 shape
% metafeatures: struct numeric or cell arrays with 1 by nfeatures shape
% OR meta
% frameperiod: scalar defining TR
classdef MriVolume < Volume
    properties
        mask % binary mask for analysis
        linind % 1 by nfeatures linear indices into mask
        xyz % 3 by nfeatures coordinates into mask
        header % mapped volume header for mask (optional)
        voxsize % size of voxels in mm
    end

    methods
        function mrivol = MriVolume(data,mask,varargin)
            if ieNotDefined('data')
                data = {[]};
            end
            if ieNotDefined('mask')
                mask = [];
            end
            % instance from Base class
            mrivol = mrivol@Volume(data,varargin{:});
            % then any non-base class inputs
            getArgs(varargin,{'header',[],'voxsize',[]},'verbose=0',...
                'suppressUnknownArgMessage=1');
            if isa(data,'MriVolume')
                mask = setifunset(mask,data.mask,true);
                header = setifunset(header,data.header);
                voxsize = setifunset(voxsize,data.voxsize);
            end
            if isa(mask,'MriVolume')
                % interesting!
                header = setifunset(header,mask.header,true);
                voxsize = setifunset(mask.voxsize,true);
                mask = mask.mask;
            end
            [mrivol.header,mrivol.voxsize,mrivol.mask] = deal(header,...
                voxsize,mask);
            % transpose because find returns a column vector and we want
            % our samples in rows
            mrivol.linind = find(mrivol.mask)';
            mrivol.xyz = linind2coord(mrivol,mrivol.linind);
            % this is inelegant but necessary for searchlights and other
            % ROI volumes where there may be only a mask and no data
            if isempty(mrivol.nfeatures)
                mrivol.nfeatures = numel(mrivol.linind);
            end
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
            datavec2nifti(datavec,self.mask,outpath,self.header);
        end

        function featind = linind2featind(self,linind)
        % return only those of the input linear indices that are inside the
        % mask.
        % featind = linind2featind(linind)
            [junk,featind] = intersect(self.linind,linind);
        end

        function smooth(self,fwhm)
        % smooth the data in place along the feature dimension (space)
        % using a gaussian filter with fwhm mm. Smoothing is done in the
        % original 3D shape of the volume and the mask is re-applied after
        % smoothing.
        %
        % smooth(self,fwhm)
        self.data = smoothdatavecs(self.data,fwhm,self.mask~=0,...
            self.voxsize);
        end

        function vol = copy(self,dat,meta)
            % may have to update mask as in subsref below
            vol = feval(class(self),dat,self.mask,meta,'header',...
                self.header,'voxsize',self.voxsize,...
                'frameperiod',self.frameperiod);
        end

        function varargout = subsref(self,s)
        % overloading of round bracket operator.
        % varargout = subsref(self,s)
            switch s(1).type
                case '()'
                    % basic parsing
                    [dat,meta] = basesubsref(self,s);
                    % MriVolume specific parsing
                    mask = self.mask;
                    % if we are also indexing in feature dimension
                    if (length(s.subs) > 1) && (~isempty(self.mask))
                        mask = false(self.header.dim);
                        % map featinds to linind
                        mask(self.linind(s.subs{2})) = 1;
                    end
                    % make a new instance
                    varargout{1} = feval(class(self),dat,mask,'meta',meta,...
                        'header',self.header,'voxsize',self.voxsize,...
                        'frameperiod',self.frameperiod);
                otherwise
                    % revert to builtin behaviour
                    try
                        varargout{1} = builtin('subsref',self,s);
                    catch
                        % maybe no output?
                        builtin('subsref',self,s);
                    end
            end
        end

        function o = horzcat(varargin)
        % o = horzcat(varargin)
            if nargin==1
                o = varargin{1};
            else
                error(['concatenation in feature dimension is not ' ...
                    'supported for MriVolume classes']);
            end
        end
    end
end
