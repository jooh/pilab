% Container class for storing the onsets and conditions for a study for
% later convolution with convolveonsets (see designmatrix method).
%
% INPUTS
% design: struct with the fields onsets, conind, n, chunk, covariates,
%   convargs and frameperiod
% frameperiod: optional global frameperiod that get plugged into each
%   design
% names: go into meta.features.names for compatibility with Volume and
%   processing code
% varargin: optional global convargs inputs that get plugged into each
%   design
%
% ds = DesignSpec(design,frameperiod,names,varargin)
classdef DesignSpec < handle
    properties
        data
        nsamples
        nfeatures
        meta
    end

    methods
        function ds = DesignSpec(design,frameperiod,names,varargin)
            if ~any(nargin)
                design = emptystruct('onsets','conind','n','chunk',...
                    'covariates','convargs','frameperiod');
            end
            if iscell(design)
                % attempt to unpack, hopefully into struct array
                design = [design{:}];
                assert(isstruct(design),...
                    'input design must be struct or cell array of struct');
            end
            ds.data = design;
            % make sure chunks are coded appropriately
            for c = 1:numel(ds.data)
                if ~isfield(ds.data,'chunk') || isempty(ds.data(c).chunk)
                    ds.data(c).chunk = c;
                end
                if isscalar(ds.data(c).chunk)
                    ds.data(c).chunk = ones(ds.data(c).n,1) * ...
                        ds.data(c).chunk;
                end
                if ~isfield(ds.data,'convargs') || isempty(ds.data(c).convargs)
                    ds.data(c).convargs = varargin;
                end
                if ~isfield(ds.data,'frameperiod') || isempty(ds.data(c).frameperiod)
                    ds.data(c).frameperiod = frameperiod;
                end
            end
            ds.meta.samples.chunks = vertcat(ds.data.chunk);
            ds.meta.features.names = names;
            ds.nfeatures = numel(unique(design(1).conind));
            ds.nsamples = sum([design.n]);
        end

        function [X,errs] = designmatrix(self)
            for d = 1:numel(self.data)
                [X{d},errs{d}] = convolveonsets(self.data(d).onsets,...
                    self.data(d).conind,self.data(d).frameperiod,...
                    self.data(d).n,self.data(d).convargs{:});
            end
        end

        function ds = selectbymeta(self,target,vals)
            assert(strcmp(target,'chunks'),...
                'only chunks are supported at present');
            chunks = arrayfun(@(x)self.data(x).chunk(1),1:numel(self.data))
            [~,ind] = intersect(chunks,vals,'stable');
            ds = DesignSpec(self.data(ind),[],self.meta.features.names);
        end
    end
end
