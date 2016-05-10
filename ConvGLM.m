% CovariateGLM sub-class where a set of onsets and condition labels are
% used to re-build a convolved design matrix online. For documentation on
% inputs, see convolveonsets. Any varargin are passed on to this function.
%
% This class has significant performance costs compared to using a
% pre-convolved design matrix but is useful for permutation tests based on
% reassigning the labels (conind). We override drawpermsample to support
% such tests here.
%
% model = ConvGLM(onsets,conind,frameperiod,data,covariates,[varargin])
classdef ConvGLM < CovariateGLM
    properties
        onsets
        conind
        frameperiod
        convargs
    end

    methods
        function gl = ConvGLM(onsets,conind,frameperiod,data,covariates,varargin)
            convargs = varargin;
            if nargin == 0 || isempty(onsets)
                [X,onsets,conind,frameperiod,data,covariates] = deal([]);
            else
                X = convolveonsets(onsets,conind,frameperiod,size(data,1),...
                    convargs{:});
            end
            gl = gl@CovariateGLM(X,data,covariates);
            [gl.onsets,gl.conind,gl.frameperiod,gl.convargs] = deal(onsets,...
                conind,frameperiod,convargs);
            % update number of independent observations to number of
            % trials, not number of samples
            gl.nrandsamp = numel(onsets);
        end

        function X = getdesign(self)
            nrun = numel(self);
            for r = 1:nrun
                % re-convolve
                self(r).X = convolveonsets(self(r).onsets,...
                    self(r).conind,self(r).frameperiod,self(r).nsamples,...
                    self(r).convargs{:});
            end
            % use super-class to do covariate filter
            X = getdesign@CovariateGLM(self);
        end

        function model = drawpermsample(self,inds)
        % return a new instance where the condition labels in the design
        % matrix have been re-ordered according to inds. Note that you must
        % supply the same number of inds as numel(self(1).conind).
        %
        % model = drawpermsample(self,inds)
            n = cellfun(@numel,{self.conind});
            assert(numel(inds)==n(1),'got %d inds for %d samples',...
                numel(inds),n(1));
            assert(all(n(1) == n),['design must have same number of ' ...
                'events in each run.'])
            model = copy(self);
            for r = 1:length(self)
                model(r).conind = model(r).conind(inds);
            end
            % flatten to GLM class for speed (prevents all of these
            % filtering operations from happening unnecessarily after the
            % permutation has been done)
            model = asglm(model);
        end

        function model = asglm(self)
            % flatten the ConvGLM to a simple GLM for faster code execution
            for r = 1:numel(self)
                model(r) = GLM(getdesign(self(r)),getdata(self(r)));
            end
        end
    end
end
