% GLM sub-class where a set of covariates are projected out of data and
% design matrix prior to fitting. 
%
% Mainly useful when the number of covariates might vary (e.g., ncovtouse
% is a free parameter tuned with crossvalidateproperty as part of
% GLMdenoise-like fits), or when the design matrix is convolved on the fly
% (see ConvGLM sub-class).
%
% model = CovariateGLM(X,data,covariates)
classdef CovariateGLM < GLM
    properties
        covariates
        ncovariates
        ncovtouse
    end

    methods
        function gl = CovariateGLM(X,data,covariates)
            if ~any(nargin)
                [X,data,covariates] = deal([]);
            end
            gl = gl@GLM(X,data);
            assert(gl.nsamples==size(covariates,1),...
                'covariates do not match nsamples');
            gl.covariates = covariates;
            gl.ncovariates = size(covariates,2);
            gl.ncovtouse = gl.ncovariates;
        end

        function data = getdatac(self)
            nrun = numel(self);
            % filter each run separately to preserve independence
            data = cell(nrun,1);
            for r = 1:nrun
                data{r} = getprojectionmatrix(self(r)) * ...
                    self(r).data;
            end
        end

        function X = getdesign(self)
            nrun = numel(self);
            Xcovcache = cell(nrun,1);
            for r = 1:nrun
                Xcovcache{r} = getprojectionmatrix(self(r)) * ...
                    self(r).X;
            end
            X = vertcat(Xcovcache{:});
        end

        function pmat = getprojectionmatrix(self)
            assert(numel(self)==1,['covariates should be processed ' ...
                'separately for each run']);
            pmat = projectionmatrix(self.covariates(:,1:self.ncovtouse));
        end
    end
end
