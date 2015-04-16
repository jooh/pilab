classdef CovariateGLM < GLM
    properties
        covariates
    end

    methods
        function gl = CovariateGLM(X,data,covariates)
            if ~any(nargin)
                X = [];
                data = [];
                covariates = [];
            end
            gl = gl@GLM(X,data);
            assert(gl.nsamples==size(covariates,1),...
                'covariates do not match nsamples');
            gl.covariates = covariates;
        end

        function data = getdata(self)
            nrun = numel(self);
            data = cell(nrun,1);
            % filter each run separately to preserve independence
            for r = 1:nrun
                data{r} = projectout(getdata@GLM(self(r)),...
                    getcovariates(self(r)));
            end
            data = vertcat(data{:});
        end

        function X = getdesign(self)
            nrun = numel(self);
            X = cell(nrun,1);
            for r = 1:nrun
                X{r} = projectout(getdesign@GLM(self(r)),...
                    getcovariates(self(r)));
            end
            X = vertcat(X{:});
        end

        function covariates = getcovariates(self)
            covariates = vertcat(self.covariates);
        end
    end
end
