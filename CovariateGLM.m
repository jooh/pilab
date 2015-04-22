classdef CovariateGLM < GLM
    properties
        covariates
        ncovariates
        ncovtouse
        Xcovcache
        Xcovcachen = NaN
        datacovcache
        datacovcachen = NaN
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
            for r = 1:nrun
                if self(r).ncovtouse ~= self(r).datacovcachen
                    % covariates changed so need to re-calculate projection
                    self(r).datacovcache = ...
                        getprojectionmatrix(self(r)) * self(r).data;
                    self(r).datacovcachen = self(r).ncovtouse;
                end
            end
            data = {self.datacovcache};
        end

        function X = getdesign(self)
            nrun = numel(self);
            for r = 1:nrun
                if self(r).ncovtouse ~= self(r).Xcovcachen
                    self(r).Xcovcache = getprojectionmatrix(self(r)) * ...
                        self(r).X;
                    self(r).Xcovcachen = self(r).ncovtouse;
                end
            end
            X = vertcat(self.Xcovcache);
        end

        function pmat = getprojectionmatrix(self)
            assert(numel(self)==1,['covariates should be processed ' ...
                'separately for each run']);
            pmat = projectionmatrix(self.covariates(:,1:self.ncovtouse));
        end

        function model = drawpermsample(self,inds)
            model = drawpermsample@GLM(self,inds);
            % clear the cache to prevent inaccurate filtering
            [model.Xcovcachen] = deal(NaN);
            [model.datacovcachen] = deal(NaN);
        end

        function model = drawpermflipsample(self,inds)
            model = drawpermflipsample(self,inds);
            % clear the cache to prevent inaccurate filtering
            [model.Xcovcachen] = deal(NaN);
            [model.datacovcachen] = deal(NaN);
        end
    end
end
