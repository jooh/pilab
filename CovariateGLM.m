classdef CovariateGLM < GLM
    properties
        covariates
        ncovtouse
        Xcache
        datacache
        Xcachencov
        datacachencov
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
            gl.ncovtouse = size(covariates,2);
            % prepare the cache
            pmat = projectionmatrix(covariates);
            gl.Xcache = pmat * X;
            gl.datacache = pmat * data;
            gl.Xcachencov = gl.ncovtouse;
            gl.datacachencov = gl.ncovtouse;
        end

        function data = getdatac(self)
            nrun = numel(self);
            % filter each run separately to preserve independence
            for r = 1:nrun
                if self(r).ncovtouse ~= self(r).datacachencov
                    % covariates changed so need to re-calculate projection
                    self(r).datacache = projectionmatrix(...
                        self(r).covariates(:,1:self(r).ncovtouse)) * ...
                        self(r).data;
                    self(r).datacachencov = self(r).ncovtouse;
                end
            end
            data = {self.datacache};
        end

        function X = getdesign(self)
            nrun = numel(self);
            for r = 1:nrun
                if self(r).ncovtouse ~= self(r).Xcachencov
                    self(r).Xcache = projectionmatrix(...
                        self(r).covariates(:,1:self(r).ncovtouse)) * ...
                        self(r).X;
                    self(r).Xcachencov = self(r).ncovtouse;
                end
            end
            X = vertcat(self.Xcache);
        end

        function model = drawpermsample(self,inds)
            model = drawpermsample@GLM(self,inds);
            % clear the cache to prevent inaccurate filtering
            [model.Xcachencov] = deal(NaN);
            [model.datacachencov] = deal(NaN);
        end

        function model = drawpermflipsample(self,inds)
            model = drawpermflipsample(self,inds);
            % clear the cache to prevent inaccurate filtering
            [model.Xcachencov] = deal(NaN);
            [model.datacachencov] = deal(NaN);
        end
    end
end
