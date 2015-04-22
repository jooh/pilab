classdef DenoiseGLM < CovariateGLM
    properties
        noisepcs
        nnoise
        nnoisetouse
        Xdenoisecache
        Xdenoisecachen = NaN;
        datadenoisecache
        datadenoisecachen = NaN;
    end

    methods
        function gl = DenoiseGLM(X,data,covariates,noisepcs)
            if ~any(nargin)
                [X,data,covariates,noisepcs] = deal([]);
            end
            gl = gl@CovariateGLM(X,data,covariates);
            assert(gl.nsamples==size(noisepcs,1),...
                'noisepcs do not match nsamples');
            gl.noisepcs = noisepcs;
            gl.nnoise = size(noisepcs,2);
            gl.nnoisetouse = gl.nnoise;
        end

        function data = getdenoisedatac(self)
            nrun = numel(self);
            % filter each run separately to preserve independence
            for r = 1:nrun
                if self(r).nnoisetouse ~= self(r).datadenoisecachen
                    % nnoise changed so need to re-calculate projection
                    self(r).datadenoisecache = getdenoiseprojectionmatrix(...
                        self(r)) * self(r).data;
                    self(r).datadenoisecachen = self(r).nnoisetouse;
                end
            end
            data = {self.datadenoisecache};
        end

        function data = getdenoisedata(self)
            datac = getdenoisedatac(self);
            data = vertcat(datac{:});
        end

        function design = getdenoisedesign(self)
            nrun = numel(self);
            for r = 1:nrun
                if self(r).nnoisetouse ~= self(r).Xdenoisecachen
                    self(r).Xdenoisecache = getdenoiseprojectionmatrix(...
                        self(r)) * self(r).X;
                    self(r).Xdenoisecachen = self(r).nnoisetouse;
                end
            end
            design = vertcat(self.Xdenoisecache);
        end

        function pmat = getdenoiseprojectionmatrix(self)
            assert(numel(self)==1,['covariates should be processed ' ...
                'separately for each run']);
            pmat = projectionmatrix([self.covariates(:,1:self.ncovtouse) ...
                self.noisepcs(:,1:self.nnoisetouse)]);
        end

        function est = fit(self)
            % by overriding fit we ensure that training methods work on
            % denoised designs and data, while prediction methods use only
            % detrended designs and data (inherited behavior from CovariateGLM)
            est = olsfit(getdenoisedesign(self),getdenoisedatac(self));
        end

        function mahdist = infoc(self,w,conmat)
            % NB we avoid the fit method here to fit without denoising
            c = conmat * olsfit(getdesign(self),getdatac(self));
            % Classification methods based on projecting the discriminant
            % time course (infomodel) do not need this since they do not
            % use the fit method over-ride.
            mahdist = sqrtsigned(diag(c * w'));
        end
    end
end
