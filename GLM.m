% General class for GLM fitting and cross-validation. 
% gl = GLM(X,Y);
classdef GLM < handle
    properties
        X % design matrix (nsamples by nregressors)
        Y % data (nsamples by nfeatures)
        betas % parameter estimates (nregressors by nfeatures)
        n % number of observations (rows in X and Y)
        nfeatures % number of features (columns in Y)
        npredictors % number of regressors in model (columns in X)
        df % degrees of freedom (scalar)
        mrss % mean residual sum of squares (only initialised if needed)
        covmat % covariance matrix of X (only initialised if needed)
    end

    methods
        function gl = GLM(X,Y)
            % support sub-classing
            if ~any(nargin)
                return
            end
            [gl.X,gl.Y] = deal(X,Y);
            gl.getdiagnostics;
            % fit model
            gl.betas = olsfit(X,Y);
        end

        function getdiagnostics(self)
            % store other model fit parameters for contrasts etc
            self.n = size(self.Y,1);
            assert(self.n == size(self.X,1),'design matrix does not match Y');
            self.nfeatures = size(self.Y,2);
            self.npredictors = size(self.X,2);
            self.df = self.n - rank(self.X);
        end

        function cmat = getcovmat(self)
        % Compute the X covariance matrix if necessary. Otherwise, return
        % the one stored in self.covmat.
        % cmat = getcovmat(self)
            if isempty(self.covmat)
                self.covmat = inv(self.X' * self.X);
            end
            cmat = self.covmat;
        end

        function mrss = getmrss(self)
        % Compute the mean residual sum of squares if necessary. Otherwise,
        % return the one stored in self.mrss.
            if isempty(self.mrss)
                resid = self.predictY(self.X) - self.Y;
                rss = sum(resid.^2);
                self.mrss = rss / self.df;
            end
            mrss = self.mrss;
        end

        function Yfit = predictY(self,X)
        % Yfit = predictY(self,X)
        % generate fitted (predicted) responses for a design matrix X.
            Yfit = X * self.betas;
        end

        function mse = mserr(self,Yfit)
        % mse = mserr(self,Yfit)
        % mean squared error between predicted responses and current Y.
            mse = mean((Yfit-self.Y).^2);
        end

        function R = rsquare(self,Yfit)
        % R = rsquare(self,Yfit)
        % R^2 between predicted responses and current Y.
            R = rsquare(Yfit,self.Y);
        end

        function con = contrast(self,cv)
        % con = contrast(cv)
        % compute a contrast for self fit based on contrast vector cv.
            assert(ndims(cv)==2 && all(size(cv)==[1 self.npredictors]),...
                'contrast vector must be 1 by npredictors');
            con = cv * self.betas;
        end

        function t = tmap(self,cv)
        % t = tmap(cv)
        % compute a T contrast for self fit based on contrast vector cv.
            % make sure we have computed a covariance matrix
            self.getcovmat;
            self.getmrss;
            t = self.contrast(cv) ./ sqrt(self.mrss * ...
                (cv * self.covmat * cv'));
        end

        function p = pmap(self,cv)
        % p = pmap(cv)
        % return one-tailed p values for self fit based on contrast vector
        % cv.
            p = 1-tcdf(abs(self.tmap(cv)),self.df);
        end
    end
end
