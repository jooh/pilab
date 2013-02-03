% General class for GLM fitting and cross-validation. Uses relatively heavy
% intialisation. For speed-dependent operations olsfit is probably quicker.
classdef GLM < handle
    properties
        X % design matrix (nsamples by nregressors)
        Y % data (nsamples by nfeatures)
        betas % parameter estimates (nregressors by nfeatures)
        Yfit % fitted responses (nsamples by nfeatures)
        resid % residual errors (nsamples by nfeatures)
        n % number of observations (rows in X and Y)
        nfeatures % number of features (columns in Y)
        npredictors % number of regressors in model (columns in X)
        df % degrees of freedom (scalar)
        rss % residual sum of square errors
        mrss % mean residual sum of squares
        covmat % covariance matrix of X
        mse % mean squared error of fit (1 by nfeatures)
        R2 % R^2 of fit (1 by nfeatures)
    end

    methods
        function gl = GLM(X,Y)
            [gl.X,gl.Y] = deal(X,Y);
            % fit model
            gl.betas = olsfit(X,Y);
            % generate stats for self fit
            gl.Yfit = gl.predictY(gl.X);
            gl.mse = gl.mserr(gl.Yfit);
            gl.R2 = gl.rsquare(gl.Yfit);
            % store other model fit parameters for contrasts etc
            gl.resid = gl.Yfit - Y;
            gl.n = size(Y,1);
            assert(gl.n == size(X,1),'design matrix does not match Y');
            gl.nfeatures = size(Y,2);
            gl.npredictors = size(X,2);
            gl.df = gl.n - rank(X);
            gl.covmat = inv(X' * X);
            gl.rss = sum(gl.resid.^2);
            gl.mrss = gl.rss / gl.df;
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
            t = self.contrast(cv) ./ sqrt(self.mrss * ...
                (cv * self.covmat * cv'));
        end

        function p = pmap(self,cv)
        % p = pmap(cv)
        % return one-tailed p values for self fit based on contrast vector
        % cv.
            p = 1-tcdf(abs(self.tmap(cv)));
        end
    end
end
