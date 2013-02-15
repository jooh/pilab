% General class for GLM fitting with OLS. Can be sub-classed for special
% cases (e.g. ridge regression in RidgeGLM).
%
% TODO: multiple inheritance! GLMVolume < MriVolume & GLM!
%
% Initialise with:
% gl = GLM(X,data);
%
% Get a fit with gl.getestimates;
%
% In general, uses lots of get methods to compute properties on an
% as-needed basis for speed (e.g., mrss and covmat are only computed if you
% are doing something like computing contrasts). This also affords
% flexibility for sub-classes, which need only replace the fit method.
classdef GLM < handle
    properties
        X % design matrix (nsamples by nregressors)
        data % data (nsamples by nfeatures)
        nsamples % number of observations (rows in X and data)
        nfeatures % number of features (columns in data)
        npredictors % number of regressors in model (columns in X)
        df % degrees of freedom (scalar)
        estimates % parameter estimates (nregressors by nfeatures)
        mrss % mean residual sum of squares (only initialised if needed)
        covmat % covariance matrix of X (only initialised if needed)
    end

    methods
        function gl = GLM(X,data)
            % support sub-classing
            if ~any(nargin)
                return
            end
            gl.nsamples = size(data,1);
            assert(gl.nsamples == size(X,1),'design matrix does not match data');
            [gl.X,gl.data] = deal(X,data);
            gl.nfeatures = size(gl.data,2);
            gl.npredictors = size(gl.X,2);
            gl.df = gl.nsamples - rank(gl.X);
        end

        function estimates = getestimates(self)
            if isempty(self.estimates)
                self.estimates = self.fit;
            end
            estimates = self.estimates;
        end

        function estimates = fit(self);
            estimates = olsfit(self.X,self.data);
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
                resid = self.predictY(self.X) - self.data;
                rss = sum(resid.^2);
                self.mrss = rss / self.df;
            end
            mrss = self.mrss;
        end

        function Yfit = predictY(self,X)
        % Yfit = predictY(self,X)
        % generate fitted (predicted) responses for a design matrix X.
            Yfit = X * self.getestimates;
        end

        function mse = mserr(self,Yfit)
        % mse = mserr(self,Yfit)
        % mean squared error between predicted responses and current data.
            mse = mean((Yfit-self.data).^2);
        end

        function R = rsquare(self,Yfit)
        % R = rsquare(self,Yfit)
        % R^2 between predicted responses and current data.
            R = rsquare(Yfit,self.data);
        end

        function con = contrast(self,cv)
        % con = contrast(cv)
        % compute a contrast for self fit based on contrast vector cv.
            assert(ndims(cv)==2 && all(size(cv)==[1 self.npredictors]),...
                'contrast vector must be 1 by npredictors');
            con = cv * self.getestimates;
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
