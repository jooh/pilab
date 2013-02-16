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
        % WOULD GO
        % df % degrees of freedom (scalar)
        % estimates % parameter estimates (nregressors by nfeatures)
        % mrss % mean residual sum of squares (only initialised if needed)
        % covmat % covariance matrix of X (only initialised if needed)
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
        end

        function estimates = fit(self);
            estimates = olsfit(vertcat(self.X),vertcat(self.data));
        end

        function cmat = covmat(self)
            X = vertcat(self.X);
            cmat = inv(X' * X);
        end

        function d = df(self)
            d = sum([self.nsamples]) - rank(vertcat(self.X));
        end

        function r = residuals(self)
            r = self.predictY(vertcat(self.X)) - vertcat(self.data);
        end

        function mr = mrss(self)
            resid = self.residuals;
            rss = sum(resid.^2);
            mr = rss / self.df;
        end

        function Yfit = predictY(self,X)
        % Yfit = predictY(self,X)
        % generate fitted (predicted) responses for a design matrix X.
            Yfit = X * self.fit;
        end

        function mse = mserr(self,Yfit)
        % mse = mserr(self,Yfit)
        % mean squared error between predicted responses and current data.
            mse = mean((Yfit-vertcat(self.data)).^2);
        end

        function R = rsquare(self,Yfit)
        % R = rsquare(self,Yfit)
        % R^2 between predicted responses and current data.
            R = rsquare(Yfit,vertcat(self.data));
        end

        function con = contrast(self,cv)
        % con = contrast(cv)
        % compute a contrast for self fit based on contrast vector cv.
            assert(ndims(cv)==2 && all(size(cv)==[1 self(1).npredictors]),...
                'contrast vector must be 1 by npredictors');
            con = cv * self.fit;
        end

        function t = tmap(self,cv)
        % t = tmap(cv)
        % compute a T contrast for self fit based on contrast vector cv.
            % make sure we have computed a covariance matrix
            t = self.contrast(cv) ./ sqrt(self.mrss * ...
                (cv * self.covmat * cv'));
        end

        function p = pmap(self,cv)
        % p = pmap(cv)
        % return one-tailed p values for self fit based on contrast vector
        % cv.
            p = 1-tcdf(abs(self.tmap(cv)),self.df);
        end

        function [winners,meds] = crossvalidateparameter(self,parameter,values,metric)
        % tune a parameter (char) over a range of values (1D array) using a
        % leave-one-out crossvalidation method with maximum median R2 or
        % min mserr (default R2) as
        % the selection criterion (set with 'metric': 'r2' or 'mse'). winners gives the winning value for each
        % feature (1 by nfeatures), and med gives the median R2/mserr
        % for each value (nvalues by nfeatures).
        %
        % [winners,meds] = crossvalidateparameter(parameter,values,[metric])
            if ieNotDefined('metric')
                metric = 'r2';
            end
            n = numel(self);
            nvalues = length(values);
            nfeatures = self(1).nfeatures;
            splits = logical(eye(n));
            meds = NaN([nvalues nfeatures]);
            orgvalue = [self.(parameter)];
            for v = 1:nvalues
                % set parameter in train and test
                [self.(parameter)] = deal(values(v));
                % cross-validate fit with chosen parameter
                splitres = NaN([n nfeatures]);
                parfor s = 1:n
                    test = splits(s,:);
                    train = ~splits(s,:);
                    switch metric
                        case 'r2'
                            splitres(s,:) = self(test).rsquare(...
                                self(train).predictY(self(test).X));
                        case 'mse'
                            splitres(s,:) = self(test).mserr(...
                                self(train).predictY(self(test).X));
                        otherwise
                            error('unknown metric: %s',metric);
                    end
                end
                % update with median across splits
                meds(v,:) = median(splitres,1);
            end
            % find winning value for each feature
            switch metric
                case 'r2'
                    funhand = @max;
                case 'mse'
                    funhand = @min;
            end
            [perf,inds] = funhand(meds,[],1);
            winners = values(inds);
        end

        function [estimates,sterrs] = bootstrapestimate(self,nboot)
        % estimate a bootstrap distribution by drawing nboot samples with
        % replacement from the available runs in self. Return the median of
        % the bootstrap estimates and half of the 68% range (ie, standard
        % error).
        %
        % [estimates,sterrs] = bootstrapestimate(nboot)
            n = numel(self);
            bootest = NaN([self(1).npredictors self(1).nfeatures nboot]);
            % TODO: consider parfor
            parfor b = 1:nboot
                % draw random sample with replacement
                inds = ceil(rand(n,1)*n);
                bootest(:,:,b) = self(inds).fit;
            end
            % use percentile method to find median and half of 68% range
            % (ie standard error) (from KK's GLMestimatemodel)
            percs = prctile(bootest,[16 50 84],3);
            estimates = percs(:,:,2);
            sterrs = diff(percs(:,:,[1 3]),1,3)/2;
        end
    end
end
