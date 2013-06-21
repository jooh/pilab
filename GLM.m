% General object array class for storing GLM-based models (one entry per
% run). This variant uses OLS, but supports sub-classing for special cases
% (e.g. ridge regression in RidgeGLM).
%
% Initialise with:
% gl = GLM(X,data);
%
classdef GLM < matlab.mixin.Copyable
    properties
        X % design matrix (nsamples by nregressors)
        data % data (nsamples by nfeatures)
        nsamples % number of observations (rows in X and data)
        nfeatures % number of features (columns in data)
        npredictors % number of regressors in model (columns in X)
        cvgroup % index for crossvalidation (scalar)
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
        % OLS fitted parameter estimates.
        %
        % estimates = fit(self);
            estimates = olsfit(vertcat(self.X),vertcat(self.data));
        end

        function cmat = covmat(self)
        % covariance matrix for the design.
        %
        % cmat = covmat(self)
            X = vertcat(self.X);
            cmat = inv(X' * X);
        end

        function d = df(self)
        % degrees of freedom.
        %
        % d = df()
            d = sum([self.nsamples]) - rank(vertcat(self.X));
        end

        function r = residuals(self)
        % residuals of (self-) fitted responses on data.
        %
        % r = residuals()
            r = predictY(self,vertcat(self.X)) - vertcat(self.data);
        end

        function mr = mrss(self)
        % return the mean residual sum of squares (calculated over df).
        %
        % mr = mrss()
            resid = residuals(self);
            rss = sum(resid.^2);
            mr = rss / df(self);
        end

        function Yfit = predictY(self,varargin)
        % generate fitted (predicted) responses for a design matrix X.
        % Uses varargin to support combining multiple design matrix inputs
        % at once.
        %
        % if no inputs are provided we assume you want a self fit.
        %
        % Yfit = predictY(X,[X],[X etc...])
            if nargin > 1                
                X = vertcat(varargin{:});
            else
                X = vertcat(self.X);
            end
            Yfit = X * fit(self);
        end

        function mse = mserr(self,Yfit)
        % mean squared error between predicted responses in Yfit and
        % current data.
        %
        % if no inputs are provided we assume you want a self fit.
        %
        % mse = mserr(Yfit)
            if nargin == 1 
                Yfit = predictY(self);
            end
            mse = mean((Yfit-vertcat(self.data)).^2);
        end

        function R = rsquare(self,Yfit)
        % R^2 between predicted responses in Yfit and current data.
        %
        % if no inputs are provided we assume you want a self fit.
        %
        % R = rsquare(Yfit)
            if nargin == 1
                Yfit = predictY(self);
            end
            R = rsquare(Yfit,vertcat(self.data));
        end

        function con = contrast(self,cv)
        % compute a contrast for self fit based on contrast vector cv.
        %
        % con = contrast(cv)
            assert(ndims(cv)==2 && all(size(cv)==[1 self(1).npredictors]),...
                'contrast vector must be 1 by npredictors');
            con = cv * fit(self);
        end

        function t = tmap(self,cv)
        % compute a T contrast for self fit based on contrast vector cv.
        %
        % t = tmap(cv)
            t = contrast(self,cv) ./ sqrt(mrss(self) * ...
                (cv * covmat(self) * cv'));
        end

        function p = pmap(self,cv)
        % return one-tailed p values for self fit based on contrast vector
        % cv.
        %
        % p = pmap(cv)
            p = 1-tcdf(abs(tmap(self,cv)),df(self));
        end

        function [winners,meds] = crossvalidateparameter(self,parameter,values,metric)
        % tune a parameter (char) over a range of values (1D array) using a
        % leave-one-out crossvalidation method with maximum median R2 or
        % min mserr (default R2) as the selection criterion (set with
        % 'metric': 'r2' or 'mse'). winners gives the winning value for
        % each feature (1 by nfeatures), and med gives the median R2/mserr
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
                splitres = crossvalidatefit(self,metric);
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

        function splitres = crossvalidatefit(self,metric)
        % return independent estimates of metric (r2 or mse) by
        % leave-one-out crossvalidation. the output splitres contains one
        % row for each split of the data.
        %
        % splitres = crossvalidatefit(metric)
            if ieNotDefined('metric')
                metric = 'r2';
            end
            if all(isempty([self.cvgroup]))
                % use straight leave one out
                splits = logical(eye(numel(self)));
            else
                % use cvgroup indices to leave one _group_ out
                cvgroups = horzcat(self.cvgroup);
                assert(~any(isempty(cvgroups) || ~isnumeric(cvgroups)),...
                    'cvgroup must be defined for all runs');
                groups = unique(cvgroups)';
                splits = cell2mat(arrayfun(@(g)cvgroups==g,groups,...
                    'uniformoutput',false));
            end
            nsplits = size(splits,1);
            assert(nsplits > 1,'can only crossvalidate if >1 run');
            nfeatures = self(1).nfeatures;
            splitres = NaN([nsplits nfeatures]);
            parfor s = 1:nsplits
                test = splits(s,:);
                train = ~splits(s,:);
                switch metric
                    case 'r2'
                        splitres(s,:) = rsquare(self(test),...
                            predictY(self(train),self(test).X));
                    case 'mse'
                        splitres(s,:) = mserr(self(test),...
                            predictY(self(train),self(test).X));
                    otherwise
                        error('unknown metric: %s',metric);
                end
            end
        end

        function [estimates,sterrs,bootest] = bootstraprunfit(self,nboot)
        % draw nboot samples with replacement from the available runs in
        % self. Return the median and standard error of the model fit. Uses
        % the bootstrapruns method internally - see GLM.bootstrapruns for
        % details.
        %
        % Return the median of the bootstrap estimates, half of the
        % 68% range (ie, an estimate of standard error), and the raw
        % distribution of bootstraps.
        % [estimates,sterrs,bootest] = bootstrapfit(self,nboot)
            bootest = bootstrapruns(self,nboot,'fit',...
                [self(1).npredictors,self(1).nfeatures]);
            % get percentile
            [estimates,sterrs] = bootprctile(bootest);
        end

        function bootest = bootstrapruns(self,nboot,bootmeth,outshape)
        % general method for computing bootstrap estimates for some
        % bootmeth (char - must refer to a method of the current instance
        % which takes no input arguments and returns at least one output).
        % outshape defines the dimensionality of the bootest
        % matrix ([row column]).
        %
        % we draw nboot samples with replacement from the available runs in
        % self.
        %
        % if you leave outshape undefined we hit bootmeth once to figure
        % out what the desired shape is. This convenience comes at a
        % performance cost.
        %
        % bootest = bootstrapruns(self,nboot,bootmeth,[outshape])
            n = numel(self);
            if ieNotDefined('outshape')
                % this is a little ugly but we can simply run bootmeth once
                % to see what the output looks like
                outshape = size(self.(bootmeth));
                assert(numel(outshape)<3,...
                    'bootmeth must return at most 2d outputs, got %s',...
                    mat2str(outshape));
            end
            bootest = NaN([outshape nboot]);
            parfor b = 1:nboot
                %fprintf('%d ',b)
                % draw random sample with replacement
                inds = ceil(rand(n,1)*n);
                bootest(:,:,b) = self(inds).(bootmeth);
            end
        end
    end
end
