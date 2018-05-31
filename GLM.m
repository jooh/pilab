% General object array class for storing GLM-based models (one entry per
% run). This variant uses OLS, but supports sub-classing for special cases
% (e.g. ridge regression in RidgeGLM).
%
% Initialise with:
% gl = GLM(X,data);
%
classdef GLM < Saveable
    properties
        X % design matrix (nsamples by nregressors)
        data % data (nsamples by nfeatures)
        nsamples % number of observations (rows in X and data)
        nfeatures % number of features (columns in data)
        npredictors % number of regressors in model (columns in X)
        cvgroup % index for crossvalidation (scalar)
        nrandsamp % number of samples for permutation / bootstrapping
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
            % default to number of samples (sub-classes may change this)
            gl.nrandsamp = gl.nsamples;
        end

        function estimates = fit(self)
        % OLS fitted parameter estimates.
        %
        % estimates = fit();
            estimates = olsfit(getdesign(self),getdatac(self));
        end

        function cmat = covmat(self)
        % covariance matrix for the design matrix.
        %
        % cmat = covmat()
            X = getdesign(self);
            cmat = X' * X \ eye(size(X,2));
            % the classic way would use the inverse but this is slower and
            % less precise in Matlab (and perhaps in general)
            %cmat = inv(X' * X);
        end

        function d = df(self)
        % degrees of freedom.
        %
        % d = df()
            d = sum([self.nsamples]) - rank(getdesign(self));
        end

        function r = residuals(self,varargin)
        % residuals of prediction on data. If undefined, we use
        % self-prediction.
        %
        % r = residuals([prediction])
            if nargin > 1
                prediction = vertcat(varargin{:});
            else
                prediction = predictY(self,getdesign(self));
            end
            r = getdata(self) - prediction;
        end

        function mr = mrss(self)
        % mean residual sum of squares (calculated over df).
        %
        % mr = mrss()
            mr = sum(residuals(self).^2,1) / df(self);
        end

        function Yfit = predictY(self,varargin)
        % generate fitted (predicted) responses for a design matrix X.
        %
        % if no inputs are provided we assume you want a self fit.
        %
        % Yfit = predictY([X])
            if nargin > 1                
                X = vertcat(varargin{:});
            else
                X = getdesign(self);
            end
            Yfit = X * fit(self);
        end

        function mse = mserr(self,Yfit)
        % mean squared error between predicted responses in Yfit and
        % current data.
        %
        % if no inputs are provided we assume you want a self fit
        % (self.predictY).
        %
        % mse = mserr([Yfit])
            if nargin == 1 
                Yfit = predictY(self);
            end
            mse = mean((Yfit-getdata(self)).^2);
        end

        function R = rsquare(self,Yfit)
        % R^2 between predicted responses in Yfit and current data.
        %
        % if no inputs are provided we assume you want a self fit
        % (self.predictY).
        %
        % R = rsquare([Yfit])
            if nargin == 1
                Yfit = predictY(self);
            end
            R = rsquare(Yfit,getdata(self));
        end

        function r = pearson(self,Yfit)
        % return pearson's r between the predicted data in Yfit and the
        % actual time series in self.data. If Yfit is undefined we use the
        % prediction from all entries in self.
        %
        % This method is used mainly as an alternative model fit metric to
        % rsquare or mserr for cases where we want to allow for gain and
        % offset parameters to vary between prediction and data.
        %
        % Note that an important limitation of this approach is that r will
        % always be positive (since we're comparing _fitted_ values with
        % data), unless you crossvalidate.
        %
        % r = pearson(self,Yfit)
            if nargin == 1
                % self prediction
                Yfit = predictY(self);
            end
            r = corrpairs(Yfit,getdata(self));
        end

        function con = contrast(self,conmat)
        % compute a contrast estimate for each contrast vector (rows) in
        % conmat.
        %
        % con = contrast(conmat)
            con = conmat * fit(self);
        end

        function se = standarderror(self,conmat)
        % compute the standard error for each contrast vector (rows) in
        % conmat.
        %
        % se = standarderror(conmat)
            se = sqrt(diag((conmat * covmat(self) * conmat')) * ...
                mrss(self));
        end

        function t = tmap(self,conmat)
        % compute a t statistic for each contrast vector (rows) in
        % conmat.
        %
        % t = tmap(conmat)
            t = contrast(self,conmat) ./ standarderror(self,conmat);
        end

        function p = pmap(self,conmat,tail)
        % return p values for each contrast vector (rows) in conmat. Tail
        % can be both (default), left or right.
        %
        % p = pmap(conmat,tail)
            if ieNotDefined('tail')
                tail = 'both';
            end
            % get the t statistics
            t = tmap(self,conmat);
            switch tail
                case 'both'
                    % two-tailed
                    p = 2 * tcdf(-abs(t),df(self));
                case 'left'
                    p = tcdf(t,df(self));
                case 'right'
                    p = tcdf(-t,df(self));
                otherwise
                    error('unknown tail: %s',tail)
            end
            % make sure p does not round to 0 for values beyond current
            % precision
            p(p==0) = realmin(class(p));
        end

        function [winners,res] = crossvalidateproperty(self,property,values,metric,selectfun)
        % tune a property (e.g., RidgeGLM's 'k') over a range of values (1D
        % array) using run-based crossvalidation (cvpredictionrun method)
        % and selectfun metric (default @max rsquare) as the selection
        % criterion. winners gives the winning value for each feature (1 by
        % nfeatures), and res gives the performance for each value
        % across splits (nvalues by nfeatures).
        %
        % This method should be applied to an independent training split of
        % the data to avoid circular inference.
        %
        % [winners,res] = crossvalidateproperty(property,values,[metric],[selectfun])
            if ieNotDefined('metric')
                metric = 'rsquare';
            end
            if ieNotDefined('selectfun')
                selectfun = @max;
            end
            nvalues = length(values);
            res = preallocate(self,[nvalues self(1).nfeatures]);
            for v = 1:nvalues
                % set property in train and test
                [self.(property)] = deal(values(v));
                % cross-validate fit with chosen property
                res(v,:) = cvpredictionrun(self,'predictY',metric);
            end
            % find winning value for each feature
            [~,inds] = selectfun(res,[],1);
            winners = values(inds);
        end

        function splits = preparerunsplits(self)
        % set up a logical array where each row provides logical indices
        % into the available runs (true for test, false for train). Uses
        % the cvgroups field.
        %
        % splits = preparerunsplits(self)
            if all(isempty([self.cvgroup]))
                % use straight leave one out
                splits = logical(eye(numel(self)));
            else
                cvgroups = horzcat(self.cvgroup);
                assert(~any(~isnumeric(cvgroups)),...
                    'cvgroup must be numeric');
                ugroups = unique(cvgroups)';
                splits = cell2mat(arrayfun(@(g)cvgroups==g,ugroups,...
                    'uniformoutput',false));
            end
        end

        function [res,splitres] = cvclassificationrun(self,trainmeth,testmeth,testself,varargin)
        % crossvalidate the performance of some classifier fit with
        % trainmeth (e.g. discriminant) and some testmeth (e.g.
        % infoc). this method is for cases where you want to do out of
        % sample decoding, ie, predict the columns of the design matrix. If
        % you want to predict the data samples, use cvpredictionrun.
        %
        % methargs is any number of extra arguments. These get passed to
        % both trainmeth and testmeth (e.g., a contrast vector).
        %
        % [res,splitres] = cvclassificationrun(self,trainmeth,testmeth,[testself],[methargs])
            if ~exist('testself','var') || isempty(testself)
                testself = self;
            end
            [res,splitres] = cvcrossclassificationrun(self,trainmeth,...
                testmeth,testself,varargin,varargin);
        end

        function [res,splitres] = cvcrossclassificationrun(self,trainmeth,testmeth,testself,trainargs,testargs)
        % crossvalidate the performance of some classifier fit with
        % trainmeth (e.g. discriminant) and some set of parameters (e.g. a
        % contrast vector) and apply testmeth (e.g.  infoc) with a
        % second set of parameters (e.g. a different contrast vector). this
        % method is for cases where you want to do out of sample
        % cross-decoding, ie, train on some columns of the design matrix.
        % and predict another set of columns. If you want to predict the
        % data samples, use cvpredictionrun.
        %
        % trainargs and testargs are cell arrays with any number of extra
        % arguments.
        %
        % res = cvcrossclassificationrun(self,trainmeth,testmeth,[testself],[trainargs],[testargs])
            if ~exist('trainargs','var') || isempty(trainargs)
                trainargs = {};
            end
            if ~iscell(trainargs)
                trainargs = {trainargs};
            end
            if ~exist('testargs','var') || isempty(testargs)
                testargs = {};
            end
            if ~exist('testself','var') || isempty(testself)
                testself = self;
            end
            if ~iscell(testargs)
                testargs = {testargs};
            end
            splits = preparerunsplits(self);
            nsplit = size(splits,1);
            assert(nsplit > 1,'can only crossvalidate if >1 run');
            splitres = cell(nsplit,1);
            for s = 1:nsplit
                prediction = feval(trainmeth,self(~splits(s,:)),...
                    trainargs{:});
                splitres{s} = feval(testmeth,testself(splits(s,:)),...
                    prediction,testargs{:});
            end
            res = matmean(splitres{:});
        end

        function varargout = validatedclassification(self,split,trainmeth,trainargs,testmeth,testargs,testself)
        % one-way classification: train some classifier on a training split
        % of self (self(~split)) using trainmeth with trainargs{:}. Then
        % validate the classifier predictions on a test split (self(split))
        % using testmeth with testargs{:}. Mainly useful for custom
        % splitting schemes where you don't want to crossvalidate (ie use
        % all data both for training and prediction). Note that we support
        % any number of returns (e.g. both t and p for infot)
        %
        % varargout = validatedclassification(self,split,trainmeth,trainargs,testmeth,testargs,[testself])
            if ~exist('testself','var') || isempty(testself)
                testself = self;
            end
            train = self(~split);
            test = testself(split);
            tfit = feval(trainmeth,train,trainargs{:});
            [varargout{1:nargout}] = feval(testmeth,test,tfit,testargs{:});
        end

        function model = drawpermruns(self,runinds)
        % return a model where the design matrix X has been randomly
        % reassigned to different runs according to runinds. This is a
        % useful way to build a null distribution if
        % a) your runs have identical nsamples 
        % b) your design matrices are independent
        %
        % model = drawpermruns(self,runinds)
            X = {self.X};
            model = copy(self);
            for r = 1:numel(self)
                model(r).X = X{runinds(r)};
            end
        end

        function res = cvpredictionrun(self,trainmeth,testmeth,testself)
        % crossvalidate the performance of some prediction trainmeth (e.g.,
        % predictY) using some testmeth (e.g., rsquare). this method is for
        % cases where you want to do out of sample regression, that is,
        % predict the data samples - that is, cases where trainmeth takes
        % the test design matrix as input. If you are doing out of sample
        % classification, that is predicting the columns of the design
        % matrix, use cvclassificationrun. 
        %
        % res = cvpredictionrun(self,trainmeth,testmeth,[testself])
            if ~exist('testself','var') || isempty(testself)
                testself = self;
            end
            splits = preparerunsplits(self);
            nsplit = size(splits,1);
            assert(nsplit > 1,'can only crossvalidate if >1 run');
            prediction = cell(nsplit,1);
            for s = 1:nsplit
                train = self(~splits(s,:));
                test = testself(splits(s,:));
                prediction{s} = feval(trainmeth,train,getdesign(test));
            end
            res = feval(testmeth,testself(testrunind(testself)),...
                getdata(testself,prediction{:}));
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

        function varargout = bootstrapruns(self,nboot,bootmeth,outshape,varargin)
        % general method for computing bootstrap estimates for some
        % bootmeth. Additional varargins get passed on to the bootmeth to
        % enable bootstrapping of e.g. particular contrasts.
        %
        % INPUTS
        % self - GLM or subclass
        % nboot - scalar OR a prepared nboot by nrun matrix (e.g.
        %   from preparerunboots). This is useful to yoke bootstraps
        %   across different analyses.
        % bootmeth - string or cell array defining a set of valid methods
        %   for self (one varargout is returned for each).
        % [outshape] - optional. Save some cycles by pre-defining,
        %   otherwise we hit bootmeth once to infer it.
        % [varargin] - any additional varargin are passed on to bootmeth.
        %
        % You can optionally specify multiple bootmeth in cell array form
        % to get yoked bootstrap estimates for multiple methods.
        %
        % varargout = bootstrapruns(self,nboot,bootmeth,[outshape],[varargin])
            if ieNotDefined('outshape')
                outshape = {};
            end
            if isscalar(nboot)
                bootinds = preparerunboots(self,nboot);
            else
                bootinds = nboot;
                nboot = size(bootinds,1);
            end
            [varargout,bootmeth] = preallocateperms(self,nboot,bootmeth,outshape,varargin{:});
            % update bootinds in case you asked for more than what's
            % possible
            nboot = size(bootinds,1);
            for m = 1:length(varargout)
                % redundant assignment to keep parfor happy
                temp = varargout{m};
                parfor b = 1:nboot
                    temp(:,:,b) = feval(bootmeth{m},self(bootinds(b,:)),...
                        varargin{:});
                end
                varargout{m} = temp;
            end
        end

        function model = drawpermsample(self,inds)
        % return a new instance where the data has been re-ordered
        % according to inds. Note that you must supply the same number of
        % inds as self.nsamples.
        %
        % model = drawpermsample(self,inds)
            assert(numel(inds)==self(1).nsamples,...
                'got %d inds for %d samples',numel(inds),self(1).nsamples);
            model = copy(self);
            for r = 1:length(self)
                model(r).data = model(r).data(inds,:);
            end
        end

        function model = drawpermflipsample(self,inds)
        % return a new instance where the data has been sign flipped
        % according to logical indices in inds. Note that you must supply
        % the same number of inds as self.nsamples.
        %
        % model = drawpermflipsample(self,inds)
            assert(numel(inds)==self(1).nsamples,...
                'got %d inds for %d samples',numel(inds),self(1).nsamples);
            model = copy(self);
            for r = 1:length(self)
                model(r).data(inds,:) = model(r).data(inds,:) * -1;
            end
        end

        function inds = preparesampleperms(self,nperms)
        % inds = preparesampleperms(self,nperms)
            inds = permuteindices(self(1).nrandsamp,nperms);
        end

        function inds = preparerunperms(self,nperms)
        % inds = preparerunperms(self,nperms)
            inds = permuteindices(numel(self),nperms);
        end

        function inds = preparesamplepermflips(self,nperms)
        % inds = preparesamplepermflips(self,nperms)
            inds = permflipindices(self(1).nrandsamp,nperms);
        end

        function varargout = permutesamples(self,nperms,permmeth,outshape,varargin)
        % generate a null distribution of some permmeth (e.g., fit) by
        % resampling the order of the samples in X without replacement.
        % Note that the same random order is inserted into each row.
        %
        % INPUTS:
        % self - GLM or subclass
        % nperms - scalar OR a prepared nperm by nsamples matrix (e.g.
        %   from preparesampleperms). This is useful to yoke permutations
        %   across different analyses.
        % permmeth - string or cell array defining a set of  valid methods
        %   for self (one varargout is returned for each).
        % [outshape] - optional. Save some cycles by pre-defining,
        %   otherwise we hit permmeth once to infer it.
        % [varargin] - any additional varargin are passed on to permmeth.
        %
        % varargout = permutesamples(self,nperms,permmeth,[outshape],[varargin])
            if ieNotDefined('outshape')
                outshape = {};
            end
            if isscalar(nperms)
                perminds = preparesampleperms(self,nperms);
            else
                perminds = nperms;
                nperms = size(perminds,1);
            end
            [varargout,permmeth] = preallocateperms(self,nperms,...
                permmeth,outshape,varargin{:});
            for m = 1:length(varargout)
                temp = varargout{m};
                parfor p = 1:size(perminds,1);
                    permd = drawpermsample(self,perminds(p,:));
                    temp(:,:,p) = feval(permmeth{m},permd,...
                        varargin{:});
                end
                varargout{m} = temp;
            end
        end

        function varargout = permflipsamples(self,nperms,permmeth,outshape,varargin)
        % generate a null distribution of some permmeth (e.g., fit) by
        % flipping the sign of the samples in the data. Note that the same
        % random resample is applied to the data in each run.
        %
        % INPUTS:
        % self - GLM or subclass
        % nperms - scalar OR a prepared nperm by nsamples matrix (e.g.
        %   from preparesampleperms). This is useful to yoke permutations
        %   across different analyses.
        % permmeth - string or cell array defining a set of  valid methods
        %   for self (one varargout is returned for each).
        % [outshape] - optional. Save some cycles by pre-defining,
        %   otherwise we hit permmeth once to infer it.
        % [varargin] - any additional varargin are passed on to permmeth.
        %
        % varargout = permflipsamples(self,nperms,permmeth,[outshape],[varargin])
            if ieNotDefined('outshape')
                outshape = {};
            end
            if isscalar(nperms)
                perminds = preparesamplepermflips(self,nperms);
            else
                perminds = nperms;
                nperms = size(perminds,1);
            end
            [varargout,permmeth] = preallocateperms(self,nperms,...
                permmeth,outshape,varargin{:});
            for m = 1:length(varargout)
                % redundant assigment to keep parfor happy
                temp = varargout{m};
                parfor p = 1:size(perminds,1);
                    permd = drawpermflipsample(self,perminds(p,:));
                    temp(:,:,p) = feval(permmeth{m},permd,...
                        varargin{:});
                end
                varargout{m} = temp;
            end
        end

        function varargout = permuteruns(self,nperms,permmeth,outshape,varargin)
        % generate a null distribution of some permmeth (e.g., fit) by
        % shuffling the assignment of design matrices (X) to runs.
        % Note that this method assumes that a) you used an independently
        % randomised trial sequence in each run, b) each run contains the
        % same number of samples, c) you have already projected out any
        % covariates that you do not want to include in the null
        % distribution.
        %
        % INPUTS:
        % self - GLM or subclass
        % nperms - scalar OR a prepared nperm by nrun matrix (e.g.
        %   from preparerunperms). This is useful to yoke permutations
        %   across different analyses.
        % permmeth - string or cell array defining a set of valid methods
        %   for self (one varargout is returned for each).
        % [outshape] - optional. Save some cycles by pre-defining,
        %   otherwise we hit permmeth once to infer it.
        % [varargin] - any additional varargin are passed on to permmeth.
        %
        % varargout = permuteruns(self,nperms,permmeth,[outshape],[varargin])
            if ieNotDefined('outshape')
                outshape = {};
            end
            if isscalar(nperms)
                perminds = preparerunperms(self,nperms);
            else
                perminds = nperms;
                nperms = size(perminds,1);
            end
            % preflight checks
            assert(numel(self) == size(perminds,2),...
                'runinds does not match nrun');
            assert(isequal(self(1).nsamples,self.nsamples),...
                'nsamples must match across runs');
            [varargout,permmeth] = preallocateperms(self,nperms,...
                permmeth,outshape,varargin{:});
            for m = 1:length(varargout)
                % redundant assigment to keep parfor happy
                temp = varargout{m};
                parfor p = 1:size(perminds,1);
                    permd = drawpermruns(self,perminds(p,:));
                    temp(:,:,p) = feval(permmeth{m},permd,varargin{:});
                end
                varargout{m} = temp;
            end
        end

        function inds = preparerunboots(self,nboot)
        % return unique indices for bootstrapping (wraps bootindices
        % function).
        %
        % inds = preparerunboots(self,nboot)
            inds = bootindices(numel(self),nboot);
        end

        function inds = preparesampleboots(self,nboot)
        % return unique indices for bootstrapping (wraps bootindices
        % function).
        %
        % inds = preparesampleboots(self,nboot)
            inds = bootindices(self(1).nrandsamp,nboot);
        end

        function model = drawbootsample(self,inds)
        % return a new instance where samples in both X and data have been
        % re-ordered according to inds. Note that you must supply the same
        % number of inds as self.nsamples.
        %
        % model = drawbootsample(self,inds)
            assert(numel(inds)==self(1).nsamples,...
                'got %d inds for %d samples',numel(inds),self(1).nsamples);
            model = copy(self);
            for r = 1:length(self)
                model(r).X = model(r).X(inds,:);
                model(r).data = model(r).data(inds,:);
            end
        end

        function varargout = bootstrapsamples(self,nboot,bootmeth,outshape,varargin)
        % bootstrap the estimate for some bootmeth (e.g. fit) by drawing
        % samples from the rows of X and data in concert. Note that the
        % same random draw is made in each run. If outshape is undefined we
        % infer it. The returned bootest is outshape by nboot.
        %
        % INPUTS
        % self - GLM or subclass
        % nboot - scalar OR a prepared nboot by nsamples matrix (e.g.
        %   from preparesampleboots). This is useful to yoke bootstraps
        %   across different analyses.
        % bootmeth - string or cell array defining a set of  valid methods
        %   for self (one varargout is returned for each).
        % [outshape] - optional. Save some cycles by pre-defining,
        %   otherwise we hit bootmeth once to infer it.
        % [varargin] - any additional varargin are passed on to bootmeth.
        %
        % varargout = bootstrapsamples(self,nboot,bootmeth,[outshape],[varargin])
            if ieNotDefined('outshape')
                outshape = {};
            end
            if isscalar(nboot)
                bootinds = self.preparesampleboots(nboot);
            else
                bootinds = nboot;
                nboot = size(bootinds,1);
            end
            % NB, you will get nans at the end of bootest if you ask for
            % more boots than is possible. prctile handles this well.
            [varargout,bootmeth] = preallocateperms(self,nboot,bootmeth,...
                outshape,varargin{:});
            % update bootinds in case you asked for more than what's
            % possible
            nboot = size(bootinds,1);
            for m = 1:length(varargout)
                % redundant assignment to keep parfor happy
                temp = varargout{m};
                parfor b = 1:nboot
                    bootd = drawbootsample(self,bootinds(b,:));
                    if rank(getdesign(bootd)) < bootd(1).npredictors
                        temp(:,:,b) = NaN;
                    else
                        temp(:,:,b) = feval(bootmeth{m},bootd,varargin{:});
                    end
                end
                varargout{m} = temp;
            end
        end

        function varargout = discriminant(self,varargin)
        % fit linear discriminant(s) w for the the contrast matrices in all
        % subsequent inputs. Uses a shrinkage covariance estimator (see
        % covdiag.m).
        %
        % In general, the slow part of this method is calculating the
        % covariance matrix. This is the same for all decoders so entering
        % multiple inputs can speed things up.
        %
        % varargout = discriminant(self,varargin)
            sa = covdiag(residuals(self));
            varargout = cell(1,nargout);
            for n = 1:length(varargout)
                % unit length transform each weight vector to ensure
                % discriminant outputs are in test data units.
                varargout{n} = unitlen(...
                    (contrast(self,varargin{n}) / sa)')';
            end
        end

        function model = infomodel(self,w)
        % return a new GLM instance where the data have been projected onto
        % the discriminant(s) w (fitted with e.g. GLM.discriminant). The
        % resulting model is one concatenated run with one feature
        % per discriminant (size(w,1)). See Kriegeskorte et al., 2007,
        % PNAS. 
        %
        % model = infomodel(self,w)
            X = getdesign(self);
            % project the data onto the discriminant matrix (so from
            % nsamples by nfeatures to nsamples by ndiscriminants)
            dataw = getdata(self) * w';
            % make a new GLM instance (NB, we ignore sub-classing here -
            % assume the original self.data has already been appropriately
            % pre-processed.
            model = GLM(X,dataw);
        end

        function [t,p] = infot(self,w,conmat)
        % return t estimates for contrasts conmat computed on the
        % infomodel(self,w). The resulting t summarises the strength of the
        % pattern information effect for each contrast (ie, the mean
        % distance from the linear discriminant decision boundary
        % normalised by the variance of said distance).
        %
        % if the second output (p) is requested we obtain a parametric p
        % value for the t statistic.
        %
        % [t,p] = infot(self,w,conmat)
            model = infomodel(self,w);
            % diag to get each contrast's estimate on its own discriminant
            % feature.
            cons = diag(contrast(model,conmat));
            errs = diag(standarderror(model,conmat));
            t = cons ./ errs;
            if nargout>1
                % 1-tailed p value
                % extra conversion here because larger single-precision t
                % stats tend to get rounded to 1 (ie, 0).
                p = 1-tcdf(double(t),df(model));
            end
        end

        function [dist,p] = infoc(self,w,conmat)
        % linear discriminant contrast estimate - ie, the distance between
        % the contrasted conditions along the discriminant dimension. 
        %
        % Note that the units of the data are only preserved if the weights
        % vector is set to unit length and the contrast vector's absolute
        % values from the positive and negative elements of each contrast
        % both sum to 1.
        %
        % If a full covariance estimate is used to calculate the linear
        % discriminant and this method is called on the same data that was
        % used to form the discriminant, the returned value is the
        % Mahalanobis distance.
        %
        % if the second output (p) is requested we obtain a parametric p
        % value using infot.
        %
        % INPUTS
        % w         weights vector, see GLM.discriminant
        % conmat    n by npredictors contrast matrix.
        %
        % OUTPUT
        % dist      n by 1 distance estimates
        % p         n by 1 one-tailed p values for T test (save compute by
        %               omitting this if not needed)
        %
        % [dist,p] = infoc(self,w,conmat)
            c = contrast(self,conmat);
            dist = diag(c * w');
            if nargout > 1
                % obtain p value by t test
                [~,p] = infot(self,w,conmat);
            end
        end

        function d = preallocate(self,shape)
        % pre-allocate a matrix of NaNs with the same type as self.data.
        % Normally this wouldn't require a custom method but unfortunately
        % Matlab's NaN function doesn't support gpuArray class inputs at
        % present. 
        %
        % d = preallocate(self,shape)
            if isa(self(1).data,'gpuArray')
                % gpu array is currently a special case
                d = gpuArray.nan(shape,classUnderlying(self(1).data));
            else
                d = NaN(shape,class(self(1).data));
            end
        end

        function [d,permmeth] = preallocateperms(self,nperms,permmeth,outshape,varargin)
        % d = preallocateperms(self,nperms,permmeth,outshape,[varargin])
            if ~iscell(permmeth)
                permmeth = {permmeth};
            end
            if ~iscell(outshape)
                outshape = {outshape};
            end
            nmeth = numel(permmeth);
            if isempty(outshape)
                for m = 1:nmeth
                    outshape{m} = size(feval(permmeth{m},self,varargin{:}));
                    assert(numel(outshape{m})<3,...
                        'permmeth must return 2d outputs, got %s',...
                        mat2str(outshape{m}));
                end
            end
            d = arrayfun(@(m)preallocate(self,...
                [outshape{m} nperms]),1:nmeth,'uniformoutput',false);
        end

        function X = getdesign(self)
        % return the vertically concatenated design matrix for the entered
        % runs.
        %
        % X = getdesign(self)
            X = vertcat(self.X);
        end

        function data = getdata(self,varargin)
        % return the vertically concatenated data for the entered runs.
        %
        % data = getdata(self)
            if nargin>1
                datac = varargin;
            else
                datac = getdatac(self);
            end
            data = vertcat(datac{:});
        end

        function data = getdatac(self)
        % return the data in cell array format.
        %
        % data = getdatac(self)
            data = ascol({self.data});
        end

        function ind = testrunind(self)
        % return the indices that sorts self according to cvgroup. This is
        % used in the run-based CV methods, where a prediction is
        % vertically concatenated across splits and compared against the
        % data. The training split is organised in numerical order (ie, low
        % cvsplits run before higher), so the test data must be sorted the
        % same way.
        %
        % ind = testrunind(self)
            [~,ind] = sort(horzcat(self.cvgroup));
            n = numel(self);
            if isempty(ind)
                ind = 1:n;
            end
            assert(n == numel(ind),'cvgroup cannot be part empty');
        end
    end
end
