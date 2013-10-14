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
            estimates = olsfit(vertcat(self.X),vertcat(self.data));
        end

        function cmat = covmat(self)
        % covariance matrix for the design matrix.
        %
        % cmat = covmat()
            X = vertcat(self.X);
            cmat = inv(X' * X);
        end

        function d = df(self)
        % degrees of freedom.
        %
        % d = df()
            d = sum([self.nsamples]) - rank(vertcat(self.X));
        end

        function r = residuals(self,varargin)
        % residuals of prediction on data. If undefined, we use
        % self-prediction.
        %
        % r = residuals([prediction])
            if nargin > 1
                prediction = vertcat(varargin{:});
            else
                prediction = predictY(self,vertcat(self.X));
            end
            r = vertcat(self.data) - prediction;
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
                X = vertcat(self.X);
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
            mse = mean((Yfit-vertcat(self.data)).^2);
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
            R = rsquare(Yfit,vertcat(self.data));
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

        function f = fmap(self,conmat)
        % compute f statistic for the contrast vector (or matrix) in
        % conmat. Uses various SPM8 functions internally (mostly lifted
        % from spm_ancova)
        %
        % f = fmap(self,conmat)
            est = fit(self);
            ress2 = sum(residuals(self).^2,1);
            % this sets up some mysterious structures upon which spm
            % depends
            xX = spm_sp('Set',self.X);
            xCon = spm_FcUtil('Set','','F','c',conmat',xX);
            % good luck unpacking what's going on in the following.
            % Hsqr - extra sum of squares of est from contrast (numerator
            % in F test)
            h = spm_FcUtil('Hsqr',xCon,xX);
            X1o = spm_FcUtil('X1o',xCon,xX);
            V = speye(self.nsamples);
            [trRV, trRVRV] = spm_SpUtil('trRV',xX,V);
            [trMV, trMVMV] = spm_SpUtil('trMV',X1o,V);
            f = sum((h*est).^2,1)./(ress2*trMV/trRV);
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
        end

        function [winners,meds] = crossvalidateproperty(self,property,values,metric,selectfun)
        % tune a property (e.g., RidgeGLM's 'k') over a range of values (1D
        % array) using run-based crossvalidation (cvpredictionrun method)
        % and selectfun metric (default @max rsquare) as the selection
        % criterion. winners gives the winning value for each feature (1 by
        % nfeatures), and med gives the median performance for each value
        % across splits (nvalues by nfeatures).
        %
        % This method should be applied to an independent training split of
        % the data to avoid circular inference.
        %
        % [winners,meds] = crossvalidateproperty(property,values,[metric],[selectfun])
            if ieNotDefined('metric')
                metric = 'rsquare';
            end
            if ieNotDefined('selectfun')
                selectfun = @max;
            end
            nvalues = length(values);
            meds = preallocate(self,[nvalues self(1).nfeatures]);
            for v = 1:nvalues
                % set property in train and test
                [self.(property)] = deal(values(v));
                % cross-validate fit with chosen property
                splitres = cvpredictionrun(self,'predictY','rsquare',...
                    [1 self(1).nfeatures]);
                % update with median across splits
                meds(v,:) = median(splitres,3);
            end
            % find winning value for each feature
            [~,inds] = selectfun(meds,[],1);
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

        function cvres = cvclassificationrun(self,trainmeth,testmeth,outshape,varargin)
        % crossvalidate the performance of some classifier fit with
        % trainmeth (e.g. discriminant) and some testmeth (e.g.
        % infotmap). this method is for cases where you want to do out of
        % sample decoding, ie, predict the columns of the design matrix. If
        % you want to predict the data samples, use cvpredictionrun.
        %
        % methargs is any number of extra arguments. These get passed to
        % both trainmeth and testmeth (e.g., a contrast vector).
        %
        % cvres = cvclassificationrun(self,trainmeth,testmeth,outshape,[methargs])
            cvres = cvcrossclassificationrun(self,trainmeth,testmeth,...
                outshape,varargin,varargin);
        end

        function cvres = cvcrossclassificationrun(self,trainmeth,testmeth,outshape,trainargs,testargs)
        % crossvalidate the performance of some classifier fit with
        % trainmeth (e.g. discriminant) and some set of parameters (e.g. a
        % contrast vector) and apply testmeth (e.g.  infotmap) with a
        % second set of parameters (e.g. a different contrast vector). this
        % method is for cases where you want to do out of sample
        % cross-decoding, ie, train on some columns of the design matrix.
        % and predict another set of columns. If you want to predict the
        % data samples, use cvpredictionrun.
        %
        % trainargs and testargs are cell arrays with any number of extra
        % arguments.         %
        %
        % cvres = cvcrossclassificationrun(self,trainmeth,testmeth,outshape,[trainargs],[testargs])
            if ieNotDefined('trainargs')
                trainargs = {};
            end
            if ieNotDefined('testargs')
                testargs = {};
            end
            if ieNotDefined('outshape')
                % do a self train/test to infer size of output
                temp = self.(trainmeth)(trainargs{:});
                outshape = size(self.(testmeth)(temp,testargs{:}));
            end
            assert(numel(outshape)<3,...
                'testmeth must return at most 2d outputs, got %s',...
                mat2str(outshape));
            splits = preparerunsplits(self);
            % TODO - if train/test are different we effectively have
            % 2*nsplit combinations to run
            nsplit = size(splits,1);
            assert(nsplit > 1,'can only crossvalidate if >1 run');
            cvres = preallocate(self,[outshape nsplit]);
            for s = 1:nsplit
                train = self(~splits(s,:));
                test = self(splits(s,:));
                tfit = train.(trainmeth)(trainargs{:});
                cvres(:,:,s) = test.(testmeth)(tfit,testargs{:});
            end
        end

        function model = drawpermrun(self,runinds)
        % return a model where the design matrix X has been randomly
        % reassigned to different runs according to runinds. This is a
        % useful way to build a null distribution if
        % a) your runs have identical nsamples 
        % b) your design matrices are independent
        %
        % model = permuterun(self,runinds)
            X = {self.X};
            model = self.copy;
            nrun = length(self);
            assert(nrun == numel(runinds),'runinds does not match nrun');
            for r = 1:length(model)
                model(r).X = X{runinds(r)};
            end
        end

        function cvres = cvpredictionrun(self,trainmeth,testmeth,outshape)
        % crossvalidate the performance of some prediction trainmeth (e.g.,
        % predictY) using some testmeth (e.g., rsquare). this method is for
        % cases where you want to do out of sample regression, that is,
        % predict the data samples - that is, cases where trainmeth takes
        % the test design matrix as input. If you are doing out of sample
        % classification, that is predicting the columns of the design
        % matrix, use cvclassificationrun. 
        %
        % cvres = cvpredictionrun(self,trainmeth,testmeth,outshape)
            if ieNotDefined('outshape')
                % do a self test to infer size of output
                prediction = self.(trainmeth)(self.X);
                outshape = size(self.(testmeth)(prediction));
            end
            assert(numel(outshape)<3,...
                'testmeth must return at most 2d outputs, got %s',...
                mat2str(outshape));
            splits = preparerunsplits(self);
            nsplit = size(splits,1);
            assert(nsplit > 1,'can only crossvalidate if >1 run');
            cvres = preallocate(self,[outshape nsplit]);
            for s = 1:nsplit
                train = self(~splits(s,:));
                test = self(splits(s,:));
                prediction = train.(trainmeth)(test.X);
                cvres(:,:,s) = test.(testmeth)(prediction);
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
            if ieNotDefined('outshape')
                % this is a little ugly but we can simply run bootmeth once
                % to see what the output looks like
                outshape = size(self.(bootmeth));
            end
            assert(numel(outshape)<3,...
                'bootmeth must return at most 2d outputs, got %s',...
                mat2str(outshape));
            bootinds = bootindices(numel(self),nboot);
            % update bootinds in case you asked for more than what's
            % possible
            nboot = size(bootinds,1);
            bootest = preallocate(self,[outshape nboot]);
            parfor b = 1:nboot
                bootest(:,:,b) = self(bootinds(b,:)).(bootmeth);
            end
        end

        function model = drawpermsample(self,inds)
        % return a new instance where the samples in X have been re-ordered
        % according to inds. Note that you must supply the same number of
        % inds as self.nsamples.
        %
        % model = drawpermsample(self,inds)
            assert(numel(inds)==self(1).nsamples,...
                'got %d inds for %d samples',numel(inds),self(1).nsamples);
            model = self.copy;
            for r = 1:length(self)
                model(r).X = model(r).X(inds,:);
            end
        end

        function inds = preparesampleperms(self,nperms)
            inds = permuteindices(self(1).nrandsamp,nperms);
        end

        function nulldist = permutesamples(self,nperms,permmeth,outshape)
        % generate a null distribution of some permmeth (e.g., fit) by
        % resampling the order of the samples in X without replacement.
        % Note that the same random resample is applied to the X in each
        % run. If outshape is undefined we infer it. The returned nulldist
        % is outshape by nperms.
        %
        % nulldist = permutesamples(self,nperms,permmeth,outshape)
            if ieNotDefined('outshape')
                outshape = size(self.(permmeth));
            end
            assert(numel(outshape)<3,...
                'permmeth must return at most 2d outputs, got %s',...
                mat2str(outshape));
            nulldist = preallocate(self,[outshape nperms]);
            perminds = self.preparesampleperms(nperms);
            for p = 1:size(perminds,1);
                permd = self.drawpermsample(perminds(p,:));
                nulldist(:,:,p) = permd.(permmeth);
            end
        end

        function inds = preparesampleboots(self,nboot)
        % return unique indices for bootstrapping (wraps bootindices
        % function).
        %
        % inds = preparesampleboots(self,nboot)
            inds = bootindices(self.nrandsamp,nboot);
        end

        function model = drawbootsample(self,inds)
        % return a new instance where samples in both X and data have been
        % re-ordered according to inds. Note that you must supply the same
        % number of inds as self.nsamples.
        %
        % model = drawbootsample(self,inds)
            assert(numel(inds)==self(1).nsamples,...
                'got %d inds for %d samples',numel(inds),self(1).nsamples);
            model = self.copy;
            for r = 1:length(self)
                model(r).X = model(r).X(inds,:);
                model(r).data = model(r).data(inds,:);
            end
        end

        function bootest = bootstrapsamples(self,nboot,bootmeth,outshape)
        % bootstrap the estimate for some bootmeth (e.g. fit) by drawing
        % samples from the rows of X and data in concert. Note that the
        % same random draw is made to each run. If outshape is undefined we
        % infer it. The returned bootest is outshape by nboot.
        %
        % bootest = bootstrapsamples(self,nboot,bootmeth,outshape)
            if ieNotDefined('outshape')
                % this is a little ugly but we can simply run bootmeth once
                % to see what the output looks like
                outshape = size(self.(bootmeth));
            end
            assert(numel(outshape)<3,...
                'bootmeth must return at most 2d outputs, got %s',...
                mat2str(outshape));
            bootinds = self.preparesampleboots(nboot);
            % update bootinds in case you asked for more than what's
            % possible
            nboot = size(bootinds,1);
            bootest = preallocate(self,[outshape nboot]);
            parfor b = 1:nboot
                bootd = self.drawbootsample(bootinds(b,:));
                bootest(:,:,b) = bootd.(bootmeth);
            end
        end

        function w = discriminant(self,conmat)
        % fit linear discriminant(s) w for the contrast(s) in conmat. Uses
        % a shrinkage covariance estimator (see covdiag.m).
        %
        % w = discriminant(self,conmat)
            sa = covdiag(residuals(self));
            w = contrast(self,conmat) / sa;
        end

        function model = infomodel(self,w)
        % return a new GLM instance where the data have been projected onto
        % the discriminant(s) w (fitted with e.g. GLM.discriminant). The
        % resulting model is one concatenated run with one feature
        % per discriminant (size(w,1)). See Kriegeskorte et al., 2007,
        % PNAS. 
        %
        % model = infomodel(self,w)
            X = vertcat(self.X);
            % project the data onto the discriminant matrix (so from
            % nsamples by nfeatures to nsamples by ndiscriminants)
            dataw = vertcat(self.data) * w';
            % make a new GLM instance (NB, we ignore sub-classing here -
            % assume the original self.data has already been appropriately
            % pre-processed, e.g. as part of constructing a CovGLM)
            model = GLM(X,dataw);
        end

        function [t,model] = infotmap(self,w,conmat)
        % return t estimates for contrasts conmat computed on the
        % infomodel(self,w). The resulting t summarises the strength of the
        % pattern information effect for each contrast.
            model = infomodel(self,w);
            % diag to get each contrast's estimate on its own discriminant
            % feature.
            cons = diag(contrast(model,conmat));
            errs = diag(standarderror(model,conmat));
            t = cons ./ errs;
        end

        function mahdist = infomahalanobis(self,w,conmat)
        % return the mahalonobis distance between the points in conmat
        % computed on the infomodel(self,w).
            model = infomodel(self,w);
            % the approach taken here is to take the square root of each
            % contrast estimate. If you estimated w on the same data this
            % is equivalent to the mahalanobis distance.
            mahdist = sqrt(diag(contrast(model,conmat)))';
        end

        function d = preallocate(self,shape)
        % pre-allocate a matrix of NaNs with the same type as self.data.
        % Normally this wouldn't require a custom method but unfortunately
        % Matlab's NaN function is doesn't support gpuArray class inputs at
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
    end
end
