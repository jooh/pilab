% GLM-derived base class for RSA.
% gl = RSA(modelrdms,datardms)
classdef RSA < GLM
    properties
        ncon % number of conditions (rows/columns in RDM)
        Xrdm % original model rdms before de-NaN and transforms
        datardm % original data rdms
        validvec % logical indices for mapping de-NaN'ed sample vector to full length
        issplitdata = false; % split data RDMs invalidate some tests
    end

    methods
        function gl = RSA(modelrdms,datardms)
            if nargin == 0 || (isempty(modelrdms) && isempty(datardms))
                Xvec = [];
                datavec = [];
            else
                % store data in matrix form (for e.g. bootstrapping)
                Xrdm = asrdmmat(modelrdms);
                issplitdata = issplitdatardm(Xrdm);
                datardm = asrdmmat(datardms);
                [ncon,~,npredictors] = size(Xrdm);
                % should start out the same size for sanity
                assert(ncon == size(datardm,1),['data rdm does not '...
                    'match model rdm']);
                % handle NaN content
                nanx = isnan(Xrdm);
                % first drop any conditions that are all NaNs (ie, all
                % row/columns). For this, it's convenient to set the diagonal
                % to true as well
                nanx(repmat(logical(eye(ncon)),[1 1 size(Xrdm,3)])) = true;
                goodcon = arrayfun(@(x)~all(nanx(:,:,x),2),...
                    1:npredictors,'uniformoutput',false);
                % nan conditions must be consistent across predictors
                assert(npredictors==1 || isequal(goodcon{:}),...
                    'inconsistent nan rows across predictors');
                % reduce to valid cons
                goodcon = goodcon{1};
                Xrdm = Xrdm(goodcon,goodcon,:);
                % update ncon 
                ncon = sum(goodcon);
                datardm = datardm(goodcon,goodcon,:);
                % convert to vector and strip NaNs
                Xvec = rdm2vec(Xrdm);
                datavec = rdm2vec(datardm);
                nanmask = any(isnan(Xvec),2) | any(isnan(datavec),2);
                Xvec(nanmask,:) = [];
                datavec(nanmask,:) = [];
                assert(~isempty(Xvec),'no valid (non-NaN) data entered');
            end
            % use super-class to initialise            
            gl = gl@GLM(Xvec,datavec);
            if ~isempty(datavec)
                % note that we set nrandsamp to be ncon because randomisation
                % of the samples is in fact randomisation on the RDM conditions
                % (rows/columns) in RSA.
                [gl.ncon,gl.Xrdm,gl.datardm,gl.validvec,gl.nrandsamp,...
                    gl.issplitdata] = deal(ncon,Xrdm,datardm,~nanmask,...
                    ncon,issplitdata);
            end
        end

        function r = residuals(self,varargin)
        % residuals of prediction on data. If undefined, we use
        % self-prediction. Because these are RDMs, we average and obtain a
        % single residual RDM (see also predictY).
        %
        % r = residuals([prediction])
            if nargin > 1
                prediction = matmean(varargin{:});
            else
                prediction = predictY(self,getdesign(self));
            end
            r = getdata(self) - prediction;
        end

        function fullr = withnans(self,r)
            fullr = NaN([size(self(1).validvec,1) size(r,2)],class(r));
            fullr(repmat(self(1).validvec,[1 size(r,2)])) = r;
        end


        function fullr = residualsfull(self,varargin)
        % The RSA class always removes NaNs before fitting, so the standard
        % residuals function may not return a valid rdvec. This function
        % re-introduces the NaNs to so that conversion to rdmat works
        % again.
        %
        % r = residualsfull(self,varargin)
            r = residuals(self,varargin{:});
            fullr = withnans(self,r);
        end

        function fully = predictYfull(self,varargin)
            y = predictY(self,varargin{:});
            fully = withnans(self,y);
        end

        function fulldata = getdatafull(self,varargin)
            data = getdata(self,varargin{:});
            fulldata = withnans(self,data);
        end

        function tau = kendall_a(self,Yfit)
        % return kendall's tau alpha between the dissimilarity vector in
        % Yfit and the mean RDM across the splits in self.data. If Yfit is
        % undefined we use the mean prediction from all entries in self.
        %
        % This method is used mainly as an alternative model fit metric to
        % rsquare or mserr for cases where we want to test whether the
        % model gets the rank order of the dissimilarities right.
        %
        % tau = kendall_a(self,Yfit)
            if nargin == 1
                % self prediction
                Yfit = predictY(self);
            end
            thisdata = getdata(self);
            tau = NaN([1 self(1).nfeatures],class(self(1).data));
            % as with PluginRSA.fit, this isn't an optimal place for
            % parfor but analyses of large datasets (e.g. searchlight
            % volumes) will be too slow if we leave this out.
            parfor f = 1:self(1).nfeatures
                tau(1,f) = kendall_a(Yfit(:,f),thisdata(:,f));
            end
        end

        function rho = spearman(self,Yfit)
        % return spearman's rho between the dissimilarity vector in Yfit
        % and the mean RDM across the splits in self.data. If Yfit is
        % undefined we use the mean prediction from all entries in self.
        %
        % This method is used mainly as an alternative model fit metric to
        % rsquare or mserr for cases where we want to test whether the
        % model gets the rank order of the dissimilarities right.
        %
        % rho = spearman(self,Yfit)
            if nargin == 1
                % self prediction
                Yfit = predictY(self);
            end
            % average the dissimilarities from separate splits
            % before comparing to the prediction
            rho = corrpairs(ranktrans(Yfit),ranktrans(getdata(self)));
        end

        function r = pearson(self,Yfit)
        % return pearson's r between the dissimilarity vector in Yfit
        % and the mean RDM across the splits in self.data. If Yfit is
        % undefined we use the mean prediction from all entries in self.
        %
        % This method is used mainly as an alternative model fit metric to
        % rsquare or mserr for cases where we want to test whether the
        % model gets the dissimilarity structure right even if the offset
        % and scale of the dissimilarities is wrong.
        %
        % r = pearson(self,Yfit)
            if nargin == 1
                % self prediction
                Yfit = predictY(self);
            end
            % average the dissimilarities from separate splits
            % before comparing to the prediction
            r = corrpairs(Yfit,getdata(self));
        end

        function cloneargs(self,oldclass)
            % does nothing for base RSA case. (in sub-classes this method
            % is used to insure that subclass properties get set properly
            % during resampling.
            %
            % cloneargs(self,oldclass)
        end

        function permglm = drawpermsample(self,inds)
        % return a new instance where the conditions in the Xrdm have been
        % re-ordered according to inds. Note that you must supply the same
        % number of inds as self.ncon. Overrides GLM base class behaviour.
        %
        % model = drawpermsample(self,inds)
            assert(~any([self.issplitdata]),['sample permutation test ' ...
                'is invalid for split data model RDMs. Use a SplitRSA-' ...
                'derived class instead']);
            permX = self(1).Xrdm(inds,inds,:);
            permglm = copy(self);
            for r = 1:length(permglm)
                permglm(r) = feval(class(self),permX,self(r).datardm);
                cloneargs(permglm(r),self(r));
            end
        end

        function [bootglm,removedprop] = drawbootsample(self,inds)
        % return a new instance where the conditions in the Xrdm and
        % datardms have been resampled with replacement in concert.
        % Overrides GLM base class behaviour.
        %
        % Unlike sampling without replacement, RDM bootstrapping can shift
        % diagonal dissimilarities to off-diagonal positions.  These should
        % be removed to avoid underestimating the variance (since data and
        % design will have matching zeros by definition). removedprop is a
        % diagnostic output of the proportion of removed dissimilarities.
        %
        % [bootglm,removedprop] = drawbootsample(self,inds)
            assert(~any([self.issplitdata]),['sample bootstrap ' ...
                'is invalid for split data model RDMs. Use a SplitRSA-' ...
                'derived class instead']);
            assert(all(all(vertcat(self.data)~=0)),...
                'cannot sample when zeros are present in data');
            % process each run in the new instance
            for r = 1:numel(self)
                datasample = resamplerdm(self(r).datardm,inds);
                modelsample = resamplerdm(self(r).Xrdm,inds);
                % create a new instance of whatever class the current
                % instance is
                bootglm(r) = feval(class(self),modelsample,datasample);
                % pass any custom properties from the current instance to
                % the boot instance (e.g., the RidgeRSA k parameter)
                cloneargs(bootglm(r),self(r));
            end
            removedprop = (self(1).nsamples - bootglm(1).nsamples) / ...
                self(1).nsamples;
            assert(removedprop < 1,['bootsample removed everything. ' ...
                'only one unique inds?']);
        end

        function model = selectconditions(self,cons)
        % return a new instance only containing the samples corresponding
        % to the condition indices cons (numerical or logical). note that
        % cons is not sorted so that the condition order can be permuted.
        %
        % model = selectconditions(self,cons)
            % just NaN out the correct parts of the Xrdm
            X = NaN(size(self(1).Xrdm));
            X(cons,:,:) = self(1).Xrdm(cons,:,:);
            X(:,cons,:) = self(1).Xrdm(:,cons,:);
            X(diagind(size(X))) = 0;
            for r = 1:numel(self)
                model(r) = feval(class(self),X,self(r).datardm);
                cloneargs(model(r),self(r));
            end
        end

        function [medr2,meanpredict,r2bycon,predictbycon] = cvpredictionconditions(self)
        % predict the dissimilarities for one condition based on fitted
        % dissimilarities for all other conditions - ie,
        % leave-one-condition out rather than leave-one-run out (see
        % cvpredictionrun).
        %
        % Note that although each train/test split is independent (unlike
        % leave-one-dissimilarity out), the test estimates are dependent
        % across folds. Use e.g. bootstrapping to account for
        % this.
        %
        % [medr2,meanpredict,r2bycon,predictbycon] = cvpredictionconditions(self)
            % store prediction performance for each condition
            r2bycon = NaN([self(1).ncon self(1).nfeatures]);
            % store the fitted response for each condition. Because the X
            % is identical across splits in RSA the prediction is the same
            % across splits. So only need to store one. (need 2d matrix
            % because the predictions are overlapping).
            predictbycon = zeros([self(1).ncon self(1).nsamples ...
                self(1).nfeatures]);
            % vectorised data matrix (one entry per iteration, later
            % collapsed)
            allinds = 1:self(1).ncon;
            for con = allinds
                % use test / ~test to explicitly divide train/test
                % dissimilarities
                test = allinds==con;
                % split data by creating separate train and test model
                % instances
                trainmodel = drawbootsample(self,~test);
                testmodel = selectconditions(self,test);
                % predict dissimilarities for test based on training fit
                prediction = trainmodel.predictY(testmodel.X);
                % and get r2 for model prediction
                r2bycon(con,:) = testmodel.rsquare(prediction);
                % because we are working with RDMs, each X in testmodel is
                % identical. So although the prediction spans all splits we
                % just need to store the first nsamples
                predictbycon(con,:,:) = datavec2mat(...
                    prediction(1:testmodel(1).nsamples,:),repmat(...
                    testmodel(1).validvec,[1 testmodel(1).nfeatures]),...
                    'NaN');
            end
            % mean of all the predictions for each dissimilarity (sort of
            % like making a weave - every dissimilarity features in 2
            % predictions. Once when you hit the right row, once when you
            % hit the right column).
            meanpredict = squeeze(nansum(predictbycon,1)) / 2;
            medr2 = median(r2bycon,1);
        end

        function Yfit = predictY(self,varargin)
        % generate a fitted (predicted) rdm vector for a design matrix X.
        %
        % if no inputs are provided we assume you want a self fit.
        %
        % unlike the super-class GLM behaviour where the prediction will
        % span all the (concatenated) design matrices across runs, here we
        % generate a prediction for a single run (ie one RDM).
        %
        % Yfit = predictY([X])
            if nargin > 1
                Xcell = varargin;
                assert(isequal(Xcell{1},Xcell{:}),['mismatched input ' ...
                    'model RDMs across runs']);
            else
                Xcell = {getdesign(self)};
            end
            % we ignore all subsequent varargin since we assume they
            % are all the same
            Yfit = predictY@GLM(self,Xcell{1});
        end

        function mr = mrss(self)
        % placeholder sub-class of GLM method - triggers an error since you
        % almost certainly don't want to be calculating this on dependent
        % data.
        %
        % mr = mrss(self)
            mr = [];
            error('RSA:noParametricStats',['MRSS and other parametric ' ...
                'inferential stats are not well defined for ' ...
                'RSA since dissimilarities are dependent']);
        end

        function X = getdesign(self)
        % get the design matrix for the first run, and test that all others
        % are the same.
        %
        % X = getdesign(self)
            X = getdesign@GLM(self(1));
            assert(isequal(self(1).X,self.X),['mismatched model RDMs ' ...
                'across runs']);
        end

        function data = getdatac(self)
            % return a cell wrapped first run
            data = {getdata(self)};
        end

        function data = getdata(self,varargin)
        % get data by averaging the RDMs across instances.
        %
        % data = getdata(self)
            if nargin>1
                indata = varargin;
            else
                indata = {self.data};
            end
            data = matmean(indata{:});
        end
    end
end
