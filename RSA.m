% GLM-derived base class for RSA.
% gl = RSA(modelrdms,datardms)
classdef RSA < GLM
    properties
        ncon % number of conditions (rows/columns in RDM)
        Xrdm % original model rdms before de-NaN and transforms
        datardm % original data rdms
        validvec % logical indices for mapping de-NaN'ed sample vector to full length
    end

    methods
        function gl = RSA(modelrdms,datardms)
            if nargin == 0
                Xvec = [];
                datavec = [];
            else
                % store data in matrix form (for e.g. bootstrapping)
                Xrdm = asrdmmat(modelrdms);
                datardm = asrdmmat(datardms);
                [ncon,ncon,npredictors] = size(Xrdm);
                % handle NaN content
                nanx = isnan(Xrdm);
                % first drop any conditions that are all NaNs (ie, all
                % row/columns). For this, it's convenient to set the diagonal
                % to true as well
                nanx(repmat(logical(eye),[1 1 size(Xrdm,3)])) = true;
                goodcon = arrayfun(@(x)~all(nanx(:,:,x),2),...
                    1:npredictors,'uniformoutput',false);
                % nan conditions must be consistent across predictors
                assert(npredictors==1 || isequal(goodcon{:}),...
                    'inconsistent nan rows across predictors');
                % reduce to valid cons
                Xrdm = Xrdm(goodcon{1},goodcon{1},:);
                datardm = datardm(goodcon{1},goodcon{1},:);
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
            if nargin > 0
                % note that we set nrandsamp to be ncon because randomisation
                % of the samples is in fact randomisation on the RDM conditions
                % (rows/columns) in RSA.
                [gl.ncon,gl.Xrdm,gl.datardm,gl.validvec,gl.nrandsamp] = ...
                    deal(ncon,Xrdm,datardm,~nanmask,ncon);
            end
        end

        function cloneargs(self,oldclass)
            % does nothing for base RSA case. (in sub-classes this method
            % is used to insure that subclass properties get set properly
            % during resampling.
        end

        function permglm = drawpermsample(self,inds)
        % return a new instance where the conditions in the Xrdm have been
        % re-ordered according to inds. Note that you must supply the same
        % number of inds as self.ncon. Overrides GLM base class behaviour.
        %
        % model = drawpermsample(self,inds)
            permX = self(1).Xrdm(inds,inds,:);
            permglm = self.copy;
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
            assert(all(all(vertcat(self.data)~=0)),...
                'cannot sample when zeros are present in data');
            % this approach assumes that zeros only happen in the data
            % dissimilarities when diagonals have been resampled to
            % off-diagonal positions
            % find zeros in the boot sampled rdm mat - ie, diagonals that
            % moved.
            matmask = self(1).datardm(inds,inds,1)==0;
            % skip true diagonals
            matmask(logical(eye(self(1).ncon))) = 0;
            % pre-allocate to biggest needed size so we can later just
            % slice for speed
            matmask = repmat(matmask,[1 1 max([self.nfeatures ...
                self.npredictors])]);
            % process each run in the new instance
            for r = 1:numel(self)
                datasample = self(r).datardm(inds,inds,:);
                modelsample = self(r).Xrdm(inds,inds,:);
                % make off-diagonal zeros NaN (these then get removed on
                % RSA class init)
                datasample(matmask(:,:,1:self(r).nfeatures)) = NaN;
                modelsample(matmask(:,:,1:self(r).npredictors)) = NaN;
                % create a new instance of whatever class the current
                % instance is
                bootglm(r) = feval(class(self),modelsample,datasample);
                % pass any custom properties from the current instance to
                % the boot instance (e.g., the RidgeRSA k parameter)
                cloneargs(bootglm(r),self(r));
            end
            removedprop = (self(1).nsamples - bootglm(1).nsamples) / ...
                self(1).nsamples;
        end

        function model = selectconditions(self,cons)
        % return a new instance only containing the samples corresponding
        % to the condition indices cons (numerical or logical). note that
        % cons is not sorted so that the condition order can be permuted.
        %
        % model = selectconditions(self,cons)
            % nchoosek here rather than nsamples since we may have stripped
            % NaNs from vector
            sampleinds = rdmvecindices(cons,nchoosek(self(1).ncon,2));
            % make numerical indices into conditions
            if islogical(cons)
                cons = find(cons);
            end
            model = self.copy;
            % maybe parfor
            for r = 1:length(model)
                % update model and data RDMs
                model(r).Xrdm = model(r).Xrdm(cons,cons);
                model(r).datardm = model(r).datardm(cons,cons,:);
                % our sampleinds vector may be invalid if NaNs have been
                % removed. the correct vector is then
                inds = sampleinds(model(r).validvec);
                % select samples
                model(r).data = model(r).data(inds,:);
                model(r).X = model(r).X(inds,:);
                % conversely, we must now update the validvec to only
                % include sampleinds
                model(r).validvec = model(r).validvec(sampleinds);
            end
            % all bets are off if this varies by run so may as well assign
            % once for all
            [model.ncon] = deal(length(cons));
        end

        function [r2bycon,meanpredict,predictbycon] = crossvalidateconditions(self)
        % predict the dissimilarities for one condition based on fitted
        % dissimilarities for all other conditions - ie,
        % leave-one-condition out rather than leave-one-run out (see
        % crossvalidatefit) or leave-one-dissimilarity out.
        %
        % Note that although each train/test split is independent (unlike
        % leave-one-dissimilarity out), the test estimates are dependent
        % across folds. Use e.g. bootstrapping to account for
        % this.
        %
        % Yfit = crossvalidateconditions()
            allinds = 1:self(1).ncon;
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
            for con = allinds
                % use train / ~train to explicitly divide train/test
                % dissimilarities
                test = allinds==con;
                % split data by creating separate train and test model
                % instances
                trainmodel = selectconditions(self,~test);
                testmodel = selectconditions(self,test);
                % predict dissimilarities for test based on training fit
                prediction = trainmodel.predictY(testmodel.X);
                % and get r2 for model prediction
                r2bycon(con,:) = testmodel.rsquare(prediction);
                % because we are working with RDMs, each X in testmodel is
                % identical. So although the prediction spans all splits we
                % just need to store the first nsamples
                predictbycon(con,testind,:) = prediction(...
                    1:testmodel(1).nsamples,:);
            end
            % mean of all the predictions for each dissimilarity (sort of
            % like making a weave - every dissimilarity features in 2
            % predictions. Once when you hit the right row, once when you
            % hit the right column).
            meanpredict = squeeze(sum(predictbycon,1)) / 2;
        end
    end
end
