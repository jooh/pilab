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
            % store unaltered data in matrix form (for e.g. bootstrapping)
            Xrdm = asrdmmat(modelrdms);
            datardm = asrdmmat(datardms);
            ncon = size(Xrdm,1);
            % convert to vector and strip NaNs
            Xvec = rdm2vec(Xrdm);
            datavec = rdm2vec(datardm);
            nanmask = any(isnan(Xvec),2) | any(isnan(datavec),2);
            Xvec(nanmask,:) = [];
            datavec(nanmask,:) = [];
            assert(~isempty(Xvec),'no valid (non-NaN) data entered');
            % use super-class to initialise            
            gl = gl@GLM(Xvec,datavec);
            [gl.ncon,gl.Xrdm,gl.datardm,gl.validvec] = deal(ncon,Xrdm,...
                datardm,~nanmask);
        end

        function cloneargs(self,oldclass)
            % does nothing for base RSA case.
        end

        function permglm = permuteconditions(self)
        % return a new instance where the rows and columns in the Xrdm
        % have been permuted without replacement.
        %
        % permglm = permuteconditions()
            xtemp = self.Xrdm;
            % set diagonal to NaN to make the check easier
            xtemp(repmat(logical(eye(self.ncon)),...
                [1 1 self.npredictors])) = NaN;
            xnan = isnan(xtemp);
            goodcons = arrayfun(@(x)~all(xnan(:,:,x),2),...
                (1:self.npredictors)','uniformoutput',false);
            keyboard;
            goodcons = cell2mat(goodcons);
            permglm = [];
        end

        function [bootglm,removedprop] = bootsampleconditions(self)
        % return a new glm instance where the Xrdm and datardm in self have
        % been sampled with replacement.         
        %
        % Unlike sampling without replacement, RDM bootstrapping can shift
        % diagonal dissimilarities to off-diagonal positions.  These should
        % be removed to avoid underestimating the variance (since data and
        % design will have matching zeros). removedprop is a diagnostic
        % output of the proportion of removed dissimilarities.
        %
        % [bootglm,removedprop] = bootsampleconditions(self)
            nrun = length(self);
            % random condition sample with replacement
            consamp = ceil(rand(self(1).ncon,1)*self(1).ncon);
            % this approach assumes that zeros only happen in the data
            % dissimilarities when diagonals have been resampled to
            % off-diagonal positions
            assert(all(all(vertcat(self.data)~=0)),...
                'cannot sample when zeros are present in data');
            % find zeros in rdm mat
            matmask = self(1).datardm(consamp,consamp,1)==0;
            % restrict to off-diagonals
            matmask(logical(eye(self(1).ncon))) = 0;
            % process each run in the new instance - maybe parfor?
            for r = 1:nrun
                datasample = self(r).datardm(consamp,consamp,:);
                modelsample = self(r).Xrdm(consamp,consamp,:);
                % make off-diagonal zeros NaN
                datasample(repmat(matmask,[1 1 self(r).nfeatures])) = NaN;
                modelsample(repmat(matmask,[1 1 self(r).npredictors])) = NaN;
                % create a new instance
                bootglm(r) = feval(class(self),modelsample,datasample);
                % pass any custom properties for the sub-class, e.g. ridge
                % k
                cloneargs(bootglm(r),self(r));
            end
            removedprop = (self(1).nsamples - bootglm(1).nsamples) / ...
                self(1).nsamples;
        end

        function [bootest,removedprops] = bootstrapconditions(self,nboot,bootmeth,outshape)
        % general method for computing bootstrap estimates for some
        % bootmeth (char) by resampling the conditions in the RDM. See
        % GLM.bootstrapruns for details of approach and arguments.
        %
        %
        % [bootest,removedprops] = bootstrapconditions(self,nboot,bootmeth,outshape)
            if ieNotDefined('outshape')
                outshape = size(self.(bootmeth));
                assert(numel(outshape)<3,...
                    'bootmeth must return at most 2d outputs, got %s',...
                    mat2str(outshape));
            end
            bootest = NaN([outshape nboot]);
            allinds = 1:self(1).ncon;
            nrun = length(self);
            % convert data back to 3D RDMs for ease of indexing
            datamats = arrayfun(@(ind)asrdmmat(self(ind).data),...
                1:nrun,'uniformoutput',false);
            modelmats = arrayfun(@(ind)asrdmmat(self(ind).X),...
                1:nrun,'uniformoutput',false);
            removedprops = NaN([1 nboot]);
            parfor b = 1:nboot
                [bootglm,removedprops(b)] = self.bootsampleconditions;
                bootest(:,:,b) = bootglm.(bootmeth);
            end
        end

        function model = selectconditions(self,cons)
        % return a new instance only containing the samples corresponding
        % to the condition indices cons (numerical or logical).
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
