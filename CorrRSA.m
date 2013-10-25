% RSA based on matlab's builtin corr function, applied to each feature in
% a parfor loop. This is slower than vectorised forms of RSA but enables
% custom RSA types such as Kendall's tau.
%
% gl = CorrRSA(modelrdms,datardms,corrtype)
classdef CorrRSA < RSA
    properties
        corrtype = []; % type input to matlab's corr function
    end

    methods
        function gl = CorrRSA(modelrdms,datardms,corrtype)
            % use super-class to initialise            
            gl = gl@RSA(modelrdms,datardms);
            assert(gl.npredictors==1,'only 1 model RDM input supported');
            % must allow undefined corrtype to support resampling
            % (cloneargs method is only called on new instance _after_ it
            % has been created)
            if ieNotDefined('corrtype')
                corrtype = [];
            end
            gl.corrtype = corrtype;
        end

        function cloneargs(self,oldinstance)
            % make sure corrtype survives resampling.
            [self.corrtype] = deal(oldinstance.corrtype);
        end

        function estimates = fit(self)
        % apply some corrtype to every feature in self.data.
        % estimates = fit(self)
            estimates = NaN([1 self(1).nfeatures]);
            X = vertcat(self.X);
            data = vertcat(self.data);
            inds = 1:self(1).nsamples;
            corrtype = self(1).corrtype;
            % parfor here means we will encounter nested parfor loops when
            % e.g. bootstrapping or permutation testing. However, if we
            % don't do this e.g searchlight mapping is likely to turn out
            % to be ridiculously slow
            parfor f = 1:self(1).nfeatures
                estimates(1,f) = corr(X,data(inds,f),'type',corrtype);
            end
        end
    end
end
