% GLM sub-class for regularised ridge fits.
% gl = RidgeGLM(X,data,k,[unscale])
classdef RidgeGLM < GLM
    properties
        k % regularisation parameter (gets scaled by k*nsamples in ridgevec)
        unscale % return data to original scale (default 0)
    end

    methods
        function gl = RidgeGLM(X,data,k,unscale)
            if nargin == 0
                % support no input case
                [X,data,k] = deal([]);
            end
            % initialise with super-class constructor
            gl = gl@GLM(X,data);
            gl.k = k;
            if ieNotDefined('unscale')
                unscale = 0;
            end
            gl.unscale = unscale;
        end

        function estimates = fit(self)
        % estimates = fit(self)
        % Fit ridge regression model. All runs must have the same k
        % parameter. k can be a scalar
            % check if we can just do one quick fit
            ks = unique(self(1).k);
            assert(all([self.unscale] == self(1).unscale),...
                'only 1 unique unscale option is permitted');
            if length(ks)==1
                estimates = ridgevec(vertcat(self.data),vertcat(self.X),...
                    self(1).k,self(1).unscale);
                return
            end
            % failing that, fit each unique k separately
            estimates = NaN([self(1).npredictors self(1).nfeatures]);
            data = vertcat(self.data);
            X = vertcat(self.X);
            for k = ks
                inds = self(1).k == k;
                estimates(:,inds) = ridgevec(data(:,inds),X,k,...
                    self(1).unscale);
            end
        end
    end
end
