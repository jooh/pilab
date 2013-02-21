% CovGLM sub-class for regularised ridge fits.
% gl = RidgeGLM(X,data,k,[covariatedegree]);
classdef RidgeGLM < CovGLM
    properties
        k % regularisation parameter (gets scaled by k*nsamples in ridgevec)
    end

    methods
        function gl = RidgeGLM(X,data,k,covariatedegree)
            % initialise with super-class constructor
            if ieNotDefined('covariatedegree')
                covariatedegree = [];
            end
            gl = gl@CovGLM(X,data,covariatedegree);
            gl.k = k;
        end

        function estimates = fit(self)
        % estimates = fit(self)
        % Fit ridge regression model. All runs must have the same k
        % parameter. k can be a scalar
            % check if we can just do one quick fit
            ks = unique(self(1).k);
            if length(ks)==1
                estimates = ridgevec(vertcat(self.data),vertcat(self.X),...
                    self(1).k);
                return
            end
            % failing that, fit each unique k separately
            estimates = NaN([self(1).npredictors self(1).nfeatures]);
            data = vertcat(self.data);
            X = vertcat(self.X);
            for k = ks
                inds = self(1).k == k;
                estimates(:,inds) = ridgevec(data(:,inds),X,k);
            end
        end
    end
end
