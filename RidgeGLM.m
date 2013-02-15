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
            estimates = ridgevec(self.data,self.X,self.k);
        end
    end
end
