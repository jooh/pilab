% GLM sub-class for covariate-based fits based on projecting out the
% covariate from X and Y. Note the the covariatedegree input is optional -
% if you provide none you get base class (GLM) behaviour.
% gl = CovGLM(X,Y,[covariatedegree]);
classdef CovGLM < GLM
    properties
        covariatedegree
    end

    methods
        function gl = CovGLM(X,Y,covariatedegree)
            if ~ieNotDefined('covariatedegree')
                if covariatedegree == 0
                    % special case - simple de-mean can be done much faster
                    % like this
                    Y = bsxfun(@minus,Y,mean(Y));
                    X = bsxfun(@minus,X,mean(X));
                else
                    covariates = constructpolynomialmatrix(size(Y,1),...
                        (0:covariatedegree));
                    % get Y/X with covariates projected out
                    projector = projectionmatrix(covariates);
                    Y = projector * Y;
                    X = projector * X;
                end
            end
            % initialise with super-class constructor
            gl = gl@GLM(X,Y);
            % can only assign to CovGLM AFTER initialise
            if ~ieNotDefined('covariatedegree')
                gl.covariatedegree = covariatedegree;
            end
        end
    end
end
