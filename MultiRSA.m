% Multiple regression RSA by squaring the predictors (and square rooting
% the estimates).
%
% gl = MultiRSA(modelrdms,datardms)
classdef MultiRSA < RSA
    methods
        function gl = MultiRSA(modelrdms,datardms)
            if nargin == 0 || (isempty(modelrdms) && isempty(datardms))
                modelrdms = [];
                datardms = [];
            end
            % use super-class to initialise            
            gl = gl@RSA(modelrdms,datardms);
            gl.X = sqsigned(gl.X);
            gl.data = sqsigned(gl.data);
        end

        function cmat = covmat(self)
            % covariance in squared model, then back to data units
            cmat = sqrtsigned(covmat@RSA(self));
        end

        function Yfit = predictY(self,varargin)
        % generate a fitted (predicted) rdm vector for a design matrix X,
        % and square root transform the prediction back to data units
        %
        % Yfit = predictY([X])
            Yfit = sqrtsigned(predictY@RSA(self,varargin{:}));
        end

        function con = contrast(self,conmat)
            con = sqrtsigned(contrast@RSA(self,conmat));
        end

        function data = getdata(self)
        % get data by reversing the square transform and averaging the RDMs
        % across instances.
        % data = getdata(self)
            data = mean(sqrtsigned(zcat(self.data)),3);
        end
    end
end
