% Multiple regression RSA by squaring the predictors (and square rooting
% the contrast estimates). Note that all other methods preserve the squared
% data.
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

        function con = contrast(self,conmat)
            con = sqrtsigned(contrast@RSA(self,conmat));
        end
    end
end
