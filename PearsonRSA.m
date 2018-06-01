% Pearson r-based representational similarity analysis. See also the RSA, GLM
% superclasses.
%
% gl = PearsonRSA(modelrdms,datardms)
classdef PearsonRSA < RSA
    methods
        function gl = PearsonRSA(modelrdms,datardms)
            % use super-class to initialise            
            gl = gl@RSA(modelrdms,datardms);
            % then Z score so that linear fits become
            % equivalent to Pearson r
            gl.X = zscore(gl.X,0,1);
            gl.data = zscore(gl.data,0,1);
        end
    end
end
