% GLM sub-class for representational similarity analysis.
% gl = PearsonRSA(modelrdms,datardms)
classdef PearsonRSA < RSA
    methods
        function gl = PearsonRSA(modelrdms,datardms)
            % use super-class to initialise            
            gl = gl@RSA(modelrdms,datardms);
            % then rank trans and Z score so that linear fits become
            % equivalent to Spearman rho
            gl.X = zscore(gl.X,0,1);
            gl.data = zscore(gl.data,0,1);
        end
    end
end
