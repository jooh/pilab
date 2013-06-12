% GLM sub-class for representational similarity analysis.
% gl = RSAGLM(modelrdms,datardms)
classdef RankRSA < RSA
    methods
        function gl = RankRSA(modelrdms,datardms)
            % use super-class to initialise            
            gl = gl@RSA(modelrdms,datardms);
            % then rank trans and Z score so that linear fits become
            % equivalent to Spearman rho
            gl.X = zscore(tiedrank(gl.X),0,1);
            gl.data = zscore(tiedrank(gl.data),0,1);
        end
    end
end
