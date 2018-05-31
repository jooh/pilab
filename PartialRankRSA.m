% Partial Spearman rho RSA.
% gl = PartialRankRSA(modelrdms,datardms,partialrdms)
classdef PartialRankRSA < PearsonRSA
    methods
        function gl = PartialRankRSA(modelrdms,datardms,partialrdms)
            % partial r is basically the r between modelrdms and datardms
            % after partialrdms has been fitted and subtracted from each.
            tempmodel = RankRSA(partialrdms,modelrdms);
            resmodel = residualsfull(tempmodel);
            tempdata = RankRSA(partialrdms,datardms);
            resdata = residualsfull(tempdata);
            % then use super-class to initialise            
            % NB, confusingly, we DON'T want to use RankRSA as superclass
            % here because then we would incur another rank transform of
            % the residuals. 
            gl = gl@PearsonRSA(resmodel,resdata);
        end
    end
end
