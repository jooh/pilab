% Partial Spearman rho RSA.
% gl = PartialRankRSA(modelrdms,datardms,partialrdms)
classdef PartialRankRSA < RSA
    methods
        function gl = PartialRankRSA(modelrdms,datardms,partialrdms)
            % partial r is basically the r between modelrdms and datardms
            % after partialrdms has been fitted and subtracted from each.
            tempmodel = RankRSA(partialrdms,modelrdms);
            resmodel = residuals(tempmodel);
            tempdata = RankRSA(partialrdms,datardms);
            resdata = residuals(tempdata);
            % then use super-class to initialise            
            % NB, confusingly, we DON'T want to use RankRSA as superclass
            % here because then we would incur another rank transform of
            % the residuals. 
            gl = gl@RSA(resmodel,resdata);
        end
    end
end
