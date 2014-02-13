% RSA sub-class for estimating mean dissimilarity. The modelrdm input is
% used as a NaN mask to optionall estimate the mean over a subset of
% dissimilarities. Assumes that you're using ldtRDMs or something
% similarly symmetrical about 0 under the null.
%
% gl = MeanRSA(modelrdms,datardms)
classdef MeanRSA < RSA
    methods
        function gl = MeanRSA(modelrdms,datardms)
            gl = gl@RSA(modelrdms,datardms);
            % check that X is of the correct type for this class
            assert(gl.npredictors == 1,...
                'only one model RDM is supported');
            uv = unique(gl.X);
            assert(length(uv)==1,...
                'model RDM must contain exactly 1 unique dissimilarity');
            % means are easy
            gl.X(:) = 1;
        end
    end
end
