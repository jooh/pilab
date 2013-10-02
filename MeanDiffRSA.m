% RSA sub-class for fitting a single, binary model RDM to data. By
% de-meaning the data and model the fitted response becomes equivalent to
% the mean difference between modelrdm==0 and modelrdm==1. Note that
% statistically this test is nearly equivalent to RankRSA...
%
% gl = MeanDiffRSA(modelrdms,datardms)
classdef MeanDiffRSA < RSA
    methods
        function gl = MeanDiffRSA(modelrdms,datardms)
            gl = gl@RSA(modelrdms,datardms);
            % check that X is of the correct type for this class
            assert(gl.npredictors == 1,...
                'only one model RDM is supported');
            uv = unique(gl.X);
            assert(length(uv)==2,...
                'model RDM must contain exactly 2 unique dissimilarities');
            % de-mean and set to unit length so that fits become equivalent
            % to mean difference
            gl.X(gl.X==uv(1)) = 0;
            gl.X(gl.X==uv(2)) = 1;
            gl.X = gl.X - mean(gl.X);
            % de-mean each feature to avoid needing a constant
            gl.data = bsxfun(@minus,gl.data,mean(gl.data,1));
        end
    end
end
