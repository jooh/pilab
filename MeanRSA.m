% RSA sub-class for estimating mean dissimilarity. We estimate the mean
% over non-zero entries in each model RDM (note that the magnitude of the
% original dissimilarity is discarded).
%
% Note that we crash if there are overlapping valid dissimilarities across
% models, since the resulting estimates are not means in this case. If you
% want to estimate such overlapping means you will need to use separate
% models.
%
% gl = MeanRSA(modelrdms,datardms)
classdef MeanRSA < RSA
    methods
        function gl = MeanRSA(modelrdms,datardms)
            if nargin==0
                modelrdms = [];
                datardms = [];
            end
            gl = gl@RSA(modelrdms,datardms);
            % so we need to do some pre-processing to make this work
            for reg = 1:gl.npredictors
                uv = unique(gl.X(gl.X(:,reg)~=0,reg));
                assert(numel(uv)==1,['each model RDM must contain ' ...
                    'exactly 1 unique non-zero dissimilarity']);
                blips = gl.X(:,reg) == uv;
                gl.X(blips,reg) = 1;
                gl.X(~blips,reg) = 0;
            end
            assert(~any(sum(gl.X,2)>1),...
                'overlapping dissimilarities across model RDMs detected');
        end
    end
end
