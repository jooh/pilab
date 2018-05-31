% RankRSA & RidgeGLM sub-class for regularised representational similarity
% analysis.
%
% gl = RidgeRSA(modelrdms,datardms,k)
classdef RidgeRSA < RankRSA & RidgeGLM
    methods
        function gl = RidgeRSA(modelrdms,datardms,k)
            gl = gl@RidgeGLM;
            gl = gl@RankRSA(modelrdms,datardms);
            if ieNotDefined('k')
                k = 0;
            end
            gl.k = k;
            gl.unscale = 1;
        end

        function cloneargs(self,oldinstance)
            % make sure k parameter survives resampling.
            [self.k] = deal(oldinstance.k);
        end
    end
end
