% GLM sub-class for regularised ridge fits.
% gl = RidgeGLM(X,Y,k);
classdef RidgeGLM < GLM
    properties
        k % regularisation parameter
    end

    methods
        function gl = RidgeGLM(X,Y,k)
            [gl.X,gl.Y,gl.k] = deal(X,Y,k);
            gl.getdiagnostics;
            stx = std(X,[],1) == 0 ;
            if ~any(stx)
                % no constant in X, so don't rescale fit
                rescale = 0;
                Xfit = X;
            else
                % make sure only one, well-placed constant
                assert(sum(stx)==1 && find(stx)==gl.npredictors,...
                    ['only one constant is supported, and it must be '...
                    'the last column in X']);
                % do rescale
                rescale = 1;
                % strip constant since ridgevec inserts its own in betas
                Xfit = X(:,1:gl.npredictors-1);
            end
            gl.betas = ridgevec(gl.Y,Xfit,gl.k,rescale);
        end
    end
end
