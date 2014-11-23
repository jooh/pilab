% Spearman rho-based representational similarity analysis. See also the RSA,
% GLM superclasses.
%
% gl = RankRSA(modelrdms,datardms)
classdef RankRSA < RSA
    methods
        function gl = RankRSA(modelrdms,datardms)
            if nargin==0
                modelrdms = [];
                datardms = [];
            end
            % use super-class to initialise            
            gl = gl@RSA(modelrdms,datardms);
            % then rank trans and Z score so that linear fits become
            % equivalent to Spearman rho
            gl.X = zscore(ranktrans(gl.X),0,1);
            gl.data = zscore(ranktrans(gl.data),0,1);
            % finally reduce precision since we are working on ranked data
            % anyway (you'd need a very large number of unique ranks to run
            % into any precision trouble with single floats)
            if isa(gl.X,'double')
                gl.X = single(gl.X);
            end
            if isa(gl.data,'double')
                gl.data = single(gl.data);
            end
        end
    end
end
