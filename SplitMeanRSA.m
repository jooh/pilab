% split RDM variant of MeanRSA
%
% gl = SplitMeanRSA(modelrdms,datardms)
classdef SplitMeanRSA < MeanRSA & SplitRSA
    methods
        function gl = SplitMeanRSA(modelrdms,datardms)
            gl = gl@SplitRSA;
            gl = gl@MeanRSA(modelrdms,datardms);
            % ideally this should be set by calling SplitRSA constructor
            gl.nrandsamp = gl.ncon / 2;
        end
    end
end
