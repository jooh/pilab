classdef SplitMeanDiffRSA < MeanDiffRSA & SplitRSA
    methods
        function gl = SplitMeanDiffRSA(modelrdms,datardms)
            gl = gl@SplitRSA;
            gl = gl@MeanDiffRSA(modelrdms,datardms);
            % ideally this should be set by calling SplitRSA constructor
            gl.nrandsamp = gl.ncon / 2;
        end
    end
end

