classdef SplitRankRSA < RankRSA & SplitRSA
    methods
        function gl = SplitRankRSA(modelrdms,datardms)
            gl = gl@SplitRSA;
            gl = gl@RankRSA(modelrdms,datardms);
            % ideally this should be set by calling SplitRSA constructor
            gl.nrandsamp = gl.ncon / 2;
        end
    end
end
