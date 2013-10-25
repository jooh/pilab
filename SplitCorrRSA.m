% Split RDM variant of CorrRSA.
%
% gl = CorrRSA(modelrdms,datardms,corrtype)
classdef SplitCorrRSA < CorrRSA & SplitRSA
    methods
        function gl = SplitCorrRSA(modelrdms,datardms,varargin)
            gl = gl@SplitRSA;
            gl = gl@CorrRSA(modelrdms,datardms,varargin{:});
            % ideally this should be set by calling SplitRSA constructor
            gl.nrandsamp = gl.ncon / 2;
        end
    end
end
