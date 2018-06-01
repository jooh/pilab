% Split RDM variant of PluginRSA.
%
% gl = SplitPluginRSA(modelrdms,datardms,corrtype)
classdef SplitPluginRSA < PluginRSA & SplitRSA
    methods
        function gl = SplitPluginRSA(modelrdms,datardms,varargin)
            gl = gl@SplitRSA;
            gl = gl@PluginRSA(modelrdms,datardms,varargin{:});
            % ideally this should be set by calling SplitRSA constructor
            gl.nrandsamp = gl.ncon / 2;
        end
    end
end
