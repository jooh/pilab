% Representational dissimilarity analysis. The key difference from
% GLMMetaProcessor is that the model constructor is called with a single
% input (rdvecs). This instance can be plugged into SecondLevelProcessor as
% the second stage.
%
% construction:
% ga = RSAMetaProcessor(constructor,processor,[combiner])
%
% call:
% result = call(rdvecs)

classdef RSAMetaProcessor < GLMMetaProcessor
    methods
        function rm = RSAMetaProcessor(varargin)
            rm = rm@GLMMetaProcessor(varargin{:});
            if ~nargin
                return
            end
        end

        function varargout = call(self,rdvecs)
            model = call(self.constructor,rdvecs);
            [varargout{1:self(1).processor(1).nreturn}] = ...
                call@MetaProcessor(self,model);
        end
    end
end
