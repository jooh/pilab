% MetaProcessor subclass. feval some function handle on the data before
% passing off to some processor (typically GLMMetaProcessor). Useful for
% e.g. de-meaning the data across features.
%
% construction:
% dpp = ROIPreProcessor(processor,operation)
%
% call:
% varargout = call(self,designmat,epimat,chunks)
classdef ROIPreProcessor < MetaProcessor
    properties
        operation
    end

    methods
        function dpp = ROIPreProcessor(processor,operation)
            if ~nargin
                processor = [];
                operation = [];
            end
            dpp = dpp@MetaProcessor(processor);
            dpp.operation = operation;
        end

        function varargout = call(self,epimat,designmat,chunks)
            assert(numel(self)==1,'only support for one entry');
            % apply some preprocessing
            epimat = feval(self.operation,epimat);
            % cart off to processor (probably a GLMMetaProcessor)
            [varargout{self(1).processor(1).nreturn}] = ...
                call(self.processor,epimat,designmat,chunks);
        end
    end
end
