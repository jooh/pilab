% Second level analysis. Calls firstprocessor with varargin, and then
% passes the output of this to processor. Used e.g. to calculate a data RDM
% (firstprocessor is GLMMetaProcessor), and compare the result to a model
% RDM (processor is RSAMetaProcessor).
%
% construction:
% sl = SecondLevelProcessor(firstprocessor,processor,combiner)
%
% call:
% varargout = call(varargin)

classdef SecondLevelProcessor < MetaProcessor
    properties
        firstprocessor
    end

    methods

        function sl = SecondLevelProcessor(firstprocessor,processor,combiner)
            if ieNotDefined('processor')
                processor = [];
            end
            if ieNotDefined('combiner')
                combiner = [];
            end
            sl = sl@MetaProcessor(processor,combiner);
            if ~nargin
                return;
            end
            sl.firstprocessor = firstprocessor;
        end

        function varargout = call(self,varargin)
            % make the first result
            firstres = call(self.firstprocessor,varargin{:});
            % 2. Plug the output into main processor to get the output
            [varargout{1:self(1).processor(1).nreturn}] = ...
                call(self.processor,firstres);
        end
    end
end
