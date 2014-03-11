% MetaProcessor subclass. Make a model with GLMConstructor, call some
% GLMProcessor and return some result (after passing through combiner,
% default handle taken from Processor).
%
% Initialise:
% ga = GLMMetaProcessor(constructor,processor,[combiner])
%
% Call:
% result = call(ga,design,data,chunks)
%
classdef GLMMetaProcessor < MetaProcessor
    properties
        constructor
    end

    methods
        function ga = GLMMetaProcessor(constructor,processor,combiner)
            if ieNotDefined('combiner')
                combiner = [];
            end
            ga = ga@MetaProcessor(processor,combiner);
            ga.constructor = constructor;
        end

        function varargout = call(self,design,data,chunks)
            model = call(self.constructor,design,data,chunks);
            [varargout{1:self(1).processor(1).nreturn}] = ...
                call@MetaProcessor(self,model);
        end
    end
end
