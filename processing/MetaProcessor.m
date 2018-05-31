% wrapper class for Processor instance. Generate the result using the
% processor instance's call and combinereturns methods as usual, then apply
% a second pass combiner to the initial result.
%
% the utility of this is to support custom cases where particular aspects
% of the output should be collapsed, e.g. the splits from
% cvcrossclassificationrun.
%
% construction:
% mp = MetaProcessor(processor,[combiner])
%
% call:
% varargout = call(varargin);
classdef MetaProcessor < Processor
    properties
        processor
    end

    methods
        function mp = MetaProcessor(processor,combiner)
            if ieNotDefined('combiner')
                combiner = [];
            end
            if isempty(processor)
                processor = [];
                nreturn = [];
            else
                nreturn = processor(1).nreturn;
                assert(isequal(nreturn,processor.nreturn),...
                    'processors must have same nreturn');
            end
            mp = mp@Processor(combiner,nreturn);
            mp.processor = processor;
        end

        function varargout = call(self,varargin)
            nreturn = self(1).nreturn;
            assert(isequal(nreturn,self.nreturn),...
            'entries in MetaProcessor instance must have same nreturn');
            % obtain (combined) result across processors
            [result{1,1:nreturn}] = call(self.processor,varargin{:});
            % apply a second pass of combining
            [varargout{1,1:nreturn}] = combinereturns(self,result);
        end
    end
end
