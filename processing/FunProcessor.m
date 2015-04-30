classdef FunProcessor < Processor
    properties
        funhand = [];
        funargs = {};
    end

    methods
        function fpr = FunProcessor(funhand,combiner,nreturn,varargin)
            if ieNotDefined('combiner')
                combiner = [];
            end
            if ieNotDefined('nreturn')
                nreturn = [];
            end
            fpr = fpr@Processor(combiner,nreturn);
            if ~any(nargin)
                return
            end
            fpr.funhand = funhand;
            fpr.funargs = varargin;
        end

        function varargout = call(self,varargin)
            nreturn = self(1).nreturn;
            assert(isequal(nreturn,self.nreturn),...
                'entries in FunProcessor instance must have same nreturn');
            assert(isequal(self(1).combiner,self.combiner),...
                'all FunProcessor entries must have the same combiner');
            nprocess = numel(self);
            result = cell(nprocess,nreturn);
            for proc = 1:nprocess
                [result{proc,1:nreturn}] = feval(self(proc).funhand,...
                    varargin{:},self(proc).funargs{:});
            end
            % combine the outputs across processes
            [varargout{1:nreturn}] = combinereturns(self,result);
        end
    end
end
