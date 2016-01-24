% Processor sub-class. Calling some GLM methods on an instance, and use
% some combiner function to combine the results from each instance in the
% array.
%
% So, for instance: cross-classification A to B and B to A could be
% implemented by having two instances in the array - one that specifies
% each direction. We would obtain 1 return on call, which is averaged
% across both cross-classification directions. 
%
% Or, permutation testing of contrast and standard error estimate could be
% implemented by having one instance where we use {contrast,standarderror}
% and methargs. We would obtain 2 returns on call.
%
% there may be multiple processes - numel(self.processor). These
% come out in cell array form and get concatenated in an unused
% dimension and passed to combiner (typically mean) to yield one
% matrix per return.
%
% Note that you will probably need to combine this instance with the
% GLMMetaProcessor class that also takes care of constructing the model
% instance.
%
% construction:
% gl = GLMProcessor(glmmethod,[combiner],nreturn,[glmmethodargs]);
%
% call:
% varargout = call(self,model)
classdef GLMProcessor < Processor
    properties
        glmmethod = '';
        varargs = {};
        selfind = [];
    end

    methods
        function gl = GLMProcessor(glmmethod,combiner,nreturn,varargin);
            if ieNotDefined('combiner')
                combiner = [];
            end
            if ieNotDefined('nreturn')
                nreturn = [];
            end
            gl = gl@Processor(combiner,nreturn);
            if ~any(nargin)
                return
            end
            gl.glmmethod = glmmethod;
            gl.varargs = varargin;
            selfposs = strcmp(varargin,'*SELF*');
            if any(selfposs)
                logstr('detected *SELF* placeholder\n')
                gl.selfind = find(selfposs);
            end
        end

        function varargout = call(self,model)
        % run an array of processes on the current model and combine all
        % the results.
        %
        % varargout = call(self,model)
            % check that our entries in self are homogeneous
            lastwarn('');
            nreturn = self(1).nreturn;
            assert(isequal(nreturn,self.nreturn),...
                'entries in GLMProcessor instance must have same nreturn');
            assert(isequal(self(1).combiner,self.combiner),...
                'all GLMProcessor entries must have the same combiner');
            nprocess = numel(self);
            result = cell(nprocess,nreturn);
            % we may be dealing with an array of processes to apply
            for proc = 1:nprocess
                % run the method for this process (it may return multiple
                % outputs)
                args = self(proc).varargs;
                args{self(proc).selfind} = model;
                [result{proc,1:nreturn}] = feval(self(proc).glmmethod,...
                    model,args{:});
            end
            % combine the outputs across processes
            [varargout{1:nreturn}] = combinereturns(self,result);
            [msg,warntype] = lastwarn;
            assert(~strcmp(warntype,'MATLAB:nearlySingularMatrix'),...
                'nearly singular result');
        end
    end
end
