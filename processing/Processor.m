% Processors perform some operation on data input(s) and return some
% result(s).  Processors can be stacked in array form, in which case a
% combiner (default @matmean) is used by the combinereturns method to
% reduce the output across processors to one per return.
%
% Because Matlab is not very flexible with varargout handling, nreturn must
% currently be specified manually by the user. The upside of this 'feature'
% is that you can elect not to capture later outputs of no interest.
%
% Processor itself is a base class that cannot be used directly.
%
classdef Processor < Saveable
    properties
        nreturn = 1;
        combiner = @matmean;
    end

    methods
        function pr = Processor(combiner,nreturn)
            if ~nargin
                return
            end
            if ~ieNotDefined('nreturn')
                pr.nreturn = nreturn;
            end
            if ~ieNotDefined('combiner')
                pr.combiner = combiner;
            end
        end

        function varargout = call(self,varargin)
            error('do not call base Processor instance directly')
        end

        function varargout = combinereturns(self,result)
            for ret = 1:self(1).nreturn
                % here we assume that you are stacking the returns in
                % columns and the processes in rows
                varargout{:,ret} = feval(self(1).combiner,...
                    result{:,ret});
            end
        end
    end
end
