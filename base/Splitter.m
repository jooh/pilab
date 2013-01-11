classdef Splitter < handle
    % base class for splitting data - do not call directly. Inherited
    % classes use a somewhat arcane syntax including the initialise method.
    % See OddEven or LeaveOneOut for example implementations.
    %
    % All splitters return object arrays with 3 fields per split:
    % train: logical indices
    % test: logical indices
    % result: cell array with as many entries as there are test datapoints
    properties
        train % ndata by 1 logical vector
        test % ndata by 1 logical vector
        result % 1 by ntest cell
    end

    methods
        function sp = Splitter(varargin)
            if ~nargin
                return
            end
        end

        function sp = initialise(self,sp,train,test)
        % initialise(sp,train,test)
        % Generate object array from train/test matrices
            % preallocate
            nsplits = size(train,2);
            self = repmat(sp,[nsplits 1]);
            for n = 1:nsplits
                sp(n).train = train(:,n);
                sp(n).test = test(:,n);
                sp(n).result = cell(sum(self(n).test),1);
            end
        end
    end
end
