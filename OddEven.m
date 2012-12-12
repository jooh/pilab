classdef OddEven < Splitter
    % odd/even split. Input in can be a Volume instance or an index vector
    % of chunks. See super class Splitter for return info.
    % sp = OddEven(in)
    methods
        function sp = OddEven(in)
            % need this for repmat below to work. For some reason.
            if ~nargin
                return
            end
            if isa(in,'Volume')
                chunks = in.chunks;
            else
                chunks = in;
            end
            assert(length(unique(chunks)) > 1,...
                'must have at least 2 unique chunks')
            train = [rem(chunks,2)>0 rem(chunks+1,2)>0];
            test = ~train;
            % need to explicitly pass sp here because matlab is very poor
            sp = sp.initialise(sp,train,test);
        end
    end
end
