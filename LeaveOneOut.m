classdef LeaveOneOut < Splitter
    % leave one out cross validation - provided either an index vector for
    % chunks or a Volume instance, return an object. See super class
    % Splitter for return info.
    % sp = LeaveOneOut(chunkorvol)
    methods
        function sp = LeaveOneOut(in)
            if ~nargin
                return
            end
            if isa(in,'Volume')
                chunks = in.chunks;
            else
                chunks = in;
            end
            ndata = length(chunks);
            [junk,junk,inds] = unique(chunks);
            nsplits = length(junk);
            assert(nsplits > 1,...
                'must have at least 2 unique chunks')
            tr = logical(eye(nsplits));
            % upcast to ndata
            test = tr(inds,:);
            train = ~test;
            sp = sp.initialise(sp,train,test);
        end
    end
end
