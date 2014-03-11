% Combine sessions.
classdef SessionProcessor < MetaProcessor
    properties
        sessionsplit
        nsplit
        usplit
    end

    methods
        function sp = SessionProcessor(sessionsplit,processor,combiner)
            if ieNotDefined('combiner')
                combiner = [];
            end
            sp = sp@MetaProcessor(processor,combiner);
            % each entry in sessionsplit defines how to split 1 run. The
            % number indicate with which other runs it should be combined.
            sp.sessionsplit = sessionsplit;
            sp.usplit = unique(sessionsplit);
            sp.usplit(isnan(sp.usplit)) = [];
            sp.nsplit = numel(sp.usplit);
        end

        function varargout = call(self,design,data,chunks)
            nreturn = self(1).processor(1).nreturn;
            for split = 1:self(1).nsplit
                % each entry in chunks defines the chunk to which a given
                % sample belongs. so need to first obtain chunkinds
                splitchunk = find(self(1).sessionsplit==...
                    self(1).usplit(split));
                % then find the volumes that have these chunk values
                thissplit = ismember(chunks,splitchunk);
                assert(any(thissplit),'empty split! bad sessionsplit?');
                [result{split,1:nreturn}] = ...
                    call@MetaProcessor(self,design(thissplit,:),...
                    data(thissplit,:),chunks(thissplit));
            end
            % and now the combiner puts the split results together
            [varargout{1:nreturn}] = combinereturns(self,result);
        end
    end
end
