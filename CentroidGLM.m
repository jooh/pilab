% GLM sub-class that removes the featurewise centroid and constant from the
% data. So basically this subtracts out a particular flavour of univariate
% effect - a single pattern that waxes and wanes (whether scaling is additive or
% multiplicative) over samples.
classdef CentroidGLM < GLM
    methods
        function gl = CentroidGLM(X,data)
            if nargin<1
                X = [];
                data = [];
            end
            % initialise with super class
            gl = gl@GLM(X,data);
        end

        function con = contrast(self,conmat)
            % first get all the parameter estimates as usual
            est = fit(self);
            relevantind = conmat~=0;
            assert(all(sum(relevantind,2)==2),['this code only '...
                'supports simple pairwise contrasts at present'])
            nfeat = size(est,2);
            ncon = size(conmat,1);
            featind = 1:nfeat;
            con = NaN([ncon nfeat]);
            parfor n = 1:ncon
                % save a bit of compute by indexing the relevant conditions
                % directly
                % subtract the mean across the voxels for each of the
                % relevant conditions ([2 1])
                est0 = bsxfun(@minus,est(relevantind(n,:),:),...
                    mean(est(relevantind(n,:),:),2));
                % remove the class centroid (across examples) from each
                % pattern ([1 nvoxels] transposed)
                est0cent = projectout(est0',mean(est0,1)')';
                % and get the contrast estimate at last
                con(n,featind) = conmat(n,relevantind(n,:)) * est0cent;
            end
        end

        function [t,model] = infot(self,w,conmat)
            t = [];
            model = [];
            error(['infot is not supported for CentroidGLM at present ' ...
                '- use infoc']);
        end
    end
end
