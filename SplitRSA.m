classdef SplitRSA < RSA
    methods
        function gl = SplitRSA(varargin)
            gl = gl@RSA(varargin{:});
            % the only property modification compared to regular RSA is
            % that the nrandsamp is half
            gl.nrandsamp = gl.ncon / 2;
        end

        function permglm = drawpermsample(self,inds)
        % return a new instance where the rows in the split-data quadrants
        % of the Xrdm have been re-ordered according to inds. Note that you
        % must supply the same number of inds as self.ncon/2 (since it's
        % split data). Overrides GLM and RSA base behaviours.
            permX = upcastsplitrdm(self(1).Xrdm(inds,...
                self(1).nrandsamp+1:end,:));
            permglm = self.copy;
            for r = 1:numel(self)
                permglm(r) = feval(class(self),permX,self(r).datardm);
                cloneargs(permglm(r),self(r));
            end
        end

        function [bootglm,removedprop] = drawbootsample(self,inds)
        % return a new instance where the rows in the split-data quadrants
        % of the Xrdm and datardms have been sampled with replacement in
        % concert. Overrides GLM and RSA base behaviours.
        %
        % Note that for consistency with RSA class, we also return the
        % diagnostic removedprop (see RSA.drawbootsample for details), even
        % though removedprop will _always_ be 0 for this sub-class since
        % the data contain no 0 diagonals that could get resampled to
        % off-diagonal positions.

            assert(all(all(vertcat(self.data)~=0)),...
                'cannot sample when zeros are present in data');
            % this approach assumes that zeros only happen in the data
            % dissimilarities when diagonals have been resampled to
            % off-diagonal positions
            % find zeros in the boot sampled rdm mat - ie, diagonals that
            % moved.
            matmask = self(1).datardm(inds,self(1).nrandsamp+1:end,1)==0;
            % skip true diagonals
            matmask(logical(eye(self(1).ncon))) = 0;
            % pre-allocate to biggest needed size so we can later just
            % slice for speed
            matmask = repmat(matmask,[1 1 max([self.nfeatures ...
                self.npredictors])]);
            % process each run in the new instance
            for r = 1:numel(self)
                datasample = self(r).datardm(inds,...
                    self(1).nrandsamp+1:end,:);
                modelsample = self(r).Xrdm(inds,...
                    self(1).nrandsamp+1:end,:);
                % make off-diagonal zeros NaN (these then get removed on
                % RSA class init)
                datasample(matmask(:,:,1:self(r).nfeatures)) = NaN;
                modelsample(matmask(:,:,1:self(r).npredictors)) = NaN;
                % upcast and create a new instance of whatever class the current
                % instance is
                bootglm(r) = feval(class(self),upcastsplitrdm(...
                    modelsample),upcastsplitrdm(datasample));
                % pass any custom properties from the current instance to
                % the boot instance (e.g., the RidgeRSA k parameter)
                cloneargs(bootglm(r),self(r));
            end
            removedprop = (self(1).nsamples - bootglm(1).nsamples) / ...
                self(1).nsamples;
            assert(removedprop==0,'uh oh, something is broken');
        end
    end
end
