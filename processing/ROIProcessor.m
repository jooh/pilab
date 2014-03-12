% MetaProcessor subclass. Iterate over samples in some rois to generate a
% result matrix. Each ROI dataset gets passed to the processor instance
% using call(analysis,design,data,chunks). Note that this input pattern
% places some constraints on what type of Processor instances can be
% supported (most likely a GLMMetaProcessor).
%
% construction:
% vl = ROIProcessor(rois,processor,[batchsize=1e4],[minn=10])
%
% call
% result = call(self,designvol,epivol)
classdef ROIProcessor < MetaProcessor
    properties
        rois
        batchsize = 10000;
        minn = 10;
        batchinds
        batchmat
        nbatch
    end

    methods
        function vl = ROIProcessor(rois,processor,batchsize,minn)
            vl = vl@MetaProcessor(processor,@(varargin)cat(2,varargin{:}));
            if ~ieNotDefined('batchsize')
                vl.batchsize = batchsize;
            end
            if ~ieNotDefined('minn')
                vl.minn = minn;
            end
            vl.rois = rois;
            vl.combiner = @(varargin)cat(2,varargin{:});
            % pad out last batch with NaNs
            vl.nbatch = ceil(vl.rois.nsamples/vl.batchsize);
            vl.batchmat = NaN([vl.batchsize vl.nbatch]);
            vl.batchmat(1:vl.rois.nsamples) = 1:vl.rois.nsamples;
        end

        function varargout = call(self,designvol,epivol)
            assert(isequal(epivol.xyz,self.rois.xyz),...
                'mismatched roi/epivol');
            assert(isequal(epivol.meta.samples.chunks,...
                designvol.meta.samples.chunks),...
                'mismatched epivol/designvol');
            designmat = designvol.data;

            % prepare the output
            nreturn = self(1).processor(1).nreturn;
            result = cell(self.rois.nsamples,nreturn);
            % parfor indexing in non-parfor'ed dimension only works if you
            % specify the vector outside the loop
            retind = 1:nreturn;
            goodres = true(1,self.rois.nsamples);

            % batch out ROIs to allow a smaller epivol. Avoids Matlab
            % memory problems when parfor involves > 2GB of data and avoids
            % passing a huge epivol around when only a small part of it
            % will actually be used.
            for batch = 1:self.nbatch
                if any(isnan(self.batchmat(:,batch)))
                    thisbatchsize = find(isnan(self.batchmat(:,batch)),...
                        1,'first')-1;
                else
                    thisbatchsize = self.batchsize;
                end
                % voxels in any roi and not nan
                batchvox = any(full(self.rois.data(self.batchmat(1:thisbatchsize,batch),:)~=0)...
                    ,1);
                % pick these for batch-specific epivol
                % may need linind2featind here
                batchepis = epivol(:,batchvox);
                batchrois = self.rois(:,batchvox);
                assert(numel(batchvox)==epivol.nfeatures,'mismatched mask and data');
                assert(numel(batchvox)==self.rois.nfeatures,'mismatched mask and data');

                % just arrays for speed
                chunks = epivol.meta.samples.chunks;
                roimat = batchrois.data;
                epimat = batchepis.data;

                parfor b = self.batchmat(1:thisbatchsize,batch)'
                    % um, can't we just use b directly?
                    % skip empty rois (these come out as NaN)
                    validvox = full(roimat(b,:)~=0);
                    if ~any(validvox) || sum(validvox)<self.minn
                        goodres(b) = false;
                        % empty or too small roi
                        continue
                    end
                    % here ideally one would use simply
                    % [result{b,retind}] = call@MetaProcessor(self,...
                    % but matlab seems unable to parse superclass method
                    % calls in parfor so instead we replicate this
                    % behaviour explicitly
                    [result{b,retind}] = call(self.processor,...
                        designmat,epimat(:,validvox),chunks);
                end % parfor b = self.batchmat...
            end % for batch = 1:self.nbatch
            % insure that we plug in appropriate nans
            result(~goodres,:) = {NaN(size(result{find(...
                goodres,1,'first')}))};
            % and off to the combiner to concatenate the ROI dim
            [varargout{retind}] = combinereturns(self,result);
        end % call method

        function roivol = result2roivol(self,result,searchvol,metasamples)
        % Write out a BaseVolume for results, preserving as much meta data
        % as possible.
        %
        % vol = result2roivol(self,result,[searchvol=false],[metasamples])
            if ieNotDefined('metasamples')
                metasamples = [];
            end
            if ieNotDefined('searchvol')
                searchvol = false;
            end
            assert(ndims(result)==2,'result must be 2d');
            assert(size(result,2)==self.rois.nsamples,...
                'result does not match self.rois.nsamples');
            % extra ROI properties for meta features
            metafeatures = self.rois.meta.samples;
            metafeatures.nfeatures = sum(self.rois.data~=0,2)';
            [metafeatures.com_x,metafeatures.com_y,metafeatures.com_z] =...
                deal(NaN([1 self.rois.nsamples]));
            for c = 1:self.rois.nsamples
                % compute centre of mass for this ROI
                coords = round(mean(self.rois.linind2coord(...
                    self.rois.linind(self.rois.data(c,:)~=0)),2));
                metafeatures.com_x(c) = coords(1);
                metafeatures.com_y(c) = coords(2);
                metafeatures.com_z(c) = coords(3);
            end
            if searchvol
                % need to update MriVolume mask to only include the
                % searchlights.
                % but pray tell how does one know? I guess this will only
                % work if nsamples==nfeatures
                assert(self.rois.nsamples==self.rois.nfeatures,...
                    ['mapping to searchvol only works if ' ...
                    'rois.nsamples==rois.nfeatures']);
                roivol = MriVolume(result,self.rois,'metafeatures',...
                    metafeatures,'metasamples',metasamples);
            else
                % make a headerless volume 
                roivol = BaseVolume(result,'metafeatures',...
                    metafeatures,'metasamples',metasamples);
            end
            % drop nan ROIs
            anynan = ~any(isnan(roivol.data),1);
            roivol = roivol(:,anynan);
        end
    end
end
