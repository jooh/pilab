% MetaProcessor subclass. Iterate over samples in some rois (instance of
% MriVolume or subclass) to generate results. Each ROI dataset gets passed
% to the processor instance using call(analysis,design,data,chunks). Note
% that this input pattern places some constraints on what type of Processor
% instances can be supported (most likely a GLMMetaProcessor).
%
% construction:
% vl = ROIProcessor(rois,processor,[minn=0],[runfun='runrois_serial'])
%
% call:
% varargout = call(self,data,varargin)
classdef ROIProcessor < MetaProcessor
    properties
        rois
        minn = 0;
        runfun
    end

    methods
        function vl = ROIProcessor(rois,processor,minn,runfun)
            if ieNotDefined('processor')
                processor = [];
            end
            vl = vl@MetaProcessor(processor,@(varargin)cat(2,varargin{:}));
            if ~nargin
                return
            end
            if ~ieNotDefined('minn')
                vl.minn = minn;
            end
            if ieNotDefined('runfun')
                vl.runfun = @runrois_serial;
            else
                vl.runfun = runfun;
            end
            vl.rois = rois;
            vl.combiner = @(varargin)cat(2,varargin{:});
        end

        function varargout = call(self,epimat,varargin)
            % prepare the output
            nreturn = self(1).processor(1).nreturn;
            % parfor indexing in non-parfor'ed dimension only works if you
            % specify the vector outside the loop
            retind = 1:nreturn;
            if ischar(self.runfun)
                self.runfun = str2func(self.runfun);
            end

            % convert to struct to prevent all kinds of exciting Matlab
            % memory problems when passing to parallel compute (basically,
            % the ROIProcessor instance will almost certainly get passed
            % too if you don't do this)
            pstruct = obj2struct(self.processor);
            result = feval(self.runfun,self.rois.data,epimat,pstruct,...
                nreturn,self.minn,varargin{:});
            [varargout{retind}] = combinereturns(self,result);
        end % call method

        function roivol = result2roivol(self,result,searchvol,metasamples)
        % Write out a Volume for results, preserving as much meta data
        % as possible.
        %
        % vol = result2roivol(self,result,[searchvol=false],[metasamples])
            if ieNotDefined('metasamples')
                metasamples = [];
            end
            if ieNotDefined('searchvol')
                searchvol = false;
            end
            assert(ismatrix(result),'result must be 2d');
            assert(size(result,2)==self.rois.nsamples,...
                'result does not match self.rois.nsamples');
            % extra ROI properties for meta features
            metafeatures = self.rois.meta.samples;
            metafeatures.nfeatures = full(sum(self.rois.data~=0,2)');
            if searchvol
                % need to update SPMVolume mask to only include the
                % searchlights.
                % but pray tell how does one know? I guess this will only
                % work if nsamples==nfeatures
                assert(self.rois.nsamples==self.rois.nfeatures,...
                    ['mapping to searchvol only works if ' ...
                    'rois.nsamples==rois.nfeatures']);
                roivol = SPMVolume(result,self.rois,'metafeatures',...
                    metafeatures,'metasamples',metasamples);
            else
                % This is slow so only do this for non-search vols
                [metafeatures.com_x,metafeatures.com_y,...
                    metafeatures.com_z] =...
                    deal(NaN([1 self.rois.nsamples]));
                for c = 1:self.rois.nsamples
                    % compute centre of mass for this ROI
                    coords = round(mean(self.rois.linind2coord(...
                        self.rois.linind(self.rois.data(c,:)~=0)),2));
                    metafeatures.com_x(c) = coords(1);
                    metafeatures.com_y(c) = coords(2);
                    metafeatures.com_z(c) = coords(3);
                end
                % make a headerless volume 
                roivol = Volume(result,'metafeatures',...
                    metafeatures,'metasamples',metasamples);
            end
            % drop nan ROIs
            anynan = ~any(isnan(roivol.data),1);
            roivol = roivol(:,anynan);
        end
    end
end
