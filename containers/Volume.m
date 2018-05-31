% Base class for data storage.
%
% This class overloads Matlab's native concatenation
% functionality so it is possible to do e.g. volboth = [vola; volb] to
% obtain a sensibly combined instance.
%
% Similarly, direct indexing is overloaded so that e.g.
% vol(1,:) returns a volume consisting only of the first data point.
%
% INPUTS:
% data              matrix in samples by columns format
%
% NAMED INPUTS:
% frameperiod   scalar defining e.g. TR.
% meta          struct with samples and features fields
%   or for shorthand, use these instead:
% metasamples   struct of numeric or cell arrays with nsamples by 1 shape
% metafeatures  struct numeric or cell arrays with 1 by nfeatures shape
%
% vol = Volume(data,varargin)
classdef Volume < handle
    properties
        data % nsamples by nfeatures matrix
        nsamples % number of sample points (rows in data)
        nfeatures % number of in-mask voxels (columns in data)
        meta = struct('samples',[],'features',[]); % meta inputs go here
        desc = struct('samples',[],'features',[]); % summaries of meta
        standardstruct = struct('labels',[],'chunks',[],'names',[],'order',[]);
        frameperiod % recording rate (s)
    end

    methods
        function vol = Volume(data,varargin)
            initialisemeta(vol,varargin{:});
            if ieNotDefined('data')
                data = {};
            end
            if ~iscell(data)
                data = {data};
            end
            for d = 1:length(data)
                thisdata = data{d};
                if isa(thisdata,'Volume')
                    datamat = thisdata.data;
                    % update meta
                    vol.meta = vol.updatemeta(vol.meta,thisdata.meta);
                    vol.frameperiod = setifunset(vol.frameperiod,...
                        thisdata.frameperiod);
                else
                    % assume raw matrix
                    datamat = thisdata;
                end
                if isempty(datamat)
                    continue
                end
                % parse the data
                tdsize = size(datamat);
                dtype = class(datamat);
                % initialise with class of first data. Helps conserve
                % memory if you are using a lower precision datatype
                % (otherwise the concatenation operation upcasts all data
                % to double)
                if isempty(vol.data)
                    vol.data = feval(dtype,[]);
                end
                % it had best be a n by nfeatures matrix (if we
                % know nfeatures)
                vol.nfeatures = setifunset(vol.nfeatures,tdsize(2),true);
                vol.data = [vol.data; datamat];
            end
            vol.nsamples = size(vol.data,1);
            % analyse the standard descriptives
            checkmeta(vol);
        end

        function self = initialisemeta(self,varargin)
        % handle the basic parsing of named meta arguments. Used as part of
        % instancing.
        %
        % self = initialisemeta(self,varargin)
            args = varargparse(varargin,...
                struct('metasamples',self.standardstruct,...
                'metafeatures',self.standardstruct,...
                'meta',[],...
                'frameperiod',[]),true);
            self.frameperiod = args.frameperiod;
            % insure meta samples and features contain mandatory
            % fields from self.standardstruct
            self.meta.samples = self.standardstruct;
            self.meta.features = self.standardstruct;
            % update with metafeatures/metasamples
            self.meta.samples = catstruct(self.meta.samples,args.metasamples);
            self.meta.features = catstruct(self.meta.features,...
                args.metafeatures);
            % update with meta
            if ~isempty(args.meta)
                self.meta.samples = catstruct(self.meta.samples,...
                    args.meta.samples);
                self.meta.features = catstruct(self.meta.features,...
                    args.meta.features);
            end
            % insure that meta samples are in rows and meta features in
            % columns
            self.meta.samples = self.imposestructfieldshape(...
                self.meta.samples,'column');
            self.meta.features = self.imposestructfieldshape(...
                self.meta.features,'row');
        end

        function filterbychunk(self,fname,varargin)
        % general method for fevaling some fname which takes sample data
        % from each chunk as a first input and anything else in the
        % following.
        %
        % filterbychunk(self,fname,[varargin])
            for c = 1:self.desc.samples.nunique.chunks
                chunkind = self.meta.samples.chunks == ...
                    self.desc.samples.unique.chunks(c);
                self.data(chunkind,:) = feval(fname,...
                    self.data(chunkind,:),varargin{:});
            end
        end

        function sgdetrend(self,k,f)
        % use a Savitzky-Golay filter to detrend the data across samples.
        % Operates separately on each chunk
        %
        % sgdetrend(k,f)
            filterbychunk(self,@(data,k,f)data-sgolayfilt(data,k,f,[],...
                1),k,f);
        end

        function sgfilter(self,k,f)
        % Apply a Savitzky-Golay filter to the data across samples
        %
        % sgfilter(k,f)
            filterbychunk(self,'sgolayfilt',k,f,[],1);
        end

        function checkmeta(self)
        % add desc data for each meta field, check that meta data is in
        % register.
        %
        % checkmeta(self)
            for fn = fieldnames(self.standardstruct)'
                fnstr = fn{1};
                % get standard fields - unique, indices, nunique
                [self.desc.samples.unique.(fnstr),...
                    junk,self.desc.samples.inds.(fnstr)] = ...
                    unique(self.meta.samples.(fnstr));
                self.desc.samples.nunique.(fnstr) = length(...
                    self.desc.samples.unique.(fnstr));
                [self.desc.features.unique.(fnstr),...
                    junk,self.desc.features.inds.(fnstr)] = ...
                    unique(self.meta.features.(fnstr));
                self.desc.features.nunique.(fnstr) = length(...
                    self.desc.features.unique.(fnstr));
            end
            % check that the meta samples/features are in register with
            % the data matrix
            assert(all(structfun(@(x)isempty(x) || ...
                (length(x)==self.nsamples),self.meta.samples)),...
                'mismatch between nsamples and meta.samples length');
            assert(all(structfun(@(x)isempty(x) || ...
                (length(x)==self.nfeatures),self.meta.features)),...
                'mismatch between nfeatures and meta.features length');
            if isempty(self.meta.samples.order)
                self.meta.samples.order = (1:self.nsamples)';
            end
        end

        function [sampind,featind] = findbymeta(self,varargin)
        % [sampind,featind] = findbymeta(self,varargin)
        % return logical indices into samples and features according to
        % some meta data. Mainly used internally to support selectbymeta
        % and removebymeta.
        %
        % [sampind,featind] = findbymeta(self,varargin)
            % get possible names
            validsamp = fieldnames(self.meta.samples);
            validfeat = fieldnames(self.meta.features);
            valid = union(validsamp,validfeat);
            [~,args] = varargparse(varargin,struct,true);
            used = intersect(valid,fieldnames(args));
            % begin by assuming you want everything
            sampind = true(self.nsamples,1);
            featind = true(1,self.nfeatures);
            for v = used'
                vstr = v{1};
                value = args.(vstr);
                % check that the field is in samples and is not empty
                insamp = isfield(self.meta.samples,vstr) && ...
                    ~isempty(self.meta.samples.(vstr));
                infeat = isfield(self.meta.features,vstr) && ...
                    ~isempty(self.meta.features.(vstr));
                % if you ask for meta data we can't find this is almost
                % certainly a problem
                assert(any([insamp infeat]),...
                        'meta data does not exist: %s',vstr)
                if insamp
                    inds = self.findbyvalue(self.meta.samples.(vstr),...
                        value);
                    sampind = sampind & inds;
                end
                if infeat
                    inds = self.findbyvalue(self.meta.features.(vstr),...
                        value);
                    featind = featind & inds;
                end
            end
        end

        function vol = selectbymeta(self,varargin)
        % convenience method for returning a volume with particular meta
        % data in samples/features dimension. Assumes that samples and
        % features are in register. You can pass any number of criteria and
        % any number of values for each criterion.
        % The resulting vol will be the intersection of all criteria.
        % vol = selectbymeta(self,varargin)
            [sampind,featind] = findbymeta(self,varargin{:});
            % convert to subsref format
            s = struct('type','()','subs',{{sampind,featind}});
            [dat,meta] = basesubsref(self,s);
            % finally, use the appropriate, potentially sub-classed copy
            % method.
            vol = copy(self,dat,meta);
        end

        function vol = removebymeta(self,varargin)
        % convenience method for returning a volume WITHOUT particular meta
        % data in samples/features dimension. Otherwise similar to
        % selectbymeta.
        % vol = removebymeta(self,varargin)
            [sampind,featind] = findbymeta(self,varargin{:});
            if ~all(sampind)
                sampind = ~sampind;
            end
            if ~all(featind)
                featind = ~featind;
            end
            % convert to subsref format
            s = struct('type','()','subs',{{sampind,featind}});
            [dat,meta] = basesubsref(self,s);
            % finally, use the appropriate, potentially sub-classed copy
            % method.
            vol = copy(self,dat,meta);
        end

        function vol = copy(self,dat,meta)
            vol = feval(class(self),dat,'meta',meta,...
                'frameperiod',self.frameperiod);
        end

        function varargout = subsref(self,s)
        % varargout = subsref(self,s)
            switch s(1).type
                case '()'
                    % basic parsing
                    [dat,meta] = basesubsref(self,s);
                    % make a new instance
                    varargout{1} = Volume(dat,'meta',meta,...
                        'frameperiod',self.frameperiod);
                otherwise
                    % revert to builtin behaviour
                    try
                        varargout{1} = builtin('subsref',self,s);
                    catch
                        % maybe no output?
                        builtin('subsref',self,s);
                    end
            end
        end

        function out = end(self,k,n)
        % override 'end' operator to get end in self.data rather than
        % self when doing e.g. v2 = vol(end,:);
        %
        % out = end(self,k,n)
            out = feval('end',self.data,k,n);
        end

        function [dat,meta] = basesubsref(self,s)
        % basic subsref behaviour. Should be called in subclass subsref.
        % It's a bit baffling, but apparently these must be dynamic
        % methods.
        %
        % [dat,meta] = basesubsref(self,s)
            assert(length(s)==1,...
                'cannot subindex instance');
            % if no column index, assume you want all features
            if length(s.subs)==1
                s.subs{2} = 1:self.nfeatures;
            end
            % pull out data field 
            dat = builtin('subsref',self.data,s);
            % handle special cases
            datind = s.subs{1};
            if strcmp(datind,':')
                datind = 1:self.nsamples;
            end
            % update meta data
            meta.samples = self.indexstructfields(self.meta.samples,datind);
            % if we are also indexing in feature dimension
            if length(s.subs) > 1 
                % update meta features
                meta.features = self.indexstructfields(self.meta.features,...
                    s.subs{2});
            end
        end

        function o = horzcat(varargin)
        % horizontal concatenation of volume instance (nb not supported for
        % all sub-classes).
        %
        % o = horzcat(varargin)
            o = [];
            if ~nargin
                return
            end
            ref = varargin{1};
            % assume same samples in each volume
            samples = ref.meta.samples;
            features = ref.meta.features;
            data = ref.data;
            frameperiod = ref.frameperiod;
            for n = 2:nargin
                % use parsemeta to check that samples are identical in the
                % to-be-concatenated data
                samples2 = ref.parsemeta(samples,...
                    varargin{n}.meta.samples,@ref.testmeta);
                assert(isequal(frameperiod,varargin{n}.frameperiod),...
                    'mismatched frameperiods');
                data = [data varargin{n}.data];
                features = ref.parsemeta(features,...
                    varargin{n}.meta.features,@horzcat);
            end
            o = feval(class(ref),data,'metasamples',samples,...
                'metafeatures',features,'frameperiod',...
                ref.frameperiod);
        end

        function o = vertcat(varargin)
        % Vertical concatenation of volume instances.
        % 
        % o = vertcat(varargin)
            % this is needed to insure that the correct sub-class is called
            cl = class(varargin{1});
            o = feval(cl,varargin);
        end

        function o = cat(dim,varargin)
        % concatenation of volume instances.
        %
        % o = cat(dim,varargin)
            switch dim
                case 1
                    o = vertcat(varargin{:});
                case 2
                    o = horzcat(varargin{:});
                otherwise
                    error('cannot concatenate in dim %d',dim);
            end
        end

    end

    methods(Static = true)
        function org = testmeta(org,new)
        % org = testmeta(org,new)
            assert(isequal(org,new),'mismatched features!')
        end

        function inds = findbyvalue(array,item)
        % inds = findbyvalue(array,item)

            if ischar(item) || (iscell(item) && ischar(item{1}))
                assert(iscell(array) && isa(array{1},'char'),...
                    'array must be cell with char entries');
                if ~iscell(item)
                    item = {item};
                end
                inds = strcmp(item{1},array);
                % support OR behaviour for multiple inputs
                for x = 2:length(item)
                    inds = strcmp(item{x},array) | inds;
                end
            elseif isnumeric(item)
                assert(isnumeric(array),'array must be numeric');
                inds = item(1) == array;
                % support OR behaviour for multiple inputs
                for x = 2:length(item)
                    inds = (item(x) == array) | inds;
                end
            else
                error('unrecognised input class: %s',class(item));
            end
        end

        function out = indexstructfields(in,inds)
        % Return a struct where each field has been indexed with the same
        % inds. Used internally to insure that any custom meta data/features
        % stay in register when the data gets indexed.
        % out = indexstructfields(in,inds)
            sfields = fieldnames(in);
            % special treatment of empties since we can't index these
            empties = structfun(@isempty,in);
            % get data for the remaining
            validdata = structfun(@(x)x(inds),...
                rmfield(in,sfields(empties)),'uniformoutput',false);
            % return the empties (with correct type whatever that may be)
            emptydata = structfun(@(x)x,rmfield(in,sfields(~empties)),...
                'uniformoutput',false);
            % recombine for the final struct
            out = catstruct(validdata, emptydata);
        end

        function in = imposestructfieldshape(in,direction)
            switch direction
                case 'column'
                    flattener = @(x)x(:);
                case 'row'
                    flattener = @(x)x(:)';
                otherwise
                    error('unknown direction: %s',direction);
            end
            for fn = fieldnames(in)'
                fnstr = fn{1};
                assert(ndims(in.(fnstr))<3,...
                    'input %s is >2d',fnstr);
                assert(~ismat(in.(fnstr)),...
                    'input must be row or column vector');
                in.(fnstr) = flattener(in.(fnstr));
            end
        end

        function org = parsemeta(org,new,metaoperate)
        % internal method for iterating recursively over data in nested
        % struct arrays. metaoperate is a function handle that gets applied
        % to data. 
        %
        % org = parsemeta(org,new,metaoperate)
            orgfields = fieldnames(org);
            newfields = fieldnames(new);
            % fields that need updating
            toupdate = intersect(orgfields,newfields)';
            % concatenate when needed
            for t = toupdate
                tstr = t{1};
                if isstruct(org.(tstr))
                    % recurse into sub-fields
                    org.(tstr) = parsemeta(org.(tstr),...
                        new.(tstr),metaoperate);
                else
                    % if we don't have it, assume we need it
                    if isempty(org.(tstr))
                        org.(tstr) = new.(tstr);
                    elseif ~isempty(new.(tstr))
                        % if there is some new data to process
                        org.(tstr) = metaoperate(org.(tstr),new.(tstr));
                    end
                end
            end
            % fields that need to be added
            trulynew = setdiff(newfields,orgfields)';
            for n = trulynew
                nstr = n{1};
                org.(nstr) = new.(nstr);
            end
        end

        function new = updatemeta(org,new)
        % new = updatemeta(org,new)
            new.samples = Volume.parsemeta(org.samples,new.samples,...
                @vertcat);
            new.features = Volume.parsemeta(org.features,new.features,...
                @Volume.testmeta);
        end
    end
end
