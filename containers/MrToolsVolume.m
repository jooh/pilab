classdef MrToolsVolume < MriVolume
    properties
        sessiondir % location from which we start a view
    end

    methods
        function mrvol = MrToolsVolume(data,mask,varargin)
            % pull out the arguments we need to extract data
            [args,remargs] = varargparse(varargin,...
                struct('group','MotionComp','scans','all','header',[]),true);
            voxsize = [];
            temp.frameperiod = [];
            temp.meta = struct('samples',struct,'features',struct);
            temp.sessiondir = [];

            % parse data input
            if ieNotDefined('data')
                data = {};
            end
            if ieNotDefined('mask')
                mask = [];
            end
            if ~iscell(data)
                data = {data};
            end
            filenames = {};
            chunks = [];
            outdata = [];
            orgdir = pwd;
            v = [];
            for d = 1:length(data)
                processview = false;
                openedview = false;
                thisdata = data{d};
                datamat = [];
                if isa(thisdata,'MrToolsVolume')
                    % concatenation
                    temp.sessiondir = setifunset(temp.sessiondir,...
                        thisdata.sessiondir,true);
                    datamat = thisdata.data;
                    temp.meta = Volume.updatemeta(temp.meta,thisdata.meta);
                    % also bring along the mask if possible and needed
                    mask = setifunset(mask,thisdata.mask,true);
                    args.header = setifunset(args.header,thisdata.header);
                    temp.frameperiod = setifunset(temp.frameperiod,...
                        thisdata.frameperiod,true);
                elseif ischar(thisdata)
                    % path to mrTools session?
                    temp.sessiondir = setifunset(temp.sessiondir,thisdata);
                    cd(temp.sessiondir);
                    v = newView;
                    processview = true;
                elseif isview(thisdata)
                    % open view
                    v = thisdata;
                    temp.sessiondir = setifunset(temp.sessiondir,...
                        viewGet(v,'sessiondirectory'),true);
                    cd(temp.sessiondir);
                    processview = true;
                else
                    % assume headerless data
                    datamat = thisdata;
                end
                if processview
                    % let's make a few assumptions to simplify this.
                    % we will assume that you have only supplied a single
                    % data input
                    assert(length(data)==1,['if mrTools view is used ' ...
                        'only a single input is supported']);
                    % and that the mask is a char pointing to an ROI in the
                    % subjects ROI directory (or an already loaded ROI by
                    % that name)
                    assert(ischar(mask),['mask input must be non-empty '...
                        'char for mrTools view input mode']);
                    maskname = mask;
                    mask = [];

                    % configure the mrTools session
                    v = viewSet(v,'currentgroup',args.group);
                    groupnum = viewGet(v,'currentgroup');
                    if strcmp(args.scans,'all')
                        args.scans = 1:viewGet(v,'numscans',groupnum);
                    end
                    nscan = numel(args.scans);
                    temp.frameperiod = viewGet(v,'frameperiod',scans(1),...
                        groupnum);
                    datamat = [];
                    % load in one go
                    roit = loadROITSeries(v,maskname,args.scans,groupnum,...
                        'straightXform=1','keepNAN=1');
                    % loadROITSeries sometimes results in messy rows
                    logstr('\n');
                    voxsize = setifunset(voxsize,roit{1}.voxelSize,1);
                    coords = roit{1}.scanCoords;
                    dims = viewGet(v,'scandims',args.scans(1),groupnum);
                    for s = 1:nscan
                        n = size(roit{s}.tSeries,2);
                        % mrTools quietly returns no data if you enter an
                        % invalid scan number (e.g. greater than the number
                        % of available scans).
                        assert(n>0,'no tSeries data for scan %d\n',...
                            args.scans(s));
                        chunks = [chunks; ones(n,1)*args.scans(s)];
                        datamat = [datamat; roit{s}.tSeries'];
                        % add volume number for spm_vol compatibility
                        filenames = [filenames; mat2strcell((1:n)',...
                            [viewGet(v,'tseriespath',args.scans(s),groupnum) ',%d'])];
                        % check that coordinates are matched across scans
                        coords = setifunset(coords,roit{s}.scanCoords,1);
                    end
                    args.header = setifunset(args.header,...
                        spm_vol(filenames{1}),false);
                    % drop nans
                    nanind = any(isnan(datamat),1);
                    datamat(:,nanind) = [];
                    coords(:,nanind) = [];
                    if any(nanind)
                        logstr('removed %d voxels (%2.0f%% of total)\n',...
                            sum(nanind),100*sum(nanind) / numel(nanind));
                    end

                    % setup mask
                    mask = coord2mat(coords,dims);
                    % an interesting challenge here is that scanCoords are
                    % not necessarily in the correct order to successfully
                    % convert to matrix form via linear indices 
                    maskind = find(mask)';
                    [maskx,masky,maskz] = ind2sub(dims,maskind);
                    maskcoords = [maskx; masky; maskz];
                    % these coordinates are actually the same as
                    % scanCoords, just in a different order
                    assert(isequal(sortrows(maskcoords'),...
                        sortrows(coords')),'mismatched coordinates');
                    % so these indices return the columns to the correct
                    % order
                    [~,~,resortind] = intersect(maskcoords',coords',...
                        'rows','stable');
                    assert(isequal(maskcoords,coords(:,resortind)),...
                        'resort problem');
                    datamat = datamat(:,resortind);
                end
                % short-circuit if there's no data
                if isempty(datamat)
                    continue
                end

                % initialise with class of first data. Helps conserve
                % memory if you are using a lower precision datatype
                % (otherwise the concatenation operation upcasts all data
                % to double)
                if isempty(outdata)
                    outdata = feval(class(datamat),[]);
                end
                outdata = [outdata; datamat];
            end % for d = 1:length(data)

            % if the mask is still undefined this is probably not good
            if ieNotDefined('mask')
                error('MrToolsVolume could not find a valid mask');
            end

            if ~isempty(filenames)
                temp.meta.samples.filenames = filenames;
            end

            if ~isempty(chunks)
                temp.meta.samples.chunks = chunks;
            end

            remargs.header = args.header;
            if ~isfield(remargs,'voxsize')
                remargs.voxsize = voxsize;
            end
            remargs.voxsize = setifunset(remargs.voxsize,voxsize);
            if ~isempty([fieldnames(temp.meta.features); fieldnames(temp.meta.samples)])
                remargs.meta = temp.meta;
            end
            remargs = structfields2varargs(remargs);
            mrvol = mrvol@MriVolume(outdata,mask,remargs{:});

            mrvol.sessiondir = temp.sessiondir;
            % ensure that if the user specified a custom frameperiod, it
            % doesn't clash with what mrTools has
            mrvol.frameperiod = setifunset(mrvol.frameperiod,...
                temp.frameperiod,true);

            % this is a little destructive but Matlab does not support
            % a second output argument for constructors, and we don't
            % really want to store the view in the object instance.
            if ~isempty(v)
                deleteView(v);
            end
            cd(orgdir);
        end
    end
end
