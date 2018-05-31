% SPM-based header volume
%
% SPMVolume < MriVolume
classdef SPMVolume < MriVolume

    methods
        function spmvol = SPMVolume(data,mask,varargin)
            % handle non-base inputs (varargin for sub-classes go into
            % remargs)
            [args,remargs] = varargparse(varargin,struct('header',[]),true);
            % we don't allow a voxsize input since there are no real cases
            % where it would be a smart idea to contradict the header
            voxsize = [];
            temp.frameperiod = [];
            temp.meta = struct('samples',struct,'features',struct);

            % parse mask input
            if ieNotDefined('mask')
                mask = [];
            end
            if isa(mask,'char')
                % assume it's a SPM-compatible volume
                args.header = setifunset(args.header,spm_vol(mask),true);
                mask = spm_read_vols(args.header) ~= 0;
            elseif isa(mask,'MriVolume')
                % interesting!
                % but here we need to be careful not to include data
                args.header = setifunset(args.header,mask.header,true);
                voxsize = mask.voxsize;
                mask = mask.mask;
            elseif ~isempty(mask)
                mask = mask ~= 0;
            end

            % parse data input
            if ieNotDefined('data')
                data = {};
            end
            if ~iscell(data)
                data = {data};
            end
            filenames = {};
            chunks = [];
            outdata = [];
            for d = 1:length(data)
                thisdata = data{d};
                if isa(thisdata,'char')
                    % assume it's a SPM-compatible volume
                    dV = spm_vol(thisdata);
                    % check that you aren't trying to combine data with
                    % a weird mask
                    if ~isempty(mask)
                        % make sure the headers match
                        assert(spm_check_orientations(...
                            [args.header; dV(:)]),'header mismatch');
                    else
                        % set header by first data volume instead of mask
                        args.header = setifunset(args.header,dV(1),true);
                    end
                    % if the mask remains undefined so far, initialise an
                    % all-on mask based on data size
                    if isempty(mask)
                        mask = true(size(spm_read_vols(dV(1))));
                    end
                    datamat = loadmaskedvolumes(dV,mask);
                    filenames = [filenames; {dV.fname}'];
                elseif isa(thisdata,'SPMVolume')
                    % interesting! 
                    if ~isempty(args.header)
                        assert(spm_check_orientations(...
                            [thisdata.header; args.header]),...
                            'data header does not match mask');
                    end
                    datamat = thisdata.data;
                    temp.meta = Volume.updatemeta(temp.meta,thisdata.meta);
                    % also bring along the mask if possible and needed
                    mask = setifunset(mask,thisdata.mask,true);
                    args.header = setifunset(args.header,thisdata.header);
                    temp.frameperiod = setifunset(temp.frameperiod,...
                        thisdata.frameperiod,true);
                elseif isstruct(thisdata) && isfield(thisdata,'SPMid')
                    % SPM.mat
                    SPM = thisdata;
                    temp.frameperiod = setifunset(temp.frameperiod,...
                        SPM.xY.RT,true);
                    args.header = setifunset(args.header,spm_vol(SPM.xY.P(1,:)));
                    if isempty(mask)
                        mask = true(args.header.dim);
                    end
                    % load in-mask data 
                    % the first pre-call is a slightly hacky way to ensure
                    % that the mask is as small as it can be before we
                    % embark on the longer load of all the (in-mask) voxel
                    % timecourses.
                    [~,mask] = loadmaskedvolumes(SPM.xY.P(1,:),mask);
                    [datamat,mask] = loadmaskedvolumes(SPM.xY.P,mask);
                    % determine chunks
                    maxchunk = max([chunks;0]);
                    thischunk = cell2mat(arrayfun(@(x)ones(...
                        SPM.nscan(x),1)*x,...
                        (maxchunk+1:maxchunk+numel(SPM.Sess))',...
                        'uniformoutput',false));
                    chunks = [chunks; thischunk];
                    filenames = [filenames; {SPM.xY.VY.fname}'];
                else
                    % assume headerless data
                    datamat = thisdata;
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
            % use SPM header (which I hope we have by now) to figure out voxel size
            voxsize = setifunset(voxsize,vox2mm(args.header),true);

            if ~isempty(filenames)
                temp.meta.samples.filenames = filenames;
            end

            if ~isempty(chunks)
                temp.meta.samples.chunks = chunks;
            end

            % note that we append the additional arguments needed for
            % MriVolume
            remargs.header = args.header;
            if ~isfield(remargs,'voxsize')
                remargs.voxsize = voxsize;
            end
            % soft set - so if remargs existed but was empty, update
            remargs.voxsize = setifunset(remargs.voxsize,voxsize);
            if ~isempty([fieldnames(temp.meta.features); fieldnames(temp.meta.samples)])
                % only add meta input if we put anything in it. (this is to
                % hack around the current lack of support for entering both
                % meta and metasamples/metafeatures)
                remargs.meta = temp.meta;
            end
            remargs = structfields2varargs(remargs);
            spmvol = spmvol@MriVolume(outdata,mask,remargs{:});

            % frameperiod is set at Volume level, but here we just check
            % that it matches our estimate
            spmvol.frameperiod = setifunset(spmvol.frameperiod,...
                temp.frameperiod,true);
            % make sure that volume header is appropriate for data
            % precision
            if isfield(spmvol.header,'dt')
                spmvol.header.dt = [matlabclass2spmnifti(spmvol.data),...
                    spm_platform('bigend')];
            end
        end
    end
end
