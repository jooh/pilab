% helper function for roidata2rdmvol_lindisc_batch to set up some thornier
% aspects of LDC. Separated out to support other functions.
%
% [glman_rdm,cons] = rdms_lindisc_configureprocess(varargin)
function [glman_rdm,cons] = rdms_lindisc_configureprocess(varargin)

arg = varargparse(varargin,struct('glmvarargs',{{}},'cvsplit',[],...
    'glmclass','GLM','sterrunits',0,'crossvalidate',true,...
    'crosscon',[],'ncon',[],'setclass','double','demean',false,...
    'nperms',0,'onewayvalidation',false,'returnp',false,...
    'ntrials',NaN,'permmeth','permuteruns'));

if ~iscell(arg.glmvarargs)
    if isempty(arg.glmvarargs)
        arg.glmvarargs = {};
    else
        arg.glmvarargs = {arg.glmvarargs};
    end
end

if ischar(arg.cvsplit)
    arg.cvsplit = feval(arg.cvsplit);
end

if ischar(arg.crosscon)
    arg.crosscon = eval(arg.crosscon);
elseif iscell(arg.crosscon)
    for c = 1:numel(arg.crosscon)
        if ischar(arg.crosscon{c})
            arg.crosscon{c} = eval(arg.crosscon{c});
        end
    end
end

nreturn = 1;
if arg.sterrunits
    testmeth = 'infot';
else
    testmeth = 'infoc';
end
if arg.returnp
    nreturn = 2;
end
logstr('testmeth: %s\n',testmeth);

% assemble processors
glmspec = GLMConstructor(arg.glmclass,arg.cvsplit,arg.glmvarargs{:});

if strcmp(arg.glmclass,'ConvGLM')
    glmspec = ConvGLMConstructor(arg.cvsplit);
    arg.permmeth = 'permutesamples';
end

nrun = numel(arg.cvsplit);
if ~isscalar(arg.nperms)
    pind = arg.nperms;
    arg.nperms = size(pind,1);
else
    if arg.nperms > 1
        switch arg.permmeth
            case 'permuteruns'
                assert(~isempty(arg.cvsplit),['cvsplit must be defined for ' ...
                    'permutation test']);
                % bit hacky but there's no telling how many runs we will have
                % at this stage.
                % Note that this solution isn't exactly bullet proof.
                pind = permuteindices(nrun,arg.nperms);
            case 'permutesamples'
                % hm.
                pind = permuteindices(arg.ntrials,arg.nperms);
            otherwise
                error('unknown permmeth: %s',arg.permmeth)
        end
    end
end

if isempty(arg.crosscon)
    % straight RDM
    % contrast vector of correct class
    cons = allpairwisecontrasts(feval(arg.setclass,...
        arg.ncon));
    if arg.crossvalidate
        if arg.onewayvalidation
            cvmeth = 'validatedclassification';
            % NB, the [false true] bit makes the strong assumption that
            % you will have only 2 runs, where the first is train and the
            % second test.
            assert(nrun==2,['only 2 runs for validatedclassification ' ...
                'first train, second test']);
            args = {[false,true],'discriminant',{cons},testmeth,{cons},...
                '*SELF*'};
        else
            cvmeth = 'cvclassificationrun';
            args = {'discriminant',testmeth,'*SELF*',cons};
        end
        if arg.nperms < 2
            rdm = GLMProcessor(cvmeth,[],nreturn,args{:});
        else
            rdm = GLMProcessor(arg.permmeth,[],nreturn,pind,cvmeth,[],...
                args{:});
        end
        % NB no longer any need for mean operation here.
        glman_rdm = GLMMetaProcessor(glmspec,rdm,[]);
    else
        error('currently unsupported');
    end
else
    % cross-class RDM - first set up contrast vectors
    assert(arg.crossvalidate==1,'must crossvalidate if crosscon are present');
    nc = numel(arg.crosscon{1});
    assert(nc==numel(arg.crosscon{2}),'mismatched crosscon');
    assert(nc<=(arg.ncon/2),...
        'bad crosscon for ncon');
    np = nchoosek(nc,2);
    cons = allpairwisecontrasts(feval(arg.setclass,nc));
    fullmat = zeros(np,arg.ncon,class(cons));
    crossconout{1} = fullmat;
    crossconout{1}(:,arg.crosscon{1}) = cons;
    crossconout{2} = fullmat;
    crossconout{2}(:,arg.crosscon{2}) = cons;
    % reassign cons to ensure labelling is approximately correct (will
    % capture one of the two cross-decoders that were actually run)
    cons = crossconout{1};
    basearg = {[],nreturn};
    if arg.onewayvalidation
        cvmeth = 'validatedclassification';
        assert(nrun==2,['only 2 runs for validatedclassification ' ...
            'first train, second test']);
        arga = {[false,true],'discriminant',crossconout(1),testmeth,...
            crossconout(2),'*SELF*'};
        argb = {[false,true],'discriminant',crossconout(2),testmeth,...
            crossconout(1),'*SELF*'};
    else
        cvmeth = 'cvcrossclassificationrun';
        arga = {'discriminant',testmeth,'*SELF*',crossconout{1},...
            crossconout{2}};
        argb = {'discriminant',testmeth,'*SELF*',crossconout{2},...
            crossconout{1}};
    end
    if arg.nperms < 2
        rdmcross(1) = GLMProcessor(cvmeth,basearg{:},arga{:});
        rdmcross(2) = GLMProcessor(cvmeth,basearg{:},argb{:});
    else
        rdmcross(1) = GLMProcessor(arg.permmeth,basearg{:},pind,...
            cvmeth,[],arga{:});
        rdmcross(2) = GLMProcessor(arg.permmeth,basearg{:},pind,...
            cvmeth,[],argb{:});
    end
    % still necessary? A little bit unsure about this.
    glman_rdm = GLMMetaProcessor(glmspec,rdmcross,[]);
end

if arg.demean
    % one more layer of wrapping...
    glman_rdm = ROIPreProcessor(glman_rdm,@(data)bsxfun(@minus,data,...
        mean(data,2)));
end
