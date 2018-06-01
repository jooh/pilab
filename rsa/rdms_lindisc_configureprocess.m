function [glman_rdm,cons] = rdms_lindisc_configureprocess(varargin)

getArgs(varargin,{'glmvarargs',{},'cvsplit',[],...
    'glmclass','GLM','sterrunits',0,'crossvalidate',true,...
    'crosscon',[],'ncon',[],'setclass','double','demean',false,...
    'nperms',0,'onewayvalidation',false,'returnp',false,...
    'ntrials',NaN,'permmeth','permuteruns'});

if ~iscell(glmvarargs)
    if isempty(glmvarargs)
        glmvarargs = {};
    else
        glmvarargs = {glmvarargs};
    end
end

if ischar(cvsplit)
    cvsplit = eval(cvsplit);
end

if ischar(crosscon)
    crosscon = eval(crosscon);
elseif iscell(crosscon)
    for c = 1:numel(crosscon)
        if ischar(crosscon{c})
            crosscon{c} = eval(crosscon{c});
        end
    end
end

nreturn = 1;
if sterrunits
    testmeth = 'infot';
else
    testmeth = 'infoc';
end
if returnp
    nreturn = 2;
end
logstr('testmeth: %s\n',testmeth);

% assemble processors
glmspec = GLMConstructor(glmclass,cvsplit,glmvarargs{:});

if strcmp(glmclass,'ConvGLM')
    glmspec = ConvGLMConstructor(cvsplit);
    permmeth = 'permutesamples';
end

nrun = numel(cvsplit);
if ~isscalar(nperms)
    pind = nperms;
    nperms = size(pind,1);
else
    if nperms > 1
        switch permmeth
            case 'permuteruns'
                assert(~isempty(cvsplit),['cvsplit must be defined for ' ...
                    'permutation test']);
                % bit hacky but there's no telling how many runs we will have
                % at this stage.
                % Note that this solution isn't exactly bullet proof.
                pind = permuteindices(nrun,nperms);
            case 'permutesamples'
                % hm.
                pind = permuteindices(ntrials,nperms);
            otherwise
                error('unknown permmeth: %s',permmeth)
        end
    end
end

if isempty(crosscon)
    % straight RDM
    % contrast vector of correct class
    cons = allpairwisecontrasts(feval(setclass,...
        ncon));
    if crossvalidate
        if onewayvalidation
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
        if nperms < 2
            rdm = GLMProcessor(cvmeth,[],nreturn,args{:});
        else
            rdm = GLMProcessor(permmeth,[],nreturn,pind,cvmeth,[],...
                args{:});
        end
        % NB no longer any need for mean operation here.
        glman_rdm = GLMMetaProcessor(glmspec,rdm,[]);
    else
        error('currently unsupported');
    end
else
    % cross-class RDM - first set up contrast vectors
    assert(crossvalidate==1,'must crossvalidate if crosscon are present');
    nc = numel(crosscon{1});
    assert(nc==numel(crosscon{2}),'mismatched crosscon');
    assert(nc<=(ncon/2),...
        'bad crosscon for ncon');
    np = nchoosek(nc,2);
    cons = allpairwisecontrasts(feval(setclass,nc));
    fullmat = zeros(np,ncon,class(cons));
    crossconout{1} = fullmat;
    crossconout{1}(:,crosscon{1}) = cons;
    crossconout{2} = fullmat;
    crossconout{2}(:,crosscon{2}) = cons;
    % reassign cons to ensure labelling is approximately correct (will
    % capture one of the two cross-decoders that were actually run)
    cons = crossconout{1};
    basearg = {[],nreturn};
    if onewayvalidation
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
    if nperms < 2
        rdmcross(1) = GLMProcessor(cvmeth,basearg{:},arga{:});
        rdmcross(2) = GLMProcessor(cvmeth,basearg{:},argb{:});
    else
        rdmcross(1) = GLMProcessor(permmeth,basearg{:},pind,...
            cvmeth,[],arga{:});
        rdmcross(2) = GLMProcessor(permmeth,basearg{:},pind,...
            cvmeth,[],argb{:});
    end
    % still necessary? A little bit unsure about this.
    glman_rdm = GLMMetaProcessor(glmspec,rdmcross,[]);
end

if demean
    % one more layer of wrapping...
    glman_rdm = ROIPreProcessor(glman_rdm,@(data)bsxfun(@minus,data,...
        mean(data,2)));
end
