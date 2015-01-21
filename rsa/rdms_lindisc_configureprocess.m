function glman_rdm = rdms_lindisc_configureprocess(varargin)

getArgs(varargin,{'glmvarargs',{},'cvsplit',[],...
    'glmclass','GLM','sterrunits',0,'crossvalidate',true,...
    'crosscon',[],'ncon',[],'setclass','double','demean',false,...
    'nperms',0});

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

if sterrunits
    testmeth = 'infot';
else
    testmeth = 'infoc';
end
logstr('testmeth: %s\n',testmeth);

% assemble processors
glmspec = GLMConstructor(glmclass,cvsplit,glmvarargs{:});


if nperms > 2
    assert(~isempty(cvsplit),['cvsplit must be defined for ' ...
        'permutation test']);
    % bit hacky but there's no telling how many runs we will have
    % at this stage.
    % Note that this solution isn't exactly bullet proof.
    nrun = numel(cvsplit);
    pind = permuteindices(nrun,nperms);
end

if isempty(crosscon)
    % straight RDM
    np = nchoosek(ncon,2);
    % contrast vector of correct class
    cons = allpairwisecontrasts(feval(setclass,...
        ncon));
    if crossvalidate
        cvmeth = 'cvclassificationrun';
        args = {'discriminant',testmeth,[np 1],cons};
        if nperms < 2
            rdm = GLMProcessor(cvmeth,[],1,args{:});
        else
            rdm = GLMProcessor('permuteruns',[],1,pind,cvmeth,[],args{:});
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
    cvmeth = 'cvcrossclassificationrun';
    basearg = {[],1};
    arga = {'discriminant',testmeth,[np 1],crossconout{1},crossconout{2}};
    argb = {'discriminant',testmeth,[np 1],crossconout{2},crossconout{1}};
    if nperms < 2
        rdmcross(1) = GLMProcessor(cvmeth,basearg{:},arga{:});
        rdmcross(2) = GLMProcessor(cvmeth,basearg{:},argb{:});
    else
        rdmcross(1) = GLMProcessor('permuteruns',basearg{:},pind,...
            cvmeth,[],arga{:});
        rdmcross(2) = GLMProcessor('permuteruns',basearg{:},pind,...
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
