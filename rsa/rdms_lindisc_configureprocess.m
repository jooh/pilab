function glman_rdm = rdms_lindisc_configureprocess(varargin)

getArgs(varargin,{'glmvarargs',{{}},'cvsplit',[],...
    'glmclass','GLM','sterrunits',1,'crossvalidate',1,...
    'crosscon',[],'ncon',[],'setclass','double'});

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
end

if sterrunits
    testmeth = 'infotmap';
else
    testmeth = 'infomahalanobis';
end

% assemble processors
glmspec = GLMConstructor(glmclass,cvsplit,glmvarargs{:});

if isempty(crosscon)
    % straight RDM
    np = nchoosek(ncon,2);
    % contrast vector of correct class
    cons = allpairwisecontrasts(feval(setclass,...
        ncon));
    if crossvalidate
        rdm = GLMProcessor('cvclassificationrun',[],1,'discriminant',...
            testmeth,[np 1],cons);
        % only necessary if CV
        glman_rdm = GLMMetaProcessor(glmspec,rdm,@(x)mean(x,3));
    else
        error('currently unsupported');
    end
else
    % cross-class RDM - first set up contrast vectors
    assert(crossvalidate,'must crossvalidate if crosscon are present');
    nc = numel(crosscon{1});
    assert(nc==numel(crosscon{2}),'mismatched crosscon');
    assert(nc<=(ncon/2),...
        'bad crosscon for ncon');
    np = nchoosek(nc,2);
    cons = allpairwisecontrasts(feval(class(setclass),nc));
    fullmat = zeros(np,ncon,class(cons));
    crosscon{1} = fullmat;
    crosscon{1}(:,crosscon{1}) = cons;
    crosscon{2} = fullmat;
    crosscon{2}(:,crosscon{2}) = cons;
    % then processors
    rdmcross(1) = GLMProcessor('cvcrossclassificationrun',[],1,...
        'discriminant',testmeth,[np 1],crosscon{1},crosscon{2});
    rdmcross(2) = GLMProcessor('cvcrossclassificationrun',[],1,...
        'discriminant',testmeth,[np 1],crosscon{2},crosscon{1});
    glman_rdm = GLMMetaProcessor(glmspec,rdmcross,@(x)mean(x,3));
end
