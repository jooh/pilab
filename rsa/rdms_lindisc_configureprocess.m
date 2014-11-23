function glman_rdm = rdms_lindisc_configureprocess(varargin)

getArgs(varargin,{'glmvarargs',{},'cvsplit',[],...
    'glmclass','GLM','sterrunits',1,'crossvalidate',true,...
    'crosscon',[],'ncon',[],'setclass','double','demean',false});

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
    % then processors
    rdmcross(1) = GLMProcessor('cvcrossclassificationrun',[],1,...
        'discriminant',testmeth,[np 1],crossconout{1},crossconout{2});
    rdmcross(2) = GLMProcessor('cvcrossclassificationrun',[],1,...
        'discriminant',testmeth,[np 1],crossconout{2},crossconout{1});
    glman_rdm = GLMMetaProcessor(glmspec,rdmcross,@(x)mean(x,3));
end

if demean
    % one more layer of wrapping...
    glman_rdm = ROIPreProcessor(glman_rdm,@(data)bsxfun(@minus,data,...
        mean(data,2)));
end
