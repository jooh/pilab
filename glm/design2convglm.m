% model = design2convglm(design,datamat,chunks)
function model = design2convglm(design,datamat,chunks)

% match up design and data by chunk
uc = unique(chunks);
dc = arrayfun(@(thisd)thisd.chunk(1),design);
assert(isequal(uc,sort(dc)),'mismatched chunks in design');
for c = 1:numel(uc)
    dind = find(dc==uc(c));
    assert(design(dind).n == sum(chunks==uc(c)),...
        'design/data chunk mismatch')
    model(c) = ConvGLM(design(dind).onsets,design(dind).conind,...
        design(dind).frameperiod,datamat(chunks==uc(c),:),...
        design(dind).covariates,design(dind).convargs{:});
end
