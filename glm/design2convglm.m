% model = design2convglm(design,datamat,chunks)
function model = design2convglm(design,datamat,chunks)

uc = unique(chunks);
for c = 1:numel(uc)
    % for now make very strong assumptions here about the chunks being 1:n
    assert(design(c).chunk(1)==c,'non-standard chunk mapping')
    assert(design(c).n == sum(chunks==uc(c)),...
        'design/data chunk mismatch')
    model(c) = ConvGLM(design(c).onsets,design(c).conind,...
        design(c).frameperiod,datamat(chunks==uc(c),:),...
        design(c).covariates,design(c).convargs{:});
end
