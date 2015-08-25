% conind = findconind(connames,targets)
function conind = findconind(connames,targets)

[~,conind] = intersect(connames,targets,'stable');
if isempty(conind)
    % maybe prefix from searchlight module
    prenames = cellfun(@(thiscon)['rsa_r_' thiscon],targets,...
        'uniformoutput',false);
    [~,conind] = intersect(connames,prenames,'stable');
    if isempty(conind)
        % maybe a hyphen got converted to underscore. This is mainly for
        % legacy support (stripbadcharacters no longer removes hyphens).
        [~,conind] = intersect(connames,strrep(prenames,'-','_'),'stable');
    end
end
assert(~isempty(conind),'did not find %s',targets);
