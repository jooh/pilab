% batch variant of rsapermtest - does inference, makes figures, and saves
% results to disk.
%
% INPUTS:
% disvol - volume with vectorised dissimilarities
% predictors - struct array of RDM predictors
%
% NAMED INPUTS:
% outputmode - searchlight (default) or roi
% nperms - default 10000
% resdir - directory where results get dumped (default pwd)
%
% OUTPUTS:
% rpaths - path to each saved r output
% ppaths - path to each saved permutation p output
% pfwepaths - path to each saved FWE-corrected permutation p output
% nulldistpaths - path to each null distribution mat (big files)
%
% [rpaths,ppaths,pfwepaths,nulldistpaths] = rsapermtest_batch(disvol,predictors,varargin)
function [rpaths,ppaths,pfwepaths,nulldistpaths] = rsapermtest_batch(disvol,predictors,varargin)

getArgs(varargin,{'outputmode','searchlight','nperms',10000,'resdir',pwd});

npredictors = length(predictors);

% check that parfor is available
if ~matlabpool('size')
    try
        matlabpool local
    catch
        warning('no matlabpool available')
    end
end

% make outputs
mkdirifneeded(resdir);
% pre-allocate depending on run mode
switch outputmode
    case 'searchlight'
        rpaths = {};
        ppaths = {};
        pfwepaths = {};
    case 'roi'
        % one output spanning all ROIs / predictors
        rpaths = fullfile(resdir,'roi_r.mat');
        ppaths = fullfile(resdir,'roi_p.mat');
        pfwepaths = fullfile(resdir,'roi_pFWE.mat');
        % pre-allocate results - struct array with one entry per
        % ROI and one field per analysis. Separate struct arrays
        % for each result type.
        % (a 2-level struct array might be more intuitive but
        % Matlab makes indexing such very painful so we use this
        % slightly more awkward syntax for ease of access)
        predres = cell2struct(cell([npredictors 1]),...
            [{predictors.name}]);
        res_r = repmat(predres,[1 disvol.nfeatures]);
        names = disvol.meta.features.names;
        res_p = res_r;
        res_pfwe = res_r;
    otherwise
        error('unrecognised outputmode setting: %s',outputmode)
end
% null dists are too big to go in one mat
nulldistpaths = {};
figdir = fullfile(resdir,'figures');
mkdirifneeded(figdir);

for pre = 1:npredictors
    fprintf('testing %s (%d of %d)...\n',predictors(pre).name,...
        pre,npredictors);
    % permutation test
    tic;
    [r,p,nulldists] = rsapermtest(predictors(pre),disvol.data,nperms);
    fprintf('finished in %s. ',seconds2str(toc));
    if isnan(p)
        pfwe = NaN(size(p));
    else
        % obtain Nichols / Holmes-style FWE-corrected p values
        pfwe = maxstatpfwe(nulldists);
        fprintf('Min p(FWE) = %.2f\n',min(pfwe(:)));
    end

    % save data depending on mode
    switch outputmode
        case 'searchlight'
            % write out niftis of maps
            % r map
            rout = fullfile(resdir,sprintf('%s_r.nii',...
                predictors(pre).name));
            disvol.data2file(r,rout);
            rpaths{end+1} = rout;
            if ~isnan(p)
                % log10 p map
                pout = fullfile(resdir,sprintf('%s_-log10p.nii',...
                    predictors(pre).name));
                disvol.data2file(-log10(p),pout);
                ppaths{end+1} = pout;
                % FWE-corrected p map
                pfweout = fullfile(resdir,sprintf('%s_-log10pFWE.nii',...
                    predictors(pre).name));
                disvol.data2file(-log10(pfwe),pfweout);
                pfwepaths{end+1} = pfweout;
            end
            % diagnostic figure
            F = figure;
            imagesc(makeimagestack(disvol.data2mat(r),[0 .5],1),...
                [0 1]);
            colormap(hot(1024));
            set(gca,'dataaspectratio',[1 1 1]);
            title(stripbadcharacters(predictors(pre).name,' '));
            C = colorbar;
            ylabel(C,'rho');
            set(C,'ytick',[0 1],'ylim',[0 1],'yticklabel',[0 .5]);
            axis off
            printstandard(fullfile(figdir,sprintf(...
                        'slices_r_%s',predictors(pre).name)));
            close(F);
        case 'roi'
            % update structs
            % (matlab makes this comically awkward but there you
            % go)
            rc = num2cell(r);
            [res_r.(predictors(pre).name)] = rc{:};
            pc = num2cell(p);
            [res_p.(predictors(pre).name)] = pc{:};
            pfwec = num2cell(pfwe);
            [res_pfwe.(predictors(pre).name)] = pfwec{:};
        otherwise
            error('unrecognised outputmode setting: %s',outputmode)
    end % switch outputmode

    % null distributions - same for ROI and searchlight
    nullout = fullfile(resdir,sprintf('%s_nulldist.mat',...
        predictors(pre).name));
    % save as volume with massive ndata
    nullvol = MriVolume(nulldists,disvol);
    % for mysterious reasons Matlab cannot save this in any older
    % version
    save(nullout,'nullvol','-v7');
    nulldistpaths{end+1} = nullout;
end % pre 1:npredictors

if strcmp(outputmode,'roi')
    % save data
    save(rpaths,'res_r');
    save(ppaths,'res_p');
    save(pfwepaths,'res_pfwe');
    % make bar chart for each ROI
    x = 1:npredictors;
    F = figure;
    for r = 1:disvol.nfeatures
        roistr = names{r};
        % plot rho in bars
        rho = structfun(@(x)x,res_r(r));
        % todo - errors, color scheme (see riken code for flexible bar
        % colors)
        F = barchart(rho,'x',x,'labels',stripbadcharacters(...
            fieldnames(res_r(r)),' '),'fighand',F);
        ylabel('spearman rho')
        titlestr = stripbadcharacters(roistr);
        if isfield(disvol.meta.features,'nfeatures')...
                && ~isempty(...
                disvol.meta.features.nfeatures)
            titlestr = sprintf('%s (%.0f voxels)',...
                titlestr,disvol.meta.features.nfeatures(r));
        end
        title(titlestr);
        % add p values on top of each bar
        p = structfun(@(x)x,res_p(r));
        y = max(double(rho),zeros(size(rho)));
        % adaptively change precision of p strs
        % depending on nperms
        T = addptext(x,y,p,gca,ceil(log10(nperms)));
        box off
        printstandard(fullfile(figdir,sprintf(...
            'barbyroi_%s',roistr)));
        clf(F);
    end
    close(F);

    % make bar chart for each analysis
    F = figure;
    xlabels = stripbadcharacters(names,' ');
    x = 1:length(xlabels);
    for pre = 1:npredictors
        prestr = predictors(pre).name;
        rho = [res_r.(prestr)];
        F = barchart(rho,'x',x,'labels',xlabels,'fighand',F);
        ylabel('spearman rho')
        title(stripbadcharacters(prestr,' '));
        p = [res_p.(prestr)];
        y = max(double(rho),zeros(size(rho)));
        T = addptext(x,y,p,gca,ceil(log10(nperms)));
        box off
        printstandard(fullfile(figdir,sprintf(...
            'barbypre_%s',prestr)));
        clf(F);
    end
    close(F);
end
