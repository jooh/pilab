% standard bar chart visualisations for results from some roidata struct
% (output from e.g. roidata_rsa). Detects particular prefixes ('l_', 'r_'
% for ROIs, 'contrast_' for conditions) and separates these into separate
% figures in sub-directories.
%
% plot_roidata(figdir,res,groupres,varargin)
function plot_roidata(figdir,res,groupres,varargin)

[m,errs,p,ylab] = plot_roidata_parseargs(res,varargin{:});

mkdirifneeded(figdir);

[ncon,nroi] = size(m);
if isfield(res,'z_subject')
    nsub = numel(res.z_subject);
else
    nsub = [];
end

% non-overlapping ROI and condition colours from glasbey colour map
roicolors = cmap_glasbey(nroi);
concolors = cmap_glasbey(ncon+nroi);
concolors(1:nroi,:) = [];

% invisible figure
fighand = figurebetter([],[30 30],1/2,false);
% make it current fig (breaks invisibility though)
figure(fighand);
% overlay singles in glasbey map
set(gcf,'defaultaxescolororder',cmap_glasbey(nsub));

conclean = stripbadcharacters(res.rows_contrast,' ');
roiclean = stripbadcharacters(res.cols_roi);

contrastind = find(cellfun(@(thisc)strcmp('contrast_',...
    thisc(1:min([9 numel(thisc)]))),res.rows_contrast));
allconind = 1:ncon;
notcontrastind = setdiff(allconind,contrastind);

congroups = struct('name','bycon_','ind',notcontrastind,'prefix','');
if any(contrastind)
    congroups(end+1) = struct('name','bycon_contrast_','ind',...
        contrastind,'prefix','\Delta');
end

leftind = find(cellfun(@(thisc)strcmp('l_',thisc(1:min([2 numel(...
    thisc)]))),res.cols_roi));
rightind = find(cellfun(@(thisc)strcmp('r_',thisc(1:min([2 numel(...
    thisc)]))),res.cols_roi));
notleftright = setdiff(1:nroi,[leftind rightind]);

roigroups = struct('name','byroi_bi_','ind',notleftright);
if any(leftind)
    roigroups(end+1) = struct('name','byroi_left_','ind',leftind);
end

if any(rightind)
    % nb we don't change the prefix because l_ and r_ are already in file
    % name
    roigroups(end+1) = struct('name','byroi_right_','ind',rightind);
end

logstr('generating figures for %d condition groups and %d roi groups\n',...
    numel(congroups),numel(roigroups));

for rgroup = 1:numel(roigroups)
    thisr = roigroups(rgroup);
    thisroicolor = roicolors(1:numel(thisr.ind),:);
    for cgroup = 1:numel(congroups)
        thisc = congroups(cgroup);
        thisylab = ylab;
        if ~isempty(thisc.prefix)
            if iscell(thisylab)
                thisylab{1} = [thisc.prefix thisylab{1}];
            else
                thisylab = [thisc.prefix thisylab];
            end
        end
        thisconcolor = concolors(1:numel(thisc.ind),:);
        % set the figdir
        thisfigdir = fullfile(figdir,[thisr.name thisc.name]);
        mkdirifneeded(thisfigdir);

        % plot each group separately
        for r = thisr.ind
            % awkwardly positioned plot to insure tons of white space so rotated x
            % labels don't get cropped out
            subplot(4,4,[6:7 10:11]);
            thisroi = roiclean{r};
            thism = m(thisc.ind,r);
            thiserr = errs(thisc.ind,r);
            thisp = p(thisc.ind,r);
            [fh,bh,eh,th] = plotsub(fighand,thism,thiserr,thisp,...
                thisconcolor,conclean(thisc.ind));
            ylabel(thisylab);
            title(thisroi);
            box off
            printstandard(fullfile(thisfigdir,sprintf('%s_%s',...
                thisr.name,thisroi)));
            % singles analysis
            if ~isempty(groupres)
                hold on;
                ph = plot(1:length(thism),squeeze(groupres.r(thisc.ind,r,:)),'.',...
                    'markersize',12);
                % set barchart underneath to outlines for clarity
                set(bh,'facecolor',[1 1 1],'edgecolor',[0 0 0]);
                % add a legend
                L = legend(ph,res.z_subject,'plotboxaspectratiomode','auto');
                lax = subplot(4,4,12);
                centerinaxis(L,lax);
                set(L,'box','off');
                delete(eh);
                delete(th);
                printstandard(fullfile(thisfigdir,sprintf(...
                    'singles_%s_%s',thisr.name,thisroi)));
            end
            clf(fighand);
        end % for r = thisr.ind(:)'

        % one chart per analysis
        for c = 1:thisc.ind
            subplot(4,4,[6:7 10:11]);
            thiscon = stripbadcharacters(conclean{c},'_');
            thism = m(c,thisr.ind);
            thiserr = errs(c,thisr.ind);
            thisp = p(c,thisr.ind);
            [fh,bh,eh,th] = plotsub(fighand,thism,thiserr,thisp,...
                thisroicolor,roiclean(thisr.ind));
            ylabel(thisylab);
            title(conclean{c});
            box off
            printstandard(fullfile(thisfigdir,sprintf('bycon_%s_%s',...
                thiscon,thisr.name)));
            % singles analysis
            if ~isempty(groupres)
                hold on;
                ph = plot(1:length(thism),squeeze(...
                    groupres.r(c,thisr.ind,:)),'.',...
                    'markersize',12);
                % set barchart underneath to outlines for clarity
                set(bh,'facecolor',[1 1 1],'edgecolor',[0 0 0]);
                % add a legend
                L = legend(ph,res.z_subject,...
                    'plotboxaspectratiomode','auto');
                lax = subplot(4,4,12);
                centerinaxis(L,lax);
                set(L,'box','off');
                delete(eh);
                delete(th);
                printstandard(fullfile(thisfigdir,...
                    sprintf('singles_bycon_%s_%s',thiscon,thisr.name)));
            end
            clf(fighand);
        end % for c = 1:thisc.ind(:)
    end % for cgroup = 1:numel(congroups)
end % for rgroup = 1:numel(roigroups)
close(fighand);

function [fh,bh,eh,th] = plotsub(fighand,thism,thiserr,thisp,concolors,conclean)
[fh,bh,eh] = barchart(thism,'errors',thiserr,'facecolor',concolors,...
    'x',1:length(thism),'labels',conclean,'fighand',fighand);
% nansum here because errs can be all NaN
th = addptext(1:length(thism),nansum(cat(3,thism,thiserr),3),thisp);
