% standard bar chart visualisations for results from some roidata struct
% (output from e.g. roidata_rsa). Detects particular prefixes ('l_', 'r_'
% for ROIs, 'contrast_' for conditions) and separates these into separate
% figures in sub-directories.
%
% plot_roidata(figdir,res,groupres,predictors,varargin)
function plot_roidata(figdir,res,groupres,predictors,varargin)

[m,errs,p,ylab,groupm,groupp] = plot_roidata_parseargs(res,groupres,varargin{:});
getArgs(varargin,{'extracongroups',{},'pthresh',.05});

mkdirifneeded(figdir);

[ncon,nroi] = size(m);
if isfield(res,'z_subject')
    nsub = numel(res.z_subject);
else
    nsub = 1;
end

% non-overlapping ROI and condition colours from glasbey colour map
allglas = cmap_glasbey;
% drop dark colours
allglas(sum(allglas,2)<.3,:) = [];
concolors = allglas(1:ncon,:);
allglas(1:ncon,:) = [];
roicolors = allglas(1:nroi,:);
allglas(1:nroi,:) = [];
notsigcolor = [.6 .6 .6];

% default to showing all in colour
singlesig = true([ncon nroi nsub]);
if ~isempty(groupp)
    % define points to gray out
    singlesig = groupp < 0.05;
end

% invisible figure
fighand = figurebetter([],[30 30],1/2,false);
% make it current fig (breaks invisibility though)
figure(fighand);
% overlay singles in glasbey map
subcolors = cmap_glasbey(nsub);

conclean = stripbadcharacters(res.rows_contrast,' ');
roiclean = stripbadcharacters(res.cols_roi);

taken = false([1 ncon]);
congroups = emptystruct('name','ind','prefix');
cons = [{'contrast_','view_left_','view_right_','view_within_',...
    'view_across_'} extracongroups];
for c = cons
    constr = c{1};
    hits = strfindcell(res.rows_contrast,constr);
    % remove any that were already taken
    hits = intersect(hits,find(~taken));
    if any(hits)
        % update the taken list
        taken(hits) = true;
        % add an entry to congroups
        congroups(end+1).name = constr;
        congroups(end).ind = hits;
        congroups(end).prefix = stripbadcharacters(constr,' ');
    end
end
if any(~taken)
    congroups(end+1).name = 'bycon_';
    congroups(end).ind = find(~taken);
    congroups(end).prefix = '';
end

% ROI groups
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
            set(th(thisp<pthresh),'fontweight','bold');
            ylabel(thisylab);
            box off
            printstandard(fullfile(thisfigdir,sprintf('%s_%s',...
                thisr.name,thisroi)),'formats=eps');
            % singles analysis
            if ~isempty(groupres)
                hold on;
                % set barchart underneath to outlines for clarity
                set(bh,'facecolor',[1 1 1],'edgecolor',[0 0 0]);
                % plot all the points in monochrome
                ph = plotarray(1:numel(thisc.ind),groupm(thisc.ind,r,:),...
                    '.','markersize',12);
                % then color according to subcolor
                for s = 1:nsub
                    if any(isnan(ph(:,:,s)))
                        assert(all(ph(:,:,s)),'should never happen');
                        continue;
                    end
                    set(ph(:,:,s),'markeredgecolor',subcolors(s,:));
                end
                nonsigind = (~singlesig(thisc.ind,r,:));
                % finally gray out non-sig
                set(ph(nonsigind),'markeredgecolor',notsigcolor);
                % add a legend
                validsub = ~isnan(squeeze(ph(1,1,:)));
                L = legend(squeeze(ph(1,1,validsub)),...
                    res.z_subject(validsub),...
                    'plotboxaspectratiomode','auto');
                lax = subplot(4,4,12);
                centerinaxis(L,lax);
                set(L,'box','off');
                delete(eh);
                delete(th);
                printstandard(fullfile(thisfigdir,sprintf(...
                    'singles_%s_%s',thisr.name,thisroi)),'formats=eps')
            end
            clf(fighand);
        end % for r = thisr.ind(:)'

        % one chart per analysis
        for c = thisc.ind(:)'
            subplot(4,4,[6:7 10:11]);
            thiscon = stripbadcharacters(conclean{c},'_');
            thism = m(c,thisr.ind);
            thiserr = errs(c,thisr.ind);
            thisp = p(c,thisr.ind);
            [fh,bh,eh,th] = plotsub(fighand,thism,thiserr,thisp,...
                thisroicolor,roiclean(thisr.ind));
            set(th(thisp<pthresh),'fontweight','bold');
            ylabel(thisylab);
            box off
            printstandard(fullfile(thisfigdir,sprintf('bycon_%s_%s',...
                thiscon,thisr.name)),'formats=eps');
            % singles analysis
            if ~isempty(groupres)
                hold on;
                % set barchart underneath to outlines for clarity
                set(bh,'facecolor',[1 1 1],'edgecolor',[0 0 0]);
                set(th,'visible','off');

                % plot all the points in monochrome
                ph = plotarray(repmat([1:numel(thisr.ind)],[1 1 nsub]),...
                    groupm(c,thisr.ind,:),'.','markersize',12);
                % then color according to subcolor
                for s = 1:nsub
                    if any(isnan(ph(:,:,s)))
                        assert(all(ph(:,:,s)),'should never happen');
                        continue;
                    end
                    set(ph(:,:,s),'markeredgecolor',subcolors(s,:));
                end
                nonsigind = ~singlesig(c,thisr.ind,:);
                % finally gray out non-sig
                set(ph(nonsigind),'markeredgecolor',notsigcolor);

                % add a legend
                validsub = ~isnan(squeeze(ph(1,1,:)));
                L = legend(squeeze(ph(1,1,validsub)),...
                    res.z_subject(validsub),...
                    'plotboxaspectratiomode','auto');
                lax = subplot(4,4,12);
                centerinaxis(L,lax);
                set(L,'box','off');
                delete(eh);
                delete(th);
                printstandard(fullfile(thisfigdir,...
                    sprintf('singles_bycon_%s_%s',thiscon,thisr.name)),...
                    'formats=eps');
            end
            clf(fighand);
        end % for c = thisc.ind
    end % for cgroup = 1:numel(congroups)
end % for rgroup = 1:numel(roigroups)
close(fighand);

function [fh,bh,eh,th] = plotsub(fighand,thism,thiserr,thisp,concolors,conclean)
[fh,bh,eh] = barchart(thism,'errors',thiserr,'facecolor',concolors,...
    'x',1:length(thism),'labels',conclean,'fighand',fighand);
yp = max(ylim) + range(ylim)*.1;
th = addptext(1:length(thism),max(ylim),thisp,3,'p=','rotation',45);
