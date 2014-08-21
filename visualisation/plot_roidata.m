% standard bar chart visualisations for results from some roidata struct
% (output from e.g. roidata_rsa).
%
% plot_roidata(figdir,res,groupres,varargin)
function plot_roidata(figdir,res,groupres,varargin)

[m,errs,p,ylab] = plot_roidata_parseargs(res,varargin{:});

mkdirifneeded(figdir);

[ncon,nroi] = size(m);
nsub = numel(res.z_subject);

roicolors = hsv(nroi) * .9;
concolors = hsv(ncon) * .9;

% invisible figure
fighand = figurebetter([],[30 30],1/2,false);
% make it current fig
figure(fighand);
% overlay singles in glasbey map
set(gcf,'defaultaxescolororder',cmap_glasbey(nsub));

conclean = stripbadcharacters(res.rows_contrast,' ');
roiclean = stripbadcharacters(res.cols_roi);

% one chart per ROI
for r = 1:nroi
    % awkwardly positioned plot to insure tons of white space so rotated x
    % labels don't get cropped out
    subplot(4,4,[6:7 10:11]);
    thisroi = roiclean{r};
    thism = m(:,r);
    thiserr = errs(:,r);
    thisp = p(:,r);
    [fh,bh,eh,th] = plotsub(fighand,thism,thiserr,thisp,concolors,conclean);
    ylabel(ylab);
    title(thisroi);
    box off
    printstandard(fullfile(figdir,sprintf('byroi_%s',thisroi)));
    % singles analysis
    if ~isempty(groupres)
        hold on;
        ph = plot(1:length(thism),squeeze(groupres.r(:,r,:)),'.',...
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
        printstandard(fullfile(figdir,sprintf('singles_byroi_%s',thisroi)));
    end
    clf(fighand);
end

% one chart per analysis
for c = 1:ncon
    subplot(4,4,[6:7 10:11]);
    thiscon = stripbadcharacters(conclean{c},'_');
    thism = m(c,:);
    thiserr = errs(c,:);
    thisp = p(c,:);
    [fh,bh,eh,th] = plotsub(fighand,thism,thiserr,thisp,roicolors,roiclean);
    ylabel(ylab);
    title(conclean{c});
    box off
    printstandard(fullfile(figdir,sprintf('bycon_%s',thiscon)));
    % singles analysis
    if ~isempty(groupres)
        hold on;
        ph = plot(1:length(thism),squeeze(groupres.r(c,:,:)),'.',...
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
        printstandard(fullfile(figdir,sprintf('singles_bycon_%s',thiscon)));
    end
    clf(fighand);
end

close(fighand);

function [fh,bh,eh,th] = plotsub(fighand,thism,thiserr,thisp,concolors,conclean)

[fh,bh,eh] = barchart(thism,'errors',thiserr,'facecolor',concolors,...
    'x',1:length(thism),'labels',conclean,'fighand',fighand);
% nansum here because errs can be all NaN
th = addptext(1:length(thism),nansum(cat(3,thism,thiserr),3),thisp);
