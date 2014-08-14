% standard bar chart visualisations for results from roidata_lindisc
%
% plot_lindisc(figdir,res,varargin)
function plot_lindisc(figdir,res,varargin)

[m,errs,p,ylab] = plot_lindisc_parseargs(res,varargin{:});

mkdirifneeded(figdir);

[ncon,nroi] = size(m);

roicolors = hsv(nroi) * .9;
concolors = hsv(ncon) * .9;

% invisible figure
fighand = figurebetter([],[30 30],1/2,false);
% make it current fig
figure(fighand);

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
    plotsub(fighand,thism,thiserr,thisp,concolors,conclean);
    ylabel(ylab);
    title(thisroi);
    box off
    printstandard(fullfile(figdir,sprintf('byroi_%s',thisroi)));
    clf(fighand);
end

% one chart per analysis
for c = 1:ncon
    subplot(4,4,[6:7 10:11]);
    thiscon = stripbadcharacters(conclean{c},'_');
    thism = m(c,:);
    thiserr = errs(c,:);
    thisp = p(c,:);
    plotsub(fighand,thism,thiserr,thisp,roicolors,roiclean);
    ylabel(ylab);
    title(conclean{c});
    box off
    printstandard(fullfile(figdir,sprintf('bycon_%s',thiscon)));
    clf(fighand);
end

close(fighand);

function plotsub(fighand,thism,thiserr,thisp,concolors,conclean)

barchart(thism,'errors',thiserr,'facecolor',concolors,...
    'x',1:length(thism),'labels',conclean,'fighand',fighand);
% nansum here because errs can be all NaN
T = addptext(1:length(thism),nansum(cat(3,thism,thiserr),3),thisp);
