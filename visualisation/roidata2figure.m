% Visualise roidata results with one EPS figure per region.
%
% INPUTS
% figdir        path to directory where figures are written
% res           output from e.g. roidata_rfx
% groupres      output from e.g. roidata_rfx
%
% NAMED INPUTS
%
% ARGUMENT          DEFAULT                     DEFINITION
% roiind            1:numel(res.cols_roi)       ROIs to make figures for
% conind            1:numel(res.rows_contrast)  conditions to plot (cell of
%                                                   char or numeric index).
%                                                   Can be nested cell
%                                                   array for grouped bars.
% restarget         ''                          if present draw as bar.
% grouptarget       ''                          if present draw as marker.
% errtarget         []                          if present draw as errorbar
% axisscale         .25                         size of axis in norm figure
%                                                   units
% xscalefactor      2.5                         scale width of axis by this
%                                                   factor. If kept
%                                                   consistent (along with
%                                                   fsize and axscale) bars
%                                                   should be of similar
%                                                   width across figures
%                                                   that vary in number of
%                                                   bars.
% barwidth          'adaptive'                  scale width of bars to
%                                                   accommodate the bigger
%                                                   bar group.
% fsize             [10,10]                     input for figurebetter.
% facecolor         [0,.3,.7]                   bar facecolor
% edgecolor         [1,1,1]                     bar edgecolor
% groupmarkerstyle  'o'                         single subject marker
% groupmarkersize   4                           single subject marker size
% groupfacecolor    [1,1,1]                     single subject facecolor
% groupedgecolor    [.5,.5,.5]                  single subject edgecolor
% noiseind          {'rep_lower','rep_upper'}   use to draw a shaded noise
%                                                   ceiling
% noisecolor        [.8,.8,.8]                  noise ceiling color
% precision         1                           scale axis to this number
%                                                   of significant digits.
% specialval        []                          roundplotlimits input
% specialtick       0                           roundplotlimits input
% docontrasts       true                        draw pairwise comparisons
%                                                   using contrastlines
% dozerotest        true                        add p values for each entry
%                                                   using addptext
% xticklabels       []                          labels for each x position
%                                                   (needs to be defined if
%                                                   using grouped bars)
% ptarget           'ppara'                     use this field for
%                                                   inference
% baseline          0                           draw this value as dashed
%                                                   line
% ylab              ''                          ylabel
% padding           .05                         use this proportion of
%                                                   range(ylim) of white
%                                                   space to separate plot
%                                                   elements
%
% handles = roidata2figure(res,groupres,varargin)
function outfile = roidata2figure(figdir,res,groupres,varargin)

getArgs(varargin,{'roiind',1:numel(res.cols_roi),...
    'conind',1:numel(res.rows_contrast),'restarget','',...
    'grouptarget','','errtarget',[],'axisscale',.25,...
    'xscalefactor',2.5,'barwidth','adaptive',...
    'fsize',[10,10],'facecolor',[0 .3 .7],'edgecolor',...
    [1,1,1],'groupmarkerstyle','o','groupmarkersize',4,...
    'groupedgecolor',[.5 .5 .5],'groupfacecolor',[1 1 1],...
    'noiseind',{'rep_lower','rep_upper'},...
    'noisecolor',[.8 .8 .8],'precision',1,'specialval',[],...
    'specialtick',0,'docontrasts',true,'dozerotest',true,...
    'xticklabels',[],'ptarget','ppara','baseline',0,'ylab','',...
    'padding',.05});

% input parsing
% add additional plot styling fields to res if they don't already exist.
stylefields = {'facecolor','edgecolor','groupmarkerstyle','groupmarkersize','groupedgecolor','groupfacecolor'};
for fn = stylefields
    fnstr = fn{1};
    if ~isfield(res,fnstr)
        res.(fnstr) = repmat({eval(fnstr)},[size(res.rows_contrast,1),1]);
    end
end

if iscell(roiind)
    % do intersection to find numerical index
    [~,~,roiind] = intersect(roiind,res.cols_roi,'stable');
end
if islogical(roiind)
    roiind = find(roiind);
end
assert(isnumeric(roiind),'roiind must be index')
nroi = numel(roiind);

% handle conind - may be nested
[conind,groupind] = findconind(conind,res.rows_contrast);
nbarall = numel(conind);
xpos = unique(groupind);
nbargroup = numel(xpos);
npergroup = arrayfun(@(x)sum(groupind==x),xpos,'uniformoutput',true);
maxn = max(npergroup(:));
xl = [.5 nbargroup+.5];

% work out contrasts
if any(docontrasts(:))
    % this is now pretty easy actually
    % we may need to temporarily sort the conind to get the contrast naming to
    % appear consistent
    [conind_s,sortind] = sort(conind);
    contrasts = roidata_allpairwisecontrasts(res.rows_contrast(conind_s));
    % and the set of contrasts is
    contrastind = findconind({contrasts.name},res.rows_contrast);
    % we have to do a bit of checking here because the validity of the
    % squareform transform hinges on all contrasts actually being present
    assert(numel(contrastind)==nchoosek(nbarall,2),...
        'missing pairwise contrasts: should be %d, got %d',...
        nchoosek(nbarall,2),numel(contrastind));
    pcon = asrdmmat(res.(ptarget)(contrastind,roiind));
    % reverse the sort above to bring the contrast rows and columns back in
    % register with conind
    pcon = pcon(sortind,sortind,:);
    if isscalar(docontrasts)
        docontrasts = zerodiagonal(repmat(docontrasts,[size(pcon,1),...
            size(pcon,2)]));
    end
end

% maybe find reproducibility stats
noiseind = findconind(noiseind,res.rows_contrast);

if strcmp(barwidth,'adaptive')
    barwidth = .6 / maxn;
end

if ~ieNotDefined('groupres')
    assert(isequal(groupres.cols_roi(roiind),res.cols_roi(roiind)),...
        'res and groupres must have matching ROIs');
    assert(isequal(groupres.rows_contrast(conind),...
        res.rows_contrast(conind)),...
        'res and groupres must have matching contrasts');
    assert(isempty(grouptarget) || isfield(groupres,grouptarget),...
        'could not find grouptarget in groupres');
end

% scale the size of the panel according to the number of bar groups. This
% will ensure that panels that differ in bar number have similar bar
% widths.
axwidth = axisscale * (nbargroup / xscalefactor);
% nb slight shift upward because there's usually more additional stuff
% underneath than above the figure
axpos = [axisscale axisscale*1.3 axwidth axisscale];

% and next we need the exact x coordinate for each individual point
if maxn == 1
    xperbar = xpos;
    if ieNotDefined('xticklabels')
        % can infer
        xticklabels = res.rows_contrast(conind);
    end
else
    for b = 1:nbargroup
        % half the total width of the arrangement, less half a bar
        off = barwidth*npergroup(b)/2 - barwidth/2;
        xperbar{b} = xpos(b) + linspace(-off,off,npergroup(b));
    end
    xperbar = cat(2,xperbar{:});
end

barfaces = vertcat(res.facecolor{conind});
baredges = vertcat(res.edgecolor{conind});

% ok, here we go
fighand = figurebetter([],fsize,'auto');
handles = struct('bars',[],'errs',[],'singles',[],'noise',[],...
    'noisetext',[]);
for r = 1:nroi
    clf(fighand);
    axhand = axes('position',axpos,'tickdir','out','box','off',...
        'xlim',xl,'xtick',xpos,'xcolor',[1 1 1],'ticklength',[.03 .03]);
    hold(axhand,'on');
    if ~isempty(baseline) && ~isnan(baseline)
        handles.base = line(xl,[baseline baseline],'linestyle',':',...
            'color','k','linewidth',.5);
    end
    if ~isempty(restarget)
        handles.bars = arrayfun(@(ind)bar(xperbar(ind),...
            res.(restarget)(conind(ind),roiind(r)),barwidth,...
            'edgecolor',baredges(ind,:),'facecolor',barfaces(ind,:),...
            'showbaseline','off'),1:nbarall);
    end
    if ~isempty(errtarget)
        handles.errs = errorbar(xperbar,res.(restarget)(conind,roiind(r)),...
            res.(errtarget)(conind,roiind(r)),'linestyle','none',...
            'color',[0 0 0],'linewidth',.5);
        arrayfun(@errorbar_tick,handles.errs(:),...
            zeros(numel(handles.errs),1));
    end
    if ~isempty(grouptarget)
        singledata = groupres.(grouptarget)(conind,roiind(r),:);
        singledata = reshape(singledata,[nbarall size(singledata,3)]);
        handles.singles = arrayfun(@(thisg)plot(...
            xperbar(groupind==thisg),singledata(groupind==thisg,:),...
            'marker',groupmarkerstyle,'markersize',groupmarkersize,...
            'markerfacecolor',groupfacecolor,...
            'markeredgecolor',groupedgecolor,...
            'color',groupedgecolor,'linewidth',.5),xpos,'uniformoutput',0);
        handles.singles = cat(1,handles.singles{:})';
    end
    noise_high = NaN;
    if ~isempty(noiseind)
        noise_low = res.(restarget)(noiseind(1),roiind(r));
        noise_high = res.(restarget)(noiseind(2),roiind(r));
        handles.noiseceil = errorshade(xl',noise_low * [1; 1],...
            noise_high * [1; 1],noisecolor,0);
        uistack(handles.noiseceil,'bottom');
        handles.noisetext = text(max(xl) + range(xl)*padding,...
            mean([noise_low noise_high]),'noise ceiling',...
            'horizontalalignment','left','verticalalignment','middle');
    end
    finaly = sety(cat(2,handles.bars,handles.errs,handles.singles),...
        noise_high,precision,specialval,specialtick);
    newy = NaN;
    if any(docontrasts(:))
        thispcon = pcon(:,:,r);
        thispcon(thispcon>.05) = NaN;
        thispcon(~docontrasts) = NaN;
        [handles.conlines,newy] = contrastlines(gca,...
            zerodiagonal(thispcon),'xpos',xperbar,...
            'ypos',finaly(2) + range(finaly)*padding);
        finaly(2) = max([finaly,newy]);
    end
    if dozerotest
        [handles.ptext,newy] = addptext(xperbar,...
            max(finaly) + range(finaly)*padding,...
            res.(ptarget)(conind,roiind(r)),3,true,'rotation',45);
        finaly(2) = max([finaly(2),newy]);
    end
    if ~isempty(xticklabels)
        handles.xticklabels = text(xpos,...
            repmat(finaly(1)-padding*range(finaly),size(xpos)),...
            xticklabels,'rotation',45,...
            'horizontalalignment','right','verticalalignment','middle');
    end
    handles.tit = text((xl(2)-xl(1))/2+xl(1),max(finaly),...
        res.cols_roi{roiind(r)},'fontweight','bold',...
        'horizontalalignment','center','verticalalignment','bottom');
    ylabel(ylab)
    outfile{r} = fullfile(figdir,['roidata2figure_' ...
        stripbadcharacters(res.cols_roi{roiind(r)})]);
    printstandard(outfile{r});
end
close(fighand);

% sub function for identifying indices for conditions
function [con,groups] = findconind(con,rows_contrast)

if ischar(con) && strcmpi(con,'nocontrasts')
    % remove contrasts
    hits = strfindcell(rows_contrast,'contrast_');
    con = setdiff(1:numel(rows_contrast),hits);
end

if iscell(con) && ischar(con{1})
    % special case
    [~,~,con] = intersect(con,rows_contrast,'stable');
end

if iscell(con)
    % we are indexing for a group so need to recurse
    for x = 1:numel(con)
        [con{x},groups{x}] = findconind(con{x},rows_contrast);
        groups{x}(:) = x;
    end
    % unpack
    con = vertcat(con{:});
    groups = vertcat(groups{:});
else
    % check that con makes sense
    if islogical(con)
        con = find(con);
    end
    assert(isnumeric(con),'conind must be index');
    groups = [1:numel(con)]';
end

function finaly = sety(hand,noise_high,precision,specialval,specialtick)

ydata = getdatalims(hand,'y');
ydata(2) = max([ydata(2),noise_high]);
ylim(ydata);
finaly = roundplotlimits(gca,'y',specialval,precision);
minimalticks(gca,'y',specialtick,precision);
