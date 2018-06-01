% Plot a typical shepard plot into the axis ax. This is used to assess the
% performance of MDS
%
% plotshepard(ax,dissimilarities,disparities,distances)
function plotshepard(ax,dissimilarities,disparities,distances)

dissimilarities=asrdmvec(dissimilarities);
distances=asrdmvec(distances);
ph = plot(ax,dissimilarities,distances,'bo');
leg = {'distances'};

if ~ieNotDefined('disparities')
    hold(ax,'on');
    disparities=asrdmvec(disparities);
    [dum,ord] = sortrows([disparities dissimilarities]);
    ph(end+1) = plot(ax,dissimilarities(ord),disparities(ord),'r.-');
    leg{end+1} = 'disparities';
end

%% labels
xlabel(ax,'dissimilarity');
ylabel(ax,'distance');
legend(ph,leg,'Location','NorthWest');
r_Pearson=corr(dissimilarities(:),distances(:),'type','Pearson');
r_Spearman=corr(dissimilarities(:),distances(:),'type','Spearman');
corrstr = sprintf('r(dissimilarity,distance)=%.2f, rho=%.2f',...
    corr(dissimilarities,distances,'type','Pearson'),...
    corr(dissimilarities,distances,'type','Spearman'));
title({'\fontsize{14}Shepard plot\fontsize{9}',corrstr});
