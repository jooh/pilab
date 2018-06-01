% Vectorised form of spearman correlation coefficient for the special case
% where one input (a) is to be correlated with a set of inputs (b), e.g. a
% predictor RDM correlated against each of a set of searchlight RDMs.
%
% Uses pearsonvec on ranktrans transformed data. The key slow-down in this
% function relative to pearsonvec is the rank transform, which iterates. So
% in some applications (e.g., fitting the same ranks to different
% predictors), further gains can be had by first rank transforming and
% then using pearsonvec on the ranks.
% r = spearmanvec(a,b)
function r = spearmanvec(a,b)

r = pearsonvec(ranktrans(a),ranktrans(b));
