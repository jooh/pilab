% parse input arguments to plot_roidata and return processed data ready for
% plotting. This function is mainly useful if you want to use a custom
% plugin plot function in e.g. aamod_pilab_rsa_visualisation.
%
% [m,errs,p,ylab,groupm,groupp] = plot_roidata_parseargs(res,groupres,varargin)
function [m,errs,p,ylab,groupm,groupp] = plot_roidata_parseargs(res,groupres,varargin)


mlabel = [];
getArgs(varargin,{'mtarget','mean','errtarget',[],'ptarget','ppara',...
    'mlabel','','errlabel','','groupmtarget','r','groupptarget',...
    'pperm'},'suppressUnknownArgMessage=1');

if ieNotDefined('groupres')
    groupm = [];
    groupp = [];
else
    % this extra test is mainly to make it easy to spot code where the call
    % syntax needs to be updated
    assert(isstruct(groupres),'groupres must be second input');
    groupm = groupres.(groupmtarget);
    groupp = groupres.(groupptarget);
end

% fill in a dummy label if nothing else is specified
if isempty(mlabel)
    mlabel = mtarget;
end

m = res.(mtarget);

if isempty(errtarget)
    errs = NaN(size(m));
    ylab = mlabel;
else 
    if isempty(errlabel)
        % dummy label
        errlabel = ['\pm1 ' errtarget];
    end
    ylab = {mlabel,errlabel};
    errs = res.(errtarget);
end


if isempty(ptarget)
    p = NaN(size(m));
else
    p = res.(ptarget);
end
