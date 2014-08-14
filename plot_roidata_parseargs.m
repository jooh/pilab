% parse input arguments to plot_lindisc and return processed data ready for
% plotting. This function is mainly useful if you want to use a custom
% plugin plot function in e.g. aamod_pilab_decode_lindisc_visualisation.
%
% [m,errs,p,ylab] = plot_lindisc_parseargs(res,varargin)
function [m,errs,p,ylab] = plot_lindisc_parseargs(res,varargin)

mlabel = [];
getArgs(varargin,{'mtarget','t','errtarget',[],'ptarget','ppara',...
    'mlabel','','errlabel','','pthresh',.05});

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

if ~isempty(pthresh)
    p(p>pthresh) = NaN;
end

