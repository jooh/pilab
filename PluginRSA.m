% RSA based on a custom function, applied to each feature in a parfor loop
% (with optional varargs applied). This is slower than vectorised forms of
% RSA but enables custom RSA types such as Kendall's tau.
%
% gl = PluginRSA(modelrdms,datardms,corrfun,[corrargs])
classdef PluginRSA < RSA
    properties
        corrfun = []; % type input to matlab's corr function
        corrargs = {};
    end

    methods
        function gl = PluginRSA(modelrdms,datardms,corrfun,varargin)
            % use super-class to initialise            
            gl = gl@RSA(modelrdms,datardms);
            assert(gl.npredictors==1,'only 1 model RDM input supported');
            % must allow undefined corrfun to support resampling
            % (cloneargs method is only called on new instance _after_ it
            % has been created)
            if ieNotDefined('corrfun')
                corrfun = [];
            end
            if ischar(corrfun)
                % insure function handle
                corrfun = str2func(corrfun);
            end
            gl.corrfun = corrfun;
            gl.corrargs = varargin;
        end

        function cloneargs(self,oldinstance)
            % make sure corrfun and corrargs survive resampling.
            [self.corrfun] = deal(oldinstance.corrfun);
            [self.corrargs] = deal(oldinstance.corrargs);
        end

        function estimates = fit(self)
        % apply some corrfun to every feature in self.data.
        % estimates = fit(self)
            estimates = NaN([1 self(1).nfeatures]);
            X = vertcat(self.X);
            data = vertcat(self.data);
            inds = 1:self(1).nsamples;
            fun = self.corrfun;
            args = self.corrargs;
            % parfor here means we will encounter nested parfor loops when
            % e.g. bootstrapping or permutation testing. However, if we
            % don't do this e.g searchlight mapping is likely to turn out
            % to be ridiculously slow
            parfor f = 1:self(1).nfeatures
                estimates(1,f) = fun(X,data(inds,f),args{:});
            end
        end
    end
end
