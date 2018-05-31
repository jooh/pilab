% Construct a RSA-based GLM instance. Unlike the standard GLMConstructor,
% this class takes model RDM inputs that get stored with the instance.
%
% construction:
% rs = RSAConstructor(glmclass,modelrdms,[glmvarargin])
%
% call:
% model = call(rdvecs)

classdef RSAConstructor < GLMConstructor
    properties
        modelrdms
    end

    methods
        function rs = RSAConstructor(glmclass,modelrdms,varargin)
            if ieNotDefined('glmclass')
                glmclass = [];
            end
            rs = rs@GLMConstructor(glmclass,{},varargin{:});
            if ~any(nargin)
                return
            end
            rs.modelrdms = modelrdms;
        end

        function model = call(self,rdvecs)
            model = feval(self.glmclass,self.modelrdms,rdvecs,...
                self.glmvarargs{:});
        end
    end
end
