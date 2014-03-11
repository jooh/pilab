% Convert data matrices to some GLM instance. This is basically a wrapper
% around array2glm.
%
% construction:
% gl = GLMConstructor(glmclass,cvgroup,[glmvarargin])
%
% call:
% model = call(gl,designmat,epimat,chunks)
classdef GLMConstructor < matlab.mixin.Copyable
    properties
        glmclass = 'GLM';
        glmvarargs = {};
        cvgroup = {};
    end

    methods
        function gl = GLMConstructor(glmclass,cvgroup,varargin)
            if ~any(nargin)
                return
            end
            if ~ieNotDefined('glmclass')
                gl.glmclass = glmclass;
            end
            if ~ieNotDefined('cvgroup')
                gl.cvgroup = cvgroup;
            end
            if ~iscell(gl.cvgroup)
                gl.cvgroup = num2cell(gl.cvgroup);
            end
            glm.varargs = varargin;
        end

        function model = call(self,designmat,epimat,chunks)
            model = array2glm(designmat,epimat,chunks,self.glmclass,...
                self.glmvarargs{:});
            if ~isempty(self.cvgroup)
                [model.cvgroup] = self.cvgroup{:};
            end
        end
    end
end
