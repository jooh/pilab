% Convert data and the data field in a DesignSpec instance to a ConvGLM
% instance. This is basically a wrapper around design2convglm.
%
% construction:
% gl = ConvGLMConstructor([cvgroup])
%
% call:
% model = call(gl,epimat,design,chunks)
classdef ConvGLMConstructor < Saveable
    properties
        glmvarargs = {};
        cvgroup = {};
    end

    methods
        function gl = ConvGLMConstructor(cvgroup)
            if ~any(nargin)
                return
            end
            if ~ieNotDefined('cvgroup')
                gl.cvgroup = cvgroup;
            end
            if ~iscell(gl.cvgroup)
                gl.cvgroup = num2cell(gl.cvgroup);
            end
        end

        function model = call(self,epimat,design,chunks)
            model = design2convglm(design,epimat,chunks);
            if ~isempty(self.cvgroup)
                [model.cvgroup] = self.cvgroup{:};
            end
        end
    end
end
