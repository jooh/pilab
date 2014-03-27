classdef Saveable < matlab.mixin.Copyable
    methods
        function outstruct = saveobj(self)
            outstruct = obj2struct(self);
        end
    end
    methods (Static=true)
        function obj = loadobj(instruct)
            obj = struct2obj(instruct);
        end
    end
end
