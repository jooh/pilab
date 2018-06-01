% Base class for supporting struct conversion as part of save and load with
% obj2struct and struct2obj. This tends to speed up saving of handle class
% instances. Note that this trick only works with object classes that have
% only public, writeable fields and support initialisation with no input
% arguments.
%
% Saveable < matlab.mixin.Copyable
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
