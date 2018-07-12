function [R] = struct_from_list(varargin)

for i=1:2:numel(varargin)
    R.(varargin{i})=varargin{i+1};
end