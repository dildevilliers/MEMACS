classdef PriorityQueue < handle
    properties
        Capacity
        Keys
        Values
        Indexes
        Size
    end
    
    methods
        function obj = PriorityQueue(capacity,Ndims)
            if nargin < 2 || isempty(Ndims), Ndims = 3; end
            obj.Capacity = capacity;
            obj.Keys = inf(capacity, 1); % Max-heap
            obj.Values = zeros(capacity, Ndims); % Store 3D points
            obj.Indexes = zeros(capacity, 1); % Store 3D points
            obj.Size = 0;
        end
        
        function insert(obj, value, key, index)
            if obj.Size < obj.Capacity
                obj.Size = obj.Size + 1;
                obj.Keys(obj.Size) = key;
                obj.Values(obj.Size, :) = value;
                obj.Indexes(obj.Size) = index;
            else
                % Replace the largest key if the new key is smaller
                [~, maxIdx] = max(obj.Keys);
                if key < obj.Keys(maxIdx)
                    obj.Keys(maxIdx) = key;
                    obj.Values(maxIdx, :) = value;
                    obj.Indexes(maxIdx) = index;
                end
            end
        end
        
        function [values, keys, indexes] = getElements(obj)
            [keys, idx] = sort(obj.Keys(1:obj.Size), 'ascend');
            values = obj.Values(idx, :);
            indexes = obj.Indexes(idx, :);
        end
        
        function result = maxKey(obj)
            result = max(obj.Keys(1:obj.Size));
        end
        
        function result = isNotFull(obj)
            result = obj.Size < obj.Capacity;
        end
    end
end
