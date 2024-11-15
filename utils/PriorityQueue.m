classdef PriorityQueue < handle
    properties
        Capacity
        Keys
        Values
        Size
    end
    
    methods
        function obj = PriorityQueue(capacity)
            obj.Capacity = capacity;
            obj.Keys = inf(capacity, 1); % Max-heap
            obj.Values = zeros(capacity, 3); % Store 3D points
            obj.Size = 0;
        end
        
        function insert(obj, value, key)
            if obj.Size < obj.Capacity
                obj.Size = obj.Size + 1;
                obj.Keys(obj.Size) = key;
                obj.Values(obj.Size, :) = value;
            else
                % Replace the largest key if the new key is smaller
                [~, maxIdx] = max(obj.Keys);
                if key < obj.Keys(maxIdx)
                    obj.Keys(maxIdx) = key;
                    obj.Values(maxIdx, :) = value;
                end
            end
        end
        
        function [values, keys] = getElements(obj)
            [keys, idx] = sort(obj.Keys(1:obj.Size), 'ascend');
            values = obj.Values(idx, :);
        end
        
        function result = maxKey(obj)
            result = max(obj.Keys(1:obj.Size));
        end
        
        function result = isNotFull(obj)
            result = obj.Size < obj.Capacity;
        end
    end
end
