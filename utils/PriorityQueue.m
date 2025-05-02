% Datastructure used by knnSearchKDTree

classdef PriorityQueue < handle
    properties
        Capacity
        Keys
        Values
        Indexes
        Size
        isFurthest % Logical flag: true for furthest points, false for nearest points
    end

    methods
        function obj = PriorityQueue(capacity, Ndims, isFurthest)
            if nargin < 2 || isempty(Ndims), Ndims = 3; end
            if nargin < 3 || isempty(isFurthest), isFurthest = false; end % Default: nearest neighbor
            obj.Capacity = capacity;
            obj.Keys = inf(capacity, 1);
            obj.Values = zeros(capacity, Ndims);
            obj.Indexes = zeros(capacity, 1);
            obj.Size = 0;
            obj.isFurthest = isFurthest;
        end

        function insert(obj, value, key, index)
            if obj.Size < obj.Capacity
                obj.Size = obj.Size + 1;
                obj.Keys(obj.Size) = key;
                obj.Values(obj.Size, :) = value;
                obj.Indexes(obj.Size) = index;
            else
                if obj.isFurthest
                    % For furthest points: replace smallest key if new key is larger
                    [minKey, insertIdx] = min(obj.Keys);
                    if key > minKey
                        obj.Keys(insertIdx) = key;
                        obj.Values(insertIdx, :) = value;
                        obj.Indexes(insertIdx) = index;
                    end
                else
                    % For nearest points: replace largest key if new key is smaller
                    [maxKey, insertIdx] = max(obj.Keys);
                    if key < maxKey
                        obj.Keys(insertIdx) = key;
                        obj.Values(insertIdx, :) = value;
                        obj.Indexes(insertIdx) = index;
                    end
                end
            end
        end

        function [values, keys, indexes] = getElements(obj)
            if obj.isFurthest
                [keys, idx] = sort(obj.Keys(1:obj.Size), 'descend'); % Furthest first
            else
                [keys, idx] = sort(obj.Keys(1:obj.Size), 'ascend'); % Nearest first
            end
            values = obj.Values(idx, :);
            indexes = obj.Indexes(idx);
        end

        function key = maxKey(obj)
            key = max(obj.Keys(1:obj.Size));
        end

        function key = minKey(obj)
            key = min(obj.Keys(1:obj.Size));
        end

        function result = isNotFull(obj)
            result = obj.Size < obj.Capacity;
        end
    end
end









% classdef PriorityQueue < handle
%     properties
%         Capacity
%         Keys
%         Values
%         Indexes
%         Size
%     end
% 
%     methods
%         function obj = PriorityQueue(capacity,Ndims)
%             if nargin < 2 || isempty(Ndims), Ndims = 3; end
%             obj.Capacity = capacity;
%             obj.Keys = inf(capacity, 1); % Max-heap
%             obj.Values = zeros(capacity, Ndims); % Store 3D points
%             obj.Indexes = zeros(capacity, 1); % Store 3D points
%             obj.Size = 0;
%         end
% 
%         function insert(obj, value, key, index)
%             if obj.Size < obj.Capacity
%                 obj.Size = obj.Size + 1;
%                 obj.Keys(obj.Size) = key;
%                 obj.Values(obj.Size, :) = value;
%                 obj.Indexes(obj.Size) = index;
%             else
% 
%                 % Replace the largest key if the new key is smaller
%                 [~, insertIdx] = max(obj.Keys);
%                 if key < obj.Keys(insertIdx)
%                     obj.Keys(insertIdx) = key;
%                     obj.Values(insertIdx, :) = value;
%                     obj.Indexes(insertIdx) = index;
%                 end
% 
%             end
%         end
% 
%         function [values, keys, indexes] = getElements(obj)
%             [keys, idx] = sort(obj.Keys(1:obj.Size), 'ascend');
%             values = obj.Values(idx, :);
%             indexes = obj.Indexes(idx, :);
%         end
% 
%         function result = maxKey(obj)
%             result = max(obj.Keys(1:obj.Size));
%         end
% 
%         function result = minKey(obj)
%             result = min(obj.Keys(1:obj.Size));
%         end
% 
%         function result = isNotFull(obj)
%             result = obj.Size < obj.Capacity;
%         end
%     end
% end
