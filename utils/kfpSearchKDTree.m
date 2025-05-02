function [furthestPoints, distances, indexes] = kfpSearchKDTree(tree, queryPoint, k)
% KFPSEARCHKDTREE finds the k-furthest neighbors using a k-d tree.
% Inputs:
%   - tree: k-d tree built using buildKDTree.
%   - queryPoint: 1x3 vector representing the query point.
%   - k: Number of furthest neighbors to find.
% Outputs:
%   - furthestPoints: k furthest points (Nx3 matrix).
%   - distances: Distances to the k furthest points.
%   - indexes: Indexes of the k furthest points in the original tree.

    % Priority queue for k-furthest neighbors (min-heap behavior with flag)
    minHeap = PriorityQueue(k, size(tree.point,2), true);

    % Recursive search
    function recursiveSearch(node)
        if isempty(node)
            return;
        end

        % Compute the distance to the current node
        distance = norm(node.point - queryPoint);

        % Add current point to the min-heap (furthest points)
        minHeap.insert(node.point, distance, node.index);

        % Determine which branch to search first (less important here, but maintain structure)
        axis = node.axis;
        if queryPoint(axis) < node.point(axis)
            primary = node.left;
            secondary = node.right;
        else
            primary = node.right;
            secondary = node.left;
        end

        % Always search primary branch first
        recursiveSearch(primary);

        % Check if we should explore the secondary branch
        % Use bounding box max distance for pruning efficiency
        if minHeap.isNotFull() || (maxPossibleDist(node.bounds, queryPoint) > minHeap.minKey())
            recursiveSearch(secondary);
        end
    end

    % Start the recursive search
    recursiveSearch(tree);

    % Extract results from the min-heap
    [furthestPoints, distances, indexes] = minHeap.getElements();
end

% Helper function to compute max possible distance to bounding box
function maxDist = maxPossibleDist(bounds, queryPoint)
    furthestCorner = zeros(size(queryPoint));
    for d = 1:length(queryPoint)
        if abs(bounds(1,d)-queryPoint(d)) > abs(bounds(2,d)-queryPoint(d))
            furthestCorner(d) = bounds(1,d);
        else
            furthestCorner(d) = bounds(2,d);
        end
    end
    maxDist = norm(furthestCorner - queryPoint);
end
