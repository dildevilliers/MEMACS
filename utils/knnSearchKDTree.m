function [nearestPoints, distances] = knnSearchKDTree(tree, queryPoint, k)
% KNNSEARCHKDTREE finds the k-nearest neighbors using a k-d tree.
% Inputs:
%   - tree: k-d tree built using buildKDTree.
%   - queryPoint: 1x3 vector representing the query point.
%   - k: Number of nearest neighbors to find.
% Outputs:
%   - nearestPoints: k nearest points (Nx3 matrix).
%   - distances: Distances to the k nearest points.

    % Priority queue for k-nearest neighbors (max-heap)
    maxHeap = PriorityQueue(k);

    % Recursive search
    function recursiveSearch(node)
        if isempty(node)
            return;
        end

        % Compute the distance to the current node
        distance = norm(node.point - queryPoint);

        % Add current point to the max-heap
        maxHeap.insert(node.point, distance);

        % Determine which branch to search first
        axis = node.axis;
        if queryPoint(axis) < node.point(axis)
            primary = node.left;
            secondary = node.right;
        else
            primary = node.right;
            secondary = node.left;
        end

        % Search the primary branch
        recursiveSearch(primary);

        % Check if we need to explore the secondary branch
        if maxHeap.isNotFull() || ...
           abs(queryPoint(axis) - node.point(axis)) < maxHeap.maxKey()
            recursiveSearch(secondary);
        end
    end

    % Start the search
    recursiveSearch(tree);

    % Extract results from the max-heap
    [nearestPoints, distances] = maxHeap.getElements();
end
