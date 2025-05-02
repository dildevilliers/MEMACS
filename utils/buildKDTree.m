function tree = buildKDTree(points, depth, indices)
% BUILDKDTREE constructs a k-d tree for a given set of points.
% Inputs:
%   - points: NxD matrix representing 3D points.
%   - depth: Current depth of the tree (start with 0).
%   - indices: Original indices of the points in the dataset.
% Outputs:
%   - tree: A structure representing the k-d tree.

    if isempty(points)
        tree = [];
        return;
    end

    if nargin < 3 || isempty(indices)
        indices = (1:size(points,1)).';
    end

    D = size(points,2);

    % Determine the axis to split on (cycle through 1, 2, 3)
    axis = mod(depth, D) + 1;

    % Sort points along the splitting axis
    [sortedPoints, sortIdx] = sortrows(points, axis);
    sortedIndices = indices(sortIdx); % Sort original indices in the same order

    % Find the median index
    medianIdx = floor(size(sortedPoints, 1) / 2) + 1;

    % Create the tree node
    tree.point = sortedPoints(medianIdx, :);
    tree.index = sortedIndices(medianIdx); % Store the original index
    tree.axis = axis;
    tree.left = buildKDTree(sortedPoints(1:medianIdx-1, :), depth + 1, sortedIndices(1:medianIdx-1));
    tree.right = buildKDTree(sortedPoints(medianIdx+1:end, :), depth + 1, sortedIndices(medianIdx+1:end));

    % Compute and store bounding box ( for furthest point searches)
    minBound = tree.point;
    maxBound = tree.point;
    
    if ~isempty(tree.left)
        minBound = min(minBound, tree.left.bounds(1,:));
        maxBound = max(maxBound, tree.left.bounds(2,:));
    end
    
    if ~isempty(tree.right)
        minBound = min(minBound, tree.right.bounds(1,:));
        maxBound = max(maxBound, tree.right.bounds(2,:));
    end

    tree.bounds = [minBound; maxBound];
end

