% function tree = buildKDTree(points, depth)
% % BUILDKDTREE constructs a k-d tree for a given set of points.
% % Inputs:
% %   - points: Nx3 matrix representing 3D points.
% %   - depth: Current depth of the tree (start with 0).
% % Outputs:
% %   - tree: A structure representing the k-d tree.
% 
%     if isempty(points)
%         tree = [];
%         return;
%     end
% 
%     % Determine the axis to split on (cycle through 1, 2, 3)
%     axis = mod(depth, 3) + 1;
% 
%     % Sort points along the splitting axis
%     [sortedPoints, indices] = sortrows(points, axis);
% 
%     % Find the median index
%     medianIdx = floor(size(sortedPoints, 1) / 2) + 1;
% 
%     % Create the tree node
%     tree.point = sortedPoints(medianIdx, :);
%     tree.index = medianIdx;
%     tree.left = buildKDTree(sortedPoints(1:medianIdx-1, :), depth + 1);
%     tree.right = buildKDTree(sortedPoints(medianIdx+1:end, :), depth + 1);
%     tree.axis = axis;
% end

function tree = buildKDTree(points, depth, indices)
% BUILDKDTREE constructs a k-d tree for a given set of points.
% Inputs:
%   - points: Nx3 matrix representing 3D points.
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

    % Determine the axis to split on (cycle through 1, 2, 3)
    axis = mod(depth, 3) + 1;

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
end

