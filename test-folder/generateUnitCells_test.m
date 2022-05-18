load('spatialDecomp_test.mat')
X = X(:,1:2);

xy = X;

%function [v,faces] = generateUnitCells(xy,unitCell,varargin)
% generate a list of patches according to spatial coordinates and the unitCell
%
% Input
%  xy       - midpoints of the cells
%  unitCell - spatial coordinates of the unit cell
%
% Ouput
%  v     - list of vertices
%  faces - list of faces

% compute the vertices
x = reshape(bsxfun(@plus,xy(:,1),unitCell(:,1).'),[],1);
y = reshape(bsxfun(@plus,xy(:,2),unitCell(:,2).'),[],1);

% remove equal points
eps = min(sqrt(diff(unitCell(:,1)).^2 + diff(unitCell(:,2)).^2))/10;
%xi = round(x-min(x)./eps);
%yi = round(y-min(y)./eps);
%[~,m,n] = unique((1+xi)*2*max(yi) + yi);
verts = round([x-min(x) y-min(y)]./eps);
% aa = [1, 2; 1, 2; 3, 5]
% [v,m,n] = unique(aa,'rows');
[v,m,n] = unique(verts,'rows');

[vnew, order] = sortrows(v);

x2 = x(m);
y2 = y(m);
v2 = [x2(order) y2(order)];


% set faces
faces = reshape(n, [], size(unitCell,1));
faces2 = reshape(order(n), [], size(unitCell,1));