clear all; close all; clc;

%% Import Script for EBSD Data
%
% This script was automatically created by the import wizard. You should
% run the whoole script or parts of it in order to import your data. There
% is no problem in making any changes to this script.

%% Specify Crystal and Specimen Symmetries

% crystal symmetry
CS = {... 
  'notIndexed',...
  crystalSymmetry('432', [2.9 2.9 2.9], 'mineral', 'Iron - Alpha', 'color', [0.53 0.81 0.98])};

% plotting convention
setMTEXpref('xAxisDirection','north');
setMTEXpref('zAxisDirection','outOfPlane');

%% Specify File Names

% path to files
pname = '/Volumes/SANDISK/untitled folder 2';

% which files to be imported
fname = [pname '/5894-500-coarsened.ang'];

%% Import the Data

% create an EBSD variable containing the data
ebsd = EBSD.load(fname,CS,'interface','ang',...
  'convertEuler2SpatialReferenceFrame');

%dbstop if error
ang_colnames = {'Euler 1','Euler 2','Euler 3','X','Y','IQ','CI','Phase','SEM_signal','Fit'}
cs = crystalSymmetry('432', 'mineral', 'Iron - Alpha');
varargin = {'cs',cs,'bunge','radiant','ColumnNames',ang_colnames,'CS',{'notIndexed', cs},'convertEuler2SpatialReferenceFrame','InfoLevel',true,'CS',{'notIndexed', cs},'interface','ang','convertEuler2SpatialReferenceFrame','header',134};

[d,options,header,c] = load_generic(char(fname),varargin{:});
varargin = options;



loader = loadHelper(d,varargin{:});

unitCell = calcUnitCell(loader.getColumnData({'x' 'y' 'z'}), varargin{:});
%writematrix(loader.getColumnData({'x' 'y' 'z'}),'calcUnitCell_input_d.csv')
%writematrix(uc,'calcUnitCell_output_unitCell.csv')

X = [ebsd.prop.x(:), ebsd.prop.y(:)];
[V,F,I_FD] = spatialDecomposition(X, unitCell);
%writematrix(X,'spatialDecomposition_input_X.csv')
%writematrix(ebsd.unitCell,'spatialDecomposition_input_unitCell.csv')
%writematrix(V, 'spatialDecomposition_output_V')
%writematrix(F, 'spatialDecomposition_output_F')
%save('spatialDecomposition_output_I_FD',I_FD)

%% 
dummyCoordinates = [];

% specify a bounding polyogn
method = 'tight';

if ischar(method)
  
  switch lower(method)
    case 'tight'
      x = X(:,1);  y = X(:,2);      

      % the value 0.95 adjusts the boundary a little bit towards the convex hull
      delta = get_option(varargin,'tight',0.95,'double');
      k = boundary(x,y,delta);
      
      % erase all linear dependend points
      angle = atan2( x(k(1:end-1))-x(k(2:end)),...
        y(k(1:end-1))-y(k(2:end)) );      
      k = k([true; abs(diff(angle))>eps; true]);
      
      boundingX = X(k,:);

    case {'hull','convexhull'}
      x = X(:,1);  y = X(:,2);
      
      k = convhull(x,y);
      
      % erase all linear dependend points
      angle = atan2( x(k(1:end-1))-x(k(2:end)),...
        y(k(1:end-1))-y(k(2:end)) );      
      k = k([true; abs(diff(angle))>eps; true]);
            
      boundingX = X(k,:);
      
    case 'cube'
      
      % set up a rectangular box
      envelopeX = [min(X); max(X)];
      boundingX = [ ...
        envelopeX(1),envelopeX(3);
        envelopeX(2),envelopeX(3);
        envelopeX(2),envelopeX(4);
        envelopeX(1),envelopeX(4);
        envelopeX(1),envelopeX(3) ];
      
    otherwise     
      
      error('uknown boundary type. Available options are ''convexhull'' ''tight '' and ''cube''.');   
      
  end
  
elseif isa(method,'double')
  
  boundingX = method;
  
end


radius = mean(sqrt(sum(unitCell.^2,2)));
edgeLength = sqrt(sum(diff(boundingX).^2,2));

% fill each line segment with nodes every 20 points (in average)
nto = fix((edgeLength>0)*4); fix(edgeLength*(2*radius));

cs = cumsum([1; 1 + nto]);
boundingX(cs,:) = boundingX;

% interpolation
for k=1:numel(nto)
  for dim=1:2
    boundingX(cs(k):cs(k+1),dim) = ...
      linspace(boundingX(cs(k),dim), boundingX(cs(k+1), dim),nto(k)+2);
  end  
end


% homogeneous coordinates
X(:,3) = 1;
boundingX(:,3) = 1;

% householder matrix
H = @(v) eye(3) - 2./(v(:)'*v(:))*(v(:)*v(:)') ;
% translation matrix
T  = @(s) [ 1 0 s(1); 0 1 s(2); 0 0 1];

% direction of the edge
edgeDirection = diff(boundingX);
edgeAngle = atan2(edgeDirection(:,2),edgeDirection(:,1));
edgeLength = sqrt(sum(edgeDirection.^2,2));

% shift the starting vertex
bX = squeeze(double(axis2quat(zvector,edgeAngle)* ...
  vector3d([0; radius; 1])));
offsetX = bX - boundingX(1:end-1,:);

%writematrix(edgeDirection(k,:),'householder_input')
%writematrix(H(edgeDirection(k,:)),'householder_output')
%writematrix(offsetX(k,:),'translation_input')
%writematrix(T(offsetX(k,:)),'translation_output')

for k=1:size(boundingX,1)-1
  
  % mirror the point set X on each edge
  pX = X * -(T(offsetX(k,:)) * H(edgeDirection(k,:)) * T(offsetX(k,:)))';
  
  % distance between original and mirrowed point
  dist = sqrt(sum((X(:,1:2)-pX(:,1:2)).^2,2));
 
  intendX = 2*radius*sign(edgeDirection(k,1:2));
  
  % now try to delete unnecessary points
  m = 2;
  while true
    
    tmpX = pX(dist < m*radius,1:2);
    
    right = (bsxfun(@minus, tmpX, boundingX(k,1:2)  - intendX ) * edgeDirection(k,1:2)') < 0;
    left  = (bsxfun(@minus, tmpX, boundingX(k+1,1:2)+ intendX ) * edgeDirection(k,1:2)') > 0;
   
    tmpX = tmpX( ~(right | left) ,:);
     
    if edgeLength(k)/size(tmpX,1) < radius/3
      break;
    elseif m < 2^7
      m = m*2;
    elseif m < 2^7+100
      m = m+10; 
    else
      break;
    end
    
  end
  
  dummyCoordinates = [dummyCoordinates; tmpX];
  
end
  
dummyCoordinates = unique(dummyCoordinates,'first','rows');

% remove those points which are inside the b

id = inpolygon(dummyCoordinates(:,1),dummyCoordinates(:,2),boundingX(:,1),boundingX(:,2));

dummyCoordinates(id,:) = [];

writematrix(X, 'calcBoundary_input_X')
writematrix(unitCell, 'calcBoundary_input_unitCell')
writematrix(dummyCoordinates, 'calcBoundary_output_dummyCoordinates')