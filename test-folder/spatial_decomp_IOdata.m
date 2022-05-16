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

uc = calcUnitCell(loader.getColumnData({'x' 'y' 'z'}), varargin{:});
writematrix(loader.getColumnData({'x' 'y' 'z'}),'calcUnitCell_input_d.csv')
writematrix(uc,'calcUnitCell_output_unitCell.csv')



% %%

% 
% 
% z = loader.getColumnData({'x' 'y' 'z'})
% calcUnitCell(loader.getColumnData({'x' 'y' 'z'}))
% %csvwrite('

