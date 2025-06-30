close all; clear; clc;

% Mask type
%type = 'singleaperture_noaug\';
%type = 'multiaperture_noaug\';
type = 'multiaperture\';

%% DEFINE PATHS

% Data paths
myFolderTrain = ['C:\Users\mherman7\Documents\Research_Code\SynthGen\2input\' type];

%name = {'sa_test.csv'; 'sa_train.csv'; 'sa_valid.csv'};
name = {'ma_test.csv'; 'ma_train.csv'; 'ma_valid.csv'};
%name = {'ca_test.csv'; 'ca_train.csv'; 'ca_valid.csv'};

num_arr=[];

for i = 1:3
    
    % Read csv file
    [num,txt,raw] = xlsread([myFolderTrain name{i}]);
    disp(length(num))
    
    num_arr = [num_arr;num];
    
    % Initialize the class vector
    out = cell(length(num),1);
    
    % Set class value for quadrant 1
    out(num(:,1)>=0 & num(:,2)>=0) = {1};
    % Set class value for quadrant 2
    out(num(:,1)<0 & num(:,2)>0) = {2};
    % Set class value for quadrant 3
    out(num(:,1)<=0 & num(:,2)<=0) = {3};
    % Set class value for quadrant 4
    out(num(:,1)>0 & num(:,2)<0) = {4};
    
    % Origin class
    out(num(:,1)==0 & num(:,2)==0) = {1};
    
    % Add class header
    out = [{'class'}; out];
    
    % Create new cell with classes added
    new = [raw(:,1) out raw(:,2:end)];
    
    % Save new csv fie
    writecell(new,[myFolderTrain 'output\' name{i}],'QuoteStrings',false)
    
end