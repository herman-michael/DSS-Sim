close all; clear; clc;

% Mask type
type = 'coded';

%% DEFINE PATHS

% Training paths
myFolderTrain = ['C:\training_data\' type '\train'];
myPathTrain = ['dataset/' type '/train/'];

% Testing paths
myFolderTest = ['C:\training_data\' type '\test'];
myPathTest = ['dataset/' type '/test/'];

% Validation paths
myFolderValid = ['C:\training_data\' type '\valid'];
myPathValid = ['dataset/' type '/valid/'];

%% TRAINING TARGET SAVE

% Check to make sure that folder actually exists
if ~isfolder(myFolderTrain)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s\nPlease specify a new folder.', myFolderTrain);
    uiwait(warndlg(errorMessage));
    myFolderTrain = uigetdir(); % Ask for a new one
    if myFolderTrain == 0
         % User clicked Cancel
         return;
    end
end

% Get a list of all files in the folder with the desired file name pattern
filePatternTrain = fullfile(myFolderTrain, '*.png');
theFilesTrain = dir(filePatternTrain);

dataCellArrTrain = cell(length(theFilesTrain)+1,1);
dataCellArrTrain{1} = 'image_path,x,y';

for k = 1 : length(theFilesTrain)
    baseFileName = theFilesTrain(k).name;
    dataCellSplit = split(baseFileName,'_');
    dataCellArrTrain{k+1} = [myPathTrain baseFileName ',' dataCellSplit{1} ','...
        dataCellSplit{2}];
end

% Save the target data as csv
writecell(dataCellArrTrain,'ca_train.csv','QuoteStrings',false)

%% TESTING TARGET SAVE

% Check to make sure that folder actually exists
if ~isfolder(myFolderTest)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s\nPlease specify a new folder.', myFolderTest);
    uiwait(warndlg(errorMessage));
    myFolderTest = uigetdir(); % Ask for a new one
    if myFolderTest == 0
         % User clicked Cancel
         return;
    end
end

% Get a list of all files in the folder with the desired file name pattern
filePatternTest = fullfile(myFolderTest, '*.png');
theFilesTest = dir(filePatternTest);

dataCellArrTest = cell(length(theFilesTest)+1,1);
dataCellArrTest{1} = 'image_path,x,y';

for k = 1 : length(theFilesTest)
    baseFileName = theFilesTest(k).name;
    dataCellSplit = split(baseFileName,'_');
    dataCellArrTest{k+1} = [myPathTest baseFileName ',' dataCellSplit{1} ','...
        dataCellSplit{2}];
end

% Save the target data as csv
writecell(dataCellArrTest,'ca_test.csv','QuoteStrings',false)

%% VALIDATION TARGET SAVE

% Check to make sure that folder actually exists
if ~isfolder(myFolderValid)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s\nPlease specify a new folder.', myFolderValid);
    uiwait(warndlg(errorMessage));
    myFolderValid = uigetdir(); % Ask for a new one
    if myFolderValid == 0
         % User clicked Cancel
         return;
    end
end

% Get a list of all files in the folder with the desired file name pattern
filePatternValid = fullfile(myFolderValid, '*.png');
theFilesValid = dir(filePatternValid);

dataCellArrValid = cell(length(theFilesValid)+1,1);
dataCellArrValid{1} = 'image_path,x,y';

for k = 1 : length(theFilesValid)
    baseFileName = theFilesValid(k).name;
    dataCellSplit = split(baseFileName,'_');
    dataCellArrValid{k+1} = [myPathValid baseFileName ',' dataCellSplit{1} ','...
        dataCellSplit{2}];
end

% Save the target data as csv
writecell(dataCellArrValid,'ca_valid.csv','QuoteStrings',false)

disp('Done');