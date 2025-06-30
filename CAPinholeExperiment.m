function [ r ] = CAPinholeExperiment( args )

if ~exist('args', 'var')
    args = [];
end

% Initialize the OpticStudio connection
TheApplication = InitConnection();
if isempty(TheApplication)
    % failed to initialize a connection
    r = [];
else
    try
        r = BeginApplication(TheApplication, args);
        CleanupConnection(TheApplication);
    catch err
        CleanupConnection(TheApplication);
        rethrow(err);
    end
end
end


function [r] = BeginApplication(TheApplication, args)

import ZOSAPI.*;

    % Pinhole geometry
    % 5mm to mask
    % 512x512 res
    % 1.55um pixel pitch
    % 117 apertures
    % 200um pinhole diameter
    % 0.03mm pinhole thickness
    
    close all; clc;
    
    % Saving paths
    save_path_training = 'C:\training_data\coded\training\';
    save_path_validation = 'C:\training_data\coded\validation\';
    save_path_testing = 'C:\training_data\coded\testing\';
    
    pix_detector=512; % pix
    pix_dim = 1.55; %microm/pix
    ap_dia=0.2; % mm (larger aperture)
    ap_thickness = 0.03; % mm
    ap_fl = 5 - ap_thickness/2; % mm
    ap_d = ap_fl + ap_thickness/2;
    dtheta=pix_dim/(ap_d*1e3);
    
    % Create the spot center vector to loop over
    [pix_locX,pix_locY] = meshgrid(-pix_detector/2:pix_detector/2,...
        -pix_detector/2:pix_detector/2);
    pix_locxy = [pix_locX(:), pix_locY(:)];
    iter_max = length(pix_locxy); % Max iterations

    % Create unshuffled seed
    setlena = round(0.8*iter_max);
    setlenb = round(0.1*iter_max);
    setlenc = iter_max - setlena - setlenb;
    seeda = 1*ones(1,setlena); % Training set unshuffled
    seedb = 2*ones(1,setlenb); % Validation set unshuffled
    seedc = 3*ones(1,setlenc); % Testing set unshuffled
    seed = [seeda seedb seedc]; % Full set unshuffled
    
    % Create random seeds
    seed0 = seed(randperm(length(seed))); % Shuffle case 0
    seed1 = seed(randperm(length(seed))); % Shuffle case 1
    seed2 = seed(randperm(length(seed))); % Shuffle case 2
    
    % Save seed incase the run stops
    seedsave.a =  seed0;
    seedsave.b =  seed1;
    seedsave.c =  seed2;
    save('CAPinholeExperimentSeed.mat','seedsave');
    
    % Un-comment if you need to load a unfinished run
    %load('SAPinholeExperimentSeed.mat');
    %seed0 = seedsave.a;
    %seed1 = seedsave.b;
    %seed2 = seedsave.c;
    
    % Calculate the pixel size
    % length = angle x focal length
    delta_calc = pix_dim/(ap_fl+ap_thickness/2); %mrad
    
    % Number of apertures on mask (ap_N x ap_N)
    % 1 aperture in this case
    ap_N = 13; % Per side
    ap_pitch = 75; % pix
    
    % If there is only one aperture, then the pitch should be zero
    if ap_N == 1
        ap_pitch = 0;
    end
    
     % Set up primary optical system
    TheSystem = TheApplication.PrimarySystem;
    
    %Define System Explorer
    % ISystemData represents the System Explorer in GUI.
    % We access options in System Explorer through ISystemData in ZOS-API
    SystExplorer = TheSystem.SystemData;
    
    % Change aperture diameter
    SystExplorer.Aperture.ApertureType =...
        ZOSAPI.SystemData.ZemaxApertureType.EntrancePupilDiameter;
    SystExplorer.Aperture.ApertureValue = ap_dia;
    
    % Switch to afocal mode
    SystExplorer.Aperture.AFocalImageSpace = true;
    
    % Create AM0 spectrum
    SystExplorer.Wavelengths.RemoveWavelength(1);
    SystExplorer.Wavelengths.AddWavelength(0.450,4.443927590000E-01);
    SystExplorer.Wavelengths.AddWavelength(0.506,1.000000000000E+00);
    SystExplorer.Wavelengths.AddWavelength(0.550,7.883974690000E-01);
    SystExplorer.Wavelengths.AddWavelength(0.600,6.884336030000E-01);
    SystExplorer.Wavelengths.AddWavelength(0.650,3.562352800000E-01);
    
    % Change units to Milliradians
    SystExplorer.Units.AfocalModeUnits =...
        ZOSAPI.SystemData.ZemaxAfocalModeUnits.Milliradians;
    
    % Get interface of Lens Data Editor and add 3 surfaces
    TheLDE = TheSystem.LDE;
    TheLDE.InsertNewSurfaceAt(2);
    TheLDE.InsertNewSurfaceAt(2);
    
    % Get the 5 LDE surfaces
    Surface0 = TheLDE.GetSurfaceAt(0); % Field
    Surface1 = TheLDE.GetSurfaceAt(1); % Dummy surface
    Surface2 = TheLDE.GetSurfaceAt(2); % Pinhole 1
    Surface3 = TheLDE.GetSurfaceAt(3); % Pinhole 2
    Surface4 = TheLDE.GetSurfaceAt(4); % Image plane
    
    % Set thickness for each surface
    Surface2.Thickness = ap_thickness; % mm
    Surface3.Thickness = ap_fl; % mm
    
    % Set diameter for each surface
    Surface2.MechanicalSemiDiameter = ap_dia/2;
    Surface3.MechanicalSemiDiameter = ap_dia/2;
    Surface4.SemiDiameter = 5.5;
    Surface4.MechanicalSemiDiameter = 5.5;
    Surface4.SemiDiameterCell.MakeSolveFixed();
    
     % Set stop surface
    Surface2.IsStop = true;
    
    % Create top aperture
    Circ_Aper_Top = Surface2.ApertureData.CreateApertureTypeSettings(...
        ZOSAPI.Editors.LDE.SurfaceApertureTypes.CircularAperture);
    Circ_Aper_Top.S_CircularAperture_.MinimumRadius = 0;
    Circ_Aper_Top.S_CircularAperture_.MaximumRadius = ap_dia/2;
    Surface2.ApertureData.ChangeApertureTypeSettings(Circ_Aper_Top);
    Surface2.ApertureData.PickupFrom = 0;
    
    % Create bottom aperture
    Circ_Aper_Bottom = Surface3.ApertureData.CreateApertureTypeSettings(...
        ZOSAPI.Editors.LDE.SurfaceApertureTypes.CircularAperture);
    Circ_Aper_Bottom.S_CircularAperture_.MinimumRadius = 0;
    Circ_Aper_Bottom.S_CircularAperture_.MaximumRadius = ap_dia/2;
    Surface3.ApertureData.ChangeApertureTypeSettings(Circ_Aper_Bottom);
    Surface3.ApertureData.PickupFrom = 2;
    
    % Retain untilted LDE configuration
    UntiltedLDE = TheLDE;
    
    for imgiter = 1:iter_max
        
    % Return configuration to untilted position
    TheLDE = UntiltedLDE;

    % Establish spot center
    pix_xc = pix_locxy(imgiter,1); % pix
    pix_yc = pix_locxy(imgiter,2); % pix
    
    % Obtain 2 sun angles from spot center (in degrees)
    alpha=atand(pix_xc*dtheta); % deg
    beta=atand(pix_yc*dtheta); % deg
    
    % Obtain 2 sun angles from spot center (in radians)
    alphar = round(atan(pix_xc*dtheta),6); % rad
    betar = round(atan(pix_yc*dtheta),6); % rad
    
    % Tilt the incident light
    TiltConfig = TheLDE.GetTool_TiltDecenterElements();
    TiltConfig.FirstSurface = 1; 
    TiltConfig.LastSurface = 1;
    TiltConfig.TiltX = alpha;
    TiltConfig.TiltY = beta;
    TheLDE.RunTool_TiltDecenterElements(TiltConfig);
    TheLDE.RemoveSurfacesAt(2,2);
    % Checking angles
    %TheLDE.GetSurfaceAt(1).GetSurfaceCell(ZOSAPI.Editors.LDE.SurfaceColumn.Par3).DoubleValue
    %TheLDE.GetSurfaceAt(1).GetSurfaceCell(ZOSAPI.Editors.LDE.SurfaceColumn.Par4).DoubleValue
        
    % Create analysis
    newWin = TheSystem.Analyses.New_HuygensPsf();
    
    % Settings
    newWin_Settings = newWin.GetSettings();
    newWin_Settings.ShowAsType =...
        ZOSAPI.Analysis.HuygensShowAsTypes.InverseGreyScale;
    newWin_Settings.PupilSampleSize =...
        ZOSAPI.Analysis.SampleSizes.S_32x32;
    newWin_Settings.ImageSampleSize =...
        ZOSAPI.Analysis.SampleSizes.S_64x64; % Smaller spot
   
    newWin_Settings.ImageDelta = delta_calc;

    % Run Analysis & gets results
    newWin.ApplyAndWaitForCompletion();
    newWin_Results = newWin.GetResults();

    newWin_Values =  newWin_Results.GetDataGrid(0).Values.double;
    
    numPix = cast(newWin_Results.GetDataGrid(0).Nx,"double");
    %minXVal = newWin_Results.GetDataGrid(0).MinX;
    %minYVal = newWin_Results.GetDataGrid(0).MinY;
    %sizePix = newWin_Results.GetDataGrid(0).Dx; % Pixel size
    
    % Close the analysis window after obtaining results
    newWin.Close();

%     figure(1);
%     subplot(2,2,1);
    %x = [-numPix/2 numPix/2]; %y = x;
    imgarr = zeros(numPix + 1, numPix + 1);
    imgarr(2:end,1:end-1)=flipud(newWin_Values);
%     imagesc(-x,y,imgarr);
%     set(gca,'YDir','normal');
%     colormap('jet');
%     colorbar;
%     axis square;
%     xlabel('X pixels');
%     ylabel('Y pixels');
%     title('Projected sun spot');
%    
%     subplot(2,2,2);
%     grid on;

    % Read and plot data series of cross-section
    %newWinX = TheSystem.Analyses.New_HuygensPsfCrossSection();
    
    % Settings for cross-section
    %newWin_SettingsX = newWinX.GetSettings();
    %newWin_SettingsX.PupilSampleSize = ZOSAPI.Analysis.SampleSizes.S_32x32;
    %newWin_SettingsX.ImageSampleSize = ZOSAPI.Analysis.SampleSizes.S_128x128;
    %newWin_SettingsX.ImageDelta = delta_calc;
    
    % Run Analysis & gets results
%     newWinX.ApplyAndWaitForCompletion();
%     newWin_ResultsX = newWinX.GetResults();
%     dataSeries = newWin_ResultsX.DataSeries;
%     cc=lines(double(newWin_ResultsX.NumberOfDataSeries));
%     for gridN=1:newWin_ResultsX.NumberOfDataSeries
%         data = dataSeries(gridN);
%         y_xsec = data.YData.Data.double;
%         x_xsec = data.XData.Data.double;
%         plot(x_xsec,y_xsec,'-','color',cc(gridN,:));
%     end
%     axis square;
%     xlabel('X pixels');
%     ylabel('Relative Intensity');
%     title('Huygens Cross Section');
    
    %xd = [-pix_detector/2 pix_detector/2]; %yd = xd;
    
    % Blank sensor initialization
    imgarr_full = zeros(pix_detector+1,pix_detector+1);
    
    % Spot template verticies
    %spot_temp = [1 1; 1 numPix + 1; numPix + 1 1; numPix + 1 numPix + 1];
    
    % Code mapping
    cg_index = 0:5:60;
    cg_rows = [fliplr(cg_index(1:2:end))'; cg_index(2:2:end)'];
    cg_cols = [15 0 30];
    codegroups = 3;
    cg_code = zeros(ap_N, 2, codegroups);
    
    % Loop over coding groups
    for cg = 1:codegroups
        cg_gap = [-300 0 300];
        cg_code(:,:,cg) = [cg_cols(cg) .* ones(length(cg_rows),1) cg_rows];
        
        % Loop over mask entries
        for i = 4:6
            ap_deltax_base = (i - 9/2 - 1/2)*ap_pitch + cg_gap(cg);
            for j = 1:ap_N
                if i == 4
                    ap_deltax = ap_deltax_base - cg_code(j,1,cg);
                elseif i == 6
                    ap_deltax = ap_deltax_base + cg_code(j,2,cg);
                else
                    ap_deltax = ap_deltax_base;
                end
                ap_deltay = (j - ap_N/2 - 1/2)*ap_pitch;
                
                % Initialize spot verticies
                spot_vert = zeros(4,2);
                
                % Relative to detector coordinates
                spot_vert(1,:) = [pix_yc+ap_deltay+numPix/2,...
                    pix_xc+ap_deltax-numPix/2]; % Top left
                spot_vert(2,:) = [pix_yc+ap_deltay+numPix/2,...
                    pix_xc+ap_deltax+numPix/2]; % Top right
                spot_vert(3,:) = [pix_yc+ap_deltay-numPix/2,...
                    pix_xc+ap_deltax-numPix/2]; % Bottom left
                spot_vert(4,:) = [pix_yc+ap_deltay-numPix/2,...
                    pix_xc+ap_deltax+numPix/2]; % Bottom right
                
                % Check spot bounds
                % If the spot is not within the detector bounds then skip
                if all(or(abs(spot_vert(:,1))>pix_detector/2,...
                        abs(spot_vert(:,2))>pix_detector/2),'all')
                    %disp(['Aperture ' num2str(i + j - 1) ' not in bounds']);
                    continue;
                elseif ~all(abs(spot_vert)>pix_detector/2,'all') &&...
                        ~any(abs(spot_vert)>pix_detector/2,'all')
                    %disp(['Aperture ' num2str(i + j - 1) ' is in bounds']);
                elseif ~all(or(abs(spot_vert(:,1))>pix_detector/2,...
                        abs(spot_vert(:,2))>pix_detector/2),'all') &&...
                        any(abs(spot_vert)>pix_detector/2,'all')
                    %disp(['Aperture ' num2str(i + j - 1) ' is partially in bounds']);
                end
                
                % Create vertex array for the cropped spot
                crop_vert = spot_vert;
                
                % Limit spot verticies to detector space boundaries
                crop_vert(abs(crop_vert)>pix_detector/2)=...
                    pix_detector/2.*sign(crop_vert(abs(...
                    crop_vert)>pix_detector/2));
                
                half_det = pix_detector/2 + 1;
                
                % Find the cropped detector frame verticies
                det_crop = cen2cor(crop_vert(:,1),crop_vert(:,2),half_det);
                
                % Find the cropped spot frame verticies
                v1 = cen2cor(spot_vert(:,1),spot_vert(:,2),half_det);
                v2 = cen2cor(crop_vert(:,1),crop_vert(:,2),half_det);
                spot_crop = [v2(:,1) - v1(1,1) + 1, v2(:,2) - v1(1,2) + 1];
                
                % Crop image spot from given verticies
                imgarr_crop=imcrop(imgarr,[spot_crop(1,2), spot_crop(1,1),...
                    spot_crop(2,2) - spot_crop(1,2),...
                    spot_crop(3,1) - spot_crop(1,1)]);
                
                imgr=det_crop(1,1):det_crop(3,1);
                imgc=det_crop(1,2):det_crop(2,2);
                
                % Add cropped image to detector space
                imgarr_full(imgr, imgc) = imgarr_crop;
            end
        end
    end
    
    % Convert to 8-bit grayscale image
    imgarr_full = uint8(255 * imgarr_full);

%     figure(1);
%     imcontour(imgarr_full,100);
    
%     figure(2);
%     plot(imgarr_full(961,1:length(imgarr_full)),'k');
    
    % Save path key
    % alpha_beta_0 = raw
    % alpha_beta_1 = noise
    % alpha_beta_2 = noise + thresholds
    % Number of raw images is 512x512 = 262,144
    % Raw 0 + Noise 1 + Noise & Threshold 2 = 262,144 x 3 = 786,432
    % Training images = 0.8 x 786,432 = 629,150
    % Validation / Testing images = 78,643
    
    % Random path selection from seed 0
    if seed0(imgiter) == 1
        save_path = save_path_training;
    elseif seed0(imgiter) == 2
        save_path = save_path_validation;
    else
        save_path = save_path_testing;
    end
    
    % Save no noise raw image
    % alpha_beta_0
    imwrite(imgarr_full,[save_path num2str(alphar) '_' num2str(betar) '_0.png']);
    
    % Add noise
    imgarr_noise = addNoise(imgarr_full);
    
    % Random path selection from seed 1
    if seed1(imgiter) == 1
        save_path = save_path_training;
    elseif seed1(imgiter) == 2
        save_path = save_path_validation;
    else
        save_path = save_path_testing;
    end
    
    % Save noisy image
    % alpha_beta_1
    imwrite(imgarr_noise,[save_path num2str(alphar) '_' num2str(betar) '_1.png']);
    
    % Threshold noisy image
    %thresh_l = floor(0.1.*255);
    thresh_u = ceil(0.9.*255);
    imgarr_noise(imgarr_noise>thresh_u) = thresh_u;
    %imgarr_noise(imgarr_noise<thresh_l) = 0;
    imgarr_noise = uint8(rescale(imgarr_noise,0,255));
    
    % Random path selection from seed 2
    if seed2(imgiter) == 1
        save_path = save_path_training;
    elseif seed2(imgiter) == 2
        save_path = save_path_validation;
    else
        save_path = save_path_testing;
    end
    
    % Save noisy threshold image
    % alpha_beta_2
    imwrite(imgarr_noise,[save_path num2str(alphar) '_' num2str(betar) '_2.png']);
    
    end
    
    %hold on;
    %plot(imgarr_noise(256-100,256+100-40:256+100+40),'b--')

%     subplot(2,2,3);
%     imagesc(xd,-yd,imgarr_noise);
%     set(gca,'YDir','normal');
%     colormap('gray');
%     colorbar;
%     caxis([0 255]);
%     axis square;
%     xlabel('X pixels');
%     ylabel('Y pixels');
%     title('Full Detector');

TheSystem.Close(false);

     r = [];
end

function app = InitConnection()

import System.Reflection.*;

% Find the installed version of OpticStudio.
zemaxData = winqueryreg('HKEY_CURRENT_USER', 'Software\Zemax', 'ZemaxRoot');
NetHelper = strcat(zemaxData, '\ZOS-API\Libraries\ZOSAPI_NetHelper.dll');
% Note -- uncomment the following line to use a custom NetHelper path
% NetHelper = 'C:\Users\Documents\Zemax\ZOS-API\Libraries\ZOSAPI_NetHelper.dll';
NET.addAssembly(NetHelper);

success = ZOSAPI_NetHelper.ZOSAPI_Initializer.Initialize();
% Note -- uncomment the following line to use a custom initialization path
% success = ZOSAPI_NetHelper.ZOSAPI_Initializer.Initialize('C:\Program Files\OpticStudio\');
if success == 1
    LogMessage(strcat('Found OpticStudio at: ', char(ZOSAPI_NetHelper.ZOSAPI_Initializer.GetZemaxDirectory())));
else
    app = [];
    return;
end

% Now load the ZOS-API assemblies
NET.addAssembly(AssemblyName('ZOSAPI_Interfaces'));
NET.addAssembly(AssemblyName('ZOSAPI'));

% Create the initial connection class
TheConnection = ZOSAPI.ZOSAPI_Connection();

% Attempt to create a Standalone connection

% NOTE - if this fails with a message like 'Unable to load one or more of
% the requested types', it is usually caused by try to connect to a 32-bit
% version of OpticStudio from a 64-bit version of MATLAB (or vice-versa).
% This is an issue with how MATLAB interfaces with .NET, and the only
% current workaround is to use 32- or 64-bit versions of both applications.
app = TheConnection.CreateNewApplication();
if isempty(app)
   HandleError('An unknown connection error occurred!');
end
if ~app.IsValidLicenseForAPI
    HandleError('License check failed!');
    app = [];
end

end

function LogMessage(msg)
disp(msg);
end

function HandleError(error)
ME = MXException(error);
throw(ME);
end

function  CleanupConnection(TheApplication)
% Note - this will close down the connection.

% If you want to keep the application open, you should skip this step
% and store the instance somewhere instead.
TheApplication.CloseApplication();
end


