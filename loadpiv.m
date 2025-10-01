function D = loadpiv(folderPIV, varargin)
% Extracts and formats data from "vc7" files contained in "folderPIV".
%
% NOTES:
% - Either function MUST be contained in the same directory as the
%   "readimx_WIN" folder, or you must have the readimx in your PATH.
% - Whenever new updates are made PLEASE UPDATE THE DOCUMENTATION AND ADD A
%   VERSION DESCRIPTION IN THE FORMAT SPECIFIED.
% -------------------------------------------------------------------------
% Sintax:
% 
% -------------------------------------------------------------------------
% >> Mandatory inputs:
%
% D = loadpiv(folderPIV) : Extracts PIV data from directory "folderPIV".
%           The directory must contain the ".vc7" files. Only extracts
%           coordinates, velocity components, and z-vorticity. Depending on
%           the dataset, the extracted variables might include the
%           following:
%               D.x : x-coordinates matrix
%               D.y : y-coordinates matrix
%               D.u : x-velocity component
%               D.v : y-velocity component
%               D.w : z-velocity component
%               D.vort : z-vorticity component
%
% >> Optional inputs:
%
% D = loadpiv(__, params, "nondim") : Extracts and non-dimensionalizes
%           data using parameters contained in structure "params".
%           "params" must have the form:
%                   params.L = characteristic_length
%                   params.U = characteristic_velocity
%
% D = loadpiv(__, params, "extractAllVariables") : Extracts all data
%           found in the "vc7" files. Depending on the dataset, the extract
%           ed variables might include the following:
%               D.x : x-coordinates matrix
%               D.y : y-coordinates matrix
%               D.u : x-velocity component
%               D.v : y-velocity component
%               D.w : z-velocity component
%               D.vort : z-vorticity component
%               D.corr : correlation values
%               D.uncU : x-velocity uncertainty
%               D.uncV : y-velocity uncertainty
%               D.uncW : z-velocity uncertainty
%
% >> Name-value arguments:
%
% D = loadpiv(__, "numCamFields", numCamera) : Number of independent fields
%               from different cameras. Used for processing.
%
% D = loadpiv(__, "frameRange", frames) : Specifies specific frames to
%               be extracted instead of the full dataset. Must be followed
%               by a 1D array containing the frame numbers.
%
% D = loadpiv(__, "fovCenter", [xcntr, ycntr]) : Specifies the location
%               of the desired origin of the coordinate axis. Must be
%               followed by a [x0,y0] array specified in meters.
%
% D = loadpiv(__, "fovRot", rotAngle) : Specifies the rotation of the
%               fov with respect to the origin. Must be followed by the
%               angle's value in radians.
%
% D = loadpiv(__, "Validate", minCorrelationValue) : sets a minimum
%               correlation value, otherwise put NaN.
%
% -------------------------------------------------------------------------
% >> Output data - all in a data structure "D":
%
% x : array containing x-coordinates with size [n,m].
%
% y : array containing y-coordinates with size [n,m].
%
% u : array containing x-component of flow velocity with size [n,m,N].
%
% v : array containing y-component of flow velocity with size [n,m,N].
%
% w : array containing z-component of flow velocity with size [n,m,N].
%
% vort : array containing z-component of vorticity with size [n,m,N].
%
% corr : array containing correlation values with size [n,m,N].
%
% uncU : array containing u-velocity uncertainty with size [n,m,N].
%
% uncV : array containing v-velocity uncertainty with size [n,m,N].
%
% uncW : array containing z-velocity uncertainty with size [n,m,N].
%
% -------------------------------------------------------------------------
%  UPDATES:
%
% Version 1.0.0
% 2024/10/31 - Eric Handy (with snippets from Alex, Siyang and Kenny)
% 
% Version 1.1.0
% 2024/11/06 - Eric Handy
% - Bug fixes.
% - Changed vorticity calculation from MATLAB's "curl()" function to the
%   algorithm specified by Raffel's PIV Handbook.
%
% Version 1.1.1
% 2025/05/12 - Eric Handy
% - Fixed bug that immediately identified any second parameter as "params".
% - Fixed a bug that did not assign all output variables specified when
%   using the "extractAllVariables" option.
% NOTE: These issues have not been robustly tested.
%
% Version 2.0.0
% 2025/08/06 - Kenny Breuer
% - Added validate option that returns data only if it meets a correlation
%   threshhold.
% - Output data is NaN if there is no value (instead of 0).
% - Returns results in a structre - solves several problems with ordering
%   output results.
% - Doesnt return "isvalid" - that doesnt seem to be present in all the PIV
%   data.
%
% Version 2.1.0
% 2025/09/30 - Pedro C Ormonde
% - Fixed hard-coded indices for field names ('U0','V0','TS:Correlation
%   value', etc). Now function reads from available field names.
% - Added optional input "n_camField" if data contains multiple monoPIV 
%   fields from multiple cameras. Default: n_camField = 1.
%
% =========================================================================
% SINTAX FOR UPDATES:
% Version N1.N2.N3
% YYYY/MM/DD - Name Last
% - N1 for large structural overhauls to the function (new data structure,
%   new input/output format, new algorithm, etc.).
% - N2 for small additions to fn (new arguments, new outputs etc.).
% - N3 for bug fixes and compatibility updates.
% =========================================================================

%% Check input parameters

% Default parameters:
CorrelationThreshold = 1.0;

% Parse inputs:
for optionNum = 1:length(varargin)
    var_option = varargin{optionNum};

    % Name-value arguments:
    if ~isnumeric(var_option) && ~isstruct(var_option)
        switch var_option
            case 'frameRange' % to extract specific frames
                frameRange = varargin{optionNum+1};
                if ~isnumeric(frameRange)
                    error('Value of "frameRange" must be a numerical integer array.');
                end

            case 'nondim' % non-dimensionalize data option
                nondim = 1;

            case 'extractAllVariables' % extract all optional variables
                eav = 1;

            case 'fovCenter' % specify coordinate center
                fovCenter = varargin{optionNum+1};
                if ~isvector(fovCenter) || ~isnumeric(fovCenter)
                    error('"fovCenter" coordinates invalid, must be a numeric vector [x, y].')
                end

            case 'fovRot' % specify coordinate rotation
                fovRot = varargin{optionNum+1};
                if ~isnumeric(fovRot)
                    error('"fovRot" value invalid, must be a numeric value in radians.')
                end

            case 'Validate' % specify minimum vector correlation
                CorrelationThreshold = varargin{optionNum+1};
                fprintf('Minimium correlation: %6.2f\n', CorrelationThreshold);
                
            case 'numCamFields' % specify which camera to extract data from
                n_camField = varargin{optionNum+1};
                fprintf('Extracting fields for camera: %1.0d\n',n_camField);

            otherwise
                error(['Invalid name-value argument: "',var_option,'" for input "options."'])
        end
    % Params structure:
    elseif isstruct(var_option)
        params = var_option;
    else
        error(['Input not recognized:',newline,var_option]);
    end
end


% Check for PIV parameters variable
if ~exist('params','Var') || isempty(params)
    params = struct;
end

% Check for "eav" (Extract All Variables) flag
if ~exist('eav','Var')
    eav = 0;
end

% Check for Field of View characteristics
if ~exist('fovCenter','Var')
    disp('WARNING: "fovCenter" not specified.')
    fovCenter = [0,0];
end

if ~exist('fovRot','Var')
    disp('WARNING: "fovRot" not specified.')
    fovRot = 0;
end

% Check for n_camFields variable
if ~exist('n_camField','Var')
    warning('"numCamField" not specified. Default to numCamField = 1.')
    n_camField = 1;
end

% Check for characteristic parameters and extract
if ~exist('nondim','Var')
    L = 1; U = 1;
elseif nondim == 1
    if ~isfield(params,'L') || ~isfield(params,'U')
        warning('"L" (length) and "U" (velocity) fields in struct "params" not passed. Output will be dimensional.')
        L = 1; U = 1;
    elseif ~isnumeric(params.L) || ~isnumeric(params.U)
        warning('Values for "params.L" and/or "params.U" invalid, values must be numeric and in SI units. Output will be dimensional.')
        L = 1; U = 1;
    else
        L = params.L; U = params.U;
    end
end

%% Load data

% locate current directory
folderMain = cd();

% Locate ".vc7" files in PIV directory
files = dir(fullfile(folderPIV,'*.vc7'));
if isempty(files)
    error('No ".vc7" files found in specified directory "folderPIV". Check path.')
end

% Determine range of frames to be processed
if ~exist('frameRange','Var') || isempty(frameRange)
    frameRange = 1:length(files);
end

%% Add path of IMX executable
if exist('readimx', 'file') ~= 3  % check if function is not yet defined
    if ismac
        % addpath(genpath('readimx_MAC'));
        cd('readimx_MAC');
    elseif ispc
        % addpath(genpath('readimx_WIN'));
        cd('readimx_WIN');
    else
        error('System not identified (if using LINUX I do not know what to do).')
    end
end

%% Read in the files
fnum = 0;
N = length(frameRange);

% Loop through each file
for file = frameRange
    fnum = fnum + 1;

    % Read DaVis file
    fname = readimx(fullfile(folderPIV,files(file).name));

    % Extract data from file - returns data structure from one PIV frame
    out1 = extractData(fname, CorrelationThreshold,n_camField);

    % Initialize variables
    if fnum == 1
        if out1.dimNum > 3
            error('Data structure constains unrecognized data structure, number of dimensions exceeds 3D data structure. Review data structure');
        end
        disp([num2str(out1.dimNum),'-dimensional data available.'])

        %% Rotate coordinate system if needed
        D.x = ( (out1.xRaw-fovCenter(1)).*cos(fovRot) - (out1.yRaw-fovCenter(2)).*sin(fovRot) );
        D.y = ( (out1.xRaw-fovCenter(1)).*sin(fovRot) + (out1.yRaw-fovCenter(2)).*cos(fovRot) );

        D.u    = nan([size(D.x),N]);
        D.v    = nan([size(D.x),N]);
        if out1.dimNum == 3
           D.w = nan([size(D.x),N]);
        end

        D.vort = nan([size(D.x),N]);

        if eav == 1
            D.corr = nan([size(D.x),N]);
            D.uncU = nan([size(D.x),N]);
            D.uncV = nan([size(D.x),N]);
            if out1.dimNum == 3
               D.uncW = nan([size(D.x),N]);
            end
        end
    end

    % Store data into arrays with rotation
    D.u(:,:,fnum) = ( out1.uRaw.*cos(fovRot) - out1.vRaw.*sin(fovRot) );
    D.v(:,:,fnum) = ( out1.uRaw.*sin(fovRot) + out1.vRaw.*cos(fovRot) );

    D.vort(:,:,fnum) = out1.vortRaw;

    if eav == 1
        D.corr(:,:,fnum) = out1.corrRaw;
        D.uncU(:,:,fnum) = out1.uncURaw;
        D.uncV(:,:,fnum) = out1.uncVRaw;
    end

    % For 3-dimensional data
    if out1.dimNum == 3
        D.w(:,:,fnum) = out1.wRaw;
        if eav == 1
            D.uncW(:,:,fnum) = out1.uncWRaw;
        end
    end

    % Progress update
    if mod(fnum,100) == 0
        disp(['processed ',num2str(fnum),'/',num2str(length(frameRange))])
    end
end

%% Non-dimesionalize
D.u = D.u/U;
D.v = D.u/U;
D.vort = D.vort*L/U;
D.x = D.x/L;
D.y = D.y/L;

if eav == 1
    D.uncU = D.uncU/U;
    D.uncV = D.uncV/U;
end

if out1.dimNum == 3
    D.w = D.w/U;
    if eav == 1
        D.uncW = D.uncW/U;
    end
end

% Return to main folder
cd(folderMain)
disp('Done')

end


%% Subroutines

% Extract data from DaVis data structure ----------------------------------
function out1 = extractData(dataStructure, CorrelationThreshold,n_camField)

% Components
C = dataStructure.Frames{n_camField, 1}.Components;

% Scales
S = dataStructure.Frames{n_camField, 1}.Scales;

% Grids
G = dataStructure.Frames{n_camField, 1}.Grids;

% Find field names by search (not hard-coded)
for i=1:length(C)
    field_names{i} = C{i,1}.Name;
end
% Assign variable names from indices for each field_name
idu = find(strcmp(field_names, 'U0'),1); % First (and only) index where matches occur
idv = find(strcmp(field_names, 'V0'),1);
idw = find(strcmp(field_names, 'W0'),1);
idcorr = find(strcmp(field_names, 'TS:Correlation value'),1);
idParticle_Size = find(strcmp(field_names, 'TS:Particle size'),1);
idPeak_Ratio = find(strcmp(field_names, 'TS:Peak ratio'),1);
iduncU = find(strcmp(field_names, 'TS:Uncertainty Vx'),1);
iduncV = find(strcmp(field_names, 'TS:Uncertainty Vy'),1);
iduncW = find(strcmp(field_names, 'TS:Uncertainty Vz'),1);

% Create 'out1.dimNum'variable
if isempty(idw)
    out1.dimNum = 2;
else
    out1.dimNum = 3;
end

% Check for ISVALID flag
if length(dataStructure.Frames{n_camField, 1}.ComponentNames) > 12
    %out1.dimNum = 3;
    if  length(dataStructure.Frames{n_camField, 1}.ComponentNames) < 16
        % isvalid does not exist
        isvalid_exists = 0;
    else
        isvalid_exists = 1;
    end
else
    %out1.dimNum = 2;
    if  length(dataStructure.Frames{n_camField, 1}.ComponentNames) < 12
        % isvalid does not exist
        isvalid_exists = 0;
    else
        isvalid_exists = 1;
    end
end

% Extract vector data
out1.uRaw = C{idu, 1}.Planes{1, 1} * S.I.Slope + C{idu, 1}.Scale.Offset;
out1.vRaw = C{idv, 1}.Planes{1, 1} * S.I.Slope + C{idv, 1}.Scale.Offset;
if out1.dimNum == 3
    out1.wRaw = C{idw, 1}.Planes{1, 1} * S.I.Slope + C{idw, 1}.Scale.Offset;
end

out1.corrRaw = C{idcorr, 1}.Planes{1, 1};

% Replace invalid vectors with NaN
if isvalid_exists
    isValid = C{end, 1}.Planes{1, 1};
    invalid_vectors = find((out1.corrRaw < CorrelationThreshold) | isnan(out1.corrRaw) | (isValid == 0) );
else
    invalid_vectors = find((out1.corrRaw < CorrelationThreshold) | isnan(out1.corrRaw));
end

% Mask out the invalid pixels
out1.uRaw(invalid_vectors) = NaN;
out1.vRaw(invalid_vectors) = NaN;
out1.corrRaw(invalid_vectors) = NaN;
if out1.dimNum == 3
    out1.wRaw(invalid_vectors) = NaN;
end

% Generate grids
xr = ( (1:size(out1.uRaw,1)) - 0.5 );
yr = ( (1:size(out1.uRaw,2)) - 0.5 );

[out1.xRaw, out1.yRaw] = ndgrid((xr * G.X * S.X.Slope + S.X.Offset)/1000,...
    (yr * G.Y * S.Y.Slope + S.Y.Offset)/1000);

% Reorient to compute vorticity
out1.xRaw =  flip(out1.xRaw');
out1.yRaw =  flip(out1.yRaw');
out1.uRaw =  flip(out1.uRaw');
out1.vRaw = -flip(out1.vRaw');
if out1.dimNum == 3
    out1.wRaw = flip(out1.wRaw');
end

% Compute vorticity
out1.vortRaw = calculateVorticity(out1.xRaw,out1.yRaw,out1.uRaw,out1.vRaw);

% Save data
if isvalid_exists
    out1.isValidRaw = single(flip(isValid'));
end

out1.corrRaw = flip(out1.corrRaw');

% Uncertainties
out1.uncURaw = C{iduncU, 1}.Planes{1, 1};
out1.uncVRaw = C{iduncV, 1}.Planes{1, 1};
if out1.dimNum == 3
    out1.uncWRaw = C{iduncW, 1}.Planes{1, 1};
end

out1.uncURaw = flip(out1.uncURaw');
out1.uncVRaw = flip(out1.uncVRaw');
if out1.dimNum == 3
    out1.uncWRaw = flip(out1.uncWRaw');
end

end

% Compute vorticity -------------------------------------------------------
function omega_z = calculateVorticity(xRaw,yRaw,uRaw,vRaw)

% A: ALGORITHM TAKEN FROM RAFFEL'S PIV HANDBOOK
omega_z = nan(size(xRaw));

dx = abs(xRaw(1,1) - xRaw(1,2));
dy = abs(yRaw(1,1) - yRaw(2,1));

% Dont compute the edge vorticity
for i = 2:size(xRaw,1)-1
    for j = 2:size(xRaw,2)-1
        omega_z(i,j) = ...
            ( ...
            - (1/2)*dx*(vRaw(i-1,j-1) + 2*vRaw(i,j-1) + vRaw(i+1,j-1)) ...
            - (1/2)*dy*(uRaw(i+1,j-1) + 2*uRaw(i+1,j) + uRaw(i+1,j+1)) ...
            + (1/2)*dx*(vRaw(i+1,j+1) + 2*vRaw(i,j+1) + vRaw(i-1,j+1)) ...
            + (1/2)*dy*(uRaw(i-1,j+1) + 2*uRaw(i-1,j) + uRaw(i-1,j-1)) ...
            ) ...
            / (4*dx*dy);
    end
end

% NOTE: ADD COMPUTATION OF EDGES (USE MATLAB'S IMPLEMENTATION IN "curl")

% B: Using Matlab's curl function:
% [omega_z, ~] = curl(xRaw, yRaw, uRaw, vRaw);

end



