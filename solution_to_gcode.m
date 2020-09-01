% Matt Ireland 2020
% Generate GCODE from TO solution data

clear;
% Load workspace variables
workspace = strcat('201904110923','.mat');
load(workspace,'elementDensity','nodalCoordinates','elementConnectivity',...
                'dhy','dhx','NYE','NXE','height','length');            
tic;

%% Notes

% Run this after successful completion of fea.m
% Imports element density solution and generates a GCODE path

% ** Currently only generates a toolpath corresponding to 90deg matl **

% ** Assumes an element size of 0.5 x 0.5 mm **

%% Inputs

mode = 0;
% 0 for MBB

render = 0;
% 0 for no plots
% 1 for plots

threshold = 0.5; % Lower density limit for binary image

scale = 10; % Upscale factor **keep at 10**

fileName = '201904110923_GCODE';

%% Print Settings

numLayers = 20;

firstLayerHeight = 0.2; % mm

layerHeight = 0.100; % mm

retractionDist = 100; % mm

xOrigin = 51.355; % mm
yOrigin = 86.150; % mm

% These values calculated from excel & slicer
layerOne_extruderMultiplier = 0.2; % extruder multiplier
layer_extruderMultiplier = 0.2;

layerOne_travelRate = 4320; % feedrate, non-extrusion
layer_travelRate = 4320;

layerOne_extrusionRate = 1800; % feedrate, extrusion

layer_retractionRate = 3600; % feedrate, retraction

feedrate_change = layer_travelRate/layerOne_travelRate; 

%% MBB 90deg Toolpath Generation
if mode == 0
%% Make a New Image
disp('New Image');

% Make the mirror image density map
elementDensityMirror = flip(elementDensity,2);

% Place the matrices side-by-side
fullMap = [elementDensityMirror elementDensity];

% Flip the matrix
fullMap = flip(fullMap);

% Convert to image
fullMapImg = mat2gray(fullMap, [max(fullMap(:)) min(fullMap(:))]);

% Upscale
fullMapImg = imresize(fullMapImg,scale,'nearest');

% Generate smoothed image
[newBounds,newImg] = sg_smooth(fullMapImg,threshold,render);

%% Generate raster base-pattern
disp('Raster Base');

% Initialize vector for pattern
rasVec = zeros(1,size(newImg,2));

% Find offset spacing value
offset = scale / 2;

% Write raster spacing into vector
rasVec(offset:offset:end) = 1;

% Write raster into matrix for pathing
rasMat = repmat(rasVec,size(newImg,1),1);

% Trim top and bottom of raster matrix
for ii = 1:offset
   rasMat(ii,:) = 0;
   rasMat( (size(rasMat,1) - (ii-1) ),:) = 0;
end

%% Trim Raster to Contours
disp('Trim Raster');

% Build logical materix corresponding to areas with holes
filledBounds = imfill(newBounds,'holes');

% Subtract holes matrix from raster matrix
trimmed_rasMat = rasMat - filledBounds;

% Reset negative values to zeroes
trimmed_rasMat(trimmed_rasMat < 0) = 0;

%% Sequence Raster Path
disp('Raster Sequence');

coords = []; % m-by-5 coordinates that will be turned into the toolpath
sign = 1; % switching var for raster direction
kk = 1; bb = 1;% coordinate index val

% Find the endpoints
endpoints = bwmorph(trimmed_rasMat,'endpoints');

for jj = offset:offset:size(trimmed_rasMat,2) % by column
    % odd columns (not by index, by raster count)
    if mod(sign,2) ~= 0 
        for ii = offset:1:size(trimmed_rasMat)-offset          
            if endpoints(ii,jj) == 1
                
                if bb == 1
                % travel move
                coords(kk,:) = [jj ii 0 nan 0];
                end
                
                if bb == 2
                % prime nozzle
                coords(kk,:) = [jj ii 0 0 0];
                kk = kk + 1;
                % extract to point
                coords(kk,:) = [jj ii 0 1 0];
                kk = kk + 1;
                % retract before moving
                coords(kk,:) = [jj ii 0 -1 0];
                bb = 0;
                end
                kk = kk + 1; bb = bb + 1;
            end                       
        end
    end
    % even columns (not by index, by raster count)
    if mod(sign,2) == 0
        for ii = size(trimmed_rasMat)-offset:-1:offset      
            if endpoints(ii,jj) == 1
                
                if bb == 1
                % travel move
                coords(kk,:) = [jj ii 0 nan 0];
                end
                
                if bb == 2
                % prime nozzle
                coords(kk,:) = [jj ii 0 0 0];
                kk = kk + 1;
                % extract to point
                coords(kk,:) = [jj ii 0 1 0];
                kk = kk + 1;
                % retract before moving
                coords(kk,:) = [jj ii 0 -1 0];
                bb = 0;
                end
                kk = kk + 1; bb = bb + 1;
            end        
        end
    end
    % Alternate path up/down
    sign = sign + 1; 
    
end

%% Build Extruder Lengths
disp('Extrusion');

% Initialize a travel distance array
extrusionDist = zeros(size(coords,1),1);

for ii = 1:size(coords,1)
    % Travel moves
    if isnan(coords(ii,4))
        extrusionDist(ii) = 0;
    end
    % Retraction moves
    if coords(ii,4) == -1
        extrusionDist(ii) = -retractionDist * scale;
    end
    % Priming moves
    if coords(ii,4) == 0
        extrusionDist(ii) = retractionDist * scale;
    end
    % Extrusion moves, distance dependent
    if coords(ii,4) == 1
        dx = coords(ii,1) - coords(ii-2,1);
        dy = coords(ii,2) - coords(ii-2,2);
        ds = (dx^2 + dy^2)^(1/2);
        extrusionDist(ii) = ds;
    end
end

% Fix scale of travel distances
extrusionDist = extrusionDist ./ scale;

% Multiply travel dist by extruder multiplier
extrusionDist = extrusionDist .* layer_extruderMultiplier;

% Add extrusion values to coordinate matrix
coords(:,4) = extrusionDist;

%% Build Feedrate Vector
disp('Feedrate');

% Initialize
feedRate = extrusionDist;

% Find retraction/prime movements
for ii = 1:size(feedRate,1)
    if mod(ii,2) == 0
        feedRate(ii) = nan;        
    end
end

% Name first layer extrude movement rates
feedRate(feedRate > 0) = layerOne_extrusionRate;

% Name first layer travel movement rates
feedRate(feedRate == 0) = layerOne_travelRate;

% Name first layer travel movement rates
feedRate(isnan(feedRate)) = layer_retractionRate;

%% Add Layer Change Toolhead Path
disp('End of Path');

% Y Movement
coords(size(coords,1)+1,:) = [coords(end,1) max(coords(:,2))+10 0 0 0];
extrusionDist(size(extrusionDist,1)+1,:) = 0;
feedRate(size(feedRate,1)+1,:) = layerOne_travelRate;

% X Movement
coords(size(coords,1)+1,:) = [min(coords(:,1))-10 coords(end,2) 0 0 0];
extrusionDist(size(extrusionDist,1)+1,:) = 0;
feedRate(size(feedRate,1)+1,:) = layerOne_travelRate;

%% Add Comment Line
disp('Layer Change Comment');

coords(size(coords,1)+1,:) = nan;
extrusionDist(size(extrusionDist,1)+1,:) = nan;
feedRate(size(feedRate,1)+1,:) = nan;

%% Add Multiple Layers, Write in Z & E Coordinates, Add Feedrate Values
disp('Build Layers');
% Write in first layer Z values
coords(1:(end-1),3) = firstLayerHeight;

% Write in first layer feed rates
coords(:,5) = feedRate;

% Initialize full coordinate matrix
fullCoords = zeros(size(coords,1)*numLayers,5);

% Repeat path over specified number of layers
ii = 0; aa = 1;
for jj = 1:numLayers
   for kk = 1:size(coords,1)
      
      if ii < 1 % first layer parameters
          coords(kk,4) = extrusionDist(kk,:) .* layerOne_extruderMultiplier; % extruder multiplier 
      end
      
      if ii > 0 % after first layer
          coords(:,3) = firstLayerHeight + layerHeight * aa; % new layer heights
          coords(kk,5) = coords(kk,5) .* feedrate_change; % new feedrates
      end
      
      % write out layer path to full materix
      fullCoords(ii+kk,:) = coords(kk,:);
   end
   aa = aa + 1;
   ii = ii + size(coords,1);
end

%% Rescale and Reposition XY Coordinates
disp('Rescale & Reposition');

fullCoords(:,1:2) = fullCoords(:,1:2) ./ scale;

fullCoords(:,1) = fullCoords(:,1) + xOrigin;
fullCoords(:,2) = fullCoords(:,2) + yOrigin;

%% Write Output File
disp('Write Out');
% Initialize Cell Array
textOut = cell(size(fullCoords,1),1);

% Generate cell data
for ii = 1:size(fullCoords,1)
    textOut{ii} = writeGCode(fullCoords(ii,:));   
end

% Convert to table
tableOut = cell2table(textOut);

% Write out table
writetable(tableOut,fileName,'FileType','text','QuoteStrings',false,'WriteVariableNames',false);

%% Plot input and output comparison
if render == 1
figure(2)
fig = gcf;
% set(fig, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
fontSize = 18;
subplot(2,1,1)
map = colormap('gray');
map = flip(map);
% Element density
imagesc(fullMap);
title('Original Element Density Map', 'FontSize', fontSize, 'Interpreter', 'None');
colormap(map);
pbaspect([size(newImg,2) size(newImg,1) 1]);

subplot(2,1,2)
map = colormap('gray');
map = flip(map);
% New Element density
imagesc(newImg);
title('Upscaled and Smoothed Density Map', 'FontSize', fontSize, 'Interpreter', 'None');
colormap(map);
pbaspect([size(newImg,2) size(newImg,1) 1]);
end %end plot render

% Write out execution time
disp(['Time to complete: ' num2str(toc)])

end %end mbb mode if

%% Savitky-Golay Edge Smoothing
function [newBounds,newImg] = sg_smooth(fullMapImg,thresholdValue,render)
clc;    % Clear the command window
close all;  % Close all figures (except those of imtool)
format long g;
format compact;
fontSize = 18;

% Resize image and apply gaussian blur
grayImage = fullMapImg;
grayImage = imgaussfilt(grayImage,10);

% Get the dimensions of the image
% numberOfColorBands should be = 1
[~, ~, numberOfColorBands] = size(grayImage);

% Convert to gray scale by taking only the green channel if not
if numberOfColorBands > 1
    grayImage = grayImage(:, :, 2);
end

% Display the original gray scale image
if render == 1
    subplot(4, 1, 1);
    imshow(grayImage, []);
    axis on;
    title('Blurred Grayscale Image', 'FontSize', fontSize, 'Interpreter', 'None');
    
    % Enlarge figure to full screen
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
end



% Threshold the image
binaryImage = grayImage > thresholdValue;

% Get rid of holes in the blobs
binaryImage = imfill(binaryImage, 'holes');

%---------------------------------------------------------------------------
% Choose filter parameters (need for min blob size)
windowWidth = 11;
polynomialOrder = 3;
% Extract the blobs
blobs = ExtractBlobs(binaryImage,windowWidth);
%---------------------------------------------------------------------------

% Display the binary image
if render == 1
    subplot(4, 1, 2);
    imshow(blobs, []);
    title('Binary Image', 'FontSize', fontSize, 'Interpreter', 'None');
    
    %% Get the boundaries
    % Display the original gray scale image
    subplot(4, 1, 3);
    imshow(grayImage, []);
    axis image;
    hold on;
end
% Find boundaries from binary image
boundaries = bwboundaries(blobs);
numberOfBoundaries = size(boundaries, 1);

% Plot
if render == 1
    for k = 1 : numberOfBoundaries
        thisBoundary = boundaries{k};
        plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 2);
    end
    hold off;
    caption = sprintf('%d Boundaries', numberOfBoundaries);
    title(caption, 'FontSize', fontSize);
    axis on;
    
    % Display the original gray scale image
    subplot(4, 1, 4);
    imshow(grayImage, []);
    axis image;
    hold on;
end

% Initialize an image drawn from filtered boundaries
newImg = zeros(size(grayImage,1),size(grayImage,2));

% Iterate through each boundary set
for ii = 1:numberOfBoundaries
    % Get the x and y coordinates.
    x{ii}(:,1) = boundaries{ii}(:, 2);
    y{ii}(:,1) = boundaries{ii}(:, 1);
    
    % Now smooth with a Savitzky-Golay sliding polynomial filter
    smoothX{ii}(:,1) = sgolayfilt(x{ii}(:,1), polynomialOrder, windowWidth);
    smoothY{ii}(:,1) = sgolayfilt(y{ii}(:,1), polynomialOrder, windowWidth);
    
    for jj = 1:length(smoothX{ii})
        % Write the smooth contour indices into a new image
        sx{ii}(jj,1) = round(smoothX{ii}(jj,1));
        sy{ii}(jj,1) = round(smoothY{ii}(jj,1));
        newImg(sy{ii}(jj,1),sx{ii}(jj,1)) = 1;
    end
    
    % Show the smooth boundaries
    if render == 1
        plot(smoothX{ii}(:), smoothY{ii}(:), 'r-', 'LineWidth', 2);
        hold on;
        grid on;
        title('Smoothed Boundaries', 'FontSize', fontSize);
    end
end

newBounds = newImg;

% Fill the holes in the new image
newImg = imfill(newImg, 'holes');

% Invert new image
newImg = (newImg - 1)./ -1;

end

%% Blob Extraction
% Extract blob properties and ignore ones with dimensions smaller than windowWidth
function binaryImage = ExtractBlobs(binaryImage,windowWidth)
try
    % Get all the blob properties
    [labeledImage, numberOfBlobs] = bwlabel(binaryImage);
    blobMeasurements = regionprops(labeledImage, 'Area','BoundingBox');
    
    % Find blobs smaller than filter window
    bm = reshape([blobMeasurements.BoundingBox],4,[]);
    smallBlobs = round(nnz(bm(bm(3:4,:) < windowWidth))/2);
    
    % Ignore blobs smaller than filter window
    numberToExtract = numberOfBlobs - smallBlobs;
    
    % Get all the areas
    allAreas = [blobMeasurements.Area];
    
    if numberToExtract > 0
        % sort in order of largest to smallest
        
        [~, sortIndexes] = sort(allAreas, 'descend');
        
    end
    
    % Extract the "numberToExtract" largest blob(a)s using ismember()
    biggestBlob = ismember(labeledImage, sortIndexes(1:numberToExtract));
    
    % Convert from integer labeled image into binary (logical) image
    binaryImage = biggestBlob > 0;
    
catch ME
    errorMessage = sprintf('Error in function ExtractBlobs().\n\nError Message:\n%s', ME.message);
    fprintf(1, '%s\n', errorMessage);
    uiwait(warndlg(errorMessage));
end
end

%% Write GCODE
function [code] = writeGCode(line)

%string for the code
code = [];

 if ~isnan(line(1))

    %convert values to strings
    x = strcat(' X', num2str(line(1), '%2.3f'));
    y = strcat(' Y', num2str(line(2), '%2.3f'));
    z = strcat(' Z', num2str(line(3), '%2.3f'));
    e = strcat(' E', num2str(line(4), '%2.4f'));
    f = strcat(' F', num2str(line(5), '%2.0f'));
    
    % G1 for extrusion movements
    if line(4) ~= 0
        code = sprintf('%sG1%s%s%s%s%s',code,x,y,z,e,f);
    end
    
    
    % G0 for extrusion movements
    if line(4) == 0
        code = sprintf('%sG0%s%s%s%s%s',code,x,y,z,f);
    end
    
 else
     code = '; Layer Change' ;
 end

end