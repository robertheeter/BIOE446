close all
clc

% DUE DECEMBER 10TH
% need image processing toolbox add-on

%% Problem 1: 1D pattern formation in the Gierer-Meinhardt model (DEMO)
disp('PROBLEM 1')
clear

%% Problem 1, Part A
% use "im2double"!!!
% imadjust to adjust contrast/intensity
% see lecture for other useful commands

% abc are easy images
% def images are hard

% segment mRNA, segment cells
% dots are single mRNAs; mRNA production happens in bursts
% count how many mRNA in cell

%addpath(genpath('data/'));
%filename = 'watershed_example.png';


%% STEP 1

addpath(genpath('data/CMV_smFISH/'));
dapi_file = 'DAPI/DAPI_B.tif';
mrna_file = 'mRNA/mRNA_B.tif';

dapi_img = im2double(imread(dapi_file));
mrna_img = im2double(imread(mrna_file));

figure
imshow(dapi_img)
title('Original DAPI Image')

figure
imshow(mrna_img)
title('Original mRNA Image')

[areas, labels, num_objects] = isolate_nuclei(dapi_img);


%% STEP 2 (currently done with part D) 

% Loop through each object and check if its a nucleus (>10000)
% If nucleus, then make a mask and isolate mRNA

total = [];

for i = 1:num_objects
    if areas(i) > 10000
        nucleus_mask = (labels == i);
    
        [num_spots, intensities] = isolate_mrna(mrna_img, nucleus_mask);
        num_spots
        total = [total; intensities];
    end
end

histogram(total, 50)



%% Functions

function [areas, labels, num_objects] = isolate_nuclei(dapi_img)
    % Reduce background noise using median filtering
    filtered_dapi = medfilt2(dapi_img, [3, 3]);
    
    % Adaptive thresholding to binarize the image
    threshold = adaptthresh(filtered_dapi, 0.5);
    binary_dapi = imbinarize(filtered_dapi, threshold);
    
    % Fill holes in the binary image
    filled_dapi = imfill(binary_dapi, 'holes');
    
    % Distance transform to find regional maxima for watershed
    D = -bwdist(~filled_dapi);
    
    % Filtering out tiny local minimum
    mask = imextendedmin(D, 2);
    D2 = imimposemin(D, mask);
    Ld2 = watershed(D2);
    
    dapi_img_2 = filled_dapi;
    dapi_img_2(Ld2 == 0) = 0;

    figure
    imshow(dapi_img_2)
    title('Segmented DAPI Image')

    % Get areas of objects
    areas = table2array(regionprops('Table', dapi_img_2, 'Area'));
    fprintf('There are %i detected nuclei.\n', sum(areas > 10000))
    
    % Get object labels and number of objects
    [labels, num_objects] = bwlabel(dapi_img_2);
end


function [num_spots, intensities] = isolate_mrna(mrna_img, nucleus_mask)
    % Apply thresholding to isolate mRNA signal
    level = multithresh(mrna_img);
    mrna_binary = imbinarize(mrna_img, level);
    
    % Isolate nucleus from mRNA image
    mrna_isolated = mrna_img.*nucleus_mask;
    mrna_binary_isolated = mrna_binary.*nucleus_mask;

    % figure
    % imshow(mrna_binary_isolated)
    % title('Segmented mRNA Image')
    
    mrna_labeled = bwlabel(mrna_binary_isolated);

    pixel_intensity = table2array(regionprops('Table', mrna_labeled, mrna_isolated, 'MeanIntensity'));
    areas = table2array(regionprops('Table', mrna_labeled, 'Area'));
    intensities = areas.*pixel_intensity;
    num_spots = length(intensities);
end




