% IMAGES: C

%% GLOBAL VARIABLES
global area_size_thresh
global single_mrna_thresh
global nascent_intensity_thresh

% threshold area size for nuclei detection
area_size_thresh = 3000;

% threshold intensity for single mRNA molecule
single_mrna_thresh  = 0.01;

% threshold intensity for nascent mRNA at transcription site
nascent_intensity_thresh = 7;


%% STEP 1
% Load in DAPI and mRNA images
% Segment DAPI image into objects
% Calculate area of each object to isolate nuclei

fprintf('IMAGES: C\n')
fprintf('RUNNING STEP 1:\n\n')

addpath(genpath('data/CMV_smFISH/'));
dapi_file = 'DAPI/DAPI_C.tif';
mrna_file = 'mRNA/mRNA_C.tif';

dapi_img = im2double(imread(dapi_file));
mrna_img = im2double(imread(mrna_file));

figure
imshow(dapi_img)
title('Original DAPI Image')

% figure
% imshow(mrna_img)
% title('Original mRNA Image')

[areas, labels, num_objects, num_nuclei, dapi_binary_img] = isolate_nuclei(dapi_img);
fprintf('There are %i detected nuclei.\n', num_nuclei)
fprintf('###############################\n\n')


%% STEP 2
% Loop through each object and check if its a nucleus (>10000)
% If nucleus, then make a mask and isolate mRNA
% Calculate total intensity of each spot

fprintf('RUNNING STEP 2:\n\n')

all_spots = [];
for i = 1:num_objects
    if areas(i) > area_size_thresh
        nucleus_mask = (labels == i);
        [num_spots, intensities] = isolate_mrna(dapi_binary_img, mrna_img, nucleus_mask);

        all_spots = [all_spots; intensities];
    end
end

fprintf('###############################\n\n')

%% STEP 3
% Plot distribution of intensities for all cells

fprintf('RUNNING STEP 3:\n\n')

figure
histogram(all_spots(all_spots <= 0.1), 50)
title('Image C')
xlabel('Intensity')
ylabel('Frequency')
% from histogram, appears that 0.02 is intensity for single mRNA
% this is equal to single_mrna_thresh

fprintf('###############################\n\n')


%% STEPS 4 and 5
% Calculate number of mRNA per cell
% Calculate total number of mRNA in image
% Check if any intensity spots for a single cell are very high (>100), and
% if so calculate the mRNA for that spot (at transcription site)

fprintf('RUNNING STEPS 4 & 5:\n\n')

% vector to store number of mRNA per cell
num_mrna = [];

% vector to store nascent mRNA (mRNA at transcription site) per cell
nascent_mrna = [];

fprintf('Checking cells for transcription sites...\n\n')
for i = 1:num_objects
    if areas(i) > area_size_thresh
        nucleus_mask = (labels == i);
        [num_spots, intensities] = isolate_mrna(dapi_binary_img, mrna_img, nucleus_mask);
        num_mrna = [num_mrna; sum(intensities)/single_mrna_thresh];
        
        if (1 <= sum(intensities > nascent_intensity_thresh)) && (sum(intensities > nascent_intensity_thresh) <= 2)
            fprintf('Cell has %i transcription site(s)\n', sum(intensities > nascent_intensity_thresh));
            nascent = sum(intensities(intensities > nascent_intensity_thresh))/single_mrna_thresh;
            nascent_mrna = [nascent_mrna; nascent];
        else
            nascent_mrna = [nascent_mrna; 0];
            fprintf('Cell does not have nascent mRNA\n')
        end
        fprintf('--------------------\n')
    end
end

% disp('Number of free mRNA per cell:')
% num_mrna
total_mrna = sum(num_mrna);
fprintf('There are %i total mRNA molecules in the image.\n\n', total_mrna)

% disp('Number of nascent mRNA per cell:')
% nascent_mrna
total_nascent_mrna = sum(nascent_mrna);
fprintf('There are %i nascent mRNA molecules in the image.\n\n', total_nascent_mrna)

total_free_mrna = total_mrna - total_nascent_mrna;
fprintf('There are %i free mRNA molecules in the image.\n\n', total_free_mrna)

fprintf('Mean nascent mRNA per cell: %i\n', mean(nascent_mrna))
fprintf('Mean free mRNA per cell: %i\n\n', mean(num_mrna - nascent_mrna))

fprintf('Variance nascent mRNA per cell: %i\n', var(nascent_mrna))
fprintf('Variance free mRNA per cell: %i\n\n', var(num_mrna - nascent_mrna))

fprintf('CV nascent mRNA per cell: %i\n', std(nascent_mrna)/mean(nascent_mrna))
fprintf('CV free mRNA per cell: %i\n\n', std(num_mrna - nascent_mrna)/mean(num_mrna - nascent_mrna))

cc = corr(nascent_mrna, (num_mrna - nascent_mrna));
fprintf('Correlation coefficient: %i\n\n', cc)

figure
scatter(nascent_mrna, (num_mrna - nascent_mrna), 'filled');
xlabel('Nascent mRNA')
ylabel('Free mRNA')
title('Image C: Nascent mRNA vs. Free mRNA Per Cell')

fprintf('###############################\n\n')



%% FUNCTIONS

function [areas, labels, num_objects, num_nuclei, dapi_img_2] = isolate_nuclei(dapi_img)
    global area_size_thresh

    % Reduce background noise using median filtering
    filtered_dapi = medfilt2(dapi_img, [3, 3]);
    
    % Adaptive thresholding to binarize the image
    threshold = adaptthresh(filtered_dapi, 0.5);
    binary_dapi = imbinarize(filtered_dapi, threshold);
    
    % Fill holes in the binary image
    filled_dapi = imfill(binary_dapi, 'holes');
    
    % Distance transform for watershed
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
    num_nuclei = sum(areas > area_size_thresh);
    
    % Get object labels and number of objects
    [labels, num_objects] = bwlabel(dapi_img_2);
end


function [num_spots, intensities] = isolate_mrna(dapi_binary_img, mrna_img, nucleus_mask)
    % Apply thresholding to isolate mRNA signal
    level = multithresh(mrna_img);
    mrna_binary = imbinarize(mrna_img, level);
    
    % Isolate nucleus from mRNA image
    mrna_isolated = mrna_img.*nucleus_mask;
    mrna_binary_isolated = mrna_binary.*nucleus_mask;

    % Isolate nucleus from DAPI image
    dapi_binary_isolated = dapi_binary_img.*nucleus_mask;
    
    % figure
    % imshow(dapi_binary_isolated)
    % title('Isolated nucleus - DAPI')

    % figure
    % imshow(mrna_binary)
    % title('Isolated nucleus - mRNA')
    
    mrna_labeled = bwlabel(mrna_binary_isolated);

    pixel_intensity = table2array(regionprops('Table', mrna_labeled, mrna_isolated, 'MeanIntensity'));
    areas = table2array(regionprops('Table', mrna_labeled, 'Area'));
    intensities = areas.*pixel_intensity;
    num_spots = length(intensities);
end




