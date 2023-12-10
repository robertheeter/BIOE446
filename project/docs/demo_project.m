%% Loading Image

%addpath(genpath('data/CMV_smFISH/DAPI'));
%filename = 'DAPI_A.tif';

filename = 'watershed_example.png';
img = imread(filename);

figure(1)
imshow(img);


%%  Applying distrance transform for watershed
img = filled_dapi;
D = -bwdist(~img);

figure(2)
imshow(D,[])

%% Getting boundaries with watershed tranform

Ld = watershed(D);

figure(3)
imshow(label2rgb(Ld));

%% Applying initial watershed
img_2 = img;
img_2(Ld == 0) = 0;

figure(4)
imshow(img_2);

%% Since initial watershed failed, filtering out tiny local minimum

mask = imextendedmin(D,2);

figure(5)
imshowpair(img,mask,'blend')

%% Using new mask to do more accurate watershed.

D2 = imimposemin(D,mask);
Ld2 = watershed(D2);
img_3 = img;

img_3(Ld2 == 0) = 0;

figure(6)
imshow(img_3);


%% Getting number and area of objects

stats = regionprops('Table',img_3,'Area');

areas = stats.Area;

n_objects = length(areas);
fprintf('There are %i detected objects.\n',n_objects)

%% Plotting distribution of areas

figure(7)
histogram(areas,50)


%% Example of multithresh v. binarize

filename = 'coins.png';

img = imread(filename);

figure(8)
imshow(img, []);

%% Using imbinarize to binarize image.

img_bn = imbinarize(img,'Adaptive');

figure(9)
imshow(img_bn,[])

%% using multithresh and imquantize

level = multithresh(img,3);
img_qt = imquantize(img,level);

figure(10)
imshow(img_qt,[]);

disp(level)