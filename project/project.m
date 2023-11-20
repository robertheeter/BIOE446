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

addpath(genpath('data'));
filename = 'watershed_example.png';

img = imread(filename);

figure(1)
imshow(img)

%% Applying the distance transform for watershedding

D = -bwdist(~img);

figure(2)
imshow(D,[]);

%% Getting boundaries with watershed transform

Ld = watershed(D);

img_2 = img;

img_2(Ld == 0) = 0;

figure(3)
imshow(img_2)

%% Filter our small local minima using imextendedmin

mask = imextendedmin(D,2);
D2 = imimposemin(D, mask);
Ld2 = watershed(D2);

img_3 = img;
img_3(Ld2 == 0) = 0;

figure(4)
imshow(img_3)

%% Counting number and area of objects in image

stats = regionprops('Table', img_3, 'Area');

%% Binarizing with multithresh and imquantize v. imbinarize

filename = 'coins.png';

img = imread(filename);
figure(5)
imshow(img)

%% Binarizing with imbinarize

img_bn = imbinarize(img, 'Adaptive');

figure(6)
imshow(img_bn, [])

%% Binarizing with multithresh and imquantize

level = multithresh(img, 2);
img_qt = imquantize(img,level);

figure(7)
imshow(img_qt,[])






%% Problem 2: 2D pattern formation in the Gray-Scott model
disp('PROBLEM 2')
% clear

%% Problem 2, Part A


%% Functions

% None
