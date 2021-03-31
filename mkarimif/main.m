image = imread('image1.jfif');

% Question 2 - Edge Detection
figure(1);
imshow(canny(image, 1.6, 100, 130));
title('edges');

% Question 3 - Corner Detection
figure(2);
% show the original image
imshow(image);
hold on;
[rows, cols] = harris(image, 5, 0.04, 5000000);
% show the corners in red dots
plot(cols, rows, 'r.');
title('corners');