% image = imread('image1.jfif');
% image = rgb2gray(image);
% imshow(canny(img, 1.4, 2.4, 2.8));
function [result_image] = canny(image, sigma, low_threshold, high_threshold)
    [image_height, image_width] = size(image);
    
    N = 2;
    x = -N : N;
    y = -transpose(x);

    % Two 1D Gaussian Filters with standard deviation sigma
    Gx = exp(-((x.^2) / (2 * sigma^2))) / (sqrt(2 * pi) * sigma);
    Gy = exp(-((y.^2) / (2 * sigma^2))) / (sqrt(2 * pi) * sigma);

    % Convoluting image with the two Gaussian filters
    % Apply 1D Gaussian in X direction
    blurred_image = conv2(image, Gx, 'same');
    % Apply 1D Gaussian in Y direction
    blurred_image = conv2(blurred_image, Gy, 'same');

    % Step 1: Compute x and y derivatives of image
    % Could also use gradient(Gx) and gradient(Gy)
    Hx = -(2^(1/2) .* x .* exp(-x.^2 / (2 * sigma^2))) / (2 * pi^(1/2) * sigma^3);
    Hy = -(2^(1/2) .* y .* exp(-y.^2 / (2 * sigma^2))) / (2 * pi^(1/2) * sigma^3);

    image_x = conv2(blurred_image, Hx, 'same');
    image_y = conv2(blurred_image, Hy, 'same'); 

    % Step 2: Compute magnitude and direction of every pixel
    magnitude = sqrt(image_x.^2 + image_y.^2);
    direction = atan2(image_y, image_x);

    % Step 3: Eliminate pixels that are not local maxima of the magnitude in the direction of the gradient
    supressed_image = thin_edges(image_height, image_width, magnitude, direction);

    % Step 4: Hysteresis Thresholding
    result_image = hystereis(supressed_image, image_height, image_width, low_threshold, high_threshold);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Helper Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function result_image = thin_edges(image_height, image_width, magnitude, direction)
        
        result_image = zeros(image_height, image_width);
        
        for i = 2 : image_height - 1
            for j = 2 : image_width - 1
                
                pixel_direction = abs(direction(i, j));
                pixel_magnitude = magnitude(i, j);
                
                if ((pixel_direction >= 0 && pixel_direction < pi / 8) || ...
                    (pixel_direction >= 15 * pi / 8 && pixel_direction <= 2 * pi) || ...
                    (pixel_direction >= 7 * pi / 8 && pixel_direction < 9 * pi / 8))
                    before = magnitude(i, j - 1);
                    after = magnitude(i, j + 1);

                elseif ((pixel_direction >= pi / 8 && pixel_direction < 3 * pi / 8) || ...
                        (pixel_direction >= 9 * pi / 8 && pixel_direction < 11 * pi / 8))
                    before = magnitude(i + 1, j - 1);
                    after = magnitude(i - 1, j + 1);

                elseif ((pixel_direction >= 3 * pi / 8 && pixel_direction < 5 * pi / 8) || ...
                        (pixel_direction >= 11 * pi / 8 && pixel_direction < 13 * pi / 8))
                    before = magnitude(i - 1, j);
                    after = magnitude(i + 1, j);

                else
                    before = magnitude(i - 1, j - 1);
                    after = magnitude(i + 1, j + 1);
                end

                if pixel_magnitude > before && pixel_magnitude > after
                    result_image(i, j) = pixel_magnitude;
                end
            end
        end
    end

    function result_image = hystereis(image, image_height, image_width, low_threshold, high_threshold)
        
        result_image = zeros(image_height, image_width);
        
        for i = 1 : image_height - 1
            for j = 1 : image_width - 1
                if (image(i, j) < low_threshold)
                    result_image(i, j) = 0;
                
                elseif (image(i, j) > high_threshold)
                    result_image(i, j) = 255;
                
                elseif (image(i, j) > low_threshold)   
                    if ((image(i + 1, j) > high_threshold) || ...
                        (image(i - 1, j) > high_threshold) || ...
                        (image(i, j + 1) > high_threshold) || ...
                        (image(i, j - 1) > high_threshold) || ...
                        (image(i + 1, j + 1) > high_threshold) || ...
                        (image(i - 1, j - 1) > high_threshold) || ...
                        (image(i + 1, j - 1) > high_threshold) || ...
                        (image(i - 1, j + 1) > high_threshold))
                        result_image(i, j) = 255;
                    end
                end
            end
        end
    end
end