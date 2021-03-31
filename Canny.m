function [result_image] = canny(image, sigma, low_threshold, high_threshold)
    
    % colour to grayscale
    image = double(rgb2gray(image));

    [image_height, image_width] = size(image);
    
    n = 2;
    x = -n : n;
    y = transpose(x);

    % Two 1D Gaussian Filters with standard deviation sigma
    Gx = exp(-((x.^2) / (2 * sigma^2)));
    Gy = exp(-((y.^2) / (2 * sigma^2)));

    % Convoluting image with the two Gaussian filters
    % Apply 1D Gaussian in X direction
    filtered_image = convolve(image, Gx);
    % Apply 1D Gaussian in Y direction
    filtered_image = convolve(filtered_image, Gy);

    % Step 1: Compute x and y derivatives of image
    % Could also use gradient(Gx) and gradient(Gy)
    Hx = -x / sigma^2 .* exp(-x.^2 / (2 * sigma^2));
    Hy = -y / sigma^2 .* exp(-y.^2 / (2 * sigma^2));

    image_x = convolve(filtered_image, Hx);
    image_y = convolve(filtered_image, Hy); 

    % Step 2: Compute magnitude and direction of every pixel
    magnitude = sqrt(image_x.^2 + image_y.^2);
    direction = atan2(image_y, image_x);

    % Step 3: Eliminate pixels that are not local maxima of the magnitude in the direction of the gradient
    supressed_image = thin_edges(image_height, image_width, magnitude, direction);

    % Step 4: Hysteresis Thresholding
    result_image = hystereis(supressed_image, image_height, image_width, low_threshold, high_threshold);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Helper Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [result] = convolve(image, kernel)
        % initialize the output result
        result = zeros(size(image));

        % find the width and heights
        [kernel_height, kernel_width] = size(kernel);
        
        % pad the image with zeros
        padded_image = padarray(image, [floor(kernel_height / 2) floor(kernel_width / 2)]);
        kernel = reshape(transpose(flipdim(flipdim(kernel, 1), 2)), [], 1);

        for i = 1 : image_height
            for j = 1 : image_width
                vectpadded_image = reshape(transpose(padded_image(i : i + kernel_height - 1, j : j + kernel_width - 1)), 1, []);
                result(i, j) = vectpadded_image * kernel;
            end
        end
    end  
  

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

                if pixel_magnitude >= before && pixel_magnitude >= after
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
                    result_image(i, j) = 1;
                
                elseif (image(i, j) > low_threshold)   
                    if ((image(i + 1, j) > high_threshold) || ...
                        (image(i - 1, j) > high_threshold) || ...
                        (image(i, j + 1) > high_threshold) || ...
                        (image(i, j - 1) > high_threshold) || ...
                        (image(i + 1, j + 1) > high_threshold) || ...
                        (image(i - 1, j - 1) > high_threshold) || ...
                        (image(i + 1, j - 1) > high_threshold) || ...
                        (image(i - 1, j + 1) > high_threshold))
                        result_image(i, j) = 1;
                    end
                end
            end
        end
    end
end