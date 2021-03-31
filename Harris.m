function [rows, cols] = harris(image, sigma, k, threshold)
    % Step 1: colour to grayscale
    image = double(rgb2gray(image));
    
    n = 2;
    x = -n : n;
    y = transpose(x);

    % Gaussian filter
    Gxy = exp(-(x.^2 + y.^2) / (2 * sigma^2));
    
    % X and Y derivatives of Gaussian
    Hx = -x / sigma^2 .* exp(-(x.^2 + y.^2) / (2 * sigma^2));
    Hy = -y / sigma^2 .* exp(-(x.^2 + y.^2) / (2 * sigma^2));

    % Step 1: Compute x and y derivatives of image
    Ix = convolve(image, Hx);
    Iy = convolve(image, Hy); 

    % Step 2: Compute products of derivatives at every pixel
    Ix2 = Ix.^2;
    Iy2 = Iy.^2;
    Ixy = Ix.*Iy;

    % Step 3: Compute the sums of the products of derivatives at each pixel
    Sx2 = convolve(Ix2, Gxy);
    Sy2 = convolve(Iy2, Gxy);
    Sxy = convolve(Ixy, Gxy);

    % Step 4: Define matrix H at each pixel
    [rows, cols] = size(image);
    R = zeros(rows, cols);
    for i = 2 : 1 : rows-1
        for j = 2 : 1 : cols-1
            % Sx2 mean
            Sx2_mean = mean(mean(Sx2(i-1:i+1, j-1:j+1)));
            
            % Sy2 mean
            Sy2_mean = mean(mean(Sy2(i-1:i+1, j-1:j+1)));
            
            % Sxy mean
            Sxy_mean = mean(mean(Sxy(i-1:i+1, j-1:j+1)));
            
            H = [Sx2_mean, Sxy_mean; Sxy_mean, Sy2_mean];

            % Step 5: Compute the response of the detector at each pixel
            R(i, j) = det(H) - (k * trace(H)^2);
        end
    end

    %Step 6: Non-maximum suppression
    window = 3;
    result = (R == ordfilt2(R, window.^2, ones(window, window))) & (R > threshold);

    % find the indices of the corners
    [rows, cols] = find(result);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Helper Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [result] = convolve(image, kernel)
        % initialize the output result
        result = zeros(size(image));

        % find the width and heights
        [image_height, image_width] = size(image);
        [kernel_height, kernel_width] = size(kernel);
        
        % pad the image with zeros
        padded_image = padarray(image, [floor(kernel_height / 2) floor(kernel_width / 2)]);
        kernel = reshape(transpose(flipdim(flipdim(kernel, 1), 2)), [], 1);

        for i = 1 : image_height
            for j = 1 : image_width
                result(i, j) = reshape(transpose(padded_image(i : i + kernel_height - 1, j : j + kernel_width - 1)), 1, []) * kernel;
            end
        end
    end 
end