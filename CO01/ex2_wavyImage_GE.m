% CO01 - Exercise 2 - Wavy Image 



function [image_out] = wavy_image(image_in)
    Y = size(image_in, 1);
    X = size(image_in, 2);
    theta = 5;
    A = 50;  % Set the Amplitude
    k = 13; % Set the frequency 
    max_shift = 2*A; % CHANGED for new wavy_image (allow for 2A extra pixels)
    image_out = zeros(Y, X+max_shift, 3);
    for y=1:Y
        
        %% previous lines of slanted_image function
        % local_shift = ceil(y*tand(theta));
        % local_x = 1:X;

        %% modify these lines for this wavy_image function
        local_shift = A + ceil(A * sin (k * y));
        local_x = 1:X;
        
        local_x = local_x + local_shift;
        image_out(y, local_x, :)=image_in(y, :, :);
    end
end


% Load the image
MyImage = imread('Hubble.jpg');

% Normalise the colours
MyImage = double(MyImage)/255;

% Call the Wavy_image function with output name MyNewImage
[MyNewImage] = wavy_image(MyImage);

% Display the new image
image(MyNewImage);