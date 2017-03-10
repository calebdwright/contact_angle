
tic %Basic Setup
clear; clc; close all;
warning('off', 'Images:initSize:adjustingMag');
[InputVideo, filepath] = uigetfile('*.tif','Please Select One Stacked Tif File');
filename = fullfile(filepath, InputVideo);
InfoImage = imfinfo(filename);
numFrames = length(InfoImage);
disp(num2str(numFrames));
width = InfoImage(1).Width;
height = InfoImage(2).Height;

vertSections = 6; %Tweakable Variables
horzSections = 4;
vidSections = 4;
numPoints = 2;
mask_size = 20;

point_array_x = [numFrames,1]; %Find strongest point in every frame
point_array_y = [numFrames,1];
for i = 1 : numFrames
    frame1 = imread(filename, 'Index', i);
    points = detectHarrisFeatures(frame1, 'FilterSize', 5, 'MinQuality', .8, 'ROI', [width/horzSections, height/vertSections, (horzSections-2)*width/horzSections, (vertSections-2)*height/vertSections]);
    strongest = points.selectStrongest(1);
    cp1 = int16(strongest.Location);
    if size(cp1,1) ~= 0
        point_array_x(i) = cp1(1);
        point_array_y(i) = cp1(2);
    else
        point_array_x(i) = width/2;
        point_array_y(i) = height/2;
    end;
end;

count_horz = zeros(horzSections, 1); %Find horizontal sections with most points
for i = 1:horzSections
    lowerHorzRange = (i-1) * width/horzSections;
    upperHorzRange = i * width/horzSections;
    for k = 1 : size(point_array_x,2)
        if point_array_x(k) <= upperHorzRange && point_array_x(k) > lowerHorzRange
            count_horz(i) = count_horz(i) + 1;
        end;
    end;
end;
most_indexed = zeros(numPoints, 1);
for i = 1:numPoints
    for k = 1:horzSections
        if count_horz(k) >= most_indexed(i)
            if i > 1 && most_indexed(i-1) ~= k
                most_indexed(i) = k;
            elseif i == 1
                most_indexed(i) = k;
            end;
        end;
    end;
end;

point_array_x = [numFrames,numPoints]; %Find contact angles for horizontal sections with most points
point_array_y = [numFrames,numPoints];
contact_angle_array = [numFrames,numPoints];
for u = 1:numPoints
    vertSection = int16(most_indexed/vertSections);
    horzSection = int16(mod(most_indexed, horzSections));
    ROI = [((most_indexed(u)-1) * width / horzSections)+1, (height / vertSections)+1, width/horzSections, (vertSections-2) * height/vertSections];
    for i = 1 : numFrames
        frame_second = imread(filename, 'Index', i);
        points = detectHarrisFeatures(frame_second, 'ROI', ROI, 'FilterSize', 5, 'MinQuality', .95);
        strongest = points.selectStrongest(1);
        cp2 = int16(strongest.Location);
        if size(cp2,1) ~= 0
            point_array_x(i,u) = cp2(1);
            point_array_y(i,u) = cp2(2);
            D = frame_second(cp2(2)-50:cp2(2)+50,cp2(1)-80:cp2(1)+80);
            D = double(imbinarize(D));
            
            Mask=zeros(2*mask_size+1,2*mask_size+1); % initilize the mask
            for k=1:mask_size+1
                Mask(k,mask_size+1-k+1:mask_size+k)=1;
            end;
            Mask=flipud(Mask); % rotate vertically
            
            conv=conv2(D,Mask,'same'); % Doing 2-D convolution operation
            convf=conv2(flipud(D),Mask,'same');
            R_square=(mask_size+1)^2; % value of R^2
            A=conv(51,81); % shaded area
            lambda=sign(0.5*R_square-A); % sign function
            alpha=pi/2*(1-lambda)+atan(lambda*(R_square*(1-lambda)-2*A)/(2*A-R_square)); % using pre-defined equation
            theta_h=(180/pi)*alpha; % converting to degree
            
            A=convf(51,81);
            lambda=sign(0.5*R_square-A); % sign function
            alpha=pi/2*(1-lambda)+atan(lambda*(R_square*(1-lambda)-2*A)/(2*A-R_square)); % using pre-defined equation
            theta_t=(180/pi)*alpha;
            
            contact_angle=180-(theta_h+theta_t);
            contact_angle_array(i,u) = contact_angle;
        else
            contact_angle_array(i,u) = 0;
            point_array_x(i,u) = height/2;
            point_array_y(i,u) = width/2;
        end;
    end;
end;

imshow(frame_second); %Plot results
hold on;
colors = ['r', 'b', 'm', 'y'];
for i = 1:numPoints
    plot(point_array_x(:,i), point_array_y(:,i), strcat(colors(i), '*'));
end;
hold off;
figure;
types = ['s', 'o', 'x', '+'];
hold on;
for i = 1:numPoints
    frameCount = [];
    k = 1;
    contact_angle_array_2 = [];
    for j = 1:size(contact_angle_array, 1)
        if contact_angle_array(j,i) > 0
            frameCount(k) = j;
            contact_angle_array_2(k) = contact_angle_array(j,i);
            k = k+1;
        end;
    end;
    contact_angle_array_3 = smooth(frameCount, contact_angle_array_2, 21, 'rloess');
    plot(frameCount, contact_angle_array_3, colors(i));
end;
hold off;
toc