
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
horzSections = 8;
vidSections = 3;
numPoints = 4;
mask_size = 20;

point_array_x = []; %Find strongest points in every frame
point_array_y = [];

for h = 1:vidSections
for i = 1 : h*numFrames/vidSections
    frame1 = imread(filename, 'Index', i);
    points = detectHarrisFeatures(frame1, 'FilterSize', 5, 'MinQuality', .8, 'ROI', [width/horzSections, height/vertSections, (horzSections-2)*width/horzSections, (vertSections-2)*height/vertSections]);
    for j = 1:length(points)
        pt = points(j).Location;
        point_array_x(i*j,h) = pt(1);
        point_array_y(i*j,h) = pt(2);
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
[count_horz_sorted, most_indexed_x] = sort(count_horz, 'descend');

count_vert = zeros(vertSections, 1); %Find vertical sections with most points (assumes all points are roughly parallel)
for i = 1:vertSections
    lowerVertRange = (i-1) * height/vertSections;
    upperVertRange = i * height/vertSections;
    for k = 1 : size(point_array_y,2)
        if point_array_y(k) <= upperVertRange && point_array_y(k) > lowerVertRange
            count_vert(i) = count_vert(i) + 1;
        end;
    end;
end;
[count_vert_sorted, most_indexed_y] = sort(count_vert, 'descend');
most_indexed_y = most_indexed_y(1);

point_array_x_2 = []; %Find contact angles for horizontal sections with most points
point_array_y_2 = [];
contact_angle_array = [];
for u = 1:numPoints
    vertSection = int16(most_indexed_x/vertSections);
    horzSection = int16(mod(most_indexed_x, horzSections));
    ROI = [((most_indexed_x(u)-1) * width / horzSections)+1, ((most_indexed_y-1) * height / vertSections)+1, width/horzSections, height/vertSections];
    for i = 1 : numFrames
        frame2 = imread(filename, 'Index', i);
        points_2 = detectHarrisFeatures(frame2, 'ROI', ROI, 'FilterSize', 5, 'MinQuality', .9);
        strongest = points_2.selectStrongest(1);
        pt = int16(strongest.Location);
        if size(pt,1) ~= 0
            point_array_x_2(i,u) = pt(1);
            point_array_y_2(i,u) = pt(2);
            D = frame2(pt(2)-(height/(vertSections*4)):pt(2)+(height/(vertSections*4)),pt(1)-(width/(horzSections*4)):pt(1)+(width/(horzSections*4)));
            D = double(imbinarize(D));
            
            Mask=zeros(2*mask_size+1,2*mask_size+1); % initilize the mask
            for k=1:mask_size+1
                Mask(k,mask_size+1-k+1:mask_size+k)=1;
            end;
            Mask=flipud(Mask); % rotate vertically
            
            conv=conv2(D,Mask,'same'); % Doing 2-D convolution operation
            convf=conv2(flipud(D),Mask,'same');
            R_square=(mask_size+1)^2; % value of R^2
            A=conv(round((height/(vertSections*4)+1)),round((width/(horzSections*4))+1)); % shaded area
            lambda=sign(0.5*R_square-A); % sign function
            alpha=pi/2*(1-lambda)+atan(lambda*(R_square*(1-lambda)-2*A)/(2*A-R_square)); % using pre-defined equation
            theta_h=(180/pi)*alpha; % converting to degree
            
            A=conv(round((height/(vertSections*4)+1)),round((width/(horzSections*4))+1));
            lambda=sign(0.5*R_square-A); % sign function
            alpha=pi/2*(1-lambda)+atan(lambda*(R_square*(1-lambda)-2*A)/(2*A-R_square)); % using pre-defined equation
            theta_t=(180/pi)*alpha;
            
            contact_angle=180-(theta_h+theta_t);
            contact_angle_array(i,u) = contact_angle;
        else
            contact_angle_array(i,u) = 0;
            point_array_x_2(i,u) = 1;
            point_array_y_2(i,u) = 1;
        end;
    end;
end;

imshow(frame2); %Plot results
hold on;
 colors = ['r', 'b', 'm', 'y'];
for i = 1:numPoints
    plot(point_array_x_2(:,i), point_array_y_2(:,i), strcat(colors(i), '*'));
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