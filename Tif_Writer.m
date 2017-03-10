clear; clc; close all;

[InputVideo, filepath] = uigetfile('*.avi','Please Select One Video File');
filename = fullfile(filepath, InputVideo);
V = VideoReader(filename);
numFrames = V.NumberOfFrames;
disp(num2str(numFrames));

tic
V = VideoReader(filename);
frame = readFrame(V);
imwrite(frame, strcat(filename, '.tif'));
for i = 2:numFrames
    frame = readFrame(V);
    imwrite(frame, strcat(filename, '.tif'), 'WriteMode', 'append');
    disp(num2str(i));
end;
toc