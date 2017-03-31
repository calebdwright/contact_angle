clear; clc; close all;

dir
FileList = dir('*.avi');
for i = 1:size(FileList, 1)
    filename = FileList(i).name;
    V = VideoReader(filename);
    numFrames = V.NumberOfFrames;
    disp(num2str(numFrames));
    V = VideoReader(filename);

    tic
    frame = readFrame(V);
    imwrite(frame, strcat(filename, '.tif'));
    for i = 2:numFrames
        frame = readFrame(V);
        imwrite(frame, strcat(filename, '.tif'), 'WriteMode', 'append');
        %disp(num2str(i));
    end;
    toc
end;