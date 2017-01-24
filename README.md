# GCaMP
% Pedro Falcon
% May 02, 2016
% ONLY FOR MUTANT MOVIES

clear; close all; clc;
%% %% Identifying file and importing movie

%Getting image and info 
[FileName,PathName] = uigetfile('*.tif','Select image to be processed/analyzed');
SplitFileName = strsplit(FileName,'.');
SplitFileName_1 = char(SplitFileName{1});
InfoImage = imfinfo(FileName);

%Getting image size and # of frames from tiff stack
nImage=InfoImage(1).Height;
mImage=InfoImage(1).Width;
NumberImages=length(InfoImage);

%Preallocating 3D matrix for tiff stack
ImportedStack=zeros(nImage,mImage,NumberImages,'uint8');
 
%Populating 3D matrix from tiff stack 
TifLink = Tiff(FileName, 'r');
for i = 1:NumberImages
   TifLink.setDirectory(i);
   ImportedStack(:,:,i )= TifLink.read();
end
TifLink.close();

%% Loading GCaMP timeseries and timepoints from Fiji preprocessing
TimePointDataPath = strcat(PathName,'/TimePointData/');
cd(TimePointDataPath);
SplitTitle = strsplit(SplitFileName_1,'_');
    %Because TimePointFileName doesn't have the appended _A
    %Timepoints in movie don't depend on analysis, and so they are connected to
    %the raw movie.  .dat file could be altered to include the analysis in the
    %timepoint label if necessary.  
TimePointFileName = strcat(char(SplitTitle{1}),'_',char(SplitTitle{2}),'_',char(SplitTitle{3}),'_TimePoints.dat');
TimePoints=csvread(TimePointFileName);

%% Geting Distances
wl = 2; %Horizontal Window Length for Envelope Function

TimePoints(2) = 231;
TimePoints(3) = 231;
lengthOfArray = TimePoints(3) - TimePoints(2) + 1; %Length of the Maximum Distance Arrays

j = 1; %Index Of Maximum Distances (notice it starts @ 1)
n= 1;
x = [1, mImage];
y = [n, n];
i = 0;

global xmaxindexTemp maxdistance yvaluePeak xmaxindex;

maxdistance = zeros(nImage,1); %Max Distances in one frame image
maximumDistancePeak = zeros(lengthOfArray,1); %Max distance using peak env
yvaluePeak = zeros(lengthOfArray,1); %Max Y values over time for boundaries
F(lengthOfArray) = struct('cdata',[],'colormap',[]);
xmaxindex = [0,0];
wdist = zeros(lengthOfArray,2);
tempUp1 = zeros(2*lengthOfArray,1);
tempUp2 = zeros(2*lengthOfArray,1);
tempDown1 = zeros(2*lengthOfArray,1);
tempDown2 = zeros(2*lengthOfArray,1);
xmaxindexTemp = zeros(1,2);
%% Calculation of Mean Values in a Movie


%% LOOP

frame = TimePoints(2); %Initial Frame

%Finding Longest Horizontal line, its Y-Index, and Pixel Intensity at its
%boundaries.
[tempYVal,pi_1,pi_2,xmaxindexTemp] = firstPI_MUTANT(ImportedStack,frame,mImage,nImage,wl);

for frame = TimePoints(2):TimePoints(3)
    
    [yvaluePeak,xmaxindex,maxdistance] = HorizontalLineLoop( ImportedStack, frame, wl,j, tempYVal,mImage, yvaluePeak, xmaxindex,xmaxindexTemp);
    maximumDistancePeak(j) = max(maxdistance);
    
    %% Finding First Distances (VERTICAL)
    n = xmaxindex(1) + 5; %First Horizontal Boundary
    %n=460;
    Q = ImportedStack(:,:,TimePoints(2)); % Getting Image Profile
    [wy1,wy2] = firstYMembraneBoundaries(Q,nImage,frame,n);


    promptt1 = '\nUpper Boundary: ';
    wy1 = input(promptt1);

    promptt2 = 'Lower Boundary: ';
    wy2 = input(promptt2);

    global tempUpperBoundary tempLowerBoundary yindeces

    tempUpperBoundary = wy1;
    tempLowerBoundary = wy2;
    close all
    
    %% Vertical Distances
    
    dataUpper = zeros(lengthOfArray,mImage); %2D Matrix 
    dataLower = zeros(lengthOfArray,mImage); %2D Matrix 
    %LengthOfArray = time. This will determine # of rows
    %mImage = columns. 
    
    yindeces = zeros(2,1);
    
    h = 1; %IMAGE COUNTER
    
    for n=(xmaxindex(1)+5):(xmaxindex(2)-5)
        J = ImportedStack(:,:,frame); 
        B = imsharpen(J, 'radius', 0.5, 'amount', 1.2); % Sharpening image
        x = [n,n];
        y = [1,nImage];
        %[e,env,lowen ] = improfile_avg1(B,nImage,n); %IMPROFILE LINE AVG (3PIXELS)
        [~, ~, e] = improfile(B,[n,n],[1,nImage]); % Getting Image Profile
        
        %FILTERING
        windowSize = 5; %Creating
        b = (1/windowSize)*ones(1,windowSize);
      
        ShifterFiltered = filter(b,1,e);
        Filler(1:(windowSize-1))=ShifterFiltered(windowSize+1);
        eFiltered=vertcat(Filler',ShifterFiltered((windowSize):length(ShifterFiltered)));
        %END OF FILTER
        
        [env,low] = envelope(eFiltered,1,'peak');
        env = (env + low)/2;

        [yheight, yindex] = findpeaks(env,'MinPeakProminence',0); %FINDING PEAKS
        
        l = length(yindex);
        secondindex = -10;
        halfLength = l/2; %Finding first half of yindeces
        halfLength = round(halfLength,0); %rounding in case of odd number
        
        global tempHeight_1 tempHeight_2;
        
        tempHeight_1 = nan;
        tempHeight_2 = nan;
        
        try
        %Finding Upper Boundary
        [tempUpperBoundary,yindeces,secondindex] = tempUpperBoundary4(yindex,tempUpperBoundary,yindeces,yheight,l);
        catch
            fprintf('TEMP_UPPER5\n');
            [tempUpperBoundary,yindeces,secondindex] = tempUpperBoundary5(yindex,tempUpperBoundary,yindeces,yheight,l );
        end
        
        try
        %Finding Lower Boundary
        [tempLowerBoundary,yindeces] = findTempLower4(secondindex, yindex,tempLowerBoundary,yindeces,yheight,l);
        catch
        end
        fprintf('ENDED: %f\n',n);
        %Saving Indeces
        datarow = frame - TimePoints(2) + 1;
        dataUpper(datarow, n) = tempUpperBoundary;
        dataLower(datarow, n) = tempLowerBoundary;
        
    %% Make Images
    set(figure(h),'visible','off');

    imagesc(ImportedStack(:,:,frame))
    hold on
    plot(xmaxindex(1),yvaluePeak(j),'r--o')
    plot(xmaxindex(2),yvaluePeak(j),'r--o')
    plot(n,tempUpperBoundary,'w--*')
    plot(n,tempLowerBoundary,'w--*')
    title(frame)
    hold off
    name = sprintf('frame_%d_index%d',j,n);
    saveas(gcf,name,'tiff');
    F(h) = getframe(gcf);
    h = h + 1;
    end 
    j = j+1;
    
end

%% PLAYING MOVIE

prompt1 = 'Press 1 and then enter to play movie: ';
playmovie = input(prompt1); %Wplaymovie

prompt2 = 'Enter FPS: ';
fps = input(prompt2);

fig = figure;
movie(fig,F,1,fps);
