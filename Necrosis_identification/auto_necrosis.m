%% Matlab code to detect Necrosis and Microvascular proliferation detection on WSI (whole slide images) of TCGA-GBM slides
%% inline comments - incorporated by Raj - 06/01/18
clear all
clc
close all
 
 

disp('Please enter a number between 1 to 8 ')
Index = input('to choose among one of 8 sample images provided:  ')    ;
disp('Please wait...')
  
% FileName = strcat(num2str(Index),'.png');
% img = imread(FileName);
img = imread('/home/raj/Dropbox/files/EmoryUni/Lab_SourceCode/Auto_discriminationof_gliomas_Necrotic_MVP/2.png);
img_Orig = img;

%% Necrosis Detection

%% Cell Segmentation
imGray = img(:,:,1); 
imGamma = imadjust(imGray,stretchlim(imGray),[]); % Contrast enhancement
bwGamma = 1- im2bw(imGamma,graythresh(imGamma));   % Thresholding
 
SE1=strel('disk',3); % Structuring - morpological elments, disk size element with radius "3"
SE2=strel('disk', 6);
SE4=strel('disk',4);
imOpened= imopen(bwGamma,SE1);      % opening
[Labels, NumOfObj] = bwlabel(imOpened,4); % Label connected components returns lable matrix

SizeConstant  = 400 ;
s = regionprops(Labels, 'Area', 'BoundingBox'); % returns struct array containing a struct for each object in the image
s(1);
area_values = [s.Area];
idx = find(( area_values> SizeConstant) );
% Adaptive opening - Image filtering, across each region, erosion followed
% by dilation
imOpened2 = imOpened;
for index = idx
    
        aaa = max(1 , floor(s(index).BoundingBox) );
        StartingPoints = aaa(1:2);
        Width = aaa(3:4);
        imOpenedTemp = imopen( imOpened( StartingPoints(2):StartingPoints(2)+Width(2)-1,StartingPoints(1):StartingPoints(1)+Width(1)-1 ) ,SE2);
        imOpened2( Labels==index ) = 0;
        imOpened2( StartingPoints(2):StartingPoints(2)+Width(2)-1,StartingPoints(1):StartingPoints(1)+Width(1)-1 ) = imOpened2( StartingPoints(2):StartingPoints(2)+Width(2)-1,StartingPoints(1):StartingPoints(1)+Width(1)-1 ) | imOpenedTemp;
        
        
   
end
        imOpened3= imopen(imOpened2,SE4);
        imOpenedFinal = imOpened3;        
         [ FinalLabels ,FinalNumOfObj ] = bwlabel(imOpenedFinal);
        
 
BlockSize1 = 50;
Step = 30;
SizeConstant1  = 1000 ;
  
%% Cell Count profile

%% Celll counting for a fixed-size window, by rolling over entire window, by fixed-size step

SizeConstant2  = 1000 ;
for i = Step:Step:size(imOpenedFinal,1)-BlockSize1-Step-1
    for j = Step:Step:size(imOpenedFinal,2)-BlockSize1-Step-1
        
        NewBlock = imOpenedFinal(i:i+BlockSize1-1,j:j+BlockSize1-1);
        [ Labels , Num ] = bwlabel(NewBlock);
            

        
        s2 = regionprops(Labels, 'Area', 'BoundingBox');
        area_values2 = [s2.Area];
        idx2 = find(( area_values2< SizeConstant2) );
        CellCount2(i/Step,j/Step) = length(idx2);
    end
end
CellCount2 = flipud(CellCount2);
H = fspecial('disk',7);

%% filters the multidimensional array Cellcount2 with the multidimensional filter 
%% on fspecial filter of disk with radius "7"
MeanCellCount2  = imfilter(CellCount2,H,'replicate'); 

MaxIm2= max(max(MeanCellCount2));
MinIm2= min(min(MeanCellCount2));
Necrosis2 = (MeanCellCount2 ) /(MaxIm2);
Threshold2 = 0.7*graythresh(Necrosis2);

NecrosisBW2 = im2bw(Necrosis2,Threshold2); %% converts grayscale image to binary image
NecProp2 = regionprops(~NecrosisBW2, 'all' ) ; %% gets measure of converted binary image

 
NecRoot1_2 = 0;
InvertedNecrosis2 = 1-Necrosis2;
Region = 20;

% Valley Detection
window  = 81;
hLocalMax = vision.LocalMaximaFinder;
hLocalMax.MaximumNumLocalMaxima = 10; 
hLocalMax.NeighborhoodSize = [window  window ]; 
jay = max(max(InvertedNecrosis2(20:size(InvertedNecrosis2,1)-20,20:size(InvertedNecrosis2,2)-20))); 
hLocalMax.Threshold = 0.98*jay; 
location = step(hLocalMax, InvertedNecrosis2);
%% Decision Tree
for i=1:size(location,1)
    ValleyValue2(i) = MeanCellCount2 (location(i,2),location(i,1) );
    for iii=1:size(NecProp2,1)
        Temp_area_indicator = [(NecProp2(iii).PixelList(:,1) == location(i,1) ) ,    ( NecProp2(iii).PixelList(:,2) == location(i,2)  )]';
        if sum( sum(Temp_area_indicator ) == 2  ) == 1
            Area_Nec = NecProp2(iii).Area;
        end
    end
    if Area_Nec > 10
        if ValleyValue2(i) > 8
            % probbably not a Nec
        else
            MaxAround = max(max( MeanCellCount2...
                ( max(location(i,2)-Region,1):min(location(i,2)+Region,size(MeanCellCount2,1) ) ...
                , max(location(i,1)-Region,1):min(location(i,1)+Region,size(MeanCellCount2,2) ) ) ));
            if ValleyValue2(i)>1
                if MaxAround > 3* ValleyValue2(i)
                    NecRoot1_2 =  NecRoot1_2 + 1 ;
                end
            elseif ValleyValue2(i)<=1
                if MaxAround > 3
                    NecRoot1_2 =  NecRoot1_2 + 1 ;
                end
            end
        end
    end
end

Decision1(Index) =  (NecRoot1_2 >= 1);
NecrosisDec=0;
if Decision1(Index) >0
    NecrosisDec = 1;
end
 
 
%% MVP Detection

Result =0;
 
 FileName  = strcat(num2str(3),'.tif');
ref = imresize(imread(FileName),1/6);
refhist = imhist(ref(:,:,2),32);


 img = imresize(img,1/6);
R = img(:,:,1);
G = img(:,:,2);
B = img(:,:,3);

Decision  = zeros(1,400);

 
G2 = histeq (G,refhist);

t=3;
s=floor(size(img,1)/t);
Ibw6 = zeros(size(G));
for iii =  0:t-1
    for jjj = 0:t-1
        G1 = G( iii*s+1: (iii+1)*s , jjj*s+1: (jjj+1)*s  );
%         G1=G;
        G1 = histeq (G1,refhist); % Color Normalization
        G1 = imadjust(G1,stretchlim(G1),[]); % Contrast enhancement
        Thresh = 0.5;
        ratio = 1;
        % Adaptive Thresholding
        while ratio >0.35
            
            b1 = im2bw(255-G1,Thresh);
            ratio  =sum(sum(b1==1))/ (size(G1,1)*size(G1,1) );
            Thresh = Thresh + 0.02;
        end
        
        
        se1=strel('rectangle',[2 2]);
        Ibwerode = imerode( b1 ,se1);
        Ibw2 = bwareaopen(Ibwerode,30);
        Ibwfilled = imfill(Ibw2,'holes');
        se2=strel('disk',3);
        Ibw3 = imclose(Ibw2,se2);
        Ibwfilled1 = imfill(Ibw3,'holes');
        se2=strel('disk',2);
        Ibw4 = bwareaopen(Ibwfilled1,100);
        se5=strel('disk',3);
        Ibw5 = imdilate(Ibw4,se5);
        Ibw6(iii*s+1: (iii+1)*s , jjj*s+1: (jjj+1)*s ) = Ibw5;
%         Ibw6 = Ibw5;
        
    end
end

 imbw = 1-Ibw6;

% median filtering
x =5;  
y =5;  
dim = size(img);
img_median_g = medfilt2(G2,[x y],'symmetric');

% Mean filtering
z = 10; 
t = fspecial('disk', z);

im1 = imfilter(img_median_g,t,'replicate');

imd = im2double(im1);

% Valley Detection
window_g = 101;
hLocalMax = vision.LocalMaximaFinder;
hLocalMax.MaximumNumLocalMaxima = 10; 
hLocalMax.NeighborhoodSize = [window_g window_g]; 
jay = max(max(1-imd(20:size(1-imd,1)-20,20:size(1-imd,2)-20))); 
hLocalMax.Threshold = 0.81*jay; 
location_g = step(hLocalMax, 1-imd);
DotSize=10;
 
[imLabeled , num]  = bwlabel(1-imbw);

imbwNew  = zeros(size(imbw));
for i=1:size(location_g,1)
    if imLabeled(location_g(i,2),location_g(i,1)) > 0 
        Temp = imLabeled(location_g(i,2),location_g(i,1));
        imbwNew ( find(imLabeled == Temp) )=  1;
    end
    
end

imbwNew = logical(imbwNew);
Properties = regionprops(imbwNew,'all');

imOut1 = img;

ROCIndex=0;
ROCIndex = ROCIndex+1;
   
Decision  = zeros(1,400);

%% Decision Tree
for i=1:length(Properties)
    if Properties(i).Area < 2500 && Properties(i).Area > 350 
        
        if Properties(i).Perimeter <250
            if Properties(i).Extent >0.3 
                if Properties(i).EquivDiameter > 13
                    if Properties(i).EquivDiameter <26
                        if Properties(i).MajorAxisLength <43
                            if Properties(i).MinorAxisLength >10
                                Decision(i) = Decision(i) + 1;
                            end

                        end

                    else
                        if Properties(i).Extent > 0.45 
                            if Properties(i).MinorAxisLength >15
                                
                                Decision(i) = Decision(i) + 1;
                            end
                        end
                    end
                end
            end
        else
            if Properties(i).MajorAxisLength > 84
                if Properties(i).MinorAxisLength >20
                    if Properties(i).Extent < 0.7
                        Decision(i) = Decision(i) + 1;
                    end
                end
            end
            
        end
    end
    
    
end

DecisionTemp(Index) = sum(Decision);

MVPDec = 0; 
if DecisionTemp(Index) >1
    MVPDec = 1;
end


%% Decision

if MVPDec == 1 || NecrosisDec ==1
    disp('Decision: This image belongs to HGG');
else
    disp('Decision: This image belongs to LGG');
end


figure(1);imshow(imresize(img_Orig,1/6));title('Original Image') % show the input image



 
      