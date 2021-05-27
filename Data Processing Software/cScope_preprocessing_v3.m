% Widefield_preprocessing

% Written by BBS 2017

% Performs the following steps on .avi data files collected with the
% miniscope camera using blue - green strobiscopic illumination.
%0) Rename the miniscope files to enable matlab to determine their
%chronological order
%1) Identification of frames dropped by the miniscope acquistion system and correction of the timestamps file
%2) Concatination of individual .avis into a single movie
%3) Downsampling of the movie in X and Y by binning. Sums the pixel values
%and converts to single format.
%4) Splitting movie into blue and green channels and saving as a .mat
%3) Motion correction (optional)
%4) Saving the blue and green channels as TIFFs (optional)
%5) Synchronization between bcontrol file and miniscope timestamps (optional, requires sessdata)

% The User is prompted for input at 4 points during the script.
% First, they are asked to identify the folder that contains the AVI files
% Second they are asked whether the Time per frame (TPF) is correct
% Third, they are asked to specify which channel is the blue illumination
% channel.
% Finaly they are asked to identify the bcontrol data file.

% Before running the script the user should specifc the following parameters
% params.savetif=0;  set to 1 if you want to save blue and green movies as tifs
% params.synchronize=1; set to 1 if you want to syncronize to session data
% params.normcorre=1; set to 1 if you want to motion correct by normcorre
% params.fps=30; 30 or 60 hz depending on your acquisition speed.
% params.binfac=4; downstampling, reduces file size, increases pixel signal
% to noise and speeds up motion correction.  binfac=4 means that you each
% pixel in the final movie corresponds to a 4by4 pixel region in the
% original movie.

%BBS Jan 2020
%fixed a bug in renumberfile to prevent behavioral movies from being
%deleted.
%fixed a bug in the indexing of the frames.

function cScope_preprocessing_v3
addpath(genpath('C:\Users\roaffan\cScope\codes\cScope-master\Data Processing Software'))

%USERS step 1 please set your parameters!
alertme=1;
params.concat=1;
params.cropme=1;      %crop out edges outside of window
    in_Xcrop=[106;469];      %2x1 double (start,end), 1:480
    in_Ycrop=[4;747];      %2x1 double (start,end), 1:752
params.normcorre=1;   %motion correct by normcorre
params.savetif=0;     %saving the tif is the longest step in this script
params.binfac=2;      %choose 2 or 4
%native movie is 752 by 480

%%OK now select the files you want to preprocess

disp('Please select a data folder')
folder_name = uigetdir; 
cd(folder_name);
disp('Processing...')
if params.concat
    %params.fps=get_fps;
    params.fps=30; %modified 12/17/19 bbs
    renumberfile; %renames the msCam .avi files with the correct number of 0s, 1.avi becomes 001.avi
    if params.fps==60
        %TPF=16.6725;
        TPF=16.6724;
        %TPF=17.08263;
    elseif params.fps==30
        
        TPF=33.35217; %works at BU 11/26/19
    end
    
    TPFcheck=1;
    while TPFcheck
        [ts.newsysClock,ts.newframeNum,ts.delete_these_frames]=miniscope_make_new_timestamps(TPF);
        prompt = 'Is the TPF correct?  Enter 1 for yes, 2 for no ';
        x = input(prompt);
        if x==1
            TPFcheck=0;
        elseif x==2
            display(TPF)
            prompt = 'Please enter new TPF';
            TPF = input(prompt);
            
        else
            disp('Please enter 1 to continue or 2 to repeat)');
        end
    end
    
    %loads each .avi;, resizes them, concatenates and saves as a .tf
    [bluemov,greenmov,params.cropX,params.cropY]=load_and_concatenate(params.binfac,...
                                                                      params.normcorre,...
                                                                      params.cropme,...
                                                                      in_Xcrop,in_Ycrop,...
                                                                      ts); 
    %cd(home_folder)
    % savefast full_mov mov crop ts %savefast can be found in the normcorre library
    save Clocks TPF ts
    save params params
else
    load full_mov
    load Clocks
end

[blue_clock,green_clock,bluemov,greenmov]=separate_and_save(ts.newsysClock,ts.newframeNum,bluemov,greenmov);

%do the final correction to account for black frames
if strcmp(folder_name, 'C:\Users\roaffan\Documents\Miniscope_DAQ_Software-master\x64\Release\data\12_14_2019_ratID36_D2S2\H15_M16_S57_WK_noBeh_D2S2')
    newblue=cat(3,bluemov(:,:,1:8670),greenmov(:,:,8672:end));
    newgreen=cat(3,greenmov(:,:,1:8670),bluemov(:,:,8671:end));
    newblue_clock=[blue_clock(1:8670);green_clock(8672:end)];
    newgreen_clock=[green_clock(1:8670);blue_clock(8671:end)];
    bluemov=newblue;
    greenmov=newgreen;
    blue_clock=newblue_clock;
    green_clock=newgreen_clock;

elseif strcmp(folder_name, 'C:\Users\roaffan\Documents\Miniscope_DAQ_Software-master\x64\Release\data\12_14_2019_ratID36_D2S2\H14_M36_S39_WK_Beh_D2S1')
    newblue=cat(3,bluemov(:,:,1:7179),greenmov(:,:,7181:end));
    newgreen=cat(3,greenmov(:,:,1:7179),bluemov(:,:,7180:end));
    newblue_clock=[blue_clock(1:7179);green_clock(7181:end)];
    newgreen_clock=[green_clock(1:7179);blue_clock(7180:end)];
    bluemov=newblue;
    greenmov=newgreen;
    blue_clock=newblue_clock;
    green_clock=newgreen_clock;

end


savefast preprocess_data_3 bluemov greenmov
% savefast preprocess_data bluemov greenmov

save Clocks_3 blue_clock green_clock TPF ts


if params.savetif
    load params
    savethetiffs(uint16(bluemov),uint16(greenmov));
end

disp('Processing completed')
if alertme
    load handel
    sound(y,Fs)
end
end

%%
function [bmov,gmov,out_Xcrop,out_Ycrop]=load_and_concatenate(binfactor,normcorre,cropme,in_Xcrop,in_Ycrop,ts)
%profile on
evs=1:2:ts.newframeNum(end);
Bfs=ismember(ts.newframeNum,evs); %even numbered frames that exist, i.e. not dropped.


warning('off','MATLAB:audiovideo:aviinfo:FunctionToBeRemoved');
curdir=pwd;
fnames=dir('msCam*.avi');
totalnumframes=0;
for i=1:length(fnames)
    AVI=aviinfo(fnames(i).name);
    totalnumframes=AVI.NumFrames+totalnumframes;
end

%% this part is important for cropping
i=round(length(fnames)/2);
V=mmread([pwd '/' fnames(i).name]);
cd(curdir);
thismov = single(accumfun(3, @(x) x.cdata(:,:,1), [V.frames]));
figure;
   
if isempty(in_Xcrop) && isempty(in_Ycrop)
    while cropme==1
        disp('Please crop the image');
        img=squeeze(mean(thismov(:,:,1:1000),3));
        subplot(2,1,1)
        imagesc(img);
        drawnow
        [Y,X] = ginput(2);
        Y=round(Y); X=round(X);

        newX=floor((X(2)-X(1))/binfactor)*binfactor-1;
        newY=floor((Y(2)-Y(1))/binfactor)*binfactor-1;
        X(2)=X(1)+newX;
        Y(2)=Y(1)+newY;
        subplot(2,1,2)
        imagesc(img(X(1):X(2),Y(1):Y(2)));
        drawnow

        prompt = 'Is the crop correct?  Enter 1 for yes, 2 for no. ';
        x = input(prompt);
        if x==1
            out_Xcrop = X;
            out_Ycrop = Y;
            cropme=0;
        elseif x==2
            cropme=1;
        else
            disp('Please enter 1 to continue or 2 to repeat the crop)');
        end
    end
else   
    img=squeeze(mean(thismov(:,:,1:1000),3));
    subplot(2,1,1)
    imagesc(img);
    drawnow 
    
    X = in_Xcrop;
    Y = in_Ycrop;
    subplot(2,1,2)
    imagesc(img(X(1):X(2),Y(1):Y(2)));
    drawnow
    disp('This is the crop defined by user-defined parameters');
    
    out_Xcrop = in_Xcrop;
    out_Ycrop = in_Ycrop;
end

%% This part initializes the movie

thismov=thismov(X(1):X(2),Y(1):Y(2),:);
off=(i-1)*1000;
theseBfs=Bfs((1:size(thismov,3))+off);

b_template=nanmean(thismov(:,:,theseBfs),3);
g_template=nanmean(thismov(:,:,~theseBfs),3);


[height,width,~]=size(thismov);
height=height/binfactor;
width=width/binfactor;

bmov=zeros(height,width,sum(Bfs),'single');
gmov=zeros(height,width,sum(~Bfs),'single');
frameoffset=0; bmovieoffset=0; gmovieoffset=0;

%% now we open and concatinate the movies
for i=1:length(fnames)
    disp(sprintf('%d of %d',[i length(fnames)]))
    V=mmread([pwd '/' fnames(i).name]);
    cd(curdir);
    thismov = single(accumfun(3, @(x) x.cdata(:,:,1), [V.frames]));
    theseframes=(1:size(thismov,3))+frameoffset;
    thismov=thismov(:,:,~ts.delete_these_frames(theseframes));
    %     try
             theseBfs=Bfs((1:size(thismov,3))+ bmovieoffset+ gmovieoffset);
    %     catch ME
    %         keyboard
    %     end
    smallmov=thismov(X(1):X(2),Y(1):Y(2),:);
    bm=smallmov(:,:,theseBfs);
    gm=smallmov(:,:,~theseBfs);
    
    if normcorre
        if i==1
            [bm_corr,~,b_template] = normcorre_batch(bm);
            [gm_corr,~,g_template] = normcorre_batch(gm);
        else
            [bm_corr,~,~] = normcorre_batch(bm,[],b_template);
            [gm_corr,~,~] = normcorre_batch(gm,[],g_template);
        end   
        bm_rebinned = rebin(bm_corr, binfactor, [1 2], @mean);
        gm_rebinned = rebin(gm_corr, binfactor, [1 2], @mean);
    else
        bm_rebinned = rebin(bm, binfactor, [1 2], @mean);
        gm_rebinned = rebin(gm, binfactor, [1 2], @mean);
        
    end
    
    %smov(:,:,(1:size(mov_rebinned,3))+frameoffset)=mov_rebinned;
    bmov(:,:,(1:size(bm_rebinned,3))+bmovieoffset)=bm_rebinned;
    gmov(:,:,(1:size(gm_rebinned,3))+gmovieoffset)=gm_rebinned;
    frameoffset=frameoffset+1000;
    bmovieoffset=bmovieoffset+size(bm_rebinned,3);
    gmovieoffset=gmovieoffset+size(gm_rebinned,3);
    
end

end

%%
function [blue_clock,green_clock,bluemov,greenmov]=separate_and_save(newsysClock,newframeNum,bluemov,greenmov)

%mov=FI(:,:,logical(~delete_these_frames));
evs=1:2:newframeNum(end);
Bfs=find(ismember(newframeNum,evs));
blue_clock=newsysClock(Bfs);
%bluemov=mov(:,:,Bfs);

ods=2:2:newframeNum(end);
Gfs=find(ismember(newframeNum,ods));
green_clock=newsysClock(Gfs);
%greenmov=mov(:,:,Gfs);

%%If the Arduino is not re-initialized between imaging sessions the frame order may be reversed
% one session you could have BGBG
% on the other you could have GBGB

%lets verifty the movie order is correct
figure;
subplot(2,1,1)
imagesc(bluemov(:,:,2))
title('1')

subplot(2,1,2)
imagesc(greenmov(:,:,2))
title('2')
check_movie=1;
while check_movie
    prompt = 'Which panel is the blue LED movie, 1 or 2? ';
    x = input(prompt);
    if x==1
        check_movie=0;
    elseif x==2
        check_movie=0;
        tmpc=green_clock; green_clock=blue_clock; blue_clock=tmpc;
        tmpm=greenmov; greenmov=bluemov; bluemov=tmpm;
    else
        dis('Please enter 1 or 2');
    end
    
end
end
%%

function savethetiffs(bluemov,greenmov)
fname='bluemov_raw.tif';
imwrite(bluemov(:,:,2), fname,'writemode', 'overwrite')
[~,~,Z]=size(bluemov);
for k = 2:Z
    imwrite(bluemov(:,:,k), fname, 'writemode', 'append');
end
% clear bluemov
%
% [~,~,Z]=size(greenmov);
% fname='greenmov.tif';
% imwrite(greenmov(:,:,2), fname,'writemode', 'overwrite')
% for k = 2:Z
%     imwrite(greenmov(:,:,k), fname, 'writemode', 'append');
% end
end



%%
function renumberfile
files=dir('msCam*.avi');
for i =1:length(files)
    n=regexp(files(i).name,'\d');
    if numel(n)==1
        newname=['msCam00' num2str(files(i).name(n)) '.avi'];
        movefile(files(i).name,newname)
    elseif numel(n)==2
        newname=['msCam0' num2str(files(i).name(n)) '.avi'];
        movefile(files(i).name,newname)
    end
end
end
%%
function [fps,TPF]=get_fps

timestamps=importdata('timestamp.dat');
sysClock=timestamps.data(:,3);
% - RA:
% sysClock=sysClock(1:2:end);
%
TPF=mode(diff(sysClock));
fps=1000/TPF;
fps=round(fps/10)*10;
end

%% Concatenate results of arrayfun() along the specified dimension.
function B = accumfun(dim, func, varargin)

if iscell(varargin{1})
    fcn = @cellfun;
else
    fcn = @arrayfun;
end
B     = fcn(func, varargin{:}, 'UniformOutput', false);
B     = cat(dim, B{:});

end

%%
function [newstack]=resize_stack(stack,n)

[X,Y,Z]=size(stack);
XX=floor(X/n);
YY=floor(Y/n);
newstack=nan(XX,YY);
for z=1:Z
    for x=1:XX
        for y=1:YY
            pixel=stack(((x-1)*n+1):x*n,((y-1)*n+1):y*n,z);
            newstack(x,y,z)=mean(pixel(:));
        end
    end
end
end
