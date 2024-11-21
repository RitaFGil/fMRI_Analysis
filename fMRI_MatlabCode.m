%% fMRI Data analysis
% Written by Rita Gil
% Modified last in 15.12.2023

% ROI Analysis:
% Load Data
% Delineate ROIs
% Detrend runs
% Average run
% Separate individual cycles
% Mean adn s.e.m. of averaged cycle

% Quantification:
% Join conditions
% Calculate the averaged cycle per run
% Interpolation of averaged cycles
% Define steady-state interval for each ROI
% Box plots for the steady state interval

% Neurometric curve:
% Interpolation of each trial
% Categorize trials as a 'pulsating' or 'continuous' light report
% Fit psichometric curve
% Boostraping with replacement for final threshold calculation (meand and standard deviation)

%% Initialization
clear all
close all
clc
pack

%Example group with three different frequency conditions (1Hz, 15Hz and 25Hz)


%Group Parameters
path='XXX'; %string where data is stored
endpath2save='XXX'; %string where data is going to be saved

% Cell with animals for each frequency condition
cond1={'Animal1','Animal2','Animal3','Animal4','Animal5'};
cond2={'Animal1','Animal2','Animal3','Animal4','Animal5'};
cond3={'Animal1','Animal2','Animal3','Animal4','Animal5'};

% Folders for each conditions (nfolderX represents which folders are from which animal)
folder1=[1,2,3,4,5,6,7,8,9,10];
nfolder1=[2,2,2,2,2];
folder2=[1,2,3,4,5,6,7,8,9,10];
nfolder2=[2,2,2,2,2];
folder3=[1,2,3,4,5,6,7,8,9,10];
nfolder3=[2,2,2,2,2];

numcond=3; %number of conditions
totalFolder=length(folder1)+length(folder2)+length(folder3); %total number of folders used

ratref='Rat1'; %Raference animal used in the fMRI preprocessing steps
search_name = 'srwtaiImage'; %Saved name for reference animal
search_name1 ='srnwtaiImage'; %Saved names of other animals in the group

nslices=8; %total number of slices
slicesout=4; %slices that will not be considered

numROI=3; %total number of ROIs
ROIstr={'VC','SC','LGN'}; %ROI string
stringcond={'1Hz','15Hz','25Hz'}; %Condition string


%fMRI Parameters
nrep=270; %total number of repetitions
TR=1.5; %TR in sec

rest=30; %Number of repetitions during resting period
stim=10; %Number of repetitions during stimulation period
lcycle=stim+rest; %Number of repetitions during 1 cycle= rest period + stimulation period
stimTime=TR*stim; %Stimulation period in seconds
restTime=TR*rest; %Resting  period in seconds

ncycle=6; %number of cycles per run
start=rest+1; %when first stim starts

%fMRI Paradigm
event=zeros(1,nrep);
event(rest+1:lcycle)=1;
event(lcycle+rest+1:2*lcycle)=1;
event(2*lcycle+rest+1:3*lcycle)=1;
event(3*lcycle+rest+1:4*lcycle)=1;
event(4*lcycle+rest+1:5*lcycle)=1;
event(5*lcycle+rest+1:6*lcycle)=1;
event1=event;

for i=1:nrep
    if event(i)==0
        event1(i)=nan;
    end
end

%Saving inputs
cd(endpath2save)
save('workspace.mat')

%% Loading data

j=1;
for i=1:length(nfolder1)
    nruns=nfolder1(i);
    rat=i;
    rat
    for pos=j:j+nruns-1
        for r=1:nrep
            clear tmp
            tmp=cond1{rat}==ratref;
            if prod(tmp)==0
                data(pos,r)=spm_vol([path,cond1{rat},filesep,num2str(folder1(pos)),'/Processed/',search_name1,'_',num2str(folder1(pos)),'_',num2str(r,'%04d'),'.nii']); %path where niftis are stored
            else
                data(pos,r)=spm_vol([path,cond1{rat},filesep,num2str(folder1(pos)),'/Processed/',search_name,'_',num2str(folder1(pos)),'_',num2str(r,'%04d'),'.nii']);
            end
        end
    end
    j
    j=j+nruns;
end

j=1;
for i=1:length(nfolder2)
    nruns=nfolder2(i);
    rat=i;
    for pos=j:j+nruns-1
        posnew=pos+length(folder1);
        
        for r=1:nrep
            clear tmp
            tmp=cond2{rat}==ratref;
            if prod(tmp)==0
                data(posnew,r)=spm_vol([path,cond2{rat},filesep,num2str(folder2(pos)),'/Processed/',search_name1,'_',num2str(folder2(pos)),'_',num2str(r,'%04d'),'.nii']);
            else
                data(posnew,r)=spm_vol([path,cond2{rat},filesep,num2str(folder2(pos)),'/Processed/',search_name,'_',num2str(folder2(pos)),'_',num2str(r,'%04d'),'.nii']);
            end
        end
    end
    j
    j=j+nruns;
end

j=1;
for i=1:length(nfolder3)
    nruns=nfolder3(i);
    rat=i;
    for pos=j:j+nruns-1
        posnew=pos+length(folder1)+length(folder2);
        
        for r=1:nrep
            clear tmp
            tmp=cond3{rat}==ratref;
            if prod(tmp)==0
                
                data(posnew,r)=spm_vol([path,cond3{rat},filesep,num2str(folder3(pos)),'/Processed/',search_name1,'_',num2str(folder3(pos)),'_',num2str(r,'%04d'),'.nii']);
                
            else
                data(posnew,r)=spm_vol([path,cond3{rat},filesep,num2str(folder3(pos)),'/Processed/',search_name,'_',num2str(folder3(pos)),'_',num2str(r,'%04d'),'.nii']);
            end
        end
    end
    j
    j=j+nruns;
end


clear temp;
nslicesEfec=nslices-slicesout; %final number of slices being used

% Selecting only slices where ROIs of interest are and rotating images
temp=spm_read_vols(data(1,1));
avData=zeros(size(temp,2),size(temp,1),nslicesEfec,nrep,totalFolder);
temp1=zeros(size(temp,1),size(temp,2),nslicesEfec,nrep,totalFolder);

clear temp
for f=1:totalFolder
    f
    for r=1:nrep
        clear temp
        temp(:,:,:,r)=spm_read_vols(data(f,r)) ;
        temp1(:,:,:,r,f)=temp(:,:,3:6,r);
    end
    avData(:,:,:,:,f)=flip(imrotate(temp1(:,:,:,:,f),-90),2);
    
end

%saving organized preprocessed raw data
save('Data.mat','data','-v7.3')
save('AvData.mat','avData','-v7.3')

%% Delineation of ROIs

close all
nslices=nslicesEfec;

%initialization of masks for the right and left hemispheres
mmr=zeros(size(avData,1),size(avData,2),(nslices),totalFolder,numROI);
mml=zeros(size(avData,1),size(avData,2),(nslices),totalFolder,numROI);


for side=1:2
    for roi=1:numROI
        for s=1:nslices
            clear pos var
            
            for f=1:totalFolder
                figure(1)
                display([ROIstr{roi} ' Side ' num2str(side)])
                close all
                imagesc(squeeze(mean(squeeze(avData(:,:,s,:,f)),3))); colormap(gray); axis off;
                
                if ~exist('pos','var')       % Check if position was loaded
                    headerMask = impoly;     % Draw ROI: h = impoly begins interactive placement of a polygon on the current axes. The function returns h, a handle to an impoly object.
                else
                    headerMask = impoly(gca, pos); shg; % Use existing ROI for the next folder
                end
                
                wait(headerMask);
                pos = headerMask.getPosition; % Saves points
                
                if side==1
                    mmr(:,:,s,f,roi) = headerMask.createMask; % Creates the mask
                else
                    mml(:,:,s,f,roi) = headerMask.createMask; % Creates the mask
                end
                
            end
        end
    end
    
    %Saving masks
    if side==1
        save('Rmask.mat','mmr');
    else
        save('Lmask.mat','mml');
    end
    
end

%Joining masks from both hemispheres
for roi=1:numROI
    mmfinal(:,:,:,:,roi)=mmr(:,:,:,:,roi)+mml(:,:,:,:,roi);
end


%%Correcting final mask for occurence of 2 (where side masks were both 1)
mmfinalcorr=zeros(size(mmfinal));

xdim=size(mmfinal,1);
ydim=size(mmfinal,2);

for roi=1:numROI
    for s=1:nslices
        for f=1:totalFolder
            temp=reshape(squeeze(mmfinal(:,:,s,f,roi)),[xdim*ydim,1]);
            for index=1:length(temp)
                if temp(index)==2
                    display('2 in mask occurred')
                    temp(index)==1
                end
            end
            mmfinalcorr(:,:,s,f,roi)=reshape(temp,[xdim,ydim]);
        end
    end
end

%Saving mask data
save('Fmask.mat','mmfinal')
save('FmaskCorrect.mat','mmfinalcorr','-v7.3')

%% Apply designed masks
clc
close all

%initialize masked data vector
Datamasked=zeros(size(avData,1),size(avData,2),size(avData,3),size(avData,4),totalFolder,numROI);

for f=1:totalFolder
    f
    for s=1:nslices
        for r=1:nrep
            for roi=1:numROI
                Datamasked(:,:,s,r,f,roi)=avData(:,:,s,r,f).*mmfinalcorr(:,:,s,f,roi);
            end
            
        end
    end
end

%saving data
save('DataMasked.mat','Datamasked','-v7.3')

%Selecting which slices to use for the each ROI (in this example the total used number of slices is 4)
clear temp
avDatamasked=zeros(nrep,totalFolder,numROI);

for roi=1:numROI
    
    for f=1:totalFolder
        f
        for r=1:nrep
            a=nonzeros(squeeze(Datamasked(:,:,1,r,f,roi))); %slice 1
            b=nonzeros(squeeze(Datamasked(:,:,2,r,f,roi))); %slice 2
            c=nonzeros(squeeze(Datamasked(:,:,3,r,f,roi))); %slice 3
            d=nonzeros(squeeze(Datamasked(:,:,4,r,f,roi))); %slice 4
            if roi==1
                clear temp
                temp(:,r,f,roi)=[a',b',c']; %First ROI uses slices 1,2 and 3
            elseif roi==2
                clear temp
                temp(:,r,f,roi)=[a',b',c'];
            elseif roi==3
                clear temp
                temp(:,r,f,roi)=[c',d'];
            end
            avDatamasked(r,f,roi)=mean(nonzeros(temp(:,r,f,roi)));
        end
    end
end

%Saving data
save('DataMasked.mat','avDatamasked','-append')

%% Detrending individual runs with a polynomial fit to the resting period

close all
xaxis=[1:nrep];

Perc_avDatamasked=zeros(size(avDatamasked));
for roi=1:numROI
    display(['ROI ' num2str(ROIstr{roi})])
    clear y3
    
    for f=1:totalFolder
        
        y3(:,f,roi) = smooth(avDatamasked(:,f,roi),10); %Smooth run for an easier fit of the polinimial
        p = polyfit(find(event==0)',y3(find(event==0),f,roi),5); %fit of a 5th degree polynomial only to resting periods
        y4(:,f,roi) = polyval(p(:),xaxis);
        
        % Calculation of percentage signal change
        Perc_avDatamasked(:,f,roi)=((avDatamasked(:,f,roi)-y4(:,f,roi))./y4(:,f,roi))*100;
        
    end
    
end

%Saving data
save('DataDetrend.mat','Perc_avDatamasked','-v7.3')

% Separate frequency conditions
Perc_avDatamasked1=Perc_avDatamasked(:,1:length(folder1),:);
Perc_avDatamasked2=Perc_avDatamasked(:,length(folder1)+1:length(folder1)+length(folder2),:);
Perc_avDatamasked3=Perc_avDatamasked(:,length(folder1)+length(folder2)+1:length(folder1)+length(folder2)+length(folder3),:);

%% Average run

for roi=1:numROI
    
    %average
    meanPerc_avDatamasked1(:,roi)=mean(squeeze(Perc_avDatamasked1(:,:,roi)),2);
    meanPerc_avDatamasked2(:,roi)=mean(squeeze(Perc_avDatamasked2(:,:,roi)),2);
    meanPerc_avDatamasked3(:,roi)=mean(squeeze(Perc_avDatamasked3(:,:,roi)),2);
    
    %standard error of the mean
    semPerc_avDatamasked1(:,roi)=std(squeeze(Perc_avDatamasked1(:,:,roi)),0,2)./sqrt(length(nfolder1));
    semPerc_avDatamasked2(:,roi)=std(squeeze(Perc_avDatamasked2(:,:,roi)),0,2)./sqrt(length(nfolder2));
    semPerc_avDatamasked3(:,roi)=std(squeeze(Perc_avDatamasked3(:,:,roi)),0,2)./sqrt(length(nfolder3));
    
end

GmeanPerc_avDatamasked(:,:,1)=meanPerc_avDatamasked1;
GmeanPerc_avDatamasked(:,:,2)=meanPerc_avDatamasked2;
GmeanPerc_avDatamasked(:,:,3)=meanPerc_avDatamasked3;

GsemPerc_avDatamasked(:,:,1)=semPerc_avDatamasked1;
GsemPerc_avDatamasked(:,:,2)=semPerc_avDatamasked2;
GsemPerc_avDatamasked(:,:,3)=semPerc_avDatamasked3;


close all
xvec=[1:nrep]; %x-axis in repetitions
xveccorr=xvec.*TR; %x-axis in sec

close all
NConditions=size(GsemPerc_avDatamasked,3); %number of conditions

%RGB colour code
r=[255,15,229,96,127];
g=[211,82,83,168,0];
b=[25,186,0,48,255];

i=1;
for cond=1:NConditions
    
    for roi=1:3
        
        limvec=[-2 2.5]; %manual defined y-axis
        
        figure(cond)
        a=shadedErrorBar(xveccorr,GmeanPerc_avDatamasked(:,roi,cond),GsemPerc_avDatamasked(:,roi,cond),'lineprops', {'color', [r(roi), g(roi), b(roi)]/255}); %plot data
        hold on
        plot(xveccorr,event*2) %plot stimulation paradigm on top of data
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        
        ylim(limvec)
        xlim([0 max(xveccorr(:))])
        
        %legend
        if cond==1
            if roi==1
                a.mainLine.DisplayName = 'V1';
                [hleg, hobj, hout, mout] = legend(findobj(gca, '-regexp', 'DisplayName', '[^'']'));
                hold on
                
            elseif roi==2
                
                a.mainLine.DisplayName = 'SC';
                [hleg, hobj, hout, mout] = legend(findobj(gca, '-regexp', 'DisplayName', '[^'']'));
            else
                
                a.mainLine.DisplayName = 'LGN';
                [hleg, hobj, hout, mout] = legend(findobj(gca, '-regexp', 'DisplayName', '[^'']'));
                
            end
        end
    end
    %title
    title(stringcond{cond})
    
    i=i+1;
end

% Saving plots
fig1=figure(1);
fig1.Renderer='Painters';
saveas(gcf,'NameXX','epsc')
saveas(gcf,'NameXX','fig')
saveas(gcf,'NameXX','png')
%% Cycle separation centered around the stimulation period

Cycle_Perc_avDatamasked1=zeros(lcycle,ncycle*length(folder1),numROI);
Cycle_Perc_avDatamasked2=zeros(lcycle,ncycle*length(folder2),numROI);
Cycle_Perc_avDatamasked3=zeros(lcycle,ncycle*length(folder3),numROI);

display(['condition ' num2str(1)])
for roi=1:numROI
    b=1;
    for f=1:length(folder1)
        m=0;
        for c=1:ncycle
            b
            display(['start ',num2str((lcycle)*m+start-rest/2),' finish ',num2str((lcycle)*m+start+lcycle-1-rest/2),]) %check if this correponds to your paradigm
            
            Cycle_Perc_avDatamasked1(:,b,roi)=Perc_avDatamasked1((lcycle)*m+start-rest/2:(lcycle)*m+start+lcycle-1-rest/2,f,roi);
            b=b+1;
            m=m+1;
        end
    end
    
end

display(['condition ' num2str(2)])
for roi=1:numROI
    b=1;
    for f=1:length(folder2)
        m=0;
        for c=1:ncycle
            b
            display(['start ',num2str((lcycle)*m+start-rest/2),' finish ',num2str((lcycle)*m+start+lcycle-1-rest/2),])
            
            Cycle_Perc_avDatamasked2(:,b,roi)=Perc_avDatamasked2((lcycle)*m+start-rest/2:(lcycle)*m+start+lcycle-1-rest/2,f,roi);
            b=b+1;
            m=m+1;
        end
    end
    
end

display(['condition ' num2str(3)])
for roi=1:numROI
    b=1;
    for f=1:length(folder3)
        m=0;
        for c=1:ncycle
            b
            display(['start ',num2str((lcycle)*m+start-rest/2),' finish ',num2str((lcycle)*m+start+lcycle-1-rest/2),])
            
            Cycle_Perc_avDatamasked3(:,b,roi)=Perc_avDatamasked3((lcycle)*m+start-rest/2:(lcycle)*m+start+lcycle-1-rest/2,f,roi);
            b=b+1;
            m=m+1;
        end
    end
    
end

%Average cycle
for roi=1:numROI
    roi
    %average
    AvCycle_Perc_avDatamasked1(:,roi)=mean(squeeze(Cycle_Perc_avDatamasked1(:,:,roi)),2);
    AvCycle_Perc_avDatamasked2(:,roi)=mean(squeeze(Cycle_Perc_avDatamasked2(:,:,roi)),2);
    AvCycle_Perc_avDatamasked3(:,roi)=mean(squeeze(Cycle_Perc_avDatamasked3(:,:,roi)),2);
    
    %standard error of the mean
    SemCycle_Perc_avDatamasked1(:,roi)=std(squeeze(Cycle_Perc_avDatamasked1(:,:,roi)),0,2)./sqrt(length(nfolder1));
    SemCycle_Perc_avDatamasked2(:,roi)=std(squeeze(Cycle_Perc_avDatamasked2(:,:,roi)),0,2)./sqrt(length(nfolder2));
    SemCycle_Perc_avDatamasked3(:,roi)=std(squeeze(Cycle_Perc_avDatamasked3(:,:,roi)),0,2)./sqrt(length(nfolder3));
end

AvCycle_Perc_avDatamasked(:,:,1)=AvCycle_Perc_avDatamasked1;
AvCycle_Perc_avDatamasked(:,:,2)=AvCycle_Perc_avDatamasked2;
AvCycle_Perc_avDatamasked(:,:,3)=AvCycle_Perc_avDatamasked3;

SemCycle_Perc_avDatamasked(:,:,1)=SemCycle_Perc_avDatamasked1;
SemCycle_Perc_avDatamasked(:,:,2)=SemCycle_Perc_avDatamasked2;
SemCycle_Perc_avDatamasked(:,:,3)=SemCycle_Perc_avDatamasked3;

%save data
save('DataDetrend.mat','Cycle_Perc_avDatamasked1','AvCycle_Perc_avDatamasked1','SemCycle_Perc_avDatamasked1',...
    'Cycle_Perc_avDatamasked2','AvCycle_Perc_avDatamasked2','SemCycle_Perc_avDatamasked2',...
    'Cycle_Perc_avDatamasked3','AvCycle_Perc_avDatamasked3','SemCycle_Perc_avDatamasked3',...
    'AvCycle_Perc_avDatamasked','SemCycle_Perc_avDatamasked','-append')
%% Plotting results

close all
xvec=[1:lcycle]; %x-axis vector in repetitions
xveccorrtmp=xvec.*TR; %x-axis vector in sec

xveccorr=[-restTime/2:TR:restTime/2+stimTime-TR];
xveccorr = round(xveccorr,2);

NConditions=size(AvCycle_Perc_avDatamasked,3);

%RGB colour code
r=[255,15,229,96,127];
g=[211,82,83,168,0];
b=[25,186,0,48,255];

i=1;
for cond=1:NConditions
    
    for roi=1:numROI
        
        limvec=[-2 2.5];
        
        figure(1)
        subplot(1,NConditions,i)
        
        a=shadedErrorBar(xveccorr,AvCycle_Perc_avDatamasked(:,roi,cond),SemCycle_Perc_avDatamasked(:,roi,cond),'lineprops', {'color', [r(roi), g(roi), b(roi)]/255}); %plot data
        hold on
        line([0 0],limvec,'Color','b','LineWidth',2) %line 1 where stimulation starts
        line([stimTime-TR  stimTime-TR],limvec,'Color','b','LineWidth',2) %line 2 where stimulation ends
        
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        ylim(limvec)
        xlim([xveccorr(1) xveccorr(end)])
        
        %legend
        if cond==1
            
            if roi==1
                a.mainLine.DisplayName = 'V1';
                [hleg, hobj, hout, mout] = legend(findobj(gca, '-regexp', 'DisplayName', '[^'']'));
                hold on
                
            elseif roi==2
                
                a.mainLine.DisplayName = 'SC';
                [hleg, hobj, hout, mout] = legend(findobj(gca, '-regexp', 'DisplayName', '[^'']'));
                
            else
                
                a.mainLine.DisplayName = 'LGN';
                [hleg, hobj, hout, mout] = legend(findobj(gca, '-regexp', 'DisplayName', '[^'']'));
                
            end
            
        end
    end
    
    %title
    title(stringcond{cond})
    
    i=i+1;
end

% Saving plots
fig1=figure(1);
fig1.Renderer='Painters';
saveas(gcf,'NameXX','epsc')
saveas(gcf,'NameXX','fig')
saveas(gcf,'NameXX','png')
%%
%%
%% Quantification (run this code for each ROI)
close all
clc

%Join all frequency conditions (example where total = 2 frequency groups)
load('path') %load data for each condition
Cycle_Perc_avDatamasked1=tCycle_Perc_avDatamasked1;
Cycle_Perc_avDatamasked2=tCycle_Perc_avDatamasked2;

numcond=2;%total number of conditions
max_numbruns=X; %maximum number of runs among all conditions
roi=1; %define roi being studied (1-VC, 2-SC and 3-LGN)

%initialize vector where data is going to be stored
Perc_avDatamasked1_runs=nan(size(Cycle_Perc_avDatamasked1,1),max_numbruns,numcond);

%Frequency condition 1
i=1;
for j=1:size(Cycle_Perc_avDatamasked1,2)/6 %/6 comes from the fact that each run has 6 cycles
    Perc_avDatamasked1_runs(:,j,1)=mean(Cycle_Perc_avDatamasked1(:,i:i+5,roi),2); % every cycle is the averaged cycles of one run
    i=i+6;
end

%Frequency condition 2
i=1;
for j=1:size(Cycle_Perc_avDatamasked2,2)/6
    Perc_avDatamasked1_runs(:,j,2)=mean(Cycle_Perc_avDatamasked2(:,i:i+5,roi),2);
    i=i+6;
end

%Saving data
save('ROIX_Data.mat','Perc_avDatamasked1_runs')

%% Calculate the steady-state region for each run

%Interpolation of averaged run
AvPerc_avDatamasked1_runs=squeeze(nanmean(Perc_avDatamasked1_runs,2)); %averaged of all runs
fs=4; %new samling frequency

%Old time vector
tx=[TR:TR:60];
tx=round(tx,2);

%New time vector
vecnew=[1/fs:1/fs:60];
vecnew=round(vecnew,2);

%New stimulation vector
stimvecResampled=zeros(size(vecnew));
stimvecResampled(find(vecnew==22.5+1/fs):find(vecnew==37.5+1/fs))=1; %interval where stimulation is ON (1s)
stimvecResampled=round(stimvecResampled,2);

%Resamling averaged data
ResampledData=[];
for freq=1:numcond
    ResampledData(:,freq) = resample(squeeze(AvPerc_avDatamasked1_runs(:,freq)),tx,fs,'spline');
end

%Go point-by-point and compare means from the 4 points around it

vecsst=[];
flagenter=1;
temp=[];
indvec=find(stimvecResampled==1);

for freq=1:numcond
    
    indmax=find(ResampledData(:,freq)==max(ResampledData(find(stimvecResampled==1),freq)));
    
    for i=indmax+4:indvec(end-4)
        
        tempB=nanmedian(ResampledData(i-4:i,freq));
        tempA=nanmedian(ResampledData(i:i+4,freq));
        
        if abs(tempB-tempA)<0.04 && flagenter==1
            vecsst(freq)=i-4;
            flagenter=2;
        end
    end
    
    flagenter=1;
end

%Saving data
save('ROIX_SST.mat','vecsst','indvec')

indvec=find(stimvecResampled==1);

%Resampling each the averaged cycles corresponding to each run separately
ResampledDataRUN=zeros(length(vecnew),size(Perc_avDatamasked1_runs,2),size(Perc_avDatamasked1_runs,3));
for f=1:size(Perc_avDatamasked1_runs,2)
    
    for cond=1:numcond
        temp=squeeze(Perc_avDatamasked1_runs(:,f,cond));
        
        if ~isnan(mean(temp))
            ResampledDataRUN(:,f,cond) = resample(temp,tx,fs,'spline');
        else
            ResampledDataRUN(:,f,cond) =nan;
        end
        
    end
    
end

%Saving Data
save('ROIX_SST_timepoints.mat','ResampledData','ResampledDataRUN','vecsst','stimvecResampled','vecnew','tx','fs','indvec','-v7.3')

%Calculate steady-state value per run
for f=1:size(ResampledDataRUN,2)
    
    for cond=1:numcond
        ResampledDataRUNCorr(:,f,cond)=ResampledDataRUN(:,f,cond)-nanmean(ResampledDataRUN(1:vecsst(cond)-1,f,cond)); %extra correction using the baseline prior to stimulation and suctract it to the calculated steady-state value
        SST_minusinitialPerc_avDatamasked1_runs(f,cond)=nanmean(ResampledDataRUNCorr(vecsst(cond):indvec(end),f,cond),1);
    end
    
end

%Steady-state interval Box plot
close all
figure('units','normalized','outerposition',[0 0 1 1])
boxplot(SST_minusinitialPerc_avDatamasked1_runs,'Labels',{'2Hz','20Hz','cont'})
hold on
notBoxPlot(SST_minusinitialPerc_avDatamasked1_runs,'style','line')
ylim([-2.5 1.8])
fig1=figure(1);
fig1.Renderer='Painters';
saveas(gcf,'NameXX','epsc')
saveas(gcf,'NameXX','fig')
saveas(gcf,'NameXX','png')

%save data
save('ROIX_SST_Boxplots.mat','SST_minusinitialPerc_avDatamasked1_runs',...
    'SST_Perc_avDatamasked1_runs','AREA_Perc_avDatamasked1_runs','-v7.3')

%%
%%
%% fMRI neurometric curve (run this code for each ROI)

cd('path') %path where quantified data is stored
load('workspace.mat')
load('ROIX_SST.mat')

%Roi we are looking at (1-VC,2-SC,3-LGN)
roi=1;

%Join different conditions (taking each individual cycle)
AvPerc_avDatamasked1_trials=nan(lcycle,maxcycles,numcond);
AvPerc_avDatamasked1_trials(:,1:length(Cycle_Perc_avDatamasked1),1)=Cycle_Perc_avDatamasked1(:,:,roi);
AvPerc_avDatamasked1_trials(:,1:length(Cycle_Perc_avDatamasked2),2)=Cycle_Perc_avDatamasked2(:,:,roi);

%Interpolation
fs=4; %new samling frequency 4HZ

%Old time vector
tx=[TR:TR:60];
tx=round(tx,2);

%New interpolated time vector
vecnew=[1/fs:1/fs:60];
vecnew=round(vecnew,2);

%Interpolated stimvec vector
stimvecResampled=zeros(size(vecnew));
stimvecResampled(find(vecnew==22.5+1/fs):find(vecnew==37.5+1/fs))=1; %stimulation period (1s)

stimvecResampled=round(stimvecResampled,2);

%Resamling each trial
ResampledDataTrial=nan(lcycle*6,maxcycles,numcond);
for c=1:maxcycles
    for freq=1:numcond
        if ~isnan(mean(AvPerc_avDatamasked1_trials(:,c,freq),1))
            ResampledDataTrial(:,c,freq) = resample(squeeze(AvPerc_avDatamasked1_trials(:,c,freq)),tx,fs,'spline');
        end
    end
end

%Calculate steady-state value per run
ResampledDataTrialCorr=nan(lcycle*6,maxcycles,numcond);
for c=1:maxcycles
    for freq=1:numcond
        
        ResampledDataTrialCorr(:,c,freq)=ResampledDataTrial(:,c,freq)-nanmean(ResampledDataTrial(1:vecsst(freq)-1,c,freq)); %extra correction using the baseline prior to stimulation and suctract it to the calculated steady-state value
        SST_minusinitialPerc_avDatamasked1_trial(c,freq)=nanmean(ResampledDataTrialCorr(vecsst(freq):indvec(end),c,freq),1);
    end
end

%save data
save('ROIX_SST_trialData.mat','ResampledDataTrialCorr','SST_minusinitialPerc_avDatamasked1_trial',...
    'SST_Perc_avDatamasked1_trial','AREA_Perc_avDatamasked1_trial')

reshplotT=SST_minusinitialPerc_avDatamasked1_trial;

%Removal of outliers (that deviate from the mean more that 3 standard deviations)
reshplotTnO=reshplotT;
c=1;
for freq=1:numcond
    
    meanval(freq)=nanmean(reshplotT(:,freq));
    deviation(freq)=nanstd(reshplotT(:,freq),0,1);
    
    for c=1:maxcycles
        if reshplotT(c,freq)>meanval(freq)+3*deviation(freq) || reshplotT(c,freq)<meanval(freq)-3*deviation(freq)
            display(['outlier corrected in Freq=' num2str(freq) ])
            reshplotTnO(c,freq)= nan;
            c=c+1;
        end
        
    end
end

%Save data
save('DataX_trialData.mat','reshplotT','reshplotTnO','-append')

%% Define thresholds for "continuous" or "pulsed" categorization

%initialize vectors where data is going to be put
DecisionTrial=nan(maxcycles,numcond);
PercPulsedLight_Trial=[];

condref=numcond-2; %take 20Hz as reference based on behavior results

%theshold is the mean value of the reference condition
thresT=round(nanmean(reshplotTnO(:,condref)));

%categorization of each trial in to "pulsating" or "continuous"
for freq=1:numcond
    for c=1:maxcycles
        
        if reshplotTnO(c,freq)>thresT && ~isnan(reshplotTnO(c,freq))
            DecisionTrial(c,freq)=1;
        elseif reshplotTnO(c,freq)<thresT && ~isnan(reshplotTnO(c,freq))
            DecisionTrial(c,freq)=2;
        end
        
    end
    
    %Calculate the percentage of pulsating reports
    PercPulsedLight_Trial(freq)=(sum(DecisionTrial(:,freq)==1)*100)./(sum(~isnan(DecisionTrial(:,freq))));
    
end

%Save data
save('DataX_trialData.mat','PercPulsedLight_Trial','-append')
%% Calculation of threshold of 50%
xvecfreq=[0.25,1,2,15,20,25,60]; %vector with tested frequencies in Hz. The continuous light is represneted by the 60Hz

%fitting a psychometric curve to the data (used in Wichmann, Felix A., and N. Jeremy Hill., Perception & psychophysics, 2001)
close all
ymeans=PercPulsedLight_Trial(1:end)/100;
[xData, yData]=prepareCurveData(xvecfreq,ymeans);
q1=fittype(@(a,b,c,d,x) d + (c-d)./(1+exp(-(2*a)*(x-b))),'independent',{'x'},'dependent',{'y'});
[f,~]=fit(xData,squeeze(yData),q1,'Lower', [0 xvecfreq(2) 0 0], 'Upper', [Inf xvecfreq(end) 1 1],'startpoint',[1,0,0.05, 0.95]);
thresh=f.b-(1/(2*f.a))*log(((f.c-f.d)/(0.50-f.d))-1) %Calculation of threshold of 50% choices to pulsating port ('chance level')

%plot results
figure('units','normalized','outerposition',[0 0 1 1])
plot(xvecfreq,ymeans,'-o')
hold on
plot(f)
ylim([0 1])
hold on
notBoxPlot(PercPulsedLight_Trial(:)/100,xvecfreq,'style','line')
legend off

%saving figure
fig1=figure(1);
fig1.Renderer='Painters';
saveas(gcf,'DataX_fMRINeurometric','epsc')
saveas(gcf,'DataX_fMRINeurometric','fig')
saveas(gcf,'DataX_fMRINeurometric','png')


%Boostraping
it=50;%iterations for bootstrapping

d=testPercPulsedLight_Trial(:,1:end);
bootstat  = bootstrp(it,@nanmean,d);

clear xData yData ymeans1 threshBOOT d
for a=1:it
    %fitting a psychometric curve to the data (used in Wichmann, Felix A., and N. Jeremy Hill., Perception & psychophysics, 2001)
    ymeans1(a,:)=bootstat(a,:)/100;
    [xData, yData(a,:)]=prepareCurveData(xvecfreq,ymeans1(a,:));
    q1=fittype(@(a,b,c,d,x) d + (c-d)./(1+exp(-(2*a)*(x-b))),'independent',{'x'},'dependent',{'y'});
    [f,~]=fit(xData,squeeze(yData(a,:))',q1,'Lower', [0 xvecfreq(1) 0 0], 'Upper', [Inf xvecfreq(end) 1 1],'startpoint',[1,0,0.05, 0.95]);
    threshBOOT(a)=f.b-(1/(2*f.a))*log(((f.c-f.d)/(0.5-f.d))-1); %Calculation of threshold of 50% choices to pulsating port ('chance level')
    
end

%Calculatr final values for the specific ROI
finalthreshhigh=mean(thresh);
finalthreshhighSTDBOOT=std(threshBOOT,0,2);

%Save data
save('DataX_trialData.mat','finalthreshhigh','finalthreshhighSTDBOOT','-append')