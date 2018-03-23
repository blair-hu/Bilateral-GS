function [fnames, IMU_data] = prepIMU(sensor_startind,Fs)
% Function: Read bilateral (shank/thigh) IMU data from raw format,
% calculate the estimated thigh/shank orientation angles using a
% complementary filter, and save output as .mat

% Input: None (to be selected using GUI)
% Output: filename_processed.mat

% Function dependencies:
% uipickfiles.m
% LPfilt.m

%%%%%
% Documented by: Blair Hu 03/21/18
%%%%%

if nargin == 0
    sensor_startind = [1 7 13 19];
    Fs = 500;
end

files = uipickfiles;
for fileind = 1:length(files)        
    [~,fname,~] = fileparts(files{fileind});
    filevar = load(fname);
    
    frinc = eval(['filevar.',fname,'.setup.DAQ_FRINC;']);
    gtruth = eval(['filevar.',fname,'.pvd.MODE;']);
    
    gtruth_daq = [];
    for i = 1:length(gtruth)
        gtruth_daq = [gtruth_daq; repmat(gtruth(i),frinc,1)];
    end
    
    % Put IMU_data in cells based on number of sensors and apply scaling
    for i = 1:length(sensor_startind)
        % Acc scaling: (UINT - 32768)/8192
        % Gyro scaling: (UINT - 32768)/65.536
        Accel = eval(['(double(filevar.',fname,'.daq.daqUINT16(:,sensor_startind(i):sensor_startind(i)+2)) - 32768)/8192;']);
        Gyro = eval(['(double(filevar.',fname,'.daq.daqUINT16(:,sensor_startind(i)+3:sensor_startind(i)+5)) - 32768)/65.536;']);
        
        % Low-pass filter IMU signals by applying 6th order Butterworth at
        % 25 Hz
        IMU_filt = LPfilt(Fs,6,25,[Accel Gyro]);      
        % Remove offset by zero-centering
        sag_vel = IMU_filt(:,4) - mean(IMU_filt(:,4)); % The fourth column refers to sagittal plane angular velocity
        % Find and threshold peaks in shank velocity
        pospeaks = findpeaks(sag_vel);
        negpeaks = findpeaks(-sag_vel);
        pospeaks = pospeaks(pospeaks > 100);
        negpeaks = negpeaks(negpeaks > 100);
        % Determine if the shank velocity signal is inverted; if so, correct it
        if mean(negpeaks) > mean(pospeaks)
            sag_vel = -sag_vel;
        end
        IMU_filt(:,4) = sag_vel;        
        IMU_data{i} = IMU_filt;      
    end
    
    % Estimate orientation angles using complementary filter
    alpha = 0.99;
    NetAcc_low = 0.99;
    NetAcc_high = 1.01;
    for i = 1:length(sensor_startind)
        Ax = IMU_data{i}(:,1);
        Ay = IMU_data{i}(:,2);
        Az = IMU_data{i}(:,3);
        Gy = IMU_data{i}(:,4);
        Gz = IMU_data{i}(:,5);
        Gx = IMU_data{i}(:,6);
        NetAcc = sqrt(Ax.^2 + Az.^2);
        pitch = atan2d(Az,Ax.^2+Ay.^2);
        % Define initial pose
        unfilt_angle(1) = pitch(1);
        filt_angle(1) = pitch(1);
        for j = 2:size(Gy,1)
            % Apply complementary filter when accel. reading is in normal range
            unfilt_angle(j) = unfilt_angle(j-1) + Gy(j)*1/Fs;
            if NetAcc(j) > NetAcc_low && NetAcc(j) < NetAcc_high
                filt_angle(j) = alpha*(filt_angle(j-1) + Gy(j)*1/Fs) + (1-alpha)*pitch(j);
            else
                filt_angle(j) = filt_angle(j-1) + Gy(j)*1/Fs;
            end
        end
        IMU_data{i} = [IMU_data{i} pitch unfilt_angle' filt_angle'];
        pitch = [];
        unfilt_angle = [];
        filt_angle = [];
    end
    
    % Identify HC, TO, and MSw gait events from shank IMU angular velocity
    % in sagittal plane
    % Parameters (somewhat arbitrary, hand-tuned) for finding peaks
    MSw_Threshold = 90;
    MSw_Width = 90;
    MSw_Separation = 400;    
    for i = 3:4 % Limit to shank IMUs (left = 3, right = 4)
        trigs_daq = zeros(length(gtruth_daq),1);
        % Get shank velocity for gait segmentation
        shank_vel = IMU_data{i}(:,4);
        % Get gait events using shank velocity
        [HC_ind,TO_ind,MSw_ind,HC_trig,TO_trig,MSw_trig] = getgaitevents(gtruth_daq,shank_vel,MSw_Threshold, MSw_Width, MSw_Separation);
        inds = [HC_ind; TO_ind; MSw_ind];
        trigs = [HC_trig; TO_trig; MSw_trig];
        for r = 1:length(inds)
            trigs_daq(inds(r)) = str2num(trigs{r});
        end
        IMU_data{i+2} = trigs_daq;
    end
    
    IMU_data = segmentgait(IMU_data);
    IMU_data_mat = cell2mat(IMU_data);
    
    disp(['Processed: ',fname])
%     save([fname,'_processed'],'IMU_data_mat');
end
end

function filt_data = LPfilt(Fs,N,Fc,data)
%%%%%%%%%%%%%%%%%%%%%%
%This function high-passfilters the EMG data to reduce the motion artifact
%Fs - sampling frequency
%N - Filter order
%Fc - cutoff frequency
%data - data to be filtered
%%%%%%%%%%%%%%%%%%%%%%%
[B,A] = butter(N, Fc/(Fs/2), 'low');
for i=1:size(data,2)
    filt_data(:,i) = filtfilt(B,A,data(:,i));
end
end

function [HC,TO,MSw,HC_trig,TO_trig,MSw_trig] = getgaitevents(MODE_DAQ,Shank_Vel,MSw_Threshold,MSw_Width,MSw_Separation)
% First identify all mid-swing peaks (largest positive velocity)
[~,MSw] = findpeaks(Shank_Vel,'MinPeakHeight',MSw_Threshold,'MinPeakWidth',MSw_Width,'MinPeakDistance',MSw_Separation);
MSw = MSw((MODE_DAQ(MSw) > 0));
MSw = MSw(and(MSw > 500,MSw < length(Shank_Vel) - 500));

HC = [];
TO = [];
for i = 1:length(MSw)
    currMSw = MSw(i);
    % Look for first large peak before each MSw event (represents toe off)
    counter = 0;
    TO_guess = [];
    while Shank_Vel(currMSw) > Shank_Vel(currMSw-1)
        currMSw = currMSw - 1;
        % Check if there is a turn in the monotonicity
        if Shank_Vel(currMSw) < Shank_Vel(currMSw-1)
            TO_guess = [TO_guess; currMSw];
            while Shank_Vel(currMSw) < Shank_Vel(currMSw-1)
                currMSw = currMSw - 1;
                counter = counter + 1;
                % Allow turns of up to 50 ms
                if counter > 50
                    currMSw = currMSw + counter;
                    counter = 0;
                    break
                end
            end
        end
    end
    Shank_Vel_TO_guess = Shank_Vel(TO_guess);
    [~,Min_guess] = min(Shank_Vel_TO_guess);
    TO = [TO; TO_guess(Min_guess)];
    % If there is another MSw event before, look for the first zero crossing 
    % after the previous MSw event (represents heel contact)
    if i > 1
        prevMSw = MSw(i-1);
        tempstart = prevMSw;
        while Shank_Vel(tempstart) > 0
            tempstart = tempstart + 1;
        end
        HC = [HC; tempstart];
        % If this is the last MSw event, look for a heel contact event
        % afterwards
        if i == length(MSw)
            lastMSw = MSw(i);
            tempstart = lastMSw;
            while Shank_Vel(tempstart) > 0
                tempstart = tempstart + 1;
            end
            HC = [HC; tempstart];
        end
    end
end

% Eliminate gait events that are too close to the beginning/end of the
% trial
HC = HC(HC > 301);
HC = HC(HC < length(MODE_DAQ)-50);
HC = intersect(HC,find(MODE_DAQ > 0));

TO = TO(TO > 301);
TO = TO(TO < length(MODE_DAQ)-50);
TO = intersect(TO,find(MODE_DAQ > 0));

MSw = MSw(MSw > 301);
MSw = MSw(MSw < length(MODE_DAQ)-50);
MSw = intersect(MSw,find(MODE_DAQ > 0));

% Determine the appropriate triggers for each gait event
if length(HC) > 0
    for i = 1:length(HC)
        prevTO = TO(max(find(TO < HC(i))));
        if ~isempty(prevTO)
            if (HC(i) - prevTO < 1000)
                HC_trig{i,1} = [num2str(MODE_DAQ(prevTO+200)),'3',num2str(MODE_DAQ(HC(i)+200)),'1'];
            else
                HC_trig{i,1} = [num2str(MODE_DAQ(HC(i)-500)),'3',num2str(MODE_DAQ(HC(i)+200)),'1'];
            end
        else % First HC after standing
            HC_trig{i,1} = [num2str(MODE_DAQ(HC(i)-300)),'3',num2str(MODE_DAQ(HC(i)+200)),'1'];
        end
    end
else
    HC_trig = {};
end

if length(TO) > 0
    for i = 1:length(TO)
        prevHC = HC(max(find(HC < TO(i))));
        if ~isempty(prevHC)
            if (TO(i) - prevHC < 1000)
                TO_trig{i,1} = [num2str(MODE_DAQ(prevHC+200)),'1',num2str(MODE_DAQ(TO(i)+200)),'2'];
            else
                TO_trig{i,1} = [num2str(MODE_DAQ(TO(i)-500)),'1',num2str(MODE_DAQ(TO(i)+200)),'2'];
            end
        else % First HC after standing
            TO_trig{i,1} = [num2str(MODE_DAQ(TO(i)-300)),'1',num2str(MODE_DAQ(TO(i)+200)),'2'];
        end
    end
else
    TO_trig = {};
end

if length(MSw) > 0
    for i = 1:length(MSw)
        MSw_trig{i,1} = [num2str(MODE_DAQ(MSw(i)-300)),'2',num2str(MODE_DAQ(MSw(i)+200)),'3'];
    end
else
    MSw_trig = {};
end
end

function IMU_data = segmentgait(IMU_data)
% _3_1 = HC
% _1_2 = TO
% _2_3 = MSW
% Identify bilateral stance (1) and swing (0) from triggers

    L_trig = IMU_data{5};
    R_trig = IMU_data{6};
    
    R_seg = ones(length(R_trig),1);
    L_seg = ones(length(L_trig),1);
    
    R_trig_ind = find(R_trig > 0);
    L_trig_ind = find(L_trig > 0);
    
    R_trig_nonzero = R_trig(R_trig_ind);
    L_trig_nonzero = L_trig(L_trig_ind);
    
    R_trig_str = num2str(R_trig_nonzero);
    L_trig_str = num2str(L_trig_nonzero);
    
    R_HC = [];
    R_TO = [];
    for i = 1:length(R_trig_str)
        try
            if strcmp(R_trig_str(i,2),'3') && strcmp(R_trig_str(i,4),'1')
                R_HC = [R_HC; R_trig_ind(i)];
            end
        end
        try
            if strcmp(R_trig_str(i,2),'1') && strcmp(R_trig_str(i,4),'2')
                R_TO = [R_TO; R_trig_ind(i)];
            end
        end
    end
    
    L_HC = [];
    L_TO = [];
    for i = 1:length(L_trig_str)
        try
            if strcmp(L_trig_str(i,2),'3') && strcmp(L_trig_str(i,4),'1')
                L_HC = [L_HC; L_trig_ind(i)];
            end
        end
        try
            if strcmp(L_trig_str(i,2),'1') && strcmp(L_trig_str(i,4),'2')
                L_TO = [L_TO; L_trig_ind(i)];
            end
        end
    end
    
    for j = 1:length(R_TO)
        try
            R_seg(R_TO(j):min(R_HC(R_HC > R_TO(j)))) = 0;
        catch
            R_seg(R_TO(j):min(L_TO(L_TO > R_TO(j)))) = 0;
        end
    end
    
    for j = 1:length(L_TO)
        try
            L_seg(L_TO(j):min(L_HC(L_HC > L_TO(j)))) = 0;
        catch
            L_seg(L_TO(j):min(R_TO(R_TO > L_TO(j)))) = 0;
        end
    end
    
    IMU_data{7} = L_seg;
    IMU_data{8} = R_seg;
end