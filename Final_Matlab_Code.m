clear; clc; close all;

%% Reading in Data from Niryo Excel and UltraSound
dataNiryo = readmatrix('GD_90.2.xlsx');
dataUltra = readmatrix('UGreg 90.xlsx');

%% File name you want to output torsion, cuvrature, and area too
% The image numbers that mark the corresponding locations of the label
% *** Make sure the size of makers = size of labels ***
file_Name = 'Test.xlsx';
markers = [29 49 64 77 96];
Labels = {"Artery", "Muscle", "Joint Line", "End CT", "Start CT"};

marker_Time = zeros(1,length(markers));

%%Finding Length of Data and find how much needs deleted
lNiryo = length(dataNiryo);
lUltra = length(dataUltra);

if lNiryo > lUltra
    l = lUltra;
else
    l = lNiryo;
end

%%Create a vector for every variable
x_Niryo = dataNiryo(:,3);
y_Niryo = dataNiryo(:,4);
z_Niryo = dataNiryo(:,5);
t_Niryo = dataNiryo(:,2);
R = dataNiryo(:,6);
P = dataNiryo(:,7);
Yaw = dataNiryo(:,8);
x_Ultra = dataUltra(2:length(dataUltra),5); % y in ultrasound data is x in robot quadrants.  This is why column 4 of ultrasound data is written into x here.
z_Ultra = dataUltra(2:length(dataUltra),4); % x in ultrasound data is y in robot quadrants.  This is why column 3 of ultrasound data is written into y here.
t_Ultra = dataUltra(2:length(dataUltra),6);
measurment_Area_Ultra = dataUltra(2:length(dataUltra),7);  % Area measured by image J trace.           This is column 7 of the ultrasound data. 
minor_Area_Ultra = dataUltra(2:length(dataUltra),8); % Area calculated using the minor dimentsion of the fitted ellipse.     This is column 8 of the ultrasound data.

for i = 1:length(markers)   % Gets the ultrasound time of the markers
    for ii = 1:length(t_Ultra)
        if ii == markers(i)
            marker_Time(i) = t_Ultra(ii);
        end
    end
end

%% CJ Interpolation code for x values in Ultra sound
% This part of the code is here because the robot takes more data points
% then the ultrasound.  This portion linearly interpolates what the
% ultrasound values are at the same times as the time stamps from the
% robot.
new_x_Ultra = zeros(1,length(t_Niryo));
new_z_Ultra = zeros(1,length(t_Niryo));
new_minor_Area_Ultra = zeros(1,length(t_Niryo));
for i = 1:length(t_Niryo)
time_1_nyrio = t_Niryo(i);
min_Diff = 100; % Set to a large initial value
compare_real = NaN;
% Finds the closest two times to the nyrio time
for ii = 1:length(t_Ultra)
compare_abs = abs(time_1_nyrio - t_Ultra(ii));
if compare_abs <= min_Diff
min_Diff = compare_abs;
compare_real = time_1_nyrio - t_Ultra(ii);
Time_Position = ii;
end
end
% Get the position in array base on time
if compare_real < 0
Position_1 = Time_Position - 1;
Position_2 = Time_Position;
else
Position_1 = Time_Position;
Position_2 = Time_Position + 1;
end
if(i == 1)
Position_1 = Time_Position;
Position_2 = Time_Position + 1;
end
if(Time_Position == length(t_Ultra))
Position_1 = length(t_Ultra) - 1;
Position_2 = length(t_Ultra);
end

% Linear interpolation for x
new_x_Ultra(i) = ((x_Ultra(Position_2) - x_Ultra(Position_1)) / ...
(t_Ultra(Position_2) - t_Ultra(Position_1))) * (time_1_nyrio - t_Ultra(Position_1)) + x_Ultra(Position_1);
% Linear interpolation for y
new_z_Ultra(i) = ((z_Ultra(Position_2) - z_Ultra(Position_1)) / ...
(t_Ultra(Position_2) - t_Ultra(Position_1))) * (time_1_nyrio - t_Ultra(Position_1)) + z_Ultra(Position_1);
% Linear interpolation of minor area
new_minor_Area_Ultra(i) = ((minor_Area_Ultra(Position_2) - minor_Area_Ultra(Position_1)) / ...
(t_Ultra(Position_2) - t_Ultra(Position_1))) * (time_1_nyrio - t_Ultra(Position_1)) + minor_Area_Ultra(Position_1);
end
%End of CJ code


%% Correcting Pitch
% Define range
lower_bound = -1.571; % -90 Degrees
upper_bound = -1.553; % -89 Degrees

% Initialize tracking flags
subtraction_active = false;
exited_range = false;

for i = 1:length(P)
    in_range = (P(i) >= lower_bound) && (P(i) <= upper_bound);

    if in_range
        if ~subtraction_active && ~exited_range
            subtraction_active = true;  % Start a new subtracting period
            exited_range = false;
        elseif subtraction_active && exited_range
            subtraction_active = false; % Stop subtracting for this cycle
            exited_range = false;
        end
    else
        if subtraction_active
            exited_range = true; % Flag that we exited range
        end
    end

    % Apply subtraction if active
    if subtraction_active
        P(i) = -1.5708 + (-1.5708 - P(i));
    else
        P(i) = P(i);
    end
end



%%Finding where to cut off first portion of data.
%The robot and ultrasound will be started at a position where there will be
%zero pitch.  The contraption will then be moved to where it can track the
%ulner nerve.  During this movement the ultrasound data will read zero in
%the z direction and -0.1651 in the x direction.

c_Beginning = 0;
tolerance = 1e-3;  % tolerance
for n = 1:length(new_x_Ultra)
    if abs(new_x_Ultra(n) - 0.1651) < tolerance
        c_Beginning = n;
    end
end

%%Cutting off the first portion of the data
x_Niryo = x_Niryo(1+c_Beginning:end,1);
y_Niryo = y_Niryo(1+c_Beginning:end,1);
z_Niryo = z_Niryo(1+c_Beginning:end,1);
t_Niryo = t_Niryo(1+c_Beginning:end,1);
R = R(1+c_Beginning:end,1);
P = P(1+c_Beginning:end,1);
Yaw = Yaw(1+c_Beginning:end,1);
new_x_Ultra = new_x_Ultra(1,1+c_Beginning:end);
new_z_Ultra = new_z_Ultra(1,1+c_Beginning:end); 
new_minor_Area_Ultra = new_minor_Area_Ultra(1,1+c_Beginning:end);
new_minor_Area_Ultra = new_minor_Area_Ultra.';

%% Dynamics Code
%%This is a dynamics portion of the code that calculates the location of
%%the ulnar nerve.  The equation this impliments is that the x-y-z postion
%%of the ulnar nerve is equal to the x-y-z position vector of the tip of
%%the ultrasound relative to the robot axes plus the position vector of the
%%ulnar nerve according to the ultrasound times the rotation matrixes of
%%the roll, pitch and yaw angles.
for i = length(new_x_Ultra):-1:1 %      1:lNiryo;
    %%Defining Rotation matrixes
    R_R = transpose([1 0 0; 0 cos(R(i)) sin(R(i)); 0 -sin(R(i)) cos(R(i))]);
    R_P = transpose([cos(P(i)) 0 -sin(P(i)); 0 1 0; sin(P(i)) 0 cos(P(i))]);
    R_Y = transpose([cos(Yaw(i)) sin(Yaw(i)) 0; -sin(Yaw(i)) cos(Yaw(i)) 0; 0, 0, 1]);
    
    %%Defining Vectors from Niryo Robot and Ultrasound machine
    P_N = [x_Niryo(i); y_Niryo(i); z_Niryo(i)];
    P_U = [new_x_Ultra(i); 0; new_z_Ultra(i)];

    %%Calculating the location of the Ulner nerve from the robot's Axis
    P_Final(:,i) = P_N + R_Y * R_P * R_R * P_U;  %This is the important equation
    
    %%Splitting answer apart to be plotted later 
    Q = P_Final(1,:);
    W = P_Final(2,:);
    E = P_Final(3,:);

end

%% This portion of the code adds ten points to the end and the beginning of each Ulnar Nerve position vector
%Finding the beginning and end points
Q_Beginning = P_Final(1,1);
Q_End = P_Final(1,length(P_Final));
W_Beginning = P_Final(2,1);
W_End = P_Final(2,length(P_Final));
E_Beginning = P_Final(3,1);
E_End = P_Final(3,length(P_Final));

%Creating a 20 point vector to go on each end
Q_B(1,1:20) = Q_Beginning;
Q_E(1,1:20) = Q_End;
W_B(1,1:20) = W_Beginning;
W_E(1,1:20) = W_End;
E_B(1,1:20) = E_Beginning;
E_E(1,1:20) = E_End;

%Resizing P Final to put extra data points on the ends
lP_Final = length(P_Final) + 40;
P_Final_Extra = rand(3, lP_Final);

%Putting ten extra points on the end of each vector
P_Final_Extra(1,:) = [Q_B,Q,Q_E];
P_Final_Extra(2,:) = [W_B,W,W_E];
P_Final_Extra(3,:) = [E_B,E,E_E];

%% Plotting Inital Results

% Finds and sets the axis scales to largest and smallest data point
upper_lim = [max(Q) max(W) max(E) max(x_Niryo) max(y_Niryo) max(z_Niryo)];
lower_lim = [min(Q) min(W) min(E) min(x_Niryo) min(y_Niryo) min(z_Niryo)];

%Plotting Position of the Ulnar Nerve without smoothing
figure;
plot3(Q, W, E, 'ro-', 'LineWidth', 3);
title("Ulnar Nerve Location")
xlabel('x');
ylabel('y');
zlabel('z');
xlim([min(lower_lim) max(upper_lim)]);
ylim([min(lower_lim) max(upper_lim)]);
zlim([min(lower_lim) max(upper_lim)]);
grid on

%Plotting just the Niryo Data.
figure
plot3(x_Niryo, y_Niryo, z_Niryo, 'bo-', 'LineWidth', 3);
title("Niryo Data")
xlabel('x');
ylabel('y');
zlabel('z');
xlim([min(lower_lim) max(upper_lim)]);
ylim([min(lower_lim) max(upper_lim)]);
zlim([min(lower_lim) max(upper_lim)]);
grid on

%Plotting just the Ultrasound Data
t = 1:1:lUltra-1;
figure
plot3(t,x_Ultra, z_Ultra, 'yo-', 'LineWidth', 3);
title("Ultrasound Data")
xlabel('Index');
ylabel('x');
zlabel('y');
xlim auto;
ylim auto;
zlim auto;
grid on

%Plotting minor area data
figure
plot(t_Niryo, new_minor_Area_Ultra, 'yo-', 'LineWidth', 3);
title("Minor area data")
xlabel('Time');
ylabel('Area mm^2');
xlim auto;
ylim auto;
grid on

%% Beginning of Ben Moyer Code
%importing data from Josh and CJ code
A = P_Final_Extra; %[x;y;z]

ulnSize = size(A);

                            %user input
                           B = A(1:3,:);

%4th order butterworth filtering on each dimension
lNiryo = length(x_Niryo);
C = transpose(B);
Fc = 0.1;
Fs = lNiryo/t_Niryo(lNiryo);
[y, x] = butter(4, Fc/(Fs/2));
% This is code to see the frequencies in our data and determine a cut off
% length
% ----- PARAMETERS -----
T = 1 / Fs; % Sampling period
L = size(C, 1); % Length of signal
t = (0:L-1) * T; % Time vector
figure;
for i = 1:3
signal = C(:,i);
Y = fft(signal);
P2 = abs(Y / L);
n = floor(L/2);
P1 = P2(1:n+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs * (0:(n)) / L;
subplot(3,1,i);
plot(f, P1, 'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('|Amplitude|');
title(['Frequency Spectrum of Axis ', char('X' + i - 1)]);
grid on;
end
inputSignalx = C(:,1);
outSignalx = filtfilt(y, x, inputSignalx);
inputSignaly = C(:,2);
outSignaly = filtfilt(y, x, inputSignaly);
inputSignalz = C(:,3);
outSignalz = filtfilt(y, x, inputSignalz);
%End of Butterworth filter

%concatenation of filtered data and removal of startup data
D = horzcat(outSignalx,outSignaly,outSignalz);
%D(1:15,:) = []; %This line cuts off the startup data

%plotting filtered data
t1 = D(:,1);
t2 = D(:,2);
t3 = D(:,3);

%intermediate variable declaration
order = 3;
pts = order + 2;
e = 0;
ihat = [double('i'),770];
jhat = [double('j'),770];
khat = [double('k'),770];
ihatstr= char(ihat);
jhatstr = char(jhat);
khatstr = char(khat);

%gets the number of data points for preallocation
size = size(D);
rows = size(1);

%variable preallocation
xpoints = zeros(rows,pts);
ypoints = zeros(rows,pts);
zpoints = zeros(rows,pts);
xto = zeros(rows,order + 1);
yto = zeros(rows,order + 1);
zto = zeros(rows,order + 1);
td = zeros(rows,1);
sf = zeros(rows-(pts-1),1);
rtp = zeros(rows-(pts-1),3);
rtdp = zeros(rows-(pts-1),3);
rttp = zeros(rows-(pts-1),3);
et = zeros(rows-(pts-1),3);
en = zeros(rows-(pts-1),3);
eb = zeros(rows-(pts-1),3);
rho = zeros(rows-(pts-1),1);
tau = zeros(rows-(pts-1),1);
dot1 = zeros(rows-(pts-1),1);
dot2 = zeros(rows-(pts-1),1);
dot3 = zeros(rows-(pts-1),1);
int1 = zeros(rows-(pts-1),1);
int5 = zeros(rows-(pts-1),3);
int6 = zeros(rows-(pts-1),3);
int7 = zeros(rows-(pts-1),1);

syms t xt yt zt xtp ytp ztp xtdp ytdp ztdp stp xttp yttp zttp;

%defining time vector
for j = 1:rows, td(j) = j/(154/19.6); end

%main data process for loop
for i = 1:(rows-(pts-1))

    %gets the current 5 point snippet of each dimension
    xpoints(i,1:pts) = D(i:i+(pts-1),1);
    ypoints(i,1:pts) = D(i:i+(pts-1),2);
    zpoints(i,1:pts) = D(i:i+(pts-1),3);

    %fits the data to the cubic spline (nth order polynomial)
    xto(i,1:order+1) = polyfit(td(i:i+(pts-1)),xpoints(i,:),order);
    yto(i,1:order+1) = polyfit(td(i:i+(pts-1)),ypoints(i,:),order);
    zto(i,1:order+1) = polyfit(td(i:i+(pts-1)),zpoints(i,:),order);

    %creates the variable form of the equations
    xt(i) = poly2sym(xto(i,1:order+1),t);
    yt(i) = poly2sym(yto(i,1:order+1),t);
    zt(i) = poly2sym(zto(i,1:order+1),t);

    %takes the first and second derivatives of the functions of time 
    xtp(i) = diff(xt(i),t);
    ytp(i) = diff(yt(i),t);
    ztp(i) = diff(zt(i),t);
    xtdp(i) = diff(xtp(i),t);
    ytdp(i) = diff(ytp(i),t);
    ztdp(i) = diff(ztp(i),t);
    xttp(i) = diff(xtdp(i),t);
    yttp(i) = diff(ytdp(i),t);
    zttp(i) = diff(ztdp(i),t);

    %gets the function ds
    stp(i) = ((xtp(i)^2)+(ytp(i)^2)+(ztp(i)^2))^0.5;

    %calculates stp, rtp, rtdp, rttp at the current point
    sf(i,1) = subs(stp(i),td(i+2));
    rtp(i,1) = subs(xtp(i),td(i+2));
    rtp(i,2) = subs(ytp(i),td(i+2));
    rtp(i,3) = subs(ztp(i),td(i+2));
    rtdp(i,1) = subs(xtdp(i),td(i+2));
    rtdp(i,2) = subs(ytdp(i),td(i+2));
    rtdp(i,3) = subs(ztdp(i),td(i+2));
    rttp(i,1) = subs(xttp(i),td(i+2));
    rttp(i,2) = subs(yttp(i),td(i+2));
    rttp(i,3) = subs(zttp(i),td(i+2));

    %calculates et, en, rho at the current point
    et(i,1) = rtp(i,1)/sf(i);
    et(i,2) = rtp(i,2)/sf(i);
    et(i,3) = rtp(i,3)/sf(i);
    dot1(i) = (rtp(i,1)*rtdp(i,1)) + (rtp(i,2)*rtdp(i,2)) + (rtp(i,3)*rtdp(i,3));
    dot2(i) = (rtdp(i,1)*rtdp(i,1)) + (rtdp(i,2)*rtdp(i,2)) + (rtdp(i,3)*rtdp(i,3));
    int1(i) = ((dot2(i)*(sf(i)^2)) - (dot1(i)^2))^0.5;
    rho(i) = (sf(i)^3)/int1(i);
    int2x = rtdp(i,1)*(sf(i)^2);
    int2y = rtdp(i,2)*(sf(i)^2);
    int2z = rtdp(i,3)*(sf(i)^2);
    int3x = rtp(i,1)*dot1(i);
    int3y = rtp(i,2)*dot1(i);
    int3z = rtp(i,3)*dot1(i);
    int4x = int2x - int3x;
    int4y = int2y - int3y;
    int4z = int2z - int3z;
    en(i,1) = int4x/(sf(i)*int1(i));
    en(i,2) = int4y/(sf(i)*int1(i));
    en(i,3) = int4z/(sf(i)*int1(i));

    % takes the cross product to find the binormal vector
    eb(i,:) = cross(et(i,:),en(i,:));
    int5(i,:) = cross(rtdp(i,:),rttp(i,:));
    dot3(i) = (rtp(i,1)*int5(i,1)) + (rtp(i,2)*int5(i,2)) + (rtp(i,3)*int5(i,3));
    int6(i,:) = cross(rtp(i,:),rtdp(i,:));
    int7(i) = (int6(i,1)^2) + (int6(i,2)^2) + (int6(i,3)^2);
    tau(i) = dot3(i)/int7(i);

end

% overall prompt while loop
while e == 0

    %asks for data point and stores it in x
    prompt = "data point of interest: ";
    x = input(prompt);

    %increments L for row searching
    for l = 1:(rows-(pts-1))

        %checks if the given point matches the current row
        if x == l

            %prints variables and names
            fprintf("\nRadius of Curvature:   %g cm\n\n",rho(l));
            fprintf("Torsion:   %g cm\n\n",tau(l));
            fprintf("Unit Tangential:   %g",et(l,1));
            fprintf(ihatstr);
            fprintf("   %g",et(l,2));
            fprintf(jhatstr);
            fprintf("   %g",et(l,3));
            fprintf(khatstr);
            fprintf("\n\n");
            fprintf("Unit Normal:   %g",en(l,1));
            fprintf(ihatstr);
            fprintf("   %g",en(l,2));
            fprintf(jhatstr);
            fprintf("   %g",en(l,3));
            fprintf(khatstr);
            fprintf("\n\n");
            fprintf("Unit Binormal:   %g",eb(l,1));
            fprintf(ihatstr);
            fprintf("   %g",eb(l,2));
            fprintf(jhatstr);
            fprintf("   %g",eb(l,3));
            fprintf(khatstr);
            fprintf("\n\n");
        end           
    end

    %prompts and stores if the user wants another point
    prompt2 = "again?: [y/n]";
    txt = input(prompt2,'s');

    %changes the while variable based on answer
    if txt == 'n', e = 1; end
    if txt == 'y', e = 0; end
end

%% End of Ben Moyer Code.  Cutting of the cushining points from the output data
%rho = transpose(rho);
%tau = transpose(tau);
trimmer = length(rho) - 18;
trimmedRho = rho(19:trimmer,1);
trimmedTau = tau(19:trimmer,1);
trimmed_et = et(19:trimmer,1);
trimmed_en = en(19:trimmer,1);
trimmed_eb = eb(19:trimmer,1);

trimmer = length(t1) - 20; 
t1 = t1(21:trimmer,1); 
t2 = t2(21:trimmer,1);
t3 = t3(21:trimmer,1);

%Making Markers
x_Marker = zeros(1,length(markers));
y_Marker = zeros(1,length(markers));
z_Marker = zeros(1,length(markers));
Q_Marker = zeros(1,length(markers));
W_Marker = zeros(1,length(markers));
E_Marker = zeros(1,length(markers));
for i = 1:length(marker_Time)
    tt = 100; %inital start point
   
    for ii = 1:length(t_Niryo)
        markerCompare = abs(t_Niryo(ii) - marker_Time(i));
        if markerCompare <= tt
            tt = markerCompare;
            markers(i) = ii;
            Q_Marker(i) = Q(ii);    %%% raw combined nyrio & ultrasound
            W_Marker(i) = W(ii);
            E_Marker(i) = E(ii);

            x_Marker(i) = t1(ii);    %%% smoothed nyrio & ultrasound
            y_Marker(i) = t2(ii);
            z_Marker(i) = t3(ii);
        end
    end
end

%% Ploting Final Results

% Unfiltered Nyrio & Ultrasound data
figure;
plot3(Q, W, E, 'ro-', 'LineWidth', 3);
hold on; % Used to put data & labels on plot
title("Ulnar Nerve Location")
xlabel('x');
ylabel('y');
zlabel('z');
xlim([min(lower_lim) max(upper_lim)]);  % manually sets the axis so they are all the same 
ylim([min(lower_lim) max(upper_lim)]);
zlim([min(lower_lim) max(upper_lim)]);
axis equal;         %auto corrects the axis so they have the same scale
offset = .02;
% inputs text labels into the grapb
text(Q_Marker, W_Marker, E_Marker-offset, Labels,'FontSize', 10, ...
        'Color', 'k', 'BackgroundColor', [1 1 1], 'EdgeColor', 'k','Margin', 2);   

for i = 1:length(Q_Marker)
    plot3([Q_Marker(i), Q_Marker(i)], [W_Marker(i), W_Marker(i)],[E_Marker(i), E_Marker(i) - offset], ...
        'Color', "r", 'LineStyle', '-');
    plot3(Q(markers(i)), W(markers(i)), E(markers(i)), 'Color', "b", 'Marker', ".", 'MarkerSize',28);
end
grid on

% String of data that is added onto the labels, includes Curvature, and
% area
data_String = "" + newline + "Curvature: " + trimmedRho(markers) + " m"...
                   + newline + "Torsion: " + trimmedTau(markers) + " m"...
                   + newline + "Area: " + new_minor_Area_Ultra(markers) + " mm^2"; 

%Plotting Filtered Data 
figure
plot3(t1,t2,t3,'go-','LineWidth',3);
hold on; %Used to put data & labels on plots
title("Ulnar Nerve Location with Butterworth Filter");
xlabel('x');
ylabel('y');
zlabel('z');
xlim([min(lower_lim) max(upper_lim)]); % manual axis that sets them all to the same, comment out axis equal to enable
ylim([min(lower_lim) max(upper_lim)]);
zlim([min(lower_lim) max(upper_lim)]); 
axis equal;  %auto corrects the axis so they have the same scale
offset = .02;
% inputs text labels into the grapb

for i = 1:length(x_Marker)
    % Creates text box with text an offset distance away from the point
    text(x_Marker(i), y_Marker(i), z_Marker(i)-offset, Labels(i) + data_String(i),'FontSize', 10, ...
        'Color', 'k', 'BackgroundColor', [1 1 1], 'EdgeColor', 'k','Margin', 2);  
    %Draws a line between marker & points of interest
    plot3([x_Marker(i), x_Marker(i)], [y_Marker(i), y_Marker(i)],[z_Marker(i), z_Marker(i)-offset], ...
        'Color', "r", 'LineStyle', '-');
    % This highlights the points in red
    plot3(t1(markers(i)), t2(markers(i)), t3(markers(i)), 'Color', "r", 'Marker', ".", 'MarkerSize',28);
end
grid on

%%Plotting Curvature
figure 
plot(trimmedRho, 'red-', 'LineWidth', 3);
title("Radius of Curvature")
xlabel('Index');
ylabel('Curvature');
grid on

%%Plotting Torsion
figure 
plot(trimmedTau, 'green-', 'LineWidth', 3);
title("Torsion");
xlabel('Index');
ylabel('Torsion');
grid on

%writes all of the path coordinate data to a file
indicators = strings(length(trimmedRho), 1);
indicators(markers) = Labels;
writematrix([indicators(:,1),trimmedRho,trimmedTau,new_minor_Area_Ultra], file_Name);



%% Potential implementation to write data into designated spot in excel sheet
% could be used to run all data at once & store
% main thing would be to function this set of code and the run a new main
% method

%{
Write array A to the 5-by-5 rectangular region, E1:I5, on the first sheet in a new spreadsheet file named testdata.xlsx.

filename = 'testdata.xlsx';
writematrix(A,filename,'Sheet',1,'Range','E1:I5')
Write cell array C to a rectangular region that starts at cell B2 on a worksheet named Temperatures. You can specify range using only the first cell.

writecell(C,filename,'Sheet','Temperatures','Range','B2');
%}