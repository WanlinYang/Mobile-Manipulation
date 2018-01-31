%% Explaination

% The YoubotSover can output the motion control result including error plot
% and .csv file. User can output error plot only for tuning by setting wantcsv = 0

% necessary functions are "youbot.m" and all ME449 functions

% Inputs are commented with "user input"

% totaltime: total time during motion
% deltat: the interval time between two calculation steps
% thetalist0: 1x8 matrix, the initial configuration of Youbot
%      elements in thetalist0:[phi, x, y, theta1, theta2, theta3, theta4, theta5]
% s(i,:): function of t(i,:), the range in s should be (0,1)
% ds(i,:): derivative of s(i,:)
% Xdlist: 4x4 matrix describing the trajectory of end-effector about {s} frame
%      each element in Xdlist is the function of s(i,:)
% DXdlist: derivative of Xdlist
%      each element in DXdlist is the function of ds(i,:)
% Kp: proportional control gain
% Ki: integral control gain
% wantcsv: =1 output .csv file for Vrep, =0 do not output for tuning


%%  initialization
 clear all
 s = zeros(1);
 ds = zeros(1);
 t = zeros(1);
 Xdlist = zeros(4,4);
 DXdlist = zeros(4,4);
 %% input totaltime and deltat
 totaltime = 5;                          %user input
 deltat = 0.01;                          %user input
 npoints = totaltime/deltat;
 for i=2:npoints
     t(i,:) = t(i-1,:)+deltat;
 end
 %% input initial position thetalist0
 thetalist0 = [0,0,0,0,0,-pi/2,pi/4,0];    %user input
 
 %% input s and derivative of s
 for i = 1:npoints                 
     s(i,:) = (3/25)*t(i,:)^2-(2/125)*t(i,:)^3; %user input
 end
 for i = 1:npoints
     ds(i,:) = 6/25*t(i,:)-6/125*t(i,:)^2;      %user input
 end
 %% input trajectory and derivative of trajectory
 for i=1:npoints
     j = 4*(i-1)+1;
     Xdlist(j,:) = [sin(pi/2*s(i,:)),0,cos(pi/2*s(i,:)),s(i,:)+1]; %user input
     Xdlist(j+1,:) = [0,1,0,s(i,:)];                               %user input
     Xdlist(j+2,:) = [-cos(pi/2*s(i,:)),0,sin(pi/2*s(i,:)),0.3+0.2*s(i,:)];%user input
     Xdlist(j+3,:) = [0,0,0,1];                                    %user input
 end
 for i=1:npoints
     j = 4*(i-1)+1;
     DXdlist(j,:) = [pi/2*cos(pi/2*s(i,:))*ds(i,:),0,-pi/2*sin(pi/2*s(i,:))*ds(i,:),ds(i,:)];%user input
     DXdlist(j+1,:) = [0,0,0,ds(i,:)];                             %user input
     DXdlist(j+2,:) = [pi/2*sin(pi/2*s(i,:))*ds(i,:),0,pi/2*cos(pi/2*s(i,:))*ds(i,:),0.2*ds(i,:)];%user input
     DXdlist(j+3,:) = [0,0,0,0];                                   %user input
 end
 %% PI control tuning and output csv file
 Kp = 5;                                 %user input
 Ki = 2;                                 %user input
 wantcsv = 0;                            %user input
 vrepthetalist = youbot(thetalist0, totaltime, deltat, Kp, Ki, Xdlist, DXdlist, wantcsv);
 