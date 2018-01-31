function vrepthetalist = youbot( thetalist0, totaltime, deltat, Kp, Ki, Xdlist, DXdlist, wantcsv)
% Receive arguments from YoubotSover
% Output (a)vrepthetalist, a 12-column matrix (b)csv file, same data with
% vrepthetalist (c) error plot, for tuning

%   Arguments

% thetalist0: 1x8 matrix, the initial configuration of Youbot
%      elements in thetalist0:[phi, x, y, theta1, theta2, theta3, theta4, theta5]
% totaltime: total time during motion
% deltat: the interval time between two calculation steps
% Kp: proportional control gain
% Ki: integral control gain
% Xdlist: 4x4 matrix describing the trajectory of end-effector about {s} frame
%      each element in Xdlist is the function of s(i,:)
% DXdlist: derivative of Xdlist
%      each element in DXdlist is the function of ds(i,:)unction
% wantcsv: =1 output .csv file for Vrep, =0 do not output for tuning

%   Example is given in YoubotSolver


 t = zeros(1);
 npoints = totaltime/deltat;
 for i=2:npoints
     t(i,:) = t(i-1,:)+deltat;
 end
 Vdlist = zeros(1,6);
 for i = 1:npoints
     j = 4*(i-1)+1;
     Vdlist(i,:) = se3ToVec(TransInv(Xdlist(j:j+3,:))*DXdlist(j:j+3,:)); 
 end
 
 Mb0 = [[1,0,0,0.1662];[0,1,0,0];[0,0,1,0.0026];[0,0,0,1]];
 M0e = [[1,0,0,0.033]; [0,1,0,0]; [0,0,1,0.6546]; [0,0,0,1]];
 Blistarm = [[0,0,1,0,0.033,0]; [0,-1,0,-0.5076,0,0];[0,-1,0,-0.3526,0,0];... 
     [0,-1,0,-0.2176,0,0]; [0,0,1,0,0,0]];
 l = 0.235; w = 0.15; r = 0.0475;
 F = r/4*[[-1/(l+w),1/(l+w),1/(l+w),-1/(l+w)];[1,1,1,1];[-1,1,-1,1]];
 F6 = [[0,0,0,0];[0,0,0,0];F;[0,0,0,0]];
 Wheel = zeros(1,4);
 thetalist = [thetalist0,Wheel];
 thetalistarm = thetalist(4:8);
 q = thetalist(1:3);
 Tsb = [[cos(q(1)),-sin(q(1)),0,q(2)];
         [sin(q(1)),cos(q(1)),0,q(3)];
         [0,0,1,0.0963];[0,0,0,1]];
 T0 = Tsb*Mb0*FKinBody(M0e,Blistarm,thetalistarm);  % Start position of end-effector
 Vlist = zeros(1,6);
 Xerr = zeros(1,6);
 Xerrint = zeros(1,6); 
 Xlist = T0;
 velocitylist = zeros(1,9);
 
 for i=1:npoints
     j = 4*(i-1)+1;
     Xerr(i,:) = MatrixLog6(TransInv(Xlist(j:j+3,:))*Xdlist(j:j+3,:))';
     if i==1
         Xerrint(i,:)=0;
     else
         Xerrint(i,:) = Xerrint(i-1,:) + Xerr(i,:)*deltat;
     end
     X = Xlist(j:j+3,:);
     Xd = Xdlist(j:j+3,:);
     AdX = Adjoint(TransInv(X)*Xd);
     Vlist(i,:) = (AdX*Vdlist(i,:)')' + Kp*Xerr(i,:) + Ki*Xerrint(i,:);   % Control law
     % odometry
     T0e = FKinBody(M0e,Blistarm,thetalistarm(i,:));
     Jbase = Adjoint(TransInv(T0e)*TransInv(Mb0))*F6;
     Jarm = JacobianBody(Blistarm, thetalistarm(i,:));
     Je = [Jbase,Jarm];
     velocitylist(i,:) = (pinv(Je)*Vlist(i,:)')';
     thetadot = velocitylist(i,5:9);    % angular velocities of arm joints 
     thetalistarm(i+1,:) = thetalistarm(i,:)+thetadot*deltat;
     u = velocitylist(i,1:4)';     % angular velocities of wheels
     Wheel(i+1,:) = Wheel(i,:) + u'*deltat;
     Vb = (F * u * deltat)';
     omegabz = Vb(1);vbx = Vb(2);vby = Vb(3);
     if  omegabz == 0
         deltaqb = [0;vbx;vby];
     else
     deltaqb = [omegabz;(vbx*sin(omegabz)+vby*(cos(omegabz)-1))/omegabz;(vby*sin(omegabz)+vbx*(1-cos(omegabz)))/omegabz];
     end
     phik = q(i,1);
     deltaq = ([[1,0,0];[0,cos(phik),-sin(phik)];[0,sin(phik),cos(phik)]]*deltaqb)';
     q(i+1,:) = q(i,:) + deltaq;
     thetalist(i+1,:) = [q(i+1,:),thetalistarm(i+1,:),Wheel(i+1,:)];
     Xlist(j+4:j+7,:) = [[cos(q(i+1,1)),-sin(q(i+1,1)),0,q(i+1,2)];
         [sin(q(i+1,1)),cos(q(i+1,1)),0,q(i+1,3)];
         [0,0,1,0.0963];[0,0,0,1]]*Mb0*FKinBody(M0e,Blistarm,thetalistarm(i+1,:));
 end
 % plot
 plottime = 0:deltat:(5-deltat);
 figure
 for i=1:6
    plot(plottime,Xerr(:,i),'-')
    hold on;
 end
 title('Error Response')
 legend('Xerr(1)','Xerr(2)','Xerr(3)','Xerr(4)','Xerr(5)','Xerr(6)')
 xlabel('Time')
 ylabel('Errors')
 
 % swtich columns for V-rep
 vrepthetalist = thetalist;
 vrepthetalist(:,1:2) = thetalist(:,2:3);
 vrepthetalist(:,3) = thetalist(:,1);
 if wantcsv
     csvwrite('youbotlist.csv',vrepthetalist);
 end

end

