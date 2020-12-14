%% Question 1

% a brick with dimention
a = 0.05; %[m]
b = 0.1;  %[m]
c = 0.2;  %[m]

% brick density
d = 2000;  %[kg/m^3]

% brick mass
M = a*b*c*d;

% Inertia matrix
Jxx = (M/12)*(b^2 + c^2);
Jyy = (M/12)*(a^2 + c^2);
Jzz = (M/12)*(a^2 + b^2);

I = [Jxx; Jyy; Jzz];
Jb = diag(I);

% plotting a brick
%patch([0,a,a,0],[0,0,b,b],[0,0,0,0],'blue')
%patch([0,a,a,0],[0,0,b,b],[c,c,c,c],'blue')
%patch([0,a,a,0],[0,0,0,0],[0,0,c,c],'red')
%patch([0,a,a,0],[b,b,b,b],[0,0,c,c],'red')
%patch([0,0,0,0],[0,0,b,b],[0,c,c,0],'green')
%patch([a,a,a,a],[0,0,b,b],[0,c,c,0],'green')

% plotting a brick and origin vector
%plotcube([0.05 0.1 0.2],[ 2  2  2],.8,[1 0 0])




%% Question 2

%Kinetic energy ellipsoids

%b3 as main-axis
w3_tilde_b3 = 2*pi;
K0 = 0.5*Jzz*w3_tilde_b3^2;
% angular velocity near b1 (or between b1 and b3)
w2_b3 = 0;
w3_b3 = w3_tilde_b3/10;
w1_tilde_b3 = sqrt(2*K0/Jxx);
w1_b3 = w1_tilde_b3*sqrt(1-(w3_b3^2)/(w3_tilde_b3^2));
w3 = [w1_b3;w2_b3;w3_b3];

%b2 as main-axis
w2_tilde_b2 = 2*pi;
% angular velocity near b3 (or between b2 and b3)
w1_b2 = 0;
w2_b2 = w2_tilde_b2/10;
w3_tilde_b2 = sqrt(2*K0/Jzz);
w3_b2 = w3_tilde_b2*sqrt(1-(w2_b2^2)/(w2_tilde_b2^2));
w2 = [w1_b2;w2_b2;w3_b2];

%b1 as main-axis
w1_tilde_b1 = 2*pi;
% angular velocity near b2 (or between b1 and b2)
w3_b1 = 0;
w1_b1 = w1_tilde_b3/10;
w2_tilde_b1 = sqrt(2*K0/Jyy);
w2_b1= w2_tilde_b1*sqrt(1-(w1_b1^2)/(w1_tilde_b1^2));
w1 = [w1_b1;w2_b1;w3_b1];

%Kinetic energy ellipsoid
K0_b1 = sqrt(2*K0/Jxx);
K0_b2 = sqrt(2*K0/Jyy);
K0_b3 = sqrt(2*K0/Jzz);

%3 different angular momentum ellipsoid from 3 different axes
H0_b1 = sqrt(Jxx^2*w1_b1^2 + Jyy^2*w2_b1^2 + Jzz^2*w3_b1^2);
H0_b2 = sqrt(Jxx^2*w1_b2^2 + Jyy^2*w2_b2^2 + Jzz^2*w3_b2^2);
H0_b3 = sqrt(Jxx^2*w1_b3^2 + Jyy^2*w2_b3^2 + Jzz^2*w3_b3^2);

%Rotating at b1 axes
figure(1)
[x,y,z] = ellipsoid(0,0,0,K0_b1,K0_b2,K0_b3);
surf(x,y,z,'FaceColor','r')
hold on
[x,y,z] = ellipsoid(0,0,0,H0_b1/Jxx, H0_b1/Jyy, H0_b1/ Jzz, 100);
surf(x,y,z,'FaceColor','b')
axis equal
xlabel('X')
ylabel('Y')
zlabel('Z')
hold off

%Rotating at b2 axes
figure(2)
[x,y,z] = ellipsoid(0,0,0,K0_b1,K0_b2,K0_b3);
surf(x,y,z,'FaceColor','r','FaceAlpha',0.5)
hold on
[x,y,z] = ellipsoid(0,0,0,H0_b2/Jxx, H0_b2/Jyy, H0_b2/ Jzz, 100);
surf(x,y,z,'FaceColor','b')
axis equal
xlabel('X')
ylabel('Y')
zlabel('Z')
hold off

%Rotating at b3 axes
figure(3)
[x,y,z] = ellipsoid(0,0,0,K0_b1,K0_b2,K0_b3);
surf(x,y,z,'FaceColor','r')
hold on
[x,y,z] = ellipsoid(0,0,0,H0_b3/Jxx, H0_b3/Jyy, H0_b3/ Jzz, 100);
surf(x,y,z,'FaceColor','b','FaceAlpha',0.5)
axis equal
xlabel('X')
ylabel('Y')
zlabel('Z')
hold off

%% Question 3

%b1 axes

syms b1 b11 b2 b22 b3 b33 
syms t1 t11 t2 t22 t3 t33

t0 = 0:0.01:10;


%b1 axes
omgb1 = @(t1, b1) [(Jyy-Jzz)*b1(2)*b1(3)/Jxx;...
                     (Jzz-Jxx)*b1(3)*b1(1)/Jyy;...
                     (Jxx-Jyy)*b1(1)*b1(2)/Jzz];
[t1, b1] = ode45(omgb1,t0,[w1(1), w1(2), w1(3)]);

omgb11 = @(t11, b11) [(Jyy-Jzz)*b11(2)*b11(3)/Jxx;...
                     (Jzz-Jxx)*b11(3)*b11(1)/Jyy;...
                     (Jxx-Jyy)*b11(1)*b11(2)/Jzz];
[t11, b11] = ode45(omgb11,t0,[-2.8745, 0.42722, -2.2108]);


%b2 axes
omgb2 = @(t2, b2) [(Jyy-Jzz)*b2(2)*b2(3)/Jxx;...
                     (Jzz-Jxx)*b2(3)*b2(1)/Jyy;...
                     (Jxx-Jyy)*b2(1)*b2(2)/Jzz];
[t2, b2] = ode45(omgb2,t0,[w2(1), w2(2), w2(3)]);

omgb22 = @(t22, b22) [(Jyy-Jzz)*b22(2)*b22(3)/Jxx;...
                     (Jzz-Jxx)*b22(3)*b22(1)/Jyy;...
                     (Jxx-Jyy)*b22(1)*b22(2)/Jzz];
[t22, b22] = ode45(omgb22,t0,[-0.53281, 0.3446, -6.1427]);

%b3 axes
omgb3 = @(t3, b3) [(Jyy-Jzz)*b3(2)*b3(3)/Jxx;...
                     (Jzz-Jxx)*b3(3)*b3(1)/Jyy;...
                     (Jxx-Jyy)*b3(1)*b3(2)/Jzz];
[t3, b3] = ode45(omgb3,t0,[w3(1), w3(2), w3(3)]);

omgb33 = @(t33, b33) [(Jyy-Jzz)*b33(2)*b33(3)/Jxx;...
                     (Jzz-Jxx)*b33(3)*b33(1)/Jyy;...
                     (Jxx-Jyy)*b33(1)*b33(2)/Jzz];
[t33, b33] = ode45(omgb33,t0,[3.1298, 0.2312, 0]);

figure(4)
[x,y,z] = ellipsoid(0,0,0,K0_b1,K0_b2,K0_b3);
surf(x,y,z,'FaceColor','r','FaceAlpha',0)
hold on

plot3(b1(:,1),b1(:,2),b1(:,3))
grid on
plot3(b11(:,1),b11(:,2),b11(:,3))
grid on

plot3(b2(:,1),b2(:,2),b2(:,3))
grid on
plot3(b22(:,1),b22(:,2),b22(:,3))
grid on

plot3(b3(:,1),b3(:,2),b3(:,3))
grid on
%plot3(b33(:,1),b33(:,2),b33(:,3))
%grid on

axis equal
xlabel('X')
ylabel('Y')
zlabel('Z')
hold off

%% Question 4

X = [a, a, -a, -a
    a, a, -a, -a
    a, a, a, a
    a, a, -a, -a
    -a, -a, -a, -a
    a, a, -a, -a] ;
Y = [b, -b, -b, b
   b, -b, -b, b
   -b, b, b, -b
     b, b, b, b
    -b, b, b, -b
     -b, -b, -b, -b] ;
Z = [c, c, c, c
    -c, -c, -c, -c
    -c, -c, c, c
    -c, c, c, -c
   -c, -c, c, c
    -c, c, c, -c] ;
C = {'green' ;'cyan' ; 'yellow' ; 'red' ; 'cyan' ; 'magenta'};


for i = 1: 200          %length(b3(:,1))
    clf
    t1 = b1(i,1);
    t2 = b1(i,2);
    t3 = b1(i,3);
    r11 = cos(t3)*cos(t2);
    r12 =-sin(t3);
    r13 = 0;
    r21 = sin(t3)*cos(t2);
    r22 = cos(t3);
    r23 = 0;
    r31 = -sin(t2);
    r32 = 0;
    r33 = 1;
    
    R = [r11,r12,r13;
        r21,r22,r23;
        r31,r32,r33]
    T = [X(:) Y(:) Z(:)]*R ;
    Xi = reshape(T(:,1),[],4) ;
    Yi = reshape(T(:,2),[],4) ;
    Zi = reshape(T(:,3),[],4) ;
    hold on
    for j = 1:6
        patch(Xi(j,:), Yi(j,:), Zi(j,:),C{j}) ;
    end
    axis([-0.5  0.5 -0.5 0.5 -0.5 0.5])
    view([30 35])
    movieVector(i) = getframe(gcf);
    pause(0.01)
    hold off
    clf
    
end
% saving the movie
writer = VideoWriter('brick_b1','MPEG-4');
writer.FrameRate = 20;

% open the VideoWriter object, write the movie and close the file
open(writer);
writeVideo(writer,movieVector);
close(writer);


