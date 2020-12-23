x_len = 0.05;
y_len = 0.10;
z_len = 0.20;
T = 1;

tetthet = 2000;
volum = x_len*y_len*z_len;
m = volum*tetthet;

%Defs for intertia
J_xx = (m/12)*(y_len^2+z_len^2); 
J_yy = (m/12)*(x_len^2+z_len^2); 
J_zz = (m/12)*(x_len^2+y_len^2); 

J = [J_xx,0,0;0,J_yy,0;0,0,J_zz];

%Defs for rotation around main axis 3
w3_ren = (2*pi)/T;

K0 = 2*pi^2*J_zz;

w1_ren = sqrt((2*K0/J_xx));
w2_ren = sqrt((2*K0/J_yy));


%[X,Y,Z] = ellipsoid(0,0,0,w1_ren,w2_ren,w3_ren);
%surf(X,Y,Z,'FaceColor',[1 0 0]);
%hold on;
%axis equal


w3 = w3_ren/10;
w1 = w1_ren*sqrt(1-(w3^2/w3_ren^2));
w2 = w2_ren*sqrt(1-(w3^2/w3_ren^2));

w1_n = [w1;0;w3];

h0 = sqrt(J_xx^2*w1_n(1)^2+J_zz^2*w1_n(3)^2);

e1 = h0/J_xx;
e2 = h0/J_yy;
e3 = h0/J_zz;


%[X1,Y1,Z1] = ellipsoid(0,0,0,e1,e2,e3);
%surf(X1,Y1,Z1,'FaceColor',[0 1 0]);
%axis equal

%figure()

w1 = w1_ren/10;
w2 = w2_ren*sqrt(1-(w1^2/w1_ren^2));
w3 = w3_ren*sqrt(1-(w1^2/w1_ren^2));


w2_n = [w1;w2;0];

h1 = sqrt(J_xx^2*w2_n(1)^2+J_yy^2*w2_n(2)^2);

y1 = h1/J_xx;
y2 = h1/J_yy;
y3 = h1/J_zz;


%[X,Y,Z] = ellipsoid(0,0,0,w1_ren,w2_ren,w3_ren);
%surf(X,Y,Z,'FaceColor',[1 0 0]);
%hold on;
%axis equal

%[X1,Y1,Z1] = ellipsoid(0,0,0,y1,y2,y3);
%surf(X1,Y1,Z1,'FaceColor',[0 1 0]);
%axis equal
%figure()

w1 = w1_ren/10;
w3 = w3_ren*sqrt(1-(w1^2/w1_ren^2));

w3_n = [w1;0;w3];

h1 = sqrt(J_xx^2*w3_n(1)^2+J_zz^2*w3_n(3)^2);

y1 = h1/J_xx;
y2 = h1/J_yy;
y3 = h1/J_zz;


%[X,Y,Z] = ellipsoid(0,0,0,w1_ren,w2_ren,w3_ren);
%surf(X,Y,Z,'FaceColor',[1 0 0]);
%hold on;
%axis equal

%[X1,Y1,Z1] = ellipsoid(0,0,0,y1,y2,y3);
%surf(X1,Y1,Z1,'FaceColor',[0 1 0]);
%axis equal
%figure()


%f = @(t,x)([(1/J_xx)*(J_yy-J_zz)*x(2)*x(3);(1/J_yy)*(J_zz-J_xx)*x(1)*x(3);(1/J_zz)*(J_xx-J_yy)*x(1)*x(2);])

%[t,xa] = ode45(f,[0 5],[w3_n(1) w3_n(2) w3_n(3)])


%[t,xa] = ode45(@(t,x)odefcn(t,x,J_xx,J_yy,J_zz),[0 5],[w3_n(1) w3_n(2) w3_n(3) 0 0 0])

[t,xa] = ode45(@(t,x)odefcn(t,x,J_xx,J_yy,J_zz),[0 100],[0 0 0 w1_n(1) w1_n(2) w1_n(3)]);
%plot3(xa(:,1),xa(:,2),xa(:,3));
%grid on
%axis equal
%hold on
[t,xb] = ode45(@(t,x)odefcn(t,x,J_xx,J_yy,J_zz),[0 100],[0 0 0 w2_n(1) w2_n(2) w2_n(3)]);
%plot3(xb(:,1),xb(:,2),xb(:,3));
%grid on
%axis equal
%hold on
[t,xc] = ode45(@(t,x)odefcn(t,x,J_xx,J_yy,J_zz),[0 100],[0 0 0 w3_n(1) w3_n(2) w3_n(3)]);
%plot3(xc(:,1),xc(:,2),xc(:,3))
%grid on
%axis equal
%hold on
%figure()


X = [x_len, x_len, -x_len, -x_len
    x_len, x_len, -x_len, -x_len
    x_len, x_len, x_len, x_len
    x_len, x_len, -x_len, -x_len
    -x_len, -x_len, -x_len, -x_len
    x_len, x_len, -x_len, -x_len] ;
Y = [y_len, -y_len, -y_len, y_len
   y_len, -y_len, -y_len, y_len
   -y_len, y_len, y_len, -y_len
     y_len, y_len, y_len, y_len
    -y_len, y_len, y_len, -y_len
     -y_len, -y_len, -y_len, -y_len] ;
Z = [z_len, z_len, z_len, z_len
    -z_len, -z_len, -z_len, -z_len
    -z_len, -z_len, z_len, z_len
    -z_len, z_len, z_len, -z_len
   -z_len, -z_len, z_len, z_len
    -z_len, z_len, z_len, -z_len] ;
C = {'blue' ;'red' ; 'green' ; 'yellow' ; 'magenta' ; 'cyan'};

for i = 1: 500%length(xc(:,1))
    t1 = xc(i,1);
    t2 = xc(i,2);
    t3 = xc(i,3);
    r11 = cos(t3)*cos(t2);
    r12 =cos(t3)*sin(t2)*sin(t1)-sin(t3)*cos(t1);
    r13 = cos(t3)*sin(t2)*cos(t1)+sin(t3)*sin(t1);
    r21 = sin(t3)*cos(t2);
    r22 = sin(t3)*sin(t2)*sin(t1)+cos(t3)*cos(t1);
    r23 = sin(t3)*sin(t2)*cos(t1)-cos(t3)*sin(t1);
    r31 = -sin(t2);
    r32 = cos(t2)*sin(t1);
    r33 = cos(t2)*cos(t1);
    
    R = [r11,r12,r13;r21,r22,r23;r31,r32,r33]
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
    flipbook(i) = getframe(gcf);
    pause(0.01)
    hold off
    clf
    
end

mywriter = VideoWriter('myFile_c','MPEG-4');
mywriter.FrameRate = 20;
open(mywriter);
writeVideo(mywriter,flipbook)
close(mywriter);

%{
f = @(t,x)([(1/J_xx)*(J_yy-J_zz)*x(2)*x(3);(1/J_yy)*(J_zz-J_xx)*x(1)*x(3);(1/J_zz)*(J_xx-J_yy)*x(1)*x(2);])

[t,xa] = ode45(f,[0 100],[w2_n(1) w2_n(2) w2_n(3)])

plot3(xa(:,1),xa(:,2),xa(:,3))
grid on
axis equal
hold on

f = @(t,x)([(1/J_xx)*(J_yy-J_zz)*x(2)*x(3);(1/J_yy)*(J_zz-J_xx)*x(1)*x(3);(1/J_zz)*(J_xx-J_yy)*x(1)*x(2);])

[t,xa] = ode45(f,[0 100],[w1_n(1) w1_n(2) w1_n(3)])

plot3(xa(:,1),xa(:,2),xa(:,3))
grid on
axis equal
hold on
%}

function dxdt = odefcn(t,x,J_xx,J_yy,J_zz)
    dxdt = zeros(6,1);
    dxdt(1) = x(4)+x(6)*cos(x(1))*tan(x(2))+x(5)*sin(x(1))*tan(x(2));
    dxdt(2) = x(5)*cos(x(1))-x(6)*sin(x(1));
    dxdt(3) = (x(6)*cos(x(1))/cos(x(2)))+(x(5)*sin(x(1))/cos(x(2)));
    dxdt(4) = (1/J_xx)*(J_yy-J_zz)*x(5)*x(6);
    dxdt(5) = (1/J_yy)*(J_zz-J_xx)*x(4)*x(6);
    dxdt(6) = (1/J_zz)*(J_xx-J_yy)*x(4)*x(5);
end 
