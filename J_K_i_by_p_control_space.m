
% �ô������ͨ������pJ��pK��Ѱ����������ϵ�swarmalator������Ϊ״̬ 
% This code can find collective behavior states of self-coupling swarmalators by adjusting pJ and pK

clc,clear,close all
% ����ȫ�ֱ��� | Define global variables
global N J K d

% ��������,����Ƶ�ʣ����ƽ��ٶ� | Number of oscillators
N=500;

% ����˫�ķֲ��ĳ�ʼ״̬ | Set initial state for double delta distribution
pJ = 1;
pK = 0.7;
JA = 0.5;
JB = -1;
KA = 1;
KB = -1;

% ����J��K��˫�ķֲ� | Generate double delta distribution for J and K
J = floor(rand(N,1)-pJ)*(-1)*(JA-JB)+JB;
K = floor(rand(N,1)-pK)*(-1)*(KA-KB)+KB;

% ʱ�䷶Χ | Time span
tspan = 0:0.1:500;

% ��ʼ�ռ�λ�ú���λ | Initial spatial positions and phases
initial_space = 2*rand(1, 2*N)-1; % �������-1~1֮��2N������� | Randomly generate 2N numbers between -1 and 1
initial_theta = 2*pi*rand(1, N)-pi; % ����-pi~pi֮��N������� | Generate N random numbers between -pi and pi
initial_space_theta = horzcat(initial_space, initial_theta); % ����������ƴ��Ϊһ������ | Concatenate two vectors

% ������뺯��������֮����빫ʽ | Define distance function between two points
d = @(x1,x2,y1,y2) sqrt((x1-x2).^2+(y1-y2).^2);

% ��ֵ�ⳣ΢�ַ��� | Numerical solution of ordinary differential equations
[t,space_theta_matrix] = ode45(@(t,space_theta) swarmalator_ode(t,space_theta),tspan,initial_space_theta);

figure;
% Ϊÿ��ʱ�̣�N������ģ��ϵͳ��Ϊ | Simulate system behavior for N oscillators at each moment
for i = 1:100:length(tspan)

    % ���ʱ��i���������ӿռ�λ�ú���λ | Get spatial positions and phases of all oscillators at moment i
    velocity = space_theta_matrix(i,:);
    space_x_velocity = velocity(1:2:2*N-1);
    space_y_velocity = velocity(2:2:2*N);  
    theta_velocity = velocity(2*N+1:1:3*N);
    
    % ���Ʀ���-pi��pi�У����Խ���������ѭ�����������Զ��Ӽ�2pi��ʹ��һֱ�ڷ�Χ�� | Restrict �� to -pi to pi, periodic cycling
    theta_velocity = mod(theta_velocity + pi,2*pi) - pi; 
    
    % ����N�������ڵ�ǰʱ�̵�λ�� | Plot positions of N oscillators at current moment
    scatter(space_x_velocity,space_y_velocity,30,theta_velocity,'filled');  
    colormap(jet); % ������ɫӳ�� | Set color mapping
    caxis([-pi pi]);   % �ֶ�������ɫ����Χ | Manually set color bar range
    axis('auto'); % �Զ����������᷶Χ | Auto set axis range
    axis equal;

    xlabel('x');
    ylabel('y');
    %grid on; % ��ʾ������ | Show grid lines
    pause(0.1); % ��ͣ0.1�� | Pause for 0.1 seconds
end

% ������ʱ�̵��������ӿռ�� | Get spatial angles of all oscillators at last moment
space_angle = atan(space_y_velocity./space_x_velocity) + ((sign(space_x_velocity)-1)/2) * pi; 
space_angle = mod(space_angle + pi,2*pi) - pi; 

% ���������S+��S-��ֵ | Calculate order parameters S+ and S-
Sup = abs((sum(exp(1i*(space_angle + theta_velocity))))/N);
Sdown  = abs((sum(exp(1i*(space_angle - theta_velocity))))/N);
S2 = (1/N)*max(abs(sum(exp(1i*(space_angle + 2.*theta_velocity)))),abs((sum(exp(1i*(space_angle - 2.*theta_velocity))))));
R  = abs((sum(exp(1i*theta_velocity)))/N);
R2  = abs((sum(exp(2*1i*theta_velocity)))/N);

% ��ʾ������ | Display calculation results
disp(['S+ = ',num2str(Sup)])
disp(['S- = ',num2str(Sdown)])
disp(['S2 = ',num2str(S2)])
disp(['R = ',num2str(R)])
disp(['R2 = ',num2str(R2)])

% swarmalatorģ�͵�ODEs | ODEs for swarmalator model
function dtheta_and_dspace_dt = swarmalator_ode(t,space_theta)
    global N J K d
    % ��ʼ���������� | Initialize return vector
    dtheta_and_dspace_dt = zeros(3*N,1);
  
    % �ֱ�õ������꣬�����꣬����λ������ | Get vectors of x-coordinates, y-coordinates, and phases
    space_x = space_theta(1:2:2*N-1);
    space_y = space_theta(2:2:2*N);
    theta = space_theta(2*N+1:1:3*N);

    for i=1:N
        % ���е����i����ľ�����ɵ����� | Vector of distances between all points and the i-th point
        distance = d(space_x,space_x(i),space_y,space_y(i));

        % �����̼������ֿ��ɵ����� | Vector of equation series part
        vq_space_x = ((space_x-space_x(i))./distance).*(1+J(i).*cos(theta-theta(i)))-(space_x-space_x(i))./(distance.^2);
        vq_space_y = ((space_y-space_y(i))./distance).*(1+J(i).*cos(theta-theta(i)))-(space_y-space_y(i))./(distance.^2);
        vq_theta = K(i).*(sin(theta-theta(i)))./distance;

        % ������������ΪNaN�����Ϊ0 | Replace NaN items in vector with 0
        vq_space_x(isnan(vq_space_x)) = 0;
        vq_space_y(isnan(vq_space_y)) = 0;
        vq_theta(isnan(vq_theta)) = 0;

        % ���㷽����ʽ | Calculate right side of equation
        dtheta_and_dspace_dt(2*i-1) = (1/N)*sum(vq_space_x);
        dtheta_and_dspace_dt(2*i) = (1/N)*sum(vq_space_y);
        dtheta_and_dspace_dt(i+2*N) = (1/N)*sum(vq_theta);

    end
    % ��ʾ��ǰʱ�� | Display current time
    disp(['t = ' ,num2str(t)]);
end
