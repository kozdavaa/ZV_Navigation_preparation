clear;
close all;
%% Preparation (load data and visualize)

% load distance from the j anchor
load measurements;

% load anchors point
% anchors(k,1) - k anchors x cordinate
% anchors(k,2) - k anchors y cordinate
load anchors;

% load real location
% real_loc(i, 1) - in i time the dev x pos
% real_loc(i, 2) - in i time the dev y pos
load real_loc;

%store the gaus newton calculated positions
nghist = [];

%kalman filter variable
kfhist = [];
m = [];
P = eye(2);
Q = eye(2)*1;

%draw the anchor points and the device path
for step=1:size(d,1)
    clf;

    % draw anchors point
    main_plot = scatter(anchors(:,1), anchors(:,2), 'r');
    axis([-12 22 -12 12]);
    axis equal;
    hold on;

    % plot real location
    plot(real_loc(1:step,1), real_loc(1:step,2),'rx-','LineWidth',2);

    % draw distance in circle
    for i=1:size(d,2)
        circle(anchors(i,1),anchors(i,2),d(step,i),'b');
    end
    %% Newton - Gauss algoritmus
    ngloc = [0;0];
    
    for iter=1:10
        J = [] ; % Jacobi matrix
        for ai = 1:size(anchors,1)
            J(ai,1) = (ngloc(1)-anchors(ai,1))/sqrt(((ngloc(1)-anchors(ai,1))^2 )+((ngloc(2)-anchors(ai,2))^2));
            J(ai,2) = (ngloc(2)-anchors(ai,2))/sqrt(((ngloc(1)-anchors(ai,1))^2 )+((ngloc(2)-anchors(ai,2))^2));

            eps(ai,1) = sqrt(((ngloc(1)-anchors(ai,1))^2 )+((ngloc(2)-anchors(ai,2))^2)) - d(step,ai);
        end;
        ngloc = ngloc - inv(J'*J)*J'*eps;
    end;
    nghist = [nghist;ngloc'];
    plot ( nghist (: ,1) ,nghist (: ,2) ,'gx-','LineWidth',2);


    %% Extended kalman filter solution
    if isempty(m)
        m = ngloc; %% Start point from Gauss Newton
    end;

    % Prediction step
    m = m;
    P = P+Q;
    H = [] ; % Jacobi matrix
    for ai = 1:size(anchors,1)
        H(ai,1) = (m(1)-anchors(ai,1))/sqrt(((m(1)-anchors(ai,1))^2 )+((m(2)-anchors(ai,2))^2));
        H(ai,2) = (m(2)-anchors(ai,2))/sqrt(((m(1)-anchors(ai,1))^2 )+((m(2)-anchors(ai,2))^2));

        h(ai,1) = sqrt(((m(1)-anchors(ai,1))^2 )+((m(2)-anchors(ai,2))^2));
    end;
    R = diag(dvar(step,:));

    % update step
    v = d(step, :)'-h;
    S = H*P*H' + R;
    K = P*H'*inv(S);
    m=m+K*v;
    P = P-K*S*K';

    % Path save and draw
    kfhist=[kfhist;m'];
    plot(kfhist(:,1),kfhist(:,2),'bx-','LineWidth',2);
    circle(kfhist(end,1),kfhist(end,2),P(1,1),'b');
    pause;
end;