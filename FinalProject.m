clear all
close all
clc
syms xi eta x y
%% Define parameters and problem domain
noElemX = 32;
noElemY = 32;

noNodeX = noElemX + 1;
noNodeY = noElemY + 1;

noNode = noNodeX * noNodeY;
noElem = noElemX * noElemY;

nodePerElem = 4;

Lx = 1;
Ly = 1;

u_b12 = 0; % Temperature boundary conditions of surface 1 and 2
h_b34 = 0; % Natural boundary conditions of surface 3 and 4

k = 50; % isotropic thermal

Q = sin(pi*x/2)*sin(pi*y/2);

[nodeCoord,elemDict] = uniformRectangularMesh4Node(Lx,Ly,noElemX, noElemY);  

% Construct shape function for each element 
N1 = shapeFunc1D(xi,1) * shapeFunc1D(eta,1);
N2 = shapeFunc1D(xi,2) * shapeFunc1D(eta,1);
N3 = shapeFunc1D(xi,2) * shapeFunc1D(eta,2);
N4 = shapeFunc1D(xi,1) * shapeFunc1D(eta,2);

N = [N1 N2 N3 N4];
dN = [diff(N,xi)
      diff(N,eta)];

u_true = 2/(50 * pi^2) * sin(x*pi/2) * sin(pi*y/2);
% Define Gaussian Integration Parameter
GaussIntegrationPt = 4;
we1 = 1;
we2 = 1;
pts = [-1/sqrt(3) -1/sqrt(3);
       1/sqrt(3) -1/sqrt(3);
       1/sqrt(3) 1/sqrt(3);
       -1/sqrt(3) 1/sqrt(3)];


%% Initialize global K and global F
globalK = zeros(noNode, noNode);
globalF = zeros(noNode, 1);

%% Iterate over each element and assembly the Ke into global K and Fe into global F

for elem=1:noElem
    % Initialize the element stiffness matrix Ke and element body force with 0
    Ke = zeros(nodePerElem,nodePerElem);
    fe_b = zeros(nodePerElem,1);

    elemNode = elemDict(elem,:);
    xe = [nodeCoord(elemNode(1),1);nodeCoord(elemNode(2),1);nodeCoord(elemNode(3),1); nodeCoord(elemNode(4),1)];
    ye = [nodeCoord(elemNode(1),2);nodeCoord(elemNode(2),2);nodeCoord(elemNode(3),2); nodeCoord(elemNode(4),2)];
    x_ = N * xe;
    y_ = N * ye;
    Jve = volumeJacobian2D(x_,y_,xi,eta);
    B = Jve\dN;

    % Calculate element stiffness Ke using Gauss Integration
    for e=1:GaussIntegrationPt
        Ke = Ke + (subs(B', [xi,eta], pts(e,:)) * subs(B,[xi,eta], pts(e,:))* subs(k,[xi,eta], pts(e,:)) * subs(det(Jve),[xi,eta], pts(e,:)) *  we1 * we2);
    end
    Ke = double(Ke);

    % Assemble Ke into global K
    globalK(elemNode,elemNode) = globalK(elemNode,elemNode) + Ke;

    % Q in parametric variables
    Q_xi_eta = subs(Q, [x,y], [x_,y_]); 
    
    % Calculate element body forces using Gauss Integration
    for e=1:GaussIntegrationPt
        fe_b = fe_b + (subs(N', [xi,eta], pts(e,:)) * subs(Q_xi_eta,[xi,eta], pts(e,:))* we1 * we2 * subs(det(Jve),[xi,eta], pts(e,:)));
    end
    fe_b = double(fe_b);
    % Assemble Fe into global F
    globalF(elemNode,:) = globalF(elemNode,:) + fe_b;
end

%% Apply Boundary conditions and solve for u 
% At any nodes that have x = 0 or y = 0, we will applied the essential BC
% u=0
x_bc = 0;
y_bc = 0;
node_BC = [];
% Extract node number that lies in the natural BC requirements
for node=1:noNode
    node_x = nodeCoord(node,1);
    node_y = nodeCoord(node,2);
    if (node_x==x_bc || node_y==y_bc || node_x==x_bc && node_y==y_bc)
        node_BC = [node_BC, node];
    end
end

% bigDes = zeros(noNode,1);
% bigDes(node_BC,1) = 0; % Natural BC imposed to all node in the chosen boundary surface
bigFBC = globalF;
bigFBC(node_BC,1) = 0;
bigKBC = globalK;

bigKBC(:, node_BC) = 0;
bigKBC(node_BC, :) = 0;
for i=1:length(node_BC)
    bigKBC(node_BC(i), node_BC(i)) = 1; % Set diagonal to 1
end

bigFBC(node_BC, 1) = 0;

% Find nodal temperature
u = bigKBC\bigFBC;

%% Plot the result of u on the line x=0.5 and compare with the exact result
x_chosen = 0.5;
y_exact = linspace(0,1,100);
u_exact = 2/(50 * pi^2) * sin(0.5*pi/2) * sin(pi*y_exact/2);
node_plot=[];


% Extract node that has x=0.5
for node=1:noNode
    node_x = nodeCoord(node,1);
    if node_x == x_chosen
        node_plot = [node_plot, node];
    end
end
figure()
u_apprx = u(node_plot);

y_mesh = nodeCoord(node_plot,2);
dudy = zeros(length(y_mesh)-1,1);
ax1 = gcf;
filename = sprintf("u_compare_%ix%i.png", noElemX, noElemY);
plot(y_mesh,u_apprx, "-o", LineWidth=1.25)
hold on
plot(y_exact, u_exact)
title("u^h on the line x = 0.5")
xlabel("y")
ylabel("u")
legend(["u^h" "u_{exact}"], Location="best")
exportgraphics(ax1, filename, "Resolution", 300)
hold off



%% Calculate dudy on the line x=0.5 and plot and compare
duy_exact = subs(diff(u_true,y),x,0.5);
duy_exact = double(subs(duy_exact,y,y_exact));
for elemNum=1:length(y_mesh)-1 
    y1 = y_mesh(elemNum);
    y2 = y_mesh(elemNum + 1);
    
    elemLength = y2 - y1;
    Be = [-1/elemLength 1/elemLength];
    de= [u_apprx(elemNum); u_apprx(elemNum+1)];
    dudy(elemNum,1) = Be * de;
end
% Plot du/dy
filename = sprintf("dudy_compare_%ix%i.png", noElemX, noElemY);
y_step = repelem(y_mesh,2);
y_step = y_step(2:end-1);
dudy_step = repelem(dudy,2);

figure()
ax2 = gcf;
plot(y_step,dudy_step, Linewidth=1.5)
hold on
plot(y_exact, duy_exact)

xlabel("y")
ylabel("u_{,y}")
title("u_{,y}^h on the line x = 0.5")
legend(["u_{,y}^h" "Exact u_{,y}"], Location="best")
exportgraphics(ax2, filename, "Resolution", 300)
hold off
%% Calculate error e over the entire domain
NoErrGaussPtsXi = 4;
NoErrGaussPtsEta = 4;

% Weights and Integration point in xi and eta direction
[GPXi, w_xi] = lgwt(4,-1,1);
[GPEta, w_eta] = lgwt(4,-1,1);
total_error = 0;
for elem=1:noElem
    element_error = 0;
    elemNode = elemDict(elem,:);
    xe = [nodeCoord(elemNode(1),1);nodeCoord(elemNode(2),1);nodeCoord(elemNode(3),1); nodeCoord(elemNode(4),1)];
    ye = [nodeCoord(elemNode(1),2);nodeCoord(elemNode(2),2);nodeCoord(elemNode(3),2); nodeCoord(elemNode(4),2)];
    x_ = N*xe;
    y_ = N*ye;
    Jve = volumeJacobian2D(x_,y_,xi,eta);
    ue_exact = sin(pi*x_/2)*sin(pi*y_/2)*1/(25*pi^2);
    ue_approx = N*u(elemNode);
    tmp = (ue_exact - ue_approx)^2 * det(Jve); % Need to integrate this over the element domain (xi=-1 to 1 and eta=-1 to 1)

    % Use Gauss Integration to integrate the error over each domain
    for i=1:NoErrGaussPtsEta
        for j=1:NoErrGaussPtsXi
            element_error = element_error + (subs(tmp,[xi,eta],[GPXi(j), GPEta(i)])*w_xi(j)*w_eta(i));
        end
    end
    total_error = total_error + (double(element_error));

end
total_error = sqrt(total_error)

%% Visualize the mesh 
% figure;
% hold on;
% axis equal;
% xlabel('X');
% ylabel('Y');
% title('Finite Element Mesh');
% 
% % Loop through each element
% for i = 1:size(elemDict, 1)
%     % Get node coordinates of the current element
%     element_nodes = elemDict(i, :); % Node indices of the element
%     coords = nodeCoord(element_nodes, :); % Get the coordinates of these nodes
% 
%     % Plot the element (close the loop by repeating the first node)
%     fill(coords(:, 1), coords(:, 2), 'white', 'EdgeColor', 'black', 'FaceAlpha', 0.5);
% 
%      % Compute the centroid of the element
%     centroid = mean(coords, 1);
% 
%     % Annotate the element number at the centroid
%     text(centroid(1), centroid(2), sprintf('%d', i), 'Color', 'black', ...
%          'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% end
% 
% % Plot the nodes
% plot(nodeCoord(:, 1), nodeCoord(:, 2), 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'black');
% 
% % Annotate the nodes
% for i = 1:size(nodeCoord, 1)
%     text(nodeCoord(i, 1), nodeCoord(i, 2), sprintf(' %d', i), 'Color', 'blue');
% end
