%% Preliminaries

% Author: Daibo Zhang (A13591601)
% UCSD MAE290B WI22 Final Project

% Clean up
clear all;
close all;
clc;
format long;

% Load parameter list
fpParam  = readtable('fpParam.csv');

% Define constants
alpha = fpParam.alpha; 
a = fpParam.a;
Omega = fpParam.Omega;
dt = fpParam.dt;
h = fpParam.h;
kappa = alpha * dt / (2 * h^2);

% Define spatial grid
xBound = [fpParam.xBound1,fpParam.xBound2];
yBound = [fpParam.yBound1,fpParam.yBound2];
[xMesh,yMesh] = meshgrid(xBound(1):h:xBound(2),yBound(1):h:yBound(2));
N = length(xMesh) - 2;

% Define discretized initial condition
icB = 0.01 * sin(pi*xMesh) .* sin(pi*yMesh);
icBi = icB(2:end-1,2:end-1);
icVecB = reshape(icBi,N*N,1);

% Define source term
spaceQ = 2.5 * sin(4*pi*xMesh) .* sin(8*pi*yMesh);
spaceQ = spaceQ(2:end-1,2:end-1);
spaceQVec = reshape(spaceQ,N*N,1);
q = @(t) (1 - exp(-a*t) * sin(Omega*t) * cos(2*Omega*t)) * spaceQVec;

% Construct matrices used in ADI scheme
plusXX = speye(N*N) + kappa * fdaMatX(N,N);
minusXX = speye(N*N) - kappa * fdaMatX(N,N);
plusYY = speye(N*N) + kappa * fdaMatY(N,N);
minusYY = speye(N*N) - kappa * fdaMatY(N,N);

%% Test Thomas algorithm implementation

% Build a random tridiagonal system
nTest = 100;
dTest = 20;
RTest = spdiags(rand(nTest,3),-dTest:dTest:dTest,nTest,nTest);
bTest = rand(nTest,1);

% Solve Rx = b using MATLAB built-in
x_exp = RTest\bTest;

% Solve using Thomas algorithm
x_tom = tridiagSolve(RTest,bTest);

% Compare the two solutions
disp('Testing Thomas algorithm implementation')
fprintf('Solving a random %d-by-%d system. Error = %1.4e \n',...
    nTest,nTest,norm(x_exp-x_tom));

%% Part b: simulate to steady state

% Define tolerance for change in solution that defines the steady state
tol = 1e-5;

% Preallocate counter
nB = 0; 

% Location of the point to plot
pLoc = [find(xBound(1):h:xBound(2)==0.55),...
    find(yBound(1):h:yBound(2)==0.45)];

% Preallocate time and solution output
tOutB = 0;
uOutB = icB;
pointOutB = uOutB(pLoc(1),pLoc(2));

% Initialize a variable that keep track of changes in solution
delU = norm(uOutB(:,:,1),'fro');

% Initialize a variable for previous step solution
uPrev = icVecB;

% Indicate start of part b computation
disp('Begin Part b: solving the reaction-diffusion equation');

% Keep solving as solution change is greater than tolerance
while delU >= tol
    
    % Start step time
    tic;
    
    % Obtain the next solution as a vector without BC
    uNext = rdeStepADI(uPrev,nB,tOutB(end),dt,...
        plusXX,minusXX,plusYY,minusYY,q);
    
    % Append solution output and time vector
    tOutB = [tOutB,tOutB(nB+1)+dt];
    uOutB(:,:,nB+2) = appendBC(reshape(uNext,N,N));
    pointOutB(nB+2) = uOutB(pLoc(1),pLoc(2),nB+2);
    
    % Update running variables
    delU = norm(uOutB(:,:,nB+2)-uOutB(:,:,nB+1),'fro');
    uPrev = uNext;
    nB = nB + 1;
    
    % Print performance tracking message
    fprintf('Part b: n=%d done, delU = %1.4e. ',nB-1,delU);
    fprintf('This step took %1.2f second \n',toc);
    
end

% Indicate start of part b computation
disp('Part b steady state solution found');

% Output result data
uB = matfile('solutionsPartB.mat');
uB = uOutB;
save('solutionsPartB.mat')

%% Part d: simulate to steady state

% Define tolerance for change in solution that defines the steady state
tol = 1e-5;

% Preallocate counter
nD = 0; 

% Location of the point to plot
pLoc = [find(xBound(1):h:xBound(2)==0.55),...
    find(yBound(1):h:yBound(2)==0.45)];

% Preallocate time and solution output
tOutD = 0;
uOutD = zeros(N+2,N+2);
pointOutD = uOutD(pLoc(1),pLoc(2));

% Initialize a variable that keep track of changes in solution
delU = 1000;

% Initialize a variable for previous step solution
uPrev = zeros(N*N,1);

% Indicate start of part b computation
disp('Begin Part b: solving the reaction-diffusion equation');

% Keep solving as solution change is greater than tolerance
while delU >= tol
    
    % Start step time
    tic;
    
    % Obtain the next solution as a vector without BC
    uNext = rdeStepADI(uPrev,nD,tOutD(end),dt,...
        plusXX,minusXX,plusYY,minusYY,q);
    
    % Append solution output and time vector
    tOutD = [tOutD,tOutD(nD+1)+dt];
    uOutD(:,:,nD+2) = appendBC(reshape(uNext,N,N));
    pointOutD(nD+2) = uOutD(pLoc(1),pLoc(2),nD+2);
    
    % Update running variables
    delU = norm(uOutD(:,:,nD+2)-uOutD(:,:,nD+1),'fro');
    uPrev = uNext;
    nD = nD + 1;
    
    % Print performance tracking message
    fprintf('Part d: n=%d done, delU = %1.4e. ',nD-1,delU);
    fprintf('This step took %1.2f second \n',toc);
    
end

% Indicate start of part d computation
disp('Part d steady state solution found');

% Output result data
uD = matfile('solutionsPartD.mat');
uD = uOutD;
save('solutionsPartD.mat')

%% Helper functions

% Constructing FDA matrix for x second partial derivative
function outMat = fdaMatX(Nx,Ny)

% A function that constructs a matrix representation of a second order
% central difference scheme to approximate a second partial derivative in
% x. The input Nx, Ny are number of interior points in the x and y grid.
% The output matrix is composed of Nx copies of Ny by Ny matrices

% Find the element of the block diagonal matrix
I = speye(Ny);

% Define the matrix used for Kronecker tensor product
P = spdiags([ones(Nx,1),-2*ones(Nx,1),ones(Nx,1)],-1:1,Nx,Nx);

% Construct the output matrix
outMat = kron(P,I);

end

% Constructing FDA matrix for x second partial derivative
function outMat = fdaMatY(Nx,Ny)

% A function that constructs a matrix representation of a second order
% central difference scheme to approximate a second partial derivative in
% y. The input Nx, Ny are number of interior points in the x and y grid.
% The output matrix is composed of Nx copies of Ny by Ny matrices

% Define the matrix used for Kronecker tensor product
P = speye(Nx);

% Find the element of the block diagonal matrix
A = spdiags([ones(Ny,1),-2*ones(Ny,1),ones(Ny,1)],-1:1,Ny,Ny);

% Construct the output matrix
outMat = kron(P,A);

end

% Append 0 Dirichlet boundary conditions to a solution matrix
function solWithBC = appendBC(solOut)

% A function that add homogeneous Dirichlet boundary to the solution, i.e.
% adds a ring of zero around the matrix representing solution values at
% interior points

% Append a ring of zero around input matrix
solWithBC = zeros(size(solOut,1)+2,size(solOut,2)+2);
solWithBC(2:end-1,2:end-1) = solOut;

end