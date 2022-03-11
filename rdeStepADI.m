function uOutVec = rdeStepADI(uIn,nIn,tIn,tStep,pXX,mXX,pYY,mYY,funcQ)

% A function that solves a 2D reaction diffusion equation using the
% alternating direction implicit (ADI) method. The equation to be solved is
% u_t = alpha (u_xx + u_yy) + Q where subscripts indicate partial
% derivatives. Q is the reaction/source term that is a function of x,y,t. 
% This function only advances the solution by a steps. The exact scheme
% used for the one step advancement differs for even and odd n. They are
% defined by helper function evenStepADI and oddStepADI. The computational
% process involves solving large systems of linear equations represented by
% tridiagonal matrices. This is accomplished efficiently using the Thomas
% Algorithm.
%
% Author: Daibo Zhang (A13591601)
% UC San Diego MAE290B WI22 Final Project
%
% Input argument
% uIn - vector representation of spatial discretized u^[n]
% nIn - the latest timepoint at which a solution is obtained
% tStep - timestep size
% pXX - the matrix I + beta * A_xx
% mXX - the matrix I - beta * A_xx
% pYY - the matrix I + beta * A_yy
% mYY - the matrix I - beta * A_yy
% funcQ - a function handle for the space discretized source term
%
% Output argiment
% uOutVec - u^[n+1] on each spatial grid points as a vector
%
% Call format
% uOutVec = rdeADI_step(uIn,nIn,tIn,tStep,pXX,mXX,pYY,mYY,funcQ)
%
% Function called
% tridiagSolve

% Perform ADI advancement based on if current step if even or odd
if mod(nIn,2) == 0
    uOutVec = evenStepADI(uIn,tIn,tStep,pXX,mXX,pYY,mYY,funcQ);
else
    uOutVec = oddStepADI(uIn,tIn,tStep,pXX,mXX,pYY,mYY,funcQ);
end % of if-statement

end % of function rdeADI_step

% Even step advancement
function uOutVec = evenStepADI(uIn,tIn,tStep,pXX,mXX,pYY,mYY,funcQ)

% Define the source term vectors
q0 = funcQ(tIn);
qH = funcQ(tIn+tStep/2);
q1 = funcQ(tIn+tStep);

% Define righthand side of ADI advancement scheme
r = pXX * pYY * uIn + ...
    (tStep/4) * pXX * (q0+qH) + ...
    (tStep/4) * mXX * (qH+q1);

% Solve for z = mYY * uOut from mXX * z = r
z = tridiagSolve(mXX,r);

% Solve for uOut from mYY * uOut = z
uOutVec = tridiagSolve(mYY,z);
    
end % of even step function

% Odd step advancement
function uOutVec = oddStepADI(uIn,tIn,tStep,pXX,mXX,pYY,mYY,funcQ)

% Define the source term vectors
q0 = funcQ(tIn);
qH = funcQ(tIn+tStep/2);
q1 = funcQ(tIn+tStep);

% Define righthand side of ADI advancement scheme
r = pXX * pYY * uIn + ...
    (tStep/4) * pYY * (q0+qH) + ...
    (tStep/4) * mYY * (qH+q1);

% Solve for z = mYY * uOut from mXX * z = r
z = tridiagSolve(mXX,r);

% Solve for uOut from mYY * uOut = z
uOutVec = tridiagSolve(mYY,z);

end % of odd step function