function x = tridiagSolve(A,y)

% A function that solves the system of equations A * x = y using the Thomas
% Algorithm. Here, A must be an nxn tridiagonal matrix with entries on its
% main diagonol and an opposing pair of sub/super diagonals. In MATLAB
% convention, matrix A has entries on its 0th, +dth, and -dth diagonols
%
% Author: Daibo Zhang (A13591601)
% UC San Diego MAE290B WI22 Final Project
%
% Input argument
% A - a nonsingular, tridiagonal sparse matrix
% y - the RHS of the equation A * x = y
%
% Output argiment
% x - the solution to A * x = y
%
% Call format
% x = tridiagSolve(A,y)
%
% Function called
% none

% Find the size of the matrix
if size(A,1) == size(A,2)
    n = length(A);
else
    disp('Error in Thomas: input matrix A is not square \n');
    return;
end

% Preallocate output
x = zeros(n,1);

% Identify the diagonal numbers of A 
[B,ds] = spdiags(A);
d = ds(end);
bs = B(:,2);
cs = B(d+1:end,3);

% Preallocate vectors used in the compuation process
gs = zeros(n-d,1); % upper diagonol of the U matrix
rs = zeros(n,1);

% Find the first elements of gs and rs
gs(1:d) = cs(1:d) ./ bs(1:d);
rs(1:d) = y(1:d) ./ bs(1:d);

% Find each element of g and r
for i = d+1:n-d
    gs(i) = A(i,i+d) / (A(i,i)-A(i,i-d)*gs(i-d));
    rs(i) = (y(i)-A(i,i-d)*rs(i-d)) / (A(i,i)-A(i,i-d)*gs(i-d));
end

% Continue building the rest of r
for i = n-d+1:n
    rs(i) = (y(i)-A(i,i-d)*rs(i-d)) / (A(i,i)-A(i,i-d)*gs(i-d));
end

% Obtain the last elements of the solution x
x(n-d+1:n) = rs(n-d+1:n);

% Find the first elements of the solution x
for i = n-d:-1:1
    x(i) = rs(i) - gs(i) * x(i+d);
end

end