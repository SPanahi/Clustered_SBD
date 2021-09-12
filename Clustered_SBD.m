%% SBD for a matrix A with clusters c
%% INPUTS:
%%    A is an n-by-n symmetric matrix
%%    c is a vector of length n with the cluster index for each element 1:n
%% OUTPUTS:
%%    P is the transformation matrix such that P'*A*P = B is block diagonal
%%    U the commuting matrix

function [P,U,B] = Clustered_SBD (A,c)

% Sort the clusters into adjacent spots
[b,perm] = sort(c);
A = A(perm,perm);

% Number of clusters
nclus = max(b);
n = size(A,1);

% Get mini-vectors to select independent clusters and store the population of each cluster
for j = 1 : nclus
    E{j} = find(b==j);
    pop(j) = length(E{j});
end

% Determine the dimensions of the matrix
nrow = 0;
ncol = 0;

for i = 1 : nclus
    for j = 1 : nclus
        nrow = nrow + pop(i)*pop(j);
    end
    ncol = ncol + pop(i)*pop(i);
end

% Allocate memory for the matrix to represent the individual pseudo-commutators
G = zeros(nrow,ncol);

% Add the true commutator portion of the matrix G
idx = 1;
for j = 1 : nclus
    G(idx:idx+pop(j)*pop(j)-1,idx:idx+pop(j)*pop(j)-1) = kron(A(E{j},E{j}),eye(pop(j))) - kron(eye(pop(j)),A(E{j},E{j}));
    idx = idx + pop(j)*pop(j);
end

% Store the starting indices of each set of columns
start = 1;
for j = 2:nclus+1
    start(j) = start(j-1) + pop(j-1)*pop(j-1);
end

% Add the pseudo-commutator portions
for i = 1 : nclus
    for j = 1:nclus
        if (i ~= j)
            G(idx:idx+pop(i)*pop(j)-1,start(i):start(i+1)-1) = kron(A(E{j},E{i}),eye(pop(i)));
            G(idx:idx+pop(i)*pop(j)-1,start(j):start(j+1)-1) = -kron(eye(pop(j)),A(E{i},E{j}));
            idx = idx + pop(i)*pop(j);
        end
    end
end

% Create the ncol-by-ncol matrix T
T = G'*G;

% Compute its eigenvalues
[V,D] = eig(T);
D = diag(D);

% Create the vector of vectorized P_j matrices
Pvec = zeros(ncol,1);

for j = 1 : ncol
    if (D(j) < 1E-12)
        Pvec = Pvec + V(:,j);
    end
end

% With the created vector of vectorized P_j matrices, recreate U
U = zeros(n);
idx = 1;
for j = 1 : nclus
    U(idx : idx+pop(j)-1, idx:idx+pop(j)-1) = reshape(Pvec(start(j):start(j+1)-1),pop(j),pop(j));
    idx = idx + pop(j);
end

U = U + U';
[P,R] = eig(U);
B = P'*A*P;
end
