%% First: Install and compile graclus 1.2
%    http://www.cs.utexas.edu/users/dml/Software/graclus.html
addpath('../Exp1-Email/graclus1.2/matlab/')

%% Load the email dataset, as well as the largest connected
% component of the triangle-motif adjacency matrix
load data/Florida_Bay_Dataset.mat

% Make it undirected
A = spones(A125+A125');
A = A- diag(diag(A));


%% Run Graclus

n = size(A,1);
ks = round(linspace(2,50,20));

GracC = zeros(n,numel(ks));

for t = 1:numel(ks)
    k = ks(t);
    [c, obj] = graclus(A, k);
    GracC(:,t) = c;
end

save('Output/Graclus_2to50.mat','GracC')
