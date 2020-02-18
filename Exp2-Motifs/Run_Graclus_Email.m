%% First: Install and compile graclus 1.2
%    http://www.cs.utexas.edu/users/dml/Software/graclus.html
addpath('graclus1.2/matlab/')

%% Load the email dataset, as well as the largest connected
% component of the triangle-motif adjacency matrix
load emailEUcore
load EmailTriCore
A = Asim;
for i = 1
    A(i,i) = 0;
end
A = A(core_inds,core_inds);

%% Run Graclus

n = size(A,1);
ks = round(linspace(2,340,20));

GracC = zeros(n,numel(ks));

for t = 1:numel(ks)
    k = ks(t);
    [c, obj] = graclus(A, k);
    GracC(:,t) = c;
end

save('Output/Graclus_clusterings_2to340.mat','GracC')
