%% First: Install/addpath to outside algorithms
%    graclus 1.2: http://www.cs.utexas.edu/users/dml/Software/graclus.html
addpath('../include/graclus1.2/matlab/')
addpath('~/software/metis-5.0.2/metismex/')
addpath('../include/hmetis-1.5-osx-i686/')
addpath('../src')

%% Load the email dataset, as well as the largest connected
% component of the triangle-motif adjacency matrix
load ../data/emailEUcore
load ../data/EmailTriCore
A = Asim;
for i = 1
    A(i,i) = 0;
end
A = A(core_inds,core_inds);

load ../data/EmailTriangle.mat
%% Run Graclus

n = size(A,1);
ks = round(linspace(2,340,20));

GracC = zeros(n,numel(ks));
GracCtri = zeros(n,numel(ks));
times = zeros(numel(ks),2);

MetC = zeros(n,numel(ks));
MetCtri = zeros(n,numel(ks));
times = zeros(numel(ks),2);

hmetC = zeros(n,numel(ks));
times = zeros(numel(ks),1);
name = 'email.hgr';
balance = 10;

for t = 1:numel(ks)
    k = ks(t);
    
    % Metis
    tic
    c = metismex('PartGraphRecursive',sparse(A),k);
    times(t,1) = toc;
    MetC(:,t) = c+1;
    
    tic
    c = metismex('PartGraphRecursive',sparse(Ata),k);
    times(t,2) = toc;
    MetCtri(:,t) = c+1;
    
    % Graclus
    tic
    [c, ~] = graclus(A, k);
    times(t,1) = toc;
    GracC(:,t) = c;
    
    tic
    [c, ~] = graclus(Ata, k);
    times(t,2) = toc;
    GracCtri(:,t) = c;
    
    % hmetis
    k = ks(t);
    [c, time] = run_hmetis(name,balance,k);
    times(t) = time;
    hmetC(:,t) = c;
end

save('Output/Metis_2to340.mat','MetC','MetCtri','times')
save('Output/Graclus_clusterings_2to340.mat','GracC','GracCtri','times')
save(strcat('Output/hmetis_2to340_',num2str(balance),'.mat'),'hmetC','times')
