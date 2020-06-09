%% First: Install/addpath to outside algorithms
%    graclus 1.2: http://www.cs.utexas.edu/users/dml/Software/graclus.html
addpath('../include/graclus1.2/matlab/')
addpath('~/software/metis-5.0.2/metismex/')
addpath('../include/hmetis-1.5-osx-i686/')
addpath('../src')

%% Load the dataset, both the original graph...
load ../data/Florida_Bay_Dataset.mat

% ...made undirected
A = spones(A125+A125');
A = A- diag(diag(A));

% and the bifan graph
load ../data/Fl_bay_Afan.mat

%% Run algorithms on both graphs

n = size(A,1);
ks = round(linspace(2,50,20));

MetC = zeros(n,numel(ks));
MetFan = zeros(n,numel(ks));
times = zeros(numel(ks),2);

GracC = zeros(n,numel(ks));
GracFan = zeros(n,numel(ks));
times = zeros(numel(ks),2);

hmetC = zeros(n,numel(ks));
times = zeros(numel(ks),1);
name = 'flbay.hgr';
balance = 10;

for t = 1:numel(ks)
    k = ks(t);
    
    % Graclus
    tic
    [c, obj] = graclus(A, k);
    times(t,1) = toc;
    GracC(:,t) = c;
    
    tic
    [c, obj] = graclus(AFan, k);
    times(t,2) = toc;
    GracFan(:,t) = c;
    
    % Metis
    tic
    c = metismex('PartGraphRecursive',sparse(A),k);
    times(t,1) = toc;
    MetC(:,t) = c+1;
    
    tic
    c = metismex('PartGraphRecursive',sparse(AFan),k);
    times(t,2) = toc;
    MetFan(:,t) = c+1;
    
    % hmetis on the hypergraph 
    tic
    [c, time] = run_hmetis(name,balance,k);
    times(t) = toc;
    hmetC(:,t) = c;

end

save(strcat('Output/hmetis_2to50_',num2str(balance),'.mat'),'hmetC','times')
save('Output/Graclus_Flbay.mat','GracC','GracFan','times')
save('Output/Metis_Flbay.mat','MetC','times','MetFan')
