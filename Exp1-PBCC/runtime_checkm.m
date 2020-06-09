% A quick peek at runtimes

datasets = {'Amazon_Amazon_Fashion','Amazon_Appliances','Cities_H','Newsgroups100','Zoo_h'};

%%
for j = 1:length(datasets)
    name = datasets{j};

    mu = 0.0;
    betas = 0.05:0.05:0.95;
    runs = zeros(length(betas),1);
    for i = 1:length(betas)
        beta = betas(i);
        load(strcat('Output/',name,'/beta_',num2str(beta),'_mu_0.0.mat'))
        runs(i) = runtime;
    end

    plot(runs)
    hold on
end

%%

figure(2)
for j = 2:length(datasets)
name = datasets{j};

mus = 0.01:0.01:0.2;
betas = 0.05:0.05:0.95;
beta = 0.5;
runs = zeros(length(betas),1);
for i = 1:length(mus)
    mu = mus(i);
    load(strcat('Output/',name,'/beta_',num2str(beta),'_mu_',num2str(mu),'.mat'))
    runs(i) = runtime;
end

plot(runs)
hold on
end