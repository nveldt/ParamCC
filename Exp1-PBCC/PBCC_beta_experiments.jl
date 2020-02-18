## Solve the PBCC objective for a range of parameter settings on a suite of small graphs
include("../src/PBCC.jl")
using MAT
using SparseArrays

graphs = ["Cities_H", "Newsgroups100", "Zoo_H", "Amazon_Gift_Cards", "Amazon_Appliances", "Amazon_Amazon_Fashion"]

## Beta experiments
for graph = graphs
    M = matread("data/"*graph*".mat")
    H = M["H"]
    # For each dataset, solve for a range of beta values, and store X
    betas = 0.05:0.05:0.95
    mu = 0.0
    for beta = betas
        outputmat = "Output/"*graph*"/beta_$beta"*"_mu_$mu.mat"
        println(outputmat)
        start = time()
        X, LPbound = PBCC_LP(H,mu,beta,false,false)
        runtime = time()-start
        delta = 2/5
        C, Labels, obj = ManyRoundX(H,X,delta,mu,beta,5)
        x,y = Biclustering(H,Labels)
        ratio = obj/LPbound
        matwrite(outputmat,Dict("X"=>X,"LPbound"=>LPbound,"C"=>C,"Labels"=>Labels,
        "row_labels"=>x, "col_labels"=>y, "obj"=>obj,"ratio"=>ratio,"runtime"=>runtime))

    end

end
