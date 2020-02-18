## Solve the PBCC objective for a range of parameter settings on a suite of small graphs
include("PBCC.jl")
using MAT
using SparseArrays

graphs = ["Cities_H", "Newsgroups100", "Zoo_H", "Amazon_Appliances", "Amazon_Amazon_Fashion","Amazon_Gift_Cards"]

## Beta experiments
for graph = graphs
    M = matread("data/"*graph*".mat")
    H = M["H"]
    # For each dataset, solve for a range of beta values, and store X
    beta = 0.5
    mus = collect(0.01:0.01:0.2)
    for mu = mus
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
