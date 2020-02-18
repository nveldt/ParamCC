## Solve the PBCC objective for a range of parameter settings on a suite of small graphs
using MAT
using SparseArrays
include("../src/PBCC.jl")
graphs = ["Cities_H", "Newsgroups100", "Zoo_H", "Amazon_Appliances", "Amazon_Amazon_Fashion"]

##
C = [1 0 0;
0 .75 0;
0 0 1;
0 .5 .5;
.5 0 .5;
1 .5 0;
0 0 0;
.75 .75 0;
0.1 0.5 0.6;
0 .75 .75;
.4 .5 .3]

mu = 0.0
s1 = 300
s2 = 250
ms = 5
lw = 2
y_label = "Ratio to LP Lower Bound"
using LaTeXStrings
x_label = L"\beta"
##
plot()
using Plots
betas = collect(0.05:0.05:0.95)
Plotdata = zeros(5,length(betas))
for ii = 1:length(graphs)
    graph = graphs[ii]
    M = matread("../data/"*graph*".mat")
    H = M["H"]
    y1 = zeros(length(beta))
    bestDs = zeros(length(betas))
    for i = 1:length(betas)
        beta = betas[i]
        outputmat = "Output/"*graph*"/beta_$beta"*"_mu_$mu.mat"

        mat = matread(outputmat)
        X = mat["X"]
        LPbound = mat["LPbound"]

        C, Labels, obj, bestD = RangeRoundX(H,X,mu,beta,10)
        x,y = Biclustering(H,Labels)
        ratio = obj/LPbound
        y1[i] = ratio
        bestDs[i] = bestD

    end
    Plotdata[ii,:] = y1
    @show bestDs
    plot!(betas,y1,linewidth = lw, markershape =:circle, grid = false,
    color = RGBA(C[ii,1],C[ii,2],C[ii,3],1), label = graph, legend = :bottomright,
    size = (s1,s2), markerstrokewidth = 0, markersize = ms,
    xlabel = x_label, ylabel = y_label)
end

##
matwrite("Plots/Plotdata_for_beta_bestround.mat", Dict("Plotdata"=>Plotdata))
plot!(legend = false)
ii = 1
annotate!(.65, 1.5, text("Cities",font(9,RGBA(C[ii,1],C[ii,2],C[ii,3],1))))
ii = 2
annotate!(.2, 2.3, text("News100",font(9,RGBA(C[ii,1],C[ii,2],C[ii,3],1))))
ii = 3
annotate!(.53, 2.5, text("Zoo",font(9,RGBA(C[ii,1],C[ii,2],C[ii,3],1))))
ii = 4
annotate!(.85, 1.9, text("Appliances",font(9,RGBA(C[ii,1],C[ii,2],C[ii,3],1))))
ii = 5
annotate!(.8, 1.1, text("Fashion",font(9,RGBA(C[ii,1],C[ii,2],C[ii,3],1))))

savefig("Plots/All_bestround.pdf")
#["Cities_H", "Newsgroups100", "Zoo_H", "Amazon_Appliances", "Amazon_Amazon_Fashion"]


##


# plot!(x,y,linewidth = lw, grid = false, markershape =:circle,
#             legend = :bottomleft, color = RGBA(C[i,1],C[i,2],C[i,3],1),
#             xaxis=:log10, label = "label $(labels[i])")
