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
x_label = L"\mu"
##

using Plots
plot()
mus = collect(0.00:0.01:0.2)
beta = 0.5
Plotdata = zeros(5,length(mus))
for ii = 1:length(graphs)
    graph = graphs[ii]
    M = matread("../data/"*graph*".mat")
    H = M["H"]
    y1 = zeros(length(mus))
    bestDs = zeros(length(mus))
    for i = 1:length(mus)
        mu = mus[i]
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
    plot!(mus,y1,linewidth = lw, markershape =:circle, grid = false,
    color = RGBA(C[ii,1],C[ii,2],C[ii,3],1), label = graph, legend = :bottomright,
    size = (s1,s2), markerstrokewidth = 0, markersize = ms,
    xlabel = x_label, ylabel = y_label)
end

##
matwrite("Plots/Plotdata_for_smaller_mu.mat", Dict("Plotdata"=>Plotdata))
plot!(legend = false)
# ii = 1
# annotate!(.65, 1.5, text("Cities",font(9,RGBA(C[ii,1],C[ii,2],C[ii,3],1))))
# ii = 2
# annotate!(.2, 2.3, text("News100",font(9,RGBA(C[ii,1],C[ii,2],C[ii,3],1))))
# ii = 3
# annotate!(.53, 2.5, text("Zoo",font(9,RGBA(C[ii,1],C[ii,2],C[ii,3],1))))
# ii = 4
# annotate!(.85, 1.9, text("Appliances",font(9,RGBA(C[ii,1],C[ii,2],C[ii,3],1))))
# ii = 5
# annotate!(.8, 1.1, text("Fashion",font(9,RGBA(C[ii,1],C[ii,2],C[ii,3],1))))

#["Cities_H", "Newsgroups100", "Zoo_H", "Amazon_Appliances", "Amazon_Amazon_Fashion"]

savefig("Plots/Mu_smaller_bestround.pdf")
