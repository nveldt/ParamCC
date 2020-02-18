## Plot results from email experiment
Val_types = ["purity", "prob", "nmi", "ari"]

mat = matread("emailEUcore.mat")
label = vec(round.(Int64,mat["labels"]))

mat3  = matread("EmailTriCore.mat")
coreinds = mat3["core_inds"]
label = label[coreinds]

mat = matread("Output/Simple_Email_clusterings.mat")
Lams = mat["Lams"]
# Get clusters from running triangle motif clustering...
triC = mat["taC"]
# And standard Louvain clustering with resolution parameter lambda
SimC = mat["SimC"]

grac = matread("Output/Graclus_clusterings_2to340.mat")
GracC = grac["GracC"]

# Check ARI scores
numpoints = length(Lams)
SimARI = zeros(numpoints)
TriARI = zeros(numpoints)
SimSize = zeros(numpoints)
TriSize = zeros(numpoints)
GracARI= zeros(numpoints)
GracSize = zeros(numpoints)

include("../src/ClusterValidation.jl")
for i = 1:length(Lams)
    lam = Lams[i]
    cSim = round.(Int64,SimC[:,i])
    SimARI[i] = ari(cSim,label)
    SimSize[i] = maximum(cSim) #length(findall(x->x>0,cSim))

    cTri = round.(Int64,triC[:,i])
    TriARI[i] = ari(cTri,label)
    TriSize[i] = maximum(cTri) #length(findall(x->x>0,cTri))

    cGrac = round.(Int64,GracC[:,i])
    GracARI[i] = ari(cGrac,label)
    GracSize[i] = maximum(cGrac) #length(findall(x->x>0,cTri))

end

##
using Plots

lw = 3
msw = 0
s1 = 300
s2 = 250
ms = 5
lw = 2

p1 = sortperm(SimSize)
p2 = sortperm(TriSize)

rightmax = minimum( [maximum(SimSize), maximum(TriSize)])
val = 4     # cluster based on ARI score
p=plot(TriSize[p2], TriARI[p2],color=:green,linewidth= lw,label = "HyperLam",
markerstrokewidth = msw, markersize = ms, markershape = :circle, legend = :topright)

plot!(SimSize[p1], SimARI[p1],color=:black,linewidth= lw,
xlim = [0, rightmax],grid = false,
xaxis = "Number of Clusters", yaxis= "ARI Scores", label = "LamLouvain",
markerstrokewidth = msw, markersize = ms, markershape = :circle)

plot!(GracSize, GracARI,color=:blue,linewidth= lw, label = "Graclus",
markerstrokewidth = msw, markersize = ms, markershape = :circle, size = (s1,s2))
savefig(p,"Plots/Email-Plot.pdf")
