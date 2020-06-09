using MAT
using LinearAlgebra
using SparseArrays
include("../src/Graph_Louvain.jl")
mat = matread("../data/emailEUcore.mat")
A = mat["A"]            # Largest connected component of email EU
ASim = mat["Asim"]      # Undirected version of same graph: ASim = spones(A+A')
m = sum(A)
label = vec(round.(Int64,mat["labels"]))
include("../src/Motif_Adj_Mat.jl")

Ata = ExtractTriMat(ASim)
Acc, p = largest_component(Ata)
core_inds = findall(x->x==true,p)
ASim = ASim[p,p]

## Performs weighted triangle clique expansion on an already undirected graph
H = ExtractTriangle_List(ASim)
n = size(ASim,1)
m = size(H,1)
open("../data/email.hgr","w") do f
    write(f, "$m $n\n")
    for i = 1:m
        write(f,"$(H[i,1]) $(H[i,2]) $(H[i,3])\n")
    end
end


## Do something similar for the floriday bay food web
mat = matread("../data/Florida_Bay_Dataset.mat")
A125 = mat["A125"]

H = ExtractFanList(A125)
n = size(A125,1)
m = size(H,1)
open("../data/flbay.hgr","w") do f
    write(f, "$m $n\n")
    for i = 1:m
        write(f,"$(H[i,1]) $(H[i,2]) $(H[i,3]) $(H[i,4])\n")
    end
end
