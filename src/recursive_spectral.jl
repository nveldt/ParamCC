using MatrixNetworks

function partition_and_map_subgraph(A,subset,minsize)
  map = collect(subset)
  A1 = A[map,map]
  sets = recursive_spectral(A1,minsize)
  for s in sets
    for i = 1:length(s)
      s[i] = map[s[i]]
    end
  end
  return sets
end

function recursive_spectral(A,minsize)
  n = size(A,1)
  if n < minsize
    sets = Vector{Vector{Int}}()
    push!(sets,Vector{Int}(1:size(A,1)))
    return sets
  else
    sc = spectral_cut(A)
    s1 = sc.set
    s2 = setdiff(1:n, s1)
    println("partitioning $n nodes into $(length(s1)) and $(length(s2)),lam $(sc.lam2) cond $(minimum(sc.sweepcut_profile.conductance))\n")
    # this is always in the largest component
    sets = partition_and_map_subgraph(A,s1,minsize)
    s2 = partition_and_map_subgraph(A,s2,minsize)
    for s in s2
      push!(sets,s)
    end
    return sets
  end
end

function set2clusvec(sets,n)
  c = zeros(n)
  next = 1
  for s in sets
    c[s] .= next
    next += 1
  end
  return c
end
## This is just a test
using MAT
mat = matread("../data/emailEUcore.mat")
A = mat["A"]            # Largest connected component of email EU
ASim = mat["Asim"]      # Undirected version of same graph: ASim = spones(A+A')

sets = recursive_spectral(ASim,100)

n = size(A,1)
c = set2clusvec(sets,n)

@show maximum(c)
