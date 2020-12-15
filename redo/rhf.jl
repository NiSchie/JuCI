include("juci.jl")

function rhf()
  sett = juci.set_config()
  mol = juci.read_mol(sett["molfile"])
  
  AOint,dims = juci.get_ints(mol,sett)
  AOops = build_AOops(AOint)

  #hcore guess
  P = calc_P(AOops,AOint,dims)
  Eguess = 0.5 * sum( P .* (AOops["H"] .+ AOops["F"] ) )
  println("Hcore guess energy: ",Eguess)

  Bmn = AOint["B"]
  C   = AOops["C"]
  H   = AOops["H"]

  @printf("\n")
  @printf("HF It.  |  HF-Energy  |abs(dE)   |abs(dD)\n")
  @printf("-----------------------------------------\n")
  for ite = 1:parse(Int,sett["hfmaxit"])

    #Build the Fock matrix
    Cocc = C[:,1:dims["nocc"]]
    tstring=@elapsed @tensor begin
      Jti[P]    := P[k,l]   * Bmn[P,k,l]
      J[m,n]    := Bmn[P,m,n] * Jti[P]
      K1[Q,m,c] := Cocc[r,c]  * Bmn[Q,m,r]
      K[m,n]    := K1[Q,m,p]  * K1[Q,n,p]
      F[m,n]    := H[m,n] + J[m,n] - K[m,n]
    end
    push!(AOops,"F" => F)
    P = calc_P(AOops,AOint,dims)
    #println("Time needed to build Fock matrix: ",tstring)
    
    E = 0.5 * sum(P .* ( AOops["H"] .+ AOops["F"] ) )
    println("HF energy: ",E)
     
  end #ite

end #function rhf


function build_AOops(AOint)
  AOops = Dict()

  push!(AOops,"Sh" => AOint["S"]^(1/2))
  push!(AOops,"H" => AOint["T"] .+ AOint["V"])
  push!(AOops,"F" => AOops["H"])

  #X matrix (canonical orthogonalization)
  s,U = eigen(AOint["S"])
  X = zeros(size(U,1),size(U,2))
  for i in 1:size(s,1)
    s[i] = 1.0/sqrt(s[i])
  end
  for i in 1:size(U,1), j in 1:size(U,2)
    X[i,j] = U[i,j] * s[j]
  end #i,j
  push!(AOops,"X" => X)

  return AOops
end #function build_AOops


function calc_P(AOops,AOint,dims)

  Ft = transpose(AOops["X"])*AOops["F"]*AOops["X"]
  eps,Ct = eigen(Hermitian(real.(Ft)))
  C = AOops["X"]*Ct
  push!(AOops,"C" => C)

  Cocc = C[:,1:dims["nocc"]]
  P    = zeros(dims["nao"],dims["nao"])
  @tensor P[m,n] = 2.0 * Cocc[m,p] * Cocc[n,p]

  return P
end #function core_guess

