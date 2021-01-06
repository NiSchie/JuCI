include("juci.jl")

using Printf
using TensorOperations
using LinearAlgebra

function rhf()
  sett = juci.set_config()
  Ethr = 1.0 * 10^(-parse(Float64,sett["econv"]))
  Dthr = 1.0 * 10^(-parse(Float64,sett["denconv"]))
  dampfac = parse(Float64,sett["hfdampfac"])
  dampstp = parse(Float64,sett["hfdampstp"])
  dampend = parse(Float64,sett["hfdampend"])

  mol = juci.read_mol(sett["molfile"])
  
  AOint,dims = juci.get_ints(mol,sett)
  AOops = build_AOops(AOint)

  #hcore guess
  P = calc_P(AOops,AOint,dims)
  Pold = zeros(dims["nao"],dims["nao"])
  Eguess = 0.5 * sum( P .* (AOops["H"] .+ AOops["F"] ) ) + mol.Vnuc
  println("Hcore guess energy: ",Eguess)

  E = zeros(Float64,0)
  append!(E,Eguess)

  @printf("\n")
  @printf("HF It.  |  HF-Energy | abs(dE)   | abs(dD)  | dampfac\n")
  @printf("-----------------------------------------\n")
  for ite = 2:parse(Int,sett["hfmaxit"])+1

    push!(AOops,"Fold" => AOops["F"])
    #Build the Fock matrix
    Cocc = AOops["C"][:,1:dims["nocc"]]
    tstring=@elapsed @tensor begin
      Jti[P]    := P[k,l]   * AOint["B"][P,k,l]
      J[m,n]    := AOint["B"][P,m,n] * Jti[P]
      K1[Q,m,c] := Cocc[r,c] * AOint["B"][Q,m,r]
      K[m,n]    := K1[Q,m,p] * K1[Q,n,p]
      F[m,n]    := AOops["H"][m,n] + J[m,n] - K[m,n]
    end

    #damping
    F = dampfac * F .+ (1.0-dampfac) * AOops["Fold"]
    dampfac = min(dampend,dampfac+dampstp)

    push!(AOops,"F" => F)
    Pold = P
    P = calc_P(AOops,AOint,dims)
    dD = P .- Pold
    dP = sqrt(sum(dD.^2))
    
    append!(E,0.5 * sum(P .* ( AOops["H"] .+ AOops["F"] ) ) + mol.Vnuc)
    dE = abs(E[ite]-E[ite-1])
    @printf("  %3d   |  %.5f | %.2e | %.2e | %.1f\n",ite-1,E[ite],abs(dE),dP,dampfac)

    if ( (abs(dE) < Ethr) & (dP < Dthr) & (ite > 5) )
      @printf("HF is converged!\n")
      break
    end #if convergence
     
  end #ite
  
end #function rhf


function build_AOops(AOint)
  AOops = Dict()

  #push!(AOops,"Sh" => AOint["S"]^(1/2))
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
  eps,Ct = eigen(Ft)
  
  C = AOops["X"]*Ct
  push!(AOops,"C" => C)

  Cocc = C[:,1:dims["nocc"]]
  P    = zeros(dims["nao"],dims["nao"])
  @tensor P[m,n] = 2.0 * Cocc[m,p] * Cocc[n,p]

  return P
end #function core_guess

