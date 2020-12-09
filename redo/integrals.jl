function integrals(sett,bas)

  @lints begin

    using Printf
  
    t_coord = Array{Array{Float64,1}}
    t_coord = [mol.coords[:,i]*PhysConsts.bohr_to_angstrom for i in 1:size(mol.coords,2)]
  
    lintmol = Lints.Molecule(mol.atchrg,t_coord)
  
    timingstring=@elapsed bas    = Lints.BasisSet(sett["bas"],lintmol)
    @printf("Time needed to construct basis:         %.4f s\n",timingstring)
    
    #build S, T and V in AO Basis (for core guess)
    S = Lints.make_S(bas)
    T = Lints.make_T(bas)
    V = Lints.make_V(bas)
  
    nao = size(S,1)
    nmo = size(S,1)
    nvir = nmo - nocc
  
    timingstring=@elapsed bas_df = Lints.BasisSet(sett["dfbas"],lintmol)
    @printf("Time needed to construct DF basis:      %.4f s\n",timingstring)
  
    #generate AO integtals and AO-DF integrals (P|Q), (P|mn), (mn|kl)
    timingstring=@elapsed PQ  = Lints.make_ERI2(bas_df)
    @printf("Time needed to construct ERI2:          %.4f s\n",timingstring)
    timingstring=@elapsed Pmn = Lints.make_ERI3(bas,bas_df)
    @printf("Time needed to construct ERI3:          %.4f s\n",timingstring)
  
    if sett["rijk"] == "false"
      timingstring=@elapsed mnkl = Lints.make_ERI4(bas)
      @printf("Time needed to construct ERI4:          %.4f s\n",timingstring)
    end #if leri4
  
    naux = size(Pmn,1)
  
    #build inverse  (P|Q)^{-1/2}
    PQh = zeros(naux,naux)
    timingstring=@elapsed PQh = PQ^(-1/2)
    PQh = real(Matrix(PQh))
    @printf("Time needed to construct (P|Q)^{-1/2}:  %.4f s\n",timingstring)
    Bmn = zeros(naux,nao,nao)
    timingstring=@elapsed @tensor Bmn[Q,m,n] = Pmn[P,m,n] * PQh[P,Q]
    @printf("Time needed to construct (P|mn):        %.4f s\n",timingstring)
  
    if lDFdebug == true
      eri4 = zeros(nao,nao,nao,nao)
      timingstring=@elapsed @tensor eri4[m,n,k,l] = Bmn[Q,m,n]*Bmn[Q,k,l]
      @printf("Time needed to construct (mn|kl) from (P|mn):  %.4f s\n",timingstring)
    end
    #
    #Fill AOintegrals Dict with Pointers to the ao integrals
    push!(AOint,"S" => S)
    push!(AOint,"T" => T)
    push!(AOint,"V" => V)
    push!(AOint,"B" => Bmn)
    if (sett["rijk"] == "false")
      push!(AOint,"ERI4" => mnkl)
    end
    

  @lints end

  return AOint