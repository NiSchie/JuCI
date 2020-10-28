module config

  using Printf

  export set_config,sett

  function set_config()

    sett_init = Dict(
      "rijk"    => "true",
      "econv"   => "6",
      "denconv" => "6",
      "molfile" => "mol.xyz",
      "bas"     => "cc-pVDZ",
      "dfbas"   => "cc-pVDZ-JKFIT"
    )

    #load settings
    sett = deepcopy(sett_init)
    if isfile("cfg.juci")
      @printf("Configuration file found.\n")
      for ln in eachline("cfg.juci")
        splitline = split(ln," ")
        filter!(e->e!="",splitline)
        println(splitline)
        sett[splitline[1]] = splitline[2]
      end
    else
      @printf("NO CONFIGURATION FILE FOUND, using standard config (HF calculation, 10^-6, cc-pVDZ).\n")
    end
    display(sett)
    @printf("\n")

    lDFdebug = false
    Ethr = 1.0 * 10^(-parse(Float64,sett["econv"]))
    Dthr = 1.0 * 10^(-parse(Float64,sett["denconv"]))

    @printf("\n")
    @printf("SCF Econv:   %.1e\n",Ethr)
    @printf("SCF Denconv: %.1e\n",Ethr)
    @printf("\n")

  end #set_config

end #module config
