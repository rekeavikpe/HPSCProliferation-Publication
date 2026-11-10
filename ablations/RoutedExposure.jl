module RoutedExposure
export ROUTES, LABELS, routed_context, response_labels
const ROUTES=Dict(
    "baseline"=>("none","none"), "shear_cell"=>("shear","none"),
    "shear_aggregate"=>("none","shear"), "shear_both"=>("shear","shear"),
    "edr_cell"=>("edr","none"), "edr_aggregate"=>("none","edr"),
    "edr_both"=>("edr","edr"), "full"=>("shear","edr"))
const LABELS=Dict("baseline"=>"Model 0", "shear_cell"=>"Shear: cell", "shear_aggregate"=>"Shear: aggregate",
    "shear_both"=>"Shear: both", "edr_cell"=>"EDR: cell", "edr_aggregate"=>"EDR: aggregate",
    "edr_both"=>"EDR: both", "full"=>"Model 6")

# Model 6 accepts a cell-growth exposure channel and an aggregate-loss channel.
# Routing changes inputs, not the archived equations or parameter storage names.
function routed_context(ctx,scenario;shared=false)
    cell,agg=ROUTES[scenario];cfg=deepcopy(ctx.config)
    original=deepcopy(cfg["bounds"])
    cfg["bounds"]["K_tau"]=copy(original[cell=="edr" ? "K_eps" : "K_tau"])
    cfg["bounds"]["K_eps"]=copy(original[agg=="shear" ? "K_tau" : "K_eps"])
    cell=="none" && (cfg["fixed"]["K_tau"]=.001)
    agg=="none" && (cfg["fixed"]["K_eps"]=.001)
    if shared && cell==agg && cell!="none"
        cfg["parameter_ties"]=Dict("K_eps"=>"K_tau")
    end
    conditions=Dict(c=>merge(q,(lookup=let lookup=q.lookup,cell=cell,agg=agg
        d->begin
            tau,edr=lookup(d)
            (cell=="shear" ? tau : cell=="edr" ? edr : 0.,
             agg=="shear" ? tau : agg=="edr" ? edr : 0.)
        end
    end,)) for (c,q) in ctx.conditions)
    merge(ctx,(config=cfg,conditions=conditions))
end
function response_labels(scenario;shared=false)
    cell,agg=ROUTES[scenario]
    shared && cell==agg && cell!="none" && return Dict("K_tau"=>cell=="edr" ? "K_{\\varepsilon}" : "K_{\\tau}")
    Dict("K_tau"=>"K_{$(cell=="edr" ? "\\varepsilon" : "\\tau"),X}",
         "K_eps"=>"K_{$(agg=="shear" ? "\\tau" : "\\varepsilon"),A}")
end
end
