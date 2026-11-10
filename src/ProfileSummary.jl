module ProfileSummary
export summarize_profile
const CUTOFF=3.841458820694124
function side_crossing(rows, center, direction)
    i=center
    while 1<=i+direction<=length(rows)
        j=i+direction; a,b=rows[i],rows[j]
        if a.twice_delta_nll<=CUTOFF && b.twice_delta_nll>CUTOFF
            good=a.stopping_met && b.stopping_met
            fraction=(CUTOFF-a.twice_delta_nll)/(b.twice_delta_nll-a.twice_delta_nll)
            value=exp(log(a.fixed_value)+fraction*log(b.fixed_value/a.fixed_value))
            return (estimate=good ? value : missing, bracket_low=min(a.fixed_value,b.fixed_value),bracket_high=max(a.fixed_value,b.fixed_value),status=good ? "approximate crossing" : "unresolved convergence at crossing")
        end
        i=j
    end
    (estimate=missing,bracket_low=missing,bracket_high=missing,status="no crossing within searched bounds")
end
function summarize_profile(input,mle)
    rows=sort(collect(input);by=r->r.fixed_value)
    center=argmin(abs.(log.([r.fixed_value/mle for r in rows])))
    lower=side_crossing(rows,center,-1); upper=side_crossing(rows,center,1)
    better=minimum(r.twice_delta_nll for r in rows)<-2e-5
    (lower=lower,upper=upper,better_minimum=better,failed_points=count(r->!r.stopping_met,rows),minimum_twice_delta=minimum(r.twice_delta_nll for r in rows))
end
end
