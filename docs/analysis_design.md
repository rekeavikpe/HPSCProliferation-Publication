# Analysis design and random splitting

“Condition” is sufficient; the proposed T prefix was only a training shorthand. Group the density series at 40 rpm first (20k, 80k, 150k cells/mL), then the additional RPMs at 80k cells/mL (30 and 60 rpm). Preserve legacy IDs in the curated data and use `data/conditions.csv` to relabel presentation. The 500-mL culture is legacy C4, not C5.

There are 31 100-mL cultures across five conditions (4, 8, 6, 7, 6), not 116 independent bioreactors. The 116 values in the existing score are condition/day/state means. A random 60/40 assignment is possible, but assigning individual daily measurements would put portions of the same culture in both subsets. Instead, randomize whole culture IDs within each condition. Keep all states and days together. If shared experimental batches are later identified, assess whether batch grouping is needed.

The included random allocation rounds the calibration count within condition: 2/2, 5/3, 4/2, 4/3 and 4/2 calibration/selection cultures. This is 19/12 overall (~61/39). It is a proposed rounding rule, not a uniquely correct one. Selection cultures are not randomly cherry-picked until a desired ranking appears; ten seeds are fixed in advance.

There is no universal evidence that all-data BIC is superior to a random 60/40 split for this exact ODE system. The recommendation to retain all 100-mL cultures for the main calibration is motivated by the small number of independent cultures, weak parameters, and the loss of estimation information in a single split. Repeated random splits can assess sensitivity to which cultures are used. They do not create new independent data or remove uncertainty about scale transfer.

Relevant evidence:

- Varoquaux (2018), *Cross-validation failure: Small sample sizes lead to large error bars*, demonstrates large and underestimated uncertainty in small-sample cross-validation. Its numerical results are from neuroimaging and must not be transferred quantitatively to these bioreactor data. https://arxiv.org/abs/1706.07581
- Roberts et al. (2017), *Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure*, explains why validation partitions should reflect dependence and the intended prediction target. This supports grouping repeated culture measurements rather than scattering time points. https://doi.org/10.1111/ecog.02881
- Schwarz (1978), *Estimating the Dimension of a Model*, defines the likelihood/complexity basis of BIC. BIC is calculated from the maximized calibration likelihood, not from a held-out score at parameters estimated elsewhere. https://doi.org/10.1214/aos/1176344136
- Picard and Berk (1990), *Data Splitting*, reviews the tradeoffs of splitting in regression. https://doi.org/10.1080/00031305.1990.10475704

Proposed primary analysis: all-100-mL working-likelihood BIC, followed by profile likelihood and repeated random whole-culture 60/40 sensitivity checks. Alternative: select using a declared predictive score on a single random 40% subset, then refit on all 100-mL cultures. These are different decision rules and should not be switched after seeing which produces the preferred model. A final choice is still pending user confirmation; the reproduction run implements the existing all-100-mL analysis only.

Scale transfer, media-exchange evaluation and other-cell-line transfer remain retrospective because the data have already been seen. Other-cell-line adaptation is a distinct fitting analysis, not independent validation.
