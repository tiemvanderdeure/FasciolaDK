### 
# Investigate if prevalence of liver fluke is within each cohort is
# correlated with the prevalence of other liver diseases.

using Turing, Mooncake

@model function binomial_correlation_lkj(
    x_a, x_b, n, slagteri_id, besid;
    nbes = maximum(besid), nslagteri = maximum(slagteri_id)
)
    N = length(x_a)

    # population-level means on logit scale (2 dims)
    μ ~ MvNormal(2, 5.0)

    # separate standard deviations for each grouping level (one per outcome dimension)
    σ_bes ~ filldist(truncated(Normal(0, 10); lower = 0), 2)
    σ_slag ~ filldist(truncated(Normal(0, 10); lower = 0), 2)

    # LKJ correlation for bes grouping (correlation only)
    L_bes ~ LKJCholesky(2, 2.0)
    Z_bes ~ filldist(Normal(0, 1), 2, nbes)      # 2 x nbes

    # LKJ correlation for slagteri grouping
    L_slagteri ~ LKJCholesky(2, 2.0)
    Z_slagteri ~ filldist(Normal(0, 1), 2, nslagteri)  # 2 x nslagteri

    # correlated unit-variance group effects (2 x n)
    U_bes = L_bes.L * Z_bes
    U_slag = L_slagteri.L * Z_slagteri

    # select per-observation columns and scale by the appropriate sigma
    # result: 2 x N matrix, column i is μ + σ_bes .* U_bes[:, besid[i]] + σ_slag .* U_slag[:, slagteri_id[i]]
    logitp = μ .+ σ_bes .* U_bes[:, besid] .+ σ_slag .* U_slag[:, slagteri_id]

    # Observations — ensure `n` has length N
    x_a ~ arraydist(BinomialLogit.(n, logitp[1, :]))
    x_b ~ arraydist(BinomialLogit.(n, logitp[2, :]))

    return
end

# read besid, slagteri_id, count and positive columns
cohorts = CSV.read(joinpath(datadir, "..", "predictors.csv"), DataFrame; 
    select = [:BES_ID, :slagteri_id, :count, :positive, :positive_other, :yr])

cohorts.besid = StatsBase.denserank(cohorts.BES_ID)
cohorts.slagteri_id = StatsBase.denserank(cohorts.slagteri_id)

# group by slagteri_id and besid to get unique cohorts
cohorts = combine(groupby(cohorts, [:slagteri_id, :besid]), 
    :count => sum => :count,
    :positive => sum => :positive,
    :positive_other => sum => :positive_other
)

# group by slagteri_id and besid to get unique cohorts
by_slagter = combine(groupby(cohorts, [:slagteri_id]), 
    :count => sum => :count,
    :positive => sum => :positive,
    :positive_other => sum => :positive_other,
    :yr => extrema => :yr_range
)
by_slagter.pct_positive_other = by_slagter.positive_other ./ by_slagter.count
by_slagter.pct_positive = by_slagter.positive ./ by_slagter.count

liverdisease_cor_model = binomial_correlation_lkj(
    cohorts.positive, cohorts.positive_other, cohorts.count, 
    cohorts.slagteri_id, cohorts.besid
)
chn = sample(liverdisease_cor_model, NUTS(; adtype = AutoMooncake(nothing)), 1000) # Takes ~1.5 hours

corcoef = chn[Symbol("L_bes.L[2, 1]")]
corslagtcoef = chn[Symbol("L_slagteri.L[2, 1]")]

quantile(corslagtcoef, [0.025, 0.975])
quantile(corcoef, [0.025, 0.975])

mean(corcoef)
mean(corslagtcoef)

open("liver_disease_cor.txt", "a") do io
    for l in corcoef
        println(io, l)
    end
end
open("liver_disease_cor_slagteri.txt", "a") do io
    for l in corslagtcoef
        println(io, l)
    end
end


corcoef = parse.(Float64,readlines("liver_disease_cor.txt"))
mean(corcoef)
quantile(corcoef, [0.05, 0.95])

using Printf
open("in_text_numbers.txt", "a") do io
    println(io, @sprintf "Liver disease correlation herds: %.2f (95%% CI: %.2f, %.2f)" mean(corcoef) quantile(corcoef, [0.05, 0.95])...)
    println(io, @sprintf "Liver disease correlation slagteri: %.2f (95%% CI: %.2f, %.2f)" mean(corslagtcoef) quantile(corslagtcoef, [0.05, 0.95])...)

end


#=
## Test this model
Sigma = [1 -0.5; -0.5 1]
n = 500

logitab = rand(filldist(MvNormal(Sigma), n))
ns = rand(50:100, n)
x_a = rand.(BinomialLogit.(ns, logitab[1,:] .* 2))
x_b = rand.(BinomialLogit.(ns, logitab[2,:]))

testmod = binomial_correlation_lkj(x_a, x_b, ns)
testchn = sample(testmod, NUTS(; adtype = AutoMooncake(nothing)), 1000)

testchn[Symbol("L.L[2, 1]")] |> mean
=#