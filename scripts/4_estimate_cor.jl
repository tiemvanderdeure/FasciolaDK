### 
# Investigate if prevalence of liver fluke is within each cohort is
# correlated with the prevalence of other liver diseases.

using Turing, Mooncake

@model function binomial_correlation_lkj(x_a, x_b, n)
    N = length(x_a)

    # Mean vector for logit(p_a), logit(p_b)
    μ ~ MvNormal(2, 5.0)

    # Standard deviations
    σ ~ filldist(truncated(Normal(0, 10); lower = 0), 2)

    # LKJ prior for correlation matrix
    L ~ LKJCholesky(2, 2.0)  # 2 dimensions, η=2 is weakly informative

    # Z-scores per sample and per dimension
    Z ~ filldist(Normal(0, 1), 2, N)

    # calculate resulting probabilities
    logitp = σ .* L.L * Z .+ μ

    # Observations
    x_a ~ arraydist(BinomialLogit.(n, logitp[1,:]))
    x_b ~ arraydist(BinomialLogit.(n, logitp[2,:]))

    return (; cor_coef = (L.L * L.U)[1,2]) # DON'T USE RETURNED!
end

liverdisease_cor_model = binomial_correlation_lkj(cohorts.positive, cohorts.positive_other, cohorts.count)
chn = sample(liverdisease_cor_model, NUTS(; adtype = AutoMooncake(nothing)), 1000) # Takes ~1.5 hours

corcoef = chn[Symbol("L.L[2, 1]")]

open("liver_disease_cor.txt", "a") do io
    for l in corcoef
        println(io, l)
    end
end

corcoef = parse.(Float64,readlines("liver_disease_cor.txt"))
mean(corcoef)
quantile(corcoef, [0.05, 0.95])


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