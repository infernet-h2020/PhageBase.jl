for φ = -10. : 1. : 10.
    p = fermi_dirac_prob(φ)
    @test p ≈ 1 / (1 + exp(-φ))
    @test fermi_dirac_logp(φ) ≈ log(p)
    @test fermi_dirac_1mp(φ) ≈ 1 - p
    @test fermi_dirac_l1mp(φ) ≈ log(1 - p)
end

for f in (fermi_dirac_prob,
          fermi_dirac_logp,
          fermi_dirac_1mp,
          fermi_dirac_l1mp)

    for μ = -10. : 1. : 10., E = -10. : 1. : 10.
        φ = μ - E
        @test f(μ, E) ≈ f(φ)
    end
end
