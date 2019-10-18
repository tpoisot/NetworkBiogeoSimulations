import Pkg; Pkg.activate(".")

@info "started"

using StatsBase
using EcologicalNetworks
using Distributions
using Plots
using StatsPlots
using DataFrames
using CSV

struct Island
	extinction::Float64   # Extinction rate
	immigration::Float64  # Immigration rate
end


extinction(island::Island) = island.extinction
immigration(island::Island) = island.immigration
α(island::Island) = island.immigration/island.extinction


function hosts_migration(mainland::T, island::Island) where {T <: AbstractEcologicalNetwork}
	return filter(sp -> rand() ≤ immigration(island), species(mainland; dims=2))
end

function sample_parasites_on_each_host(parasites, island::Island)
	return filter(sp -> rand() ≤ immigration(island), parasites)
end

function parasites_migration(mainland::T, island::Island, hosts) where {T <: AbstractEcologicalNetwork}
	d = degree(simplify(mainland[:,hosts]); dims=1)
	return collect(keys(filter(p -> rand() ≥ (1-immigration(island))^p.second, d)))
end

function simulate(mainland::T, island::Island; steps::Integer=200) where {T <: AbstractEcologicalNetwork}

	networks = typeof(mainland)[]
	for i in 1:steps
		# Pick mainland hosts and parasites
		migrating_hosts = hosts_migration(mainland, island)
		migrating_parasites = parasites_migration(mainland, island, migrating_hosts)
		if (@isdefined island_network)
			global island_network
			resident_hosts = species(island_network; dims=2)
			resident_parasites = species(island_network; dims=1)
			island_network = mainland[union(migrating_parasites, resident_parasites), union(migrating_hosts, resident_hosts)]
		else
			island_network = mainland[migrating_parasites, migrating_hosts]
		end
		# Extinctions
		remaining_hosts = filter(sp -> rand() ≥ extinction(island), species(island_network; dims=2))
		island_network = island_network[:,remaining_hosts]
		# Note that some parasites will have a degree of 0, but that's all right
		d_parasites = degree(island_network; dims=1)
		remaining_parasites = filter(sp -> rand() ≥ extinction(island)^d_parasites[sp], species(island_network; dims=1))
		island_network = island_network[remaining_parasites,:]

		push!(networks, copy(island_network))
	end
	return networks
end

function simulation(n::BipartiteNetwork, i::Island)
	S = simulate(n, i; steps=200)
	return i => [last(S)]
end

@info "generating networks"
ids = getfield.(filter(x -> occursin("Hadfield", x.Reference), web_of_life()), :ID);
networks = convert.(BinaryNetwork, web_of_life.(ids));
mainland = reduce(union, networks)

@info "generating islands"
islands = Island[]
immigrationrates = 10.0.^rand(Distributions.Uniform(-1.3, -0.3), 1000)
extinctionrates = 10.0.^rand(Distributions.Uniform(-1.3, -0.3), 1000)
for i in eachindex(extinctionrates)
	push!(islands, Island(extinctionrates[i], immigrationrates[i]))
end

@info "subsampling islands"
filter!(i -> -1.0 ≤ log10(α(i)) ≤ 1.0, islands)
n_islands = min(200, length(islands))
islands = StatsBase.sample(islands, n_islands, replace=false)

@info "simulations"
results = [simulation(mainland, i) for i in islands]

function results_to_table(results)
    mainland_spe = specificity(mainland)
    df = DataFrame(
		immigration = Float64[],
		extinction = Float64[],
		hosts = Float64[],
		parasites = Float64[],
		η = Float64[],
		ρ = Float64[],
		connectance = Float64[],
		links = Float64[]
	)
    for r in results
        for n in r.second
            if richness(n; dims=1) > 0
				push!(df, (
					r.first.immigration,
					r.first.extinction,
					richness(n; dims=2)./richness(mainland; dims=2),
					richness(n; dims=1)./richness(mainland; dims=1),
					η(n),
					ρ(n),
					connectance(n),
					links(n)/links(mainland)
				))
            end
        end
    end
    return df
end

@info "analysis"
df = results_to_table(results)
df[:α] = df[:immigration]./df[:extinction]

@info "plotting"
StatsPlots.@df df Plots.scatter(:α, :hosts, legend=:bottomright, c=:black, lab="Hosts", frame=:box)
StatsPlots.@df df Plots.scatter!(:α, :parasites, c=:grey, m=:square, msw=0, lab="Parasites")
StatsPlots.@df df Plots.scatter!(:α, :links, c=:darkgrey, m=:utriangle, msw=0, lab="Links")
Plots.xaxis!(:log10, (0.1,10.0), "\\alpha")
Plots.yaxis!((0,1), "Relative number")
Plots.savefig(joinpath("figures", "richness.png"))

StatsPlots.@df df Plots.scatter(:α, :η, legend=:topleft, lab="", mc=:white, frame=:box)
Plots.hline!([η(mainland)], c=:grey, lw=2, ls=:dash, lab="")
Plots.xaxis!(:log10, (0.1,10.0), "\\alpha")
Plots.yaxis!((0,1), "Nestedness\\eta")
Plots.savefig(joinpath("figures", "nestedness.png"))

StatsPlots.@df df Plots.scatter(:α, :ρ, legend=:topleft, lab="", mc=:white, frame=:box)
Plots.hline!([ρ(mainland)], c=:grey, lw=2, ls=:dash, lab="")
Plots.xaxis!(:log10, (0.1,10.0), "\\alpha")
Plots.yaxis!((0,1), "Spectral radius\\rho")
Plots.savefig(joinpath("figures", "radius.png"))

StatsPlots.@df df Plots.scatter(:α, :connectance, legend=:topleft, lab="", mc=:white, frame=:box)
Plots.hline!([connectance(mainland)], c=:grey, lw=2, ls=:dash, lab="")
Plots.xaxis!(:log10, (0.1,10.0), "\\alpha")
Plots.yaxis!((0,1), "Connectance")
Plots.savefig(joinpath("figures", "connectance.png"))

@info "ended"
