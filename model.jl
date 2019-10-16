import Pkg; Pkg.activate(".")

using StatsBase
using EcologicalNetworks
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


ids = getfield.(filter(x -> occursin("Hadfield", x.Reference), web_of_life()), :ID);
networks = convert.(BinaryNetwork, web_of_life.(ids));
mainland = reduce(union, networks)


islands = Island[]
immigrationrates = 10.0.^rand(Distributions.Uniform(-1.3, -0.3), 1000)
extinctionrates = 10.0.^rand(Distributions.Uniform(-1.3, -0.3), 1000)
for i in eachindex(extinctionrates)
    push!(islands, Island(extinctionrates[i], immigrationrates[i]))
end

filter!(i -> 0.0 < α(i) ≤ 2.0, islands)
islands = StatsBase.sample(islands, 200, replace=false)

function simulation(n::BipartiteNetwork, i::Island)
    S = simulate(n, i; steps=50)
    return i => [last(S)]
end

results = [simulation(mainland, i) for i in islands]


function results_to_table(results)
    mainland_spe = specificity(mainland)
    rows = NamedTuple[]
    for r in results
        for n in r.second
            if richness(n; dims=1) > 0
                local_spe = collect(values(specificity(n)))
                global_spe = [mainland_spe[s] for s in species(n; dims=1)]
                push!(rows, (
                    immigration = r.first.immigration,
                    extinction = r.first.extinction,
                    α = α(r.first),
                    hosts = richness(n; dims=2)./richness(mainland; dims=2),
                    parasites = richness(n; dims=1)./richness(mainland; dims=1),
                    nestedness = η(n),
                    connectance = connectance(n),
                    realized_specificity = local_spe,
                    potential_specificity = global_spe,
                    network = n,
                    links = links(n)/links(mainland)
                ))
            end
        end
    end
    return JuliaDB.table(rows, pkey=(:α))
end

df = results_to_table(results)


StatsPlots.@df df Plots.scatter(:α, :hosts, ylim=(0,1),
	legend=:bottomright, lw=3, c=:black, ls=:dash, lab="Hosts")
StatsPlots.@df df Plots.scatter!(:α, :parasites, c=:white, lab="Parasites")
StatsPlots.@df df Plots.scatter!(:α, :links, c=:grey, lab="Links")


StatsPlots.@df df Plots.scatter(:α, :nestedness, ylim=(0,1),
	xlim=(0,2), legend=:topleft, lab="Nestedness", c=:grey)
Plots.hline!([η(mainland)], c=:black, ls=:dot, lab="Mainland \\eta")


StatsPlots.@df df Plots.scatter(:α, :connectance, ylim=(0,1), xlim=(0,2), legend=:topleft, lab="Connectance", c=:grey)
Plots.hline!([connectance(mainland)], c=:black, ls=:dot, lab="Mainland Co")


StatsPlots.@df df Plots.scatter(:hosts, :nestedness, xlim=(0,1), ylim=(0,1), legend=:topleft, lab="Hosts", c=:white)
StatsPlots.@df df Plots.scatter!(:parasites, :nestedness, legend=:topleft, lab="Parasites", c=:black)
StatsPlots.@df df Plots.scatter!(:links, :nestedness, legend=:topleft, lab="Links", c=:grey)



StatsPlots.@df df Plots.scatter(:hosts, :parasites, legend=:topleft,
  lab="Links", c=:grey)
