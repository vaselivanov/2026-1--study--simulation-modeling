using DrWatson
@quickactivate "project"
include(srcdir("SIRPetri.jl"))
using .SIRPetri
using DataFrames, CSV, Plots

β = 0.3
γ = 0.1
tmax = 100.0

net, u0, states = build_sir_network(β, γ)

group_names = string.(states)  # ["S", "I", "R"]

df_det = simulate_deterministic(net, u0, (0.0, tmax), saveat = 0.2, rates = [β, γ])


anim = @animate for row in eachrow(df_det)
    
    u = [row[col] for col in propertynames(row) if col != :time]
    
    
    plot(
        1:length(u),
        u,
        legend = false,
        ylims = (0, 1000),                   
        xlabel = "Compartment",
        ylabel = "Population",
        title = "SIR dynamics at t = $(round(row.time, digits=2))",
        markersize = 4,                       
        markerstrokewidth = 0,
        linewidth = 2,                        
        color = [:blue :red :green],
        marker = :circle,                     
        linestyle = :solid
    )
    
    xticks!(1:length(group_names), group_names, rotation = 0)
    
    plot!(grid = true, gridalpha = 0.3)
end

gif(anim, plotsdir("sir_animation.gif"), fps = 15)
println("Анимация сохранена в plots/sir_animation.gif")