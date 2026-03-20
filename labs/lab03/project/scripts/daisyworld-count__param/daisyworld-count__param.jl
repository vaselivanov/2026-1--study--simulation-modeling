using DrWatson
@quickactivate "project"
using Agents
using DataFrames
using Plots
using CairoMakie
include(srcdir("daisyworld.jl"))
black(a) = a.breed == :black
white(a) = a.breed == :white
adata = [(black, count), (white, count)]
# Параметры эксперимента
param_dict = Dict(
:griddims => (30, 30),
:max_age => [25, 40],
:init_white => [0.2, 0.8],
:init_black => 0.2,
:albedo_white => 0.75,
:albedo_black => 0.25,
:surface_albedo => 0.4,
:solar_change => 0.005,
:solar_luminosity => 1.0,
:scenario => :default,
:seed => 165,
)
# Создаём список всех комбинаций
params_list = dict_list(param_dict)
counter = 0
for params in params_list
model = daisyworld(;params...)
agent_df, model_df = run!(model, 1000; adata)
figure = Figure(size = (600, 400));
ax = figure[1, 1] = Axis(figure, xlabel = "tick", ylabel = "daisy count")
blackl = lines!(ax, agent_df[!, :time], agent_df[!, :count_black], color = :black)
whitel = lines!(ax, agent_df[!, :time], agent_df[!, :count_white], color = :orange)
Legend(figure[1, 2], [blackl, whitel], ["black", "white"], labelsize = 12)
global counter
counter += 1
plt_name = "daisy-count_$(counter).png"
save(plotsdir(plt_name), figure)
end
