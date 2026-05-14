using StableRNGs
using Distributions
using ConcurrentSim
using ResumableFunctions
using Plots
using Statistics

#set simulation parameters
rng = StableRNG(123)
num_customers = 10 # total number of customers generated

# set queue parameters
num_servers = 2 # number of servers
mu = 1.0 / 2 # service rate
lam = 0.9 # arrival rate
arrival_dist = Exponential(1 / lam) # interarrival time distribution
service_dist = Exponential(1 / mu) # service time distribution

# Массивы для сбора данных
arrival_times = Vector{Float64}(undef, num_customers)
service_start_times = Vector{Float64}(undef, num_customers)
service_end_times = Vector{Float64}(undef, num_customers)

# define customer behavior
@resumable function customer(
    env::Environment,
    server::Resource,
    id::Integer,
    t_a::Float64,
    d_s::Distribution,
)
    @yield timeout(env, t_a) # customer arrives
    arrival_times[id] = now(env)
    println("Customer $id arrived: ", now(env))
    @yield request(server) # customer starts service
    service_start_times[id] = now(env)
    println("Customer $id entered service: ", now(env))
    service_time = rand(rng, d_s)
    @yield timeout(env, service_time) # server is busy
    service_end_times[id] = now(env)
    println("Customer $id exited service: ", now(env))
    @yield unlock(server) # customer exits service
end

# setup and run simulation
function setup_and_run()
    sim = Simulation() # initialize simulation environment
    server = Resource(sim, num_servers) # initialize servers
    arrival_time = 0.0
    for i = 1:num_customers # initialize customers
        arrival_time += rand(rng, arrival_dist)
        @process customer(sim, server, i, arrival_time, service_dist)
    end
    run(sim) # run simulation
    
    # Создание графиков
    waiting_times = service_start_times .- arrival_times
    
    # 1. Диаграмма Ганта
    p1 = plot(title="Временная диаграмма обслуживания",
              xlabel="Время", ylabel="Клиент", legend=:topright, size=(800, 600))
    for i in 1:num_customers
        plot!([arrival_times[i], service_start_times[i], service_end_times[i]], 
              [i, i, i], marker=:circle, linewidth=2, label="Customer $i")
    end
    savefig(p1, "plots/gantt_chart.png")
    println("✓ plots/gantt_chart.png")
    
    # 2. Гистограмма времени ожидания
    p2 = bar(1:length(waiting_times), waiting_times, 
             title="Время ожидания по клиентам",
             xlabel="Номер клиента", ylabel="Время ожидания",
             label="Время ожидания", color=:blue, alpha=0.7)
    hline!([mean(waiting_times)], label="Среднее = $(round(mean(waiting_times), digits=3))", 
           color=:red, linewidth=2)
    savefig(p2, "plots/waiting_time_histogram.png")
    println("✓ plots/waiting_time_histogram.png")
    
    # 3. Динамика количества клиентов
    time_points = 0:0.1:maximum(service_end_times)
    n_clients = zeros(length(time_points))
    
    for j in 1:length(time_points)
        t = time_points[j]
        count = 0
        for i in 1:num_customers
            if arrival_times[i] <= t && service_end_times[i] > t
                count += 1
            end
        end
        n_clients[j] = count
    end
    
    p3 = plot(time_points, n_clients, 
              title="Количество клиентов в системе",
              xlabel="Время", ylabel="Количество клиентов",
              label="Клиенты в системе", linewidth=2, color=:blue)
    hline!([num_servers], linestyle=:dash, linewidth=2, 
           label="Количество серверов ($num_servers)", color=:red)
    savefig(p3, "plots/system_load.png")
    println("✓ plots/system_load.png")
    
    # Краткая статистика
    println("\n=== СТАТИСТИКА СИСТЕМЫ ===")
    println("Среднее время ожидания: $(round(mean(waiting_times), digits=3))")
    println("Макс. время ожидания: $(round(maximum(waiting_times), digits=3))")
    println("Общее время симуляции: $(round(maximum(service_end_times), digits=3))")
end

setup_and_run()