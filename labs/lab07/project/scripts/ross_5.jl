using ResumableFunctions
using ConcurrentSim
using Distributions
using Random
using StableRNGs
using Plots
using Statistics

const RUNS = 5
const N = 10
const S = 3
const SEED = 150
const LAMBDA = 100
const MU = 1

const rng = StableRNG(42) # setting a random seed for reproducibility
const F = Exponential(LAMBDA)
const G = Exponential(MU)

# Глобальные переменные для сбора данных
all_working_machines = Vector{Vector{Float64}}()
all_working_times = Vector{Vector{Float64}}()
all_queue_lengths = Vector{Vector{Float64}}()
all_queue_times = Vector{Vector{Float64}}()

# Счетчик машин в очереди
queue_counter = 0

@resumable function machine(
    env::Environment,
    repair_facility::Resource,
    spares::Store{Process},
    working_machines::Vector{Float64},
    working_times::Vector{Float64},
    queue_lengths::Vector{Float64},
    queue_times::Vector{Float64},
    total_machines::Int,
)
    global queue_counter
    
    while true
        try
            @yield timeout(env, Inf)
        catch
        end
        
        # Машина работает
        @yield timeout(env, rand(rng, F))
        
        # Поломка - уменьшаем исправные машины
        push!(working_times, now(env))
        push!(working_machines, working_machines[end] - 1)
        
        get_spare = take!(spares)
        @yield get_spare | timeout(env)
        if state(get_spare) != ConcurrentSim.idle
            @yield interrupt(value(get_spare))
        else
            throw(StopSimulation("No more spares!"))
        end
        
        # Запись очереди перед запросом ремонтника
        push!(queue_times, now(env))
        push!(queue_lengths, queue_counter)
        
        queue_counter += 1  # Входим в очередь
        @yield request(repair_facility)
        queue_counter -= 1  # Выходим из очереди
        
        @yield timeout(env, rand(rng, G))
        @yield unlock(repair_facility)
        
        # Ремонт завершен - увеличиваем исправные машины
        push!(working_times, now(env))
        push!(working_machines, working_machines[end] + 1)
        
        @yield put!(spares, active_process(env))
    end
end

@resumable function monitor_system(
    env::Environment,
    repair_facility::Resource,
    working_machines::Vector{Float64},
    working_times::Vector{Float64},
    queue_lengths::Vector{Float64},
    queue_times::Vector{Float64},
    total_machines::Int,
)
    global queue_counter
    
    while true
        @yield timeout(env, 0.1)
        
        # Запись текущей длины очереди
        push!(queue_times, now(env))
        push!(queue_lengths, queue_counter)
        
        # Запись исправных машин
        push!(working_times, now(env))
        push!(working_machines, working_machines[end])
    end
end

@resumable function start_sim(
    env::Environment,
    repair_facility::Resource,
    spares::Store{Process},
    working_machines::Vector{Float64},
    working_times::Vector{Float64},
    queue_lengths::Vector{Float64},
    queue_times::Vector{Float64},
    total_machines::Int,
)
    # Запускаем мониторинг
    @process monitor_system(env, repair_facility, working_machines, working_times, 
                           queue_lengths, queue_times, total_machines)
    
    for i = 1:N
        proc = @process machine(env, repair_facility, spares, working_machines, 
                               working_times, queue_lengths, queue_times, total_machines)
        @yield interrupt(proc)
    end
    for i = 1:S
        proc = @process machine(env, repair_facility, spares, working_machines, 
                               working_times, queue_lengths, queue_times, total_machines)
        @yield put!(spares, proc)
    end
end

function sim_repair()
    global queue_counter
    queue_counter = 0
    
    total_machines = N + S
    working_machines = [Float64(total_machines)]
    working_times = [0.0]
    queue_lengths = Float64[]
    queue_times = Float64[]
    
    sim = Simulation()
    repair_facility = Resource(sim)  # 1 ремонтник
    spares = Store{Process}(sim)
    @process start_sim(sim, repair_facility, spares, working_machines, working_times, 
                      queue_lengths, queue_times, total_machines)
    msg = run(sim)
    stop_time = now(sim)
    println("At time $stop_time: $msg")
    
    # Сохраняем данные для графиков
    push!(all_working_machines, copy(working_machines))
    push!(all_working_times, copy(working_times))
    push!(all_queue_lengths, copy(queue_lengths))
    push!(all_queue_times, copy(queue_times))
    
    return stop_time, working_machines, working_times, queue_lengths, queue_times
end

# Запуск симуляций и сбор данных
results = Float64[]
all_stop_times = Float64[]

println("Запуск $(RUNS) симуляций...\n")

for i = 1:RUNS
    println("--- Симуляция $i ---")
    stop_time, _, _, _, _ = sim_repair()
    push!(results, stop_time)
    push!(all_stop_times, stop_time)
end

println("\nAverage crash time: ", sum(results)/RUNS)

# Построение графиков по последней симуляции
println("\n=== ПОСТРОЕНИЕ ГРАФИКОВ ===")

# Создаем папку для графиков
if !isdir("plots")
    mkdir("plots")
end

# Берем данные последней симуляции
working_machines = all_working_machines[end]
working_times = all_working_times[end]
queue_lengths = all_queue_lengths[end]
queue_times = all_queue_times[end]
total_machines = N + S

# График 1: Количество исправных машин во времени
p1 = plot(title="Количество исправных машин во времени",
          xlabel="Время", ylabel="Количество исправных машин",
          legend=:topright, linewidth=2)

# Подготовка данных для ступенчатого графика
plot_times = []
plot_values = []
for i in 1:length(working_times)-1
    push!(plot_times, working_times[i])
    push!(plot_values, working_machines[i])
    push!(plot_times, working_times[i+1])
    push!(plot_values, working_machines[i])
end

plot!(plot_times, plot_values, label="Исправные машины", color=:blue, step=:post)
hline!([total_machines], linestyle=:dash, label="Всего машин ($total_machines)", color=:green)
hline!([N], linestyle=:dot, label="Рабочих машин (N=$N)", color=:orange)
savefig(p1, "plots/working_machines_time.png")
println("✓ plots/working_machines_time.png")

# График 2: Длина очереди на ремонт
p2 = plot(title="Длина очереди на ремонт",
          xlabel="Время", ylabel="Длина очереди",
          legend=:topright, linewidth=2)

# Подготовка данных для ступенчатого графика очереди
if length(queue_times) > 0
    queue_times_plot = []
    queue_vals_plot = []
    for i in 1:length(queue_times)-1
        push!(queue_times_plot, queue_times[i])
        push!(queue_vals_plot, queue_lengths[i])
        push!(queue_times_plot, queue_times[i+1])
        push!(queue_vals_plot, queue_lengths[i])
    end
    
    plot!(queue_times_plot, queue_vals_plot, label="Длина очереди", color=:red, step=:post)
    if length(queue_lengths) > 0
        hline!([mean(queue_lengths)], linestyle=:dash, 
               label="Средняя = $(round(mean(queue_lengths), digits=2))", color=:purple)
    end
end
savefig(p2, "plots/queue_length_time.png")
println("✓ plots/queue_length_time.png")

# График 3: Совмещенный график
p3 = plot(title="Состояние системы во времени",
          xlabel="Время", ylabel="Количество",
          legend=:topright)

# Исправные машины
plot!(plot_times, plot_values, label="Исправные машины", color=:blue, step=:post, linewidth=2)

# Очередь (масштабируем для наглядности)
if length(queue_lengths) > 0 && maximum(queue_lengths) > 0
    scale_factor = total_machines / (maximum(queue_lengths) + 1)
    queue_scaled = queue_vals_plot .* scale_factor
    plot!(queue_times_plot, queue_scaled, label="Очередь (масштаб)", color=:red, step=:post, linewidth=1.5, linestyle=:dot)
end

hline!([total_machines], linestyle=:dash, label="Всего машин", color=:green)
hline!([N], linestyle=:dot, label="Рабочих машин (N)", color=:orange)
savefig(p3, "plots/combined_state.png")
println("✓ plots/combined_state.png")

# Дополнительная статистика
println("\n=== СТАТИСТИКА ПОСЛЕДНЕЙ СИМУЛЯЦИИ ===")
println("Время остановки: $(round(all_stop_times[end], digits=3))")
println("Средняя длина очереди: $(round(mean(queue_lengths), digits=3))")
println("Макс. длина очереди: $(round(maximum(queue_lengths), digits=0))")
println("Минимум исправных машин: $(round(minimum(working_machines), digits=0))")
println("Среднее исправных машин: $(round(mean(working_machines), digits=2))")

println("\nВсе графики сохранены в папку 'plots'")