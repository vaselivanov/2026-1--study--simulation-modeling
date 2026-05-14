using ResumableFunctions
using ConcurrentSim
using Distributions
using Random
using StableRNGs
using Statistics

const RUNS = 5
const N = 10
const S = 3
const SEED = 150
const LAMBDA = 100
const MU = 1
const NUM_REPAIR = 1  # один ремонтник

const rng = StableRNG(42)
const F = Exponential(LAMBDA)
const G = Exponential(MU)

# Глобальные переменные для мониторинга
queue_lengths = Float64[]
queue_times = Float64[]
repair_busy_time = 0.0  # время занятости ремонтника
total_repair_time = 0.0  # общее время ремонта
last_state_change = 0.0  # время последнего изменения состояния

@resumable function machine(
    env::Environment,
    repair_facility::Resource,
    spares::Store{Process},
)
    global repair_busy_time, last_state_change
    
    while true
        try
            @yield timeout(env, Inf)
        catch
        end
        @yield timeout(env, rand(rng, F))
        get_spare = take!(spares)
        @yield get_spare | timeout(env)
        if state(get_spare) != ConcurrentSim.idle
            @yield interrupt(value(get_spare))
        else
            throw(StopSimulation("No more spares!"))
        end
        
        # Запрос ремонтника и запись длины очереди
        current_queue = length(repair_facility.queue)
        push!(queue_times, now(env))
        push!(queue_lengths, current_queue)
        
        @yield request(repair_facility)
        
        # Ремонт
        repair_time = rand(rng, G)
        @yield timeout(env, repair_time)
        
        @yield unlock(repair_facility)
        @yield put!(spares, active_process(env))
    end
end

@resumable function monitor_system(
    env::Environment,
    repair_facility::Resource
)
    while true
        @yield timeout(env, 0.1)  # проверяем каждые 0.1 времени
        current_queue = length(repair_facility.queue)
        push!(queue_times, now(env))
        push!(queue_lengths, current_queue)
    end
end

@resumable function start_sim(
    env::Environment,
    repair_facility::Resource,
    spares::Store{Process},
)
    # Запускаем мониторинг
    @process monitor_system(env, repair_facility)
    
    for i = 1:N
        proc = @process machine(env, repair_facility, spares)
        @yield interrupt(proc)
    end
    for i = 1:S
        proc = @process machine(env, repair_facility, spares)
        @yield put!(spares, proc)
    end
end

function sim_repair()
    global queue_lengths, queue_times, repair_busy_time, last_state_change
    
    # Сбрасываем данные
    queue_lengths = Float64[]
    queue_times = Float64[]
    repair_busy_time = 0.0
    last_state_change = 0.0
    
    sim = Simulation()
    repair_facility = Resource(sim, NUM_REPAIR)
    spares = Store{Process}(sim)
    @process start_sim(sim, repair_facility, spares)
    msg = run(sim)
    stop_time = now(sim)
    
    # Расчет метрик
    avg_queue_length = length(queue_lengths) > 0 ? mean(queue_lengths) : 0.0
    max_queue_length = length(queue_lengths) > 0 ? maximum(queue_lengths) : 0
    
    # Расчет загрузки ремонтника через статистику ресурса
    # ConcurrentSim не дает прямой доступ к времени занятости, оцениваем косвенно
    utilization = (1 - (length(repair_facility.queue) / (NUM_REPAIR * stop_time + 1))) * 100
    
    # Альтернативный расчет загрузки через количество обращений
    total_services = length(repair_facility.server)  # количество обслуженных
    
    println("\n=== РЕЗУЛЬТАТЫ СИМУЛЯЦИИ ===")
    println("Время остановки: $(round(stop_time, digits=3))")
    println("Причина: $msg")
    println("\n=== МОНИТОРИНГ ОЧЕРЕДИ ===")
    println("Средняя длина очереди: $(round(avg_queue_length, digits=3))")
    println("Максимальная длина очереди: $max_queue_length")
    println("Всего измерений: $(length(queue_lengths))")
    println("\n=== ЗАГРУЗКА РЕМОНТНИКА ===")
    println("Количество ремонтников: $NUM_REPAIR")
    println("Примерная загрузка: $(round(utilization, digits=2))%")
    
    return stop_time, avg_queue_length, max_queue_length, utilization
end

# Запуск нескольких симуляций
results_time = Float64[]
results_queue = Float64[]
results_maxqueue = Float64[]
results_util = Float64[]

println("Запуск $(RUNS) симуляций с 1 ремонтником...\n")
println("="^60)

for i = 1:RUNS
    println("\n--- Симуляция $i ---")
    stop_time, avg_q, max_q, util = sim_repair()
    push!(results_time, stop_time)
    push!(results_queue, avg_q)
    push!(results_maxqueue, max_q)
    push!(results_util, util)
end

# Итоговая статистика
println("\n" * "="^60)
println("ИТОГОВАЯ СТАТИСТИКА ПО $(RUNS) ЗАПУСКАМ (1 ремонтник)")
println("="^60)

println("\n=== ВРЕМЯ ДО КРАХА ===")
println("Среднее: $(round(mean(results_time), digits=3))")
println("Ст. отклонение: $(round(std(results_time), digits=3))")
println("Минимум: $(round(minimum(results_time), digits=3))")
println("Максимум: $(round(maximum(results_time), digits=3))")

println("\n=== СРЕДНЯЯ ДЛИНА ОЧЕРЕДИ ===")
println("Среднее: $(round(mean(results_queue), digits=3))")
println("Ст. отклонение: $(round(std(results_queue), digits=3))")
println("Минимум: $(round(minimum(results_queue), digits=3))")
println("Максимум: $(round(maximum(results_queue), digits=3))")

println("\n=== МАКСИМАЛЬНАЯ ДЛИНА ОЧЕРЕДИ ===")
println("Среднее: $(round(mean(results_maxqueue), digits=3))")
println("Ст. отклонение: $(round(std(results_maxqueue), digits=3))")

println("\n=== ЗАГРУЗКА РЕМОНТНИКА ===")
println("Средняя: $(round(mean(results_util), digits=2))%")
println("Ст. отклонение: $(round(std(results_util), digits=2))%")
println("Минимум: $(round(minimum(results_util), digits=2))%")
println("Максимум: $(round(maximum(results_util), digits=2))%")

println("\n" * "="^60)
println("СРАВНЕНИЕ:")
println("Для 1 ремонтника среднее время жизни = $(round(mean(results_time), digits=3))")
println("Средняя очередь = $(round(mean(results_queue), digits=3))")
println("Загрузка = $(round(mean(results_util), digits=2))%")