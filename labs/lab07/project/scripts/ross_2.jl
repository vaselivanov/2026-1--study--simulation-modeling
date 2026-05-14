using ResumableFunctions
using ConcurrentSim
using Distributions
using Random
using StableRNGs
using Statistics

const RUNS = 5
const N = 10  # количество машин
const S = 3   # количество запасных частей
const SEED = 150
const LAMBDA = 100  # интенсивность отказов
const MU = 1        # интенсивность ремонта
const NUM_REPAIR = 5  # количество ремонтников

const rng = StableRNG(42) # setting a random seed for reproducibility
const F = Exponential(LAMBDA)
const G = Exponential(MU)

@resumable function machine(
    env::Environment,
    repair_facility::Resource,
    spares::Store{Process},
)
    while true
        try
            @yield timeout(env, Inf)
        catch
        end
        @yield timeout(env, rand(rng, F))  # время до отказа
        get_spare = take!(spares)
        @yield get_spare | timeout(env)
        if state(get_spare) != ConcurrentSim.idle
            @yield interrupt(value(get_spare))
        else
            throw(StopSimulation("No more spares!"))
        end
        @yield request(repair_facility)  # запрос ремонтника
        @yield timeout(env, rand(rng, G))  # время ремонта
        @yield unlock(repair_facility)     # освобождение ремонтника
        @yield put!(spares, active_process(env))  # возврат в запас
    end
end

@resumable function start_sim(
    env::Environment,
    repair_facility::Resource,
    spares::Store{Process},
)
    # Запускаем N машин
    for i = 1:N
        proc = @process machine(env, repair_facility, spares)
        @yield interrupt(proc)
    end
    # Добавляем S запасных частей
    for i = 1:S
        proc = @process machine(env, repair_facility, spares)
        @yield put!(spares, proc)
    end
end

function sim_repair()
    sim = Simulation()
    repair_facility = Resource(sim, NUM_REPAIR)  # NUM_REPAIR ремонтников
    spares = Store{Process}(sim)
    @process start_sim(sim, repair_facility, spares)
    msg = run(sim)
    stop_time = now(sim)
    println("At time $stop_time: $msg")
    return stop_time
end

results = Float64[]
for i = 1:RUNS
    push!(results, sim_repair())
end

println("Количество ремонтников: $NUM_REPAIR")
println("Среднее время до краха: ", sum(results)/RUNS)