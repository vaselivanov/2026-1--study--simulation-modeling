using Plots, Random, LinearAlgebra

# ============================================================
# 1. Решение уравнения Лапласа (метод SOR)
# ============================================================
function solve_laplace!(ϕ, mask; ω=1.9, tol=1e-5, max_iter=20000)
    nx, ny = size(ϕ)
    for _ in 1:max_iter
        max_diff = 0.0
        for i in 2:nx-1, j in 2:ny-1
            mask[i,j] && continue
            ϕ_new = 0.25 * (ϕ[i-1,j] + ϕ[i+1,j] + ϕ[i,j-1] + ϕ[i,j+1])
            diff = ϕ_new - ϕ[i,j]
            ϕ[i,j] += ω * diff
            max_diff = max(max_diff, abs(diff))
        end
        if max_diff < tol
            break
        end
    end
    return ϕ
end

# ============================================================
# 2. Выбор граничной ячейки для роста (флуктуационный критерий)
# ============================================================
function select_boundary_cell(ϕ, structure, eta)
    nx, ny = size(ϕ)
    boundary = Tuple{Int,Int}[]
    E_vals = Float64[]
    
    for i in 2:nx-1, j in 2:ny-1
        structure[i,j] && continue
        # есть ли сосед в структуре?
        has_neighbor = false
        for di in -1:1, dj in -1:1
            ni, nj = i+di, j+dj
            if 1<=ni<=nx && 1<=nj<=ny && structure[ni,nj]
                has_neighbor = true
                break
            end
        end
        if has_neighbor
            push!(boundary, (i,j))
            Ex = (ϕ[i-1,j] - ϕ[i+1,j]) / 2
            Ey = (ϕ[i,j-1] - ϕ[i,j+1]) / 2
            push!(E_vals, sqrt(Ex^2 + Ey^2))
        end
    end
    
    if isempty(boundary)
        return nothing
    end
    
    probs = E_vals .^ eta
    probs ./= sum(probs)
    idx = sample(1:length(boundary), Weights(probs))
    return boundary[idx]
end

# ============================================================
# 3. Геометрия "остриё-плоскость"
# ============================================================
function init_needle_plane(nx, ny, needle_x, needle_y; V=1.0)
    ϕ = zeros(nx, ny)
    mask = falses(nx, ny)
    # плоскость (нижняя граница)
    mask[:, 1] .= true
    # остриё
    mask[needle_x, needle_y] = true
    ϕ[mask] .= V
    return ϕ, mask
end

# ============================================================
# 4. Геометрия "точка-окружность"
# ============================================================
function init_point_circle(n, R; V=1.0)
    ϕ = zeros(n, n)
    mask = falses(n, n)
    center = (n÷2 + 1, n÷2 + 1)
    for i in 1:n, j in 1:n
        if (i-center[1])^2 + (j-center[2])^2 >= R^2
            mask[i,j] = true
        end
    end
    ϕ[center[1], center[2]] = V
    ϕ[mask] .= 0.0
    return ϕ, mask
end

# ============================================================
# 5. Основной цикл роста стримера
# ============================================================
function grow_structure!(ϕ, mask, structure, eta; max_steps=1000)
    for step in 1:max_steps
        solve_laplace!(ϕ, mask)
        # проверка достижения границы
        if all(structure[mask]) 
            break
        end
        cell = select_boundary_cell(ϕ, structure, eta)
        if cell === nothing
            break
        end
        i, j = cell
        structure[i,j] = true
        mask[i,j] = true
        ϕ[i,j] = maximum(ϕ[mask])
    end
    return structure
end

# ============================================================
# 6. Расчёт густоты ветвей
# ============================================================
function branching_density(structure, center, r)
    cnt = 0
    n = size(structure, 1)
    for i in 1:n, j in 1:n
        if structure[i,j]
            dist2 = (i-center[1])^2 + (j-center[2])^2
            if abs(dist2 - r^2) < 1
                cnt += 1
            end
        end
    end
    return cnt / (2π * r)
end

# ============================================================
# 7. Запуск расчётов и визуализация
# ============================================================

println("=== Моделирование электрического пробоя ===\n")

# --- 7.1. Остриё-плоскость, η = 2 ---
println("1. Геометрия: остриё-плоскость (η = 2)")
nx, ny = 100, 100
ϕ, mask = init_needle_plane(nx, ny, nx÷2, ny-5)
structure = copy(mask)
grow_structure!(ϕ, mask, structure, 2.0, max_steps=500)

p1 = heatmap(structure', 
             color=:greys, 
             axis=false, 
             framestyle=:none,
             title="Стример, η = 2, остриё-плоскость",
             size=(500,500))
savefig(p1, "C:/Users/Пользователь/work/study/2026-1/2026-1==study--simulation-modeling/2026-1--study--simulation-modeling/labs/lab06/project/needle_plane_eta2.png")
println("   -> Сохранён: needle_plane_eta2.png")

# --- 7.2. Точка-окружность, R = 30, η = 2 ---
println("2. Геометрия: точка-окружность (R = 30, η = 2)")
n = 150
R = 30
ϕ, mask = init_point_circle(n, R)
structure = copy(mask)
grow_structure!(ϕ, mask, structure, 2.0, max_steps=800)

p2 = heatmap(structure', 
             color=:greys, 
             axis=false, 
             framestyle=:none,
             title="Стример, R = 30, η = 2",
             size=(500,500))
savefig(p2, "circle_R30_eta2.png")
println("   -> Сохранён: circle_R30_eta2.png")

# --- 7.3. Точка-окружность, R = 60, η = 2 ---
println("3. Геометрия: точка-окружность (R = 60, η = 2)")
R = 60
ϕ, mask = init_point_circle(n, R)
structure = copy(mask)
grow_structure!(ϕ, mask, structure, 2.0, max_steps=800)

p3 = heatmap(structure', 
             color=:greys, 
             axis=false, 
             framestyle=:none,
             title="Стример, R = 60, η = 2",
             size=(500,500))
savefig(p3, "circle_R60_eta2.png")
println("   -> Сохранён: circle_R60_eta2.png")

# --- 7.4. Зависимость густоты ветвей от R ---
println("4. Исследование густоты ветвей от радиуса окружности")
radii = 20:10:70
density_vals = Float64[]

for R_test in radii
    println("   Расчёт для R = $R_test ...")
    ϕ, mask = init_point_circle(n, R_test)
    structure = copy(mask)
    grow_structure!(ϕ, mask, structure, 2.0, max_steps=500)
    center = (n÷2+1, n÷2+1)
    r_measure = max(1, R_test - 5)
    d = branching_density(structure, center, r_measure)
    push!(density_vals, d)
end

p4 = plot(radii, density_vals, 
          marker=:circle, 
          linewidth=2,
          xlabel="Радиус окружности R", 
          ylabel="Густота ветвей",
          title="Зависимость ветвистости от радиуса (η=2)",
          legend=false,
          size=(600,400))
savefig(p4, "density_vs_R.png")
println("   -> Сохранён: density_vs_R.png")

# --- 7.5. Сравнение разных η (R = 40) ---
println("5. Сравнение разных η при R = 40")
etas = [1.0, 2.0, 5.0]
R_fixed = 40
n = 150

plots_array = []  # массив для хранения графиков

for (idx, eta_test) in enumerate(etas)
    println("   η = $eta_test ...")
    ϕ_local, mask_local = init_point_circle(n, R_fixed)
    structure_local = copy(mask_local)
    grow_structure!(ϕ_local, mask_local, structure_local, eta_test, max_steps=600)
    
    p = heatmap(structure_local', 
                color=:greys, 
                axis=false, 
                framestyle=:none,
                title="η = $eta_test",
                size=(300,300))
    push!(plots_array, p)
end

# Создаём компоновку из трёх графиков
p5 = plot(plots_array..., layout=(1,3), size=(900,400))

savefig(p5, "compare_eta_R40.png")
println("   -> Сохранён: compare_eta_R40.png")