using Pkg
const DOCS_DIR = @__DIR__
const REPO_ROOT = dirname(DOCS_DIR)

Pkg.activate(DOCS_DIR)
Pkg.develop(PackageSpec(path = REPO_ROOT))
Pkg.instantiate()

using Documenter
using TSCOPF

const SRC_DIR = joinpath(DOCS_DIR, "src")

# Long manuals at docs/*.md are copied at build time.
for f in (
    "parameter_reference.md",
    "running_a_case.md",
    "configuration_map.md",
    "dynamic_controls_avr_governor.md",
)
    cp(joinpath(DOCS_DIR, f), joinpath(SRC_DIR, f); force = true)
end

makedocs(
    sitename = "TSCOPF.jl",
    authors = "Alex Junior da Cunha Coelho",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://alexupm95.github.io/tscopf-economic-perspective",
        assets = ["assets/fonts.css", "assets/custom.css"],
        mathengine = Documenter.MathJax3(Dict(
            :loader => Dict("load" => ["[tex]/require"]),
            :tex => Dict(
                "inlineMath" => [["\$", "\$"], ["\\(", "\\)"]],
                "tags" => "ams",
                "packages" => ["base", "ams", "autoload"],
            ),
        )),
    ),
    modules = [TSCOPF],
    pages = [
        "Home" => "index.md",
        "Part I — The model" => [
            "model/02_steady_state_opf.md",
            "model/03_generator_models.md",
            "model/04_network_kron_fullbus.md",
            "model/05_swing_dynamics.md",
            "model/06_tsc_opf_assembled.md",
            "model/07_duals_economics.md",
            "model/08_controls_avr_governor.md",
            "model/09_grid_forming.md",
        ],
        "Part II — Using the code" => [
            "install.md",
            "quickstart.md",
            "running_a_case.md",
            "parameter_reference.md",
            "complete_example.md",
            "configuration_map.md",
            "dynamic_controls_avr_governor.md",
            "developer_architecture.md",
        ],
        "API" => "api.md",
    ],
    checkdocs = :none,
)

deploydocs(
    repo = "github.com/alexupm95/tscopf-economic-perspective.git",
    devbranch = "main",
    devurl = "dev",
)
