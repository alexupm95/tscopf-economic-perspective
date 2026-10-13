using Pkg

const DOCS_DIR = @__DIR__
const REPO_ROOT = dirname(DOCS_DIR)
const DOCUMENTER_UUID = Base.UUID("e30172f5-a6a5-5a46-863b-614d45cd2de4")

Pkg.activate(DOCS_DIR)
Pkg.develop(PackageSpec(path = REPO_ROOT))
if !haskey(Pkg.project().dependencies, DOCUMENTER_UUID)
    Pkg.add("Documenter")
end
Pkg.instantiate()
println("docs env ready")
