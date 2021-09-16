


using Documenter, ACEatoms 

makedocs(sitename="ACEatoms.jl Documentation",
         pages = [
        "Home" => "index.md",
        "Getting Started" => "gettingstarted.md",
        "Linear Site E Models" => "linear.md",
        "Developer Docs" => "devel.md",
        "Types & Functions" => "docs.md"
         ])

deploydocs(
    repo = "github.com/ACEsuit/ACEatoms.jl.git",
    devbranch = "main"
)
