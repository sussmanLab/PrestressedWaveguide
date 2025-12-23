using Test, StaticArrays, LinearAlgebra
using PrestressedWaveguide


function findTests()
    allFiles = readdir(@__DIR__)
    testFiles = filter(f -> endswith(f, ".jl") && f != "runtests.jl", allFiles)
    return [replace(f, ".jl" => "") for f in testFiles]
end

const allTests = findTests()
const testsToRun = isempty(ARGS) ? allTests : ARGS

@testset verbose = true "PrestressedWaveguide.jl" begin
    for test in testsToRun
        include("$(test).jl")
    end
end


