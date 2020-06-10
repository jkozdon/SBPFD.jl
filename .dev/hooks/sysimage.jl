#!/usr/bin/env julia

using Pkg
Pkg.add("PackageCompiler")
using PackageCompiler
Pkg.activate(joinpath(@__DIR__, ".."))
using PackageCompiler
PackageCompiler.create_sysimage(
  :JuliaFormatter;
  precompile_execution_file = joinpath(@__DIR__, "precompile.jl"),
  sysimage_path = joinpath(
    @__DIR__,
    "..",
    "..",
    ".git",
    "hooks",
    "JuliaFormatterSysimage.so",
  ),
)
