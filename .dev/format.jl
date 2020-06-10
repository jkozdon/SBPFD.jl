#!/usr/bin/env julia
#
# This is an adapted version of format.jl from JuliaFormatter with the
# following license:
#
#    MIT License
#    Copyright (c) 2019 Dominique Luna
#
#    Permission is hereby granted, free of charge, to any person obtaining a
#    copy of this software and associated documentation files (the
#    "Software"), to deal in the Software without restriction, including
#    without limitation the rights to use, copy, modify, merge, publish,
#    distribute, sublicense, and/or sell copies of the Software, and to permit
#    persons to whom the Software is furnished to do so, subject to the
#    following conditions:
#
#    The above copyright notice and this permission notice shall be included
#    in all copies or substantial portions of the Software.
#
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#    OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
#    NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
#    DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
#    OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
#    USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# and from climaformat.jl from ClimateMachine.jl with the following license:
#
#    Copyright 2019 Climate Modeling Alliance
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using JuliaFormatter

include("formatter_options.jl")

help = """
Usage: format.jl [flags] [FILE/PATH]...

Formats the given julia files using the CLIMA formatting options.  If paths
are given it will format the julia files in the paths. Otherwise, it will
format all changed julia files.

    -v, --verbose
        Print the name of the files being formatted with relevant details.

    -h, --help
        Print this message.
"""

function parse_opts!(args::Vector{String})
  i = 1
  opts = Dict{Symbol, Union{Int, Bool}}()
  while i ≤ length(args)
    arg = args[i]
    if arg[1] != '-'
      i += 1
      continue
    end
    if arg == "-v" || arg == "--verbose"
      opt = :verbose
    elseif arg == "-h" || arg == "--help"
      opt = :help
    else
      error("invalid option $arg")
    end
    if opt in (:verbose, :help)
      opts[opt] = true
      deleteat!(args, i)
    end
  end
  return opts
end

opts = parse_opts!(ARGS)
if haskey(opts, :help)
  write(stdout, help)
  exit(0)
end
filenames = ARGS

format(filenames; formatter_options..., opts...)
