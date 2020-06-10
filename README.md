SBPFD.jl
========

An summation by parts finite difference library built on top of
[KernelAbstractions.jl].

Developer Tools
---------------

[JuliaFormatter.jl] is used to consistently format the code base. To make
formatting more seamless a git pre-commit is available. This can be enabled with
the following symbolic link executed from a top level clone of SBPFD.jl:

```sh
ln -s ../../.dev/hooks/pre-commit .git/hooks/pre-commit
```

Alternatively, a custom sysmem image can be used which will make the run time of
the formatter much more reasonable. To use this execute the following from a
top level clone of SBPFD.jl:

```sh
julia .dev/hooks/sysimage.jl
ln -s ../../.dev/hooks/pre-commit .git/hooks/pre-commit
```

License
-------

SBPFD.jl is licensed under [CC0 license](LICENSE.md).

[KernelAbstractions.jl]: https://github.com/JuliaGPU/KernelAbstractions.jl
[JuliaFormatter.jl]: https://github.com/domluna/JuliaFormatter.jl
