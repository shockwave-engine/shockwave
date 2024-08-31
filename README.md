# Shockwave

Shockwave is a 2D physics engine for rigid bodies using sequential impulses.

It improves over the [simu](https://github.com/Matthieu-Beauchamp/Simu) engine, notably by implementing a [shockwave collision algorithm](https://www.cs.columbia.edu/cg/rosi/rosi.pdf) to correctly simulate simultaneous collisions. 

# Development

## Tools

Ensure you have the following installed
- python (3.12+)
- clang-format (~18)
- cmake (3.30+)

## Building

To build locally, clone the repository and use the `tools.py` script to build or test:
```sh
./tools.py build
```
```sh
./tools.py run test
```

Use the following to see other options
```sh
./tools.py -h
```

If using the `clangd` LSP, make sure to build once for the `compile_commands.json` to be generated.

You can also use CMake with Visual Studio Code and other editors. You will need to set the `SHOCKWAVE_BUILD_TESTS=ON` cmake option to build the `shockwave-test` target.
