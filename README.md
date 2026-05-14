# CynoMathTool

Cross-platform branch for the native math library and .NET interop tests.

Linux flow:
- Native library sources now build from `cynoMathTool/*.c`.
- Rebuild the native library with `make -C build` when the `build/` tree is already configured.
- Run the managed tests with `dotnet test cynoMathToolTest/cynoMathToolTest.csproj --no-build` after the native library is built.

If the `build/` directory does not exist yet, configure it first with CMake against `cynoMathTool/`.

Examples:
- Cross-language usage samples live in `examples/`.
- The current examples cover Java, C#, and VBA interop against the native library.
