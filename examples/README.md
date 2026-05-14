# Examples

This folder shows how to call the native CynoMath library from:
- Java via JNI
- C# via `DllImport`
- VBA via `Declare PtrSafe`

The examples focus on the plain-array API because it is the simplest cross-language surface for:
- scalar pricing functions
- random number generation
- Sobol quasi-random generation
- antithetic normal pairs
- matrix operations

Library names:
- Linux: `libcynomath.so`
- Windows: `cynoMathTools.dll`

Files:
- `java/CynoMathExample.java`
- `java/CynoMathJNI.c`
- `csharp/Program.cs`
- `csharp/CynoMathExample.csproj`
- `vba/CynoMathExample.bas`

Notes:
- The Java example uses a small JNI bridge and loads `build/libcynomath.so` plus `examples/java/build/libcynomathjni.so`.
- The C# example targets the current native export style used in this repo.
- The VBA example is Windows-only and uses `LongPtr`-safe declarations for 64-bit Office.
- Sobol generation is exact for the built-in first 8 dimensions and uses a deterministic low-discrepancy fallback beyond that.
