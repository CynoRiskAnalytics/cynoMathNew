#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$ROOT_DIR/build"
SRC_DIR="$ROOT_DIR/cynoMathTool"
JAVA_SRC_DIR="$ROOT_DIR/examples/java"
JAVA_BUILD_DIR="$JAVA_SRC_DIR/build"
JAVA_CLASSES_DIR="$JAVA_BUILD_DIR/classes"
JAVA_HOME_DIR="$(dirname "$(dirname "$(readlink -f "$(command -v javac)")")")"

mkdir -p "$BUILD_DIR"
mkdir -p "$JAVA_BUILD_DIR" "$JAVA_CLASSES_DIR"

gcc -shared -fPIC -O2 \
  -I"$SRC_DIR" \
  -o "$BUILD_DIR/libcynomath.so" \
  "$SRC_DIR/BlackScholes.c" \
  "$SRC_DIR/BlackScholesBarrierOption.c" \
  "$SRC_DIR/BlackSholesDoubleBarrier.c" \
  "$SRC_DIR/Integration_fin.c" \
  "$SRC_DIR/MTRandomNumber.c" \
  "$SRC_DIR/cynoMathUtility.c" \
  "$SRC_DIR/integration_GH.c" \
  "$SRC_DIR/integration_GL.c" \
  "$SRC_DIR/matrix.c" \
  "$SRC_DIR/optimizer.c" \
  "$SRC_DIR/pch.c" \
  "$SRC_DIR/randomNumber.c" \
  "$SRC_DIR/rootFinder.c" \
  "$SRC_DIR/sobol.c" \
  "$SRC_DIR/sorting.c" \
  "$SRC_DIR/specialFunctions.c" \
  -lm

dotnet build "$ROOT_DIR/examples/csharp/CynoMathExample.csproj"

javac -d "$JAVA_CLASSES_DIR" "$JAVA_SRC_DIR/CynoMathExample.java"

gcc -shared -fPIC -O2 \
  -I"$JAVA_HOME_DIR/include" \
  -I"$JAVA_HOME_DIR/include/linux" \
  -I"$SRC_DIR" \
  -L"$BUILD_DIR" \
  -Wl,-rpath,'$ORIGIN/../../../build' \
  -o "$JAVA_BUILD_DIR/libcynomathjni.so" \
  "$JAVA_SRC_DIR/CynoMathJNI.c" \
  -lcynomath

echo "Build complete."
echo "Native library: $BUILD_DIR/libcynomath.so"
echo "C# example: $ROOT_DIR/examples/csharp/bin/Debug/net8.0/CynoMathExample.dll"
echo "Java classes: $JAVA_CLASSES_DIR"
echo "Java JNI library: $JAVA_BUILD_DIR/libcynomathjni.so"
