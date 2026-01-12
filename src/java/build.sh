#!/bin/bash
# Build FastAddVar Java plugin for Stata

# Find Stata's sfi-api.jar
STATA_PATHS=(
    "/Applications/Stata/sfi-api.jar"
    "/Applications/Stata18/sfi-api.jar"
    "/Applications/StataSE.app/Contents/Java/sfi-api.jar"
    "/Applications/StataMP.app/Contents/Java/sfi-api.jar"
    "/Applications/Stata/StataSE.app/Contents/Java/sfi-api.jar"
    "/Applications/Stata/StataMP.app/Contents/Java/sfi-api.jar"
    "$HOME/Applications/Stata/sfi-api.jar"
)

SFI_JAR=""
for path in "${STATA_PATHS[@]}"; do
    if [ -f "$path" ]; then
        SFI_JAR="$path"
        break
    fi
done

if [ -z "$SFI_JAR" ]; then
    echo "Error: Could not find sfi-api.jar"
    echo "Please set SFI_JAR environment variable to the path of sfi-api.jar"
    exit 1
fi

echo "Using SFI API: $SFI_JAR"

# Compile
javac -cp "$SFI_JAR" -d ../../build FastAddVar.java

if [ $? -eq 0 ]; then
    echo "Success! FastAddVar.class written to build/"
else
    echo "Compilation failed"
    exit 1
fi
