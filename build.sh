#!/bin/bash
mkdir -p build
cd build
cmake .. 2>/dev/null || true   # игнорируем, если уже есть
make -j$(nproc) && ./fem "$@"