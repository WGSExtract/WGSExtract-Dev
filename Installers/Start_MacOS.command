#!/bin/bash
echo "Starting WGSExtract..."
WGSEDIR=$(dirname "$0")
WGSEABS=$(cd "$WGSEDIR"; pwd)
WGSEESC=${WGSEABS/ /\\}
/usr/local/bin/python3 "$WGSEABS"/program/wgsextract.py
