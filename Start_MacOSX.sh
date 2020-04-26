#!/bin/sh
echo "Starting WGSExtract..."
WGSEDIR=`dirname "$0"`
WGSEABS=`cd "$WGSEDIR"; pwd`
/usr/local/bin/python3 $WGSEABS/programs/wgsextract/wgsextract.py
