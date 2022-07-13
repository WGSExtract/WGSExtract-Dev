#!/bin/bash
WGSEDIR=$(dirname "$0")
WGSEABS=$(cd "$WGSEDIR" ; pwd)
WGSEESC=${WGSEABS/ /\\}
/bin/bash "$WGSEESC"/Upgrade_v2tov3.sh
