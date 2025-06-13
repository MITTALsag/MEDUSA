#!/bin/bash
for i in {1..1000}; do
gnome-terminal -- bash -c "./LzAsync/DispatchQ_Test; if [ \$? -eq 0 ]; then exit 0; else echo 'Program crashed (Exit Code: \$?)'; exec bash; fi"
done