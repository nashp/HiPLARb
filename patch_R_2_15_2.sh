#!/bin/sh

echo "1392c1392
< SEXP attribute_hidden
---
> SEXP" > patchfile

patch $1/src/main/objects.c <patchfile
