#!/bin/bash

# I copy the rendered book from the Rbookdown folder "_book" to the deploy repo
# Andrea Del Cortona
# 2024/09/27

DEPLOY_DIR=$1

~/bin/quarto-1.5.57/bin/quarto render --output-dir ./docs

#rsync -ah --remove-source-files --info=progress2 docs/ "${DEPLOY_DIR}"/
#rm -rf docs
