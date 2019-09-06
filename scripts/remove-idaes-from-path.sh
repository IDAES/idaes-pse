#!/usr/bin/env bash
# Remove any existing references to idaes egg files
# from the easy-install.pth . Called from the CircleCI config.yml
# -dang 8/30/2019
f=$(dirname $(which python))/../lib/python*/site-packages/easy-install.pth
grep -v "idaes.*-.*egg" ${f}  > /tmp/easy-install.pth.tmp && \
 mv /tmp/easy-install.pth.tmp ${f}
