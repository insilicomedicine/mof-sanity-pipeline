#!/bin/bash
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Saudi Aramco -- MIT License.
# See repository LICENSE for full terms.

USER_ID=${HOST_USER_ID:-1000}
GROUP_ID=${HOST_GROUP_ID:-1000}

groupadd -f -g $GROUP_ID sanity_group

id -u sanity_user &>/dev/null || useradd -u $USER_ID -g $GROUP_ID -m sanity_user

chown -R $USER_ID:$GROUP_ID /workspace 2>/dev/null || true
chown -R $USER_ID:$GROUP_ID /opt/oxichecker 2>/dev/null || true
find /opt/libcif -writable -exec chown $USER_ID:$GROUP_ID {} \; 2>/dev/null || true

exec sudo -u sanity_user -E "$@"