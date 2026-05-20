#!/bin/bash
set -e

if [ "$(id -u)" != "0" ]; then
    exec /entrypoint.sh "$@"
fi

TARGET_UID=""
TARGET_GID=""

if [ -d /data ]; then
    MOUNT_UID=$(stat -c "%u" /data)
    MOUNT_GID=$(stat -c "%g" /data)
    if [ "$MOUNT_UID" != "0" ]; then
        TARGET_UID="$MOUNT_UID"
        TARGET_GID="$MOUNT_GID"
    fi
fi

if [ -z "$TARGET_UID" ]; then
    TARGET_UID=1000
    TARGET_GID=1000
fi

if ! getent group "$TARGET_GID" > /dev/null 2>&1; then
    groupadd -g "$TARGET_GID" runtimeuser
fi

if ! id -u "$TARGET_UID" > /dev/null 2>&1; then
    useradd -u "$TARGET_UID" -g "$TARGET_GID" -M -s /bin/bash runtimeuser 2>/dev/null || \
    useradd -u "$TARGET_UID" -g "$TARGET_GID" -M -s /bin/bash -o runtimeuser
fi

exec gosu "$TARGET_UID:$TARGET_GID" /entrypoint.sh "$@"
