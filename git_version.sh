#!/bin/sh

# Directory to put git_version.h
DIR="$(pwd)"  # Or adjust if running from out-of-tree build
TMPFILE="$(mktemp)"

# Default values in case Git is not available
BUILDNO="unknown"
BRANCH="unknown"

if git rev-parse --git-dir >/dev/null 2>&1; then
    BUILDNO=$(git describe --always --dirty)
    BRANCH=$(git rev-parse --abbrev-ref HEAD)
fi

# Generate temporary header
cat > "$TMPFILE" <<EOF
#ifndef _GIT_VERSION_H
#define _GIT_VERSION_H
#define PACKAGE_BUILD "$BUILDNO"
#define PACKAGE_BRANCH "$BRANCH"
#endif
EOF

# Only replace git_version.h if content changed
if [ ! -f "$DIR/git_version.h" ] || ! cmp -s "$TMPFILE" "$DIR/git_version.h"; then
    mv "$TMPFILE" "$DIR/git_version.h"
    echo "Updated git_version.h (Build: $BUILDNO, Branch: $BRANCH)"
else
    rm "$TMPFILE"
fi

