#!/bin/sh

# post-commit hook to automatically track version info
# cp post-commit.sh .git/hooks/post-commit

git describe --tags > version_info
