#!/usr/bin/env python
import sys

raise RuntimeError("Many things have changed in the Python 3 port. " +
                   "bump-version.py might require updates. " +
                   "Only remove this error after careful checking!")

if len(sys.argv) < 2:
  print("USAGE: python scripts/bump-version.py OST_VERSION")
  print("-> OST_VERSION is MAJOR.MINOR.PATCH (e.g. 1.9.1)")
  print("-> assumption is that a git tag OST_VERSION will exist")
  sys.exit(1)

# split up version number
version_string = sys.argv[1]
version = version_string.split('.')
major, minor, patch = (int(version[0]), int(version[1]), int(version[2]))

# fix CMakeLists
lines = open("CMakeLists.txt").readlines()
for i, line in enumerate(lines):
  if line.startswith("set (OST_VERSION_MAJOR"):
    lines[i] = "set (OST_VERSION_MAJOR %d)\n" % major
  elif line.startswith("set (OST_VERSION_MINOR"):
    lines[i] = "set (OST_VERSION_MINOR %d)\n" % minor
  elif line.startswith("set (OST_VERSION_PATCH"):
    lines[i] = "set (OST_VERSION_PATCH %d)\n" % patch
open("CMakeLists.txt", "w").writelines(lines)

# fix CHANGELOG
lines = open("CHANGELOG.txt").readlines()
for i, line in enumerate(lines):
  if line.startswith("Changes in Release") and "X" in line.upper():
    lines[i] = "Changes in Release %s\n" % version_string
open("CHANGELOG.txt", "w").writelines(lines)

# fix Docker recipe
lines = open("docker/Dockerfile").readlines()
for i, line in enumerate(lines):
  if line.startswith("ARG OPENSTRUCTURE_VERSION"):
    lines[i] = 'ARG OPENSTRUCTURE_VERSION="%s"\n' % version_string
open("docker/Dockerfile", "w").writelines(lines)

# fix Singularity recipe
lines = open("singularity/Singularity").readlines()
for i, line in enumerate(lines):
  if line.startswith("export OPENSTRUCTURE_VERSION="):
    lines[i] = 'export OPENSTRUCTURE_VERSION="%s"\n' % version_string
open("singularity/Singularity", "w").writelines(lines)
