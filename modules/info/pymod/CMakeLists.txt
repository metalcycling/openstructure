set(OST_INFO_PYMOD_SOURCES
wrap_info.cc
)

pymod(NAME info CPP ${OST_INFO_PYMOD_SOURCES} PY __init__.py)
