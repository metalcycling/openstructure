set(OST_DB_PYMOD_SOURCES
export_linear_db.cc
wrap_db.cc
)

pymod(NAME db CPP ${OST_DB_PYMOD_SOURCES} PY __init__.py)