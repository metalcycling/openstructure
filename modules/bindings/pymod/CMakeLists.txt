set(OST_BINDINGS
__init__.py
lga.py
hbplus.py
msms.py
tmtools.py 
dssp.py
clustalw.py
utils.py
naccess.py
hhblits.py
hhblits2.py
hhblits3.py
blast.py
cadscore.py
kclust.py
ialign.py
)

set(OST_BINDINGS_PYMOD_SOURCES
  export_tmalign.cc
  wrap_bindings.cc
)

pymod(NAME bindings OUTPUT_DIR ost/bindings CPP ${OST_BINDINGS_PYMOD_SOURCES} 
      PY ${OST_BINDINGS})

