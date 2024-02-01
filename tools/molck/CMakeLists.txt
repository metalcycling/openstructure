if (NOT WIN32)
executable(NAME molck SOURCES main.cc
           DEPENDS_ON ost_io ost_mol_alg STATIC)
endif(NOT WIN32)
