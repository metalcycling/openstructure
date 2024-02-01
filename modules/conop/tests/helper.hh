#ifndef OST_CONOP_TEST_HELPERS_HH
#define OST_CONOP_TEST_HELPERS_HH

#include <ost/mol/entity_handle.hh>
#include <ost/mol/residue_handle.hh>

namespace ost { namespace conop {

mol::ResidueHandle make_uracyl(mol::ChainHandle chain);
void verify_1zk_connectivity(const mol::ResidueHandle& r1);
mol::ResidueHandle make_cytosine(mol::ChainHandle chain); 
mol::ResidueHandle make_uracil1(mol::ChainHandle chain); 
mol::ResidueHandle make_uracil2(mol::ChainHandle chain); 
mol::ResidueHandle make_defective_uracil2(mol::ChainHandle chain);
mol::ResidueHandle make_leu(mol::ChainHandle chain);
mol::ResidueHandle make_arg(mol::ChainHandle chain);
mol::ResidueHandle make_defective_leu(mol::ChainHandle chain);
mol::ResidueHandle make_1zk(mol::ChainHandle chain);
void verify_1zk_connectivity(const mol::ResidueHandle& r1);
void verify_nucleotide_connectivity(const mol::ResidueHandle& res);
void verify_nucleotide_link(const mol::ResidueHandle& p3, 
                            const mol::ResidueHandle& p5);
void verify_nucleotide_nolink(const mol::ResidueHandle& p3, 
                              const mol::ResidueHandle& p5);
}}

#endif
