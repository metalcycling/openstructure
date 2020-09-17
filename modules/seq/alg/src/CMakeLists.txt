set(OST_SEQ_ALG_HEADERS
alignment_opts.hh
clip_alignment.hh
conservation.hh
contact_prediction_score.hh
contact_weight_matrix.hh
distance_map.hh
entropy.hh
global_align.hh
ins_del.hh
local_align.hh
merge_pairwise_alignments.hh
module_config.hh
pair_subst_weight_matrix.hh
semiglobal_align.hh
sequence_identity.hh
sequence_similarity.hh
subst_weight_matrix.hh
variance_map.hh
data/default_contact_weight_matrix.hh
data/default_pair_subst_weight_matrix.hh
hmm_pseudo_counts.hh
hmm_score.hh
)

set(OST_SEQ_ALG_SOURCES
clip_alignment.cc
conservation.cc
contact_prediction_score.cc
contact_weight_matrix.cc
distance_map.cc
entropy.cc
global_align.cc
ins_del.cc
local_align.cc
merge_pairwise_alignments.cc
pair_subst_weight_matrix.cc
semiglobal_align.cc
sequence_identity.cc
sequence_similarity.cc
subst_weight_matrix.cc
variance_map.cc
hmm_pseudo_counts.cc
hmm_score.cc
)

module(NAME seq_alg HEADER_OUTPUT_DIR ost/seq/alg SOURCES ${OST_SEQ_ALG_SOURCES}
       HEADERS ${OST_SEQ_ALG_HEADERS} DEPENDS_ON ost_seq)
