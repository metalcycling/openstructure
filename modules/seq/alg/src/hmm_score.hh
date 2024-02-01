//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 3.0 of the License, or (at your option)
// any later version.
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//------------------------------------------------------------------------------

#ifndef OST_SEQ_ALG_HMM_SCORE_HH
#define OST_SEQ_ALG_HMM_SCORE_HH

#include <ost/seq/profile_handle.hh>
#include <ost/seq/alignment_handle.hh>

namespace ost{ namespace seq{ namespace alg{

Real HMMScore(const ost::seq::ProfileHandle& profile_0, 
	          const ost::seq::ProfileHandle& profile_1,
	          const ost::seq::AlignmentHandle& aln,
	          int seq_0_idx, int seq_1_idx, 
	          Real match_score_offset = -0.03,
	          Real correl_score_weight = 0.1,
	          Real del_start_penalty_factor=0.6,
	          Real del_extend_penalty_factor=0.6,
	          Real ins_start_penalty_factor=0.6,
	          Real ins_extend_penalty_factor=0.6);

}}} // ns

#endif