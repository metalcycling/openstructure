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

#include <ost/db/extract_data_helper.hh>

namespace ost { namespace db {

void ExtractValidPositions(const String& entry_name, const String& chain_name,
                           LinearIndexer& indexer,
                           LinearCharacterContainer& atomseq_container,
                           LinearPositionContainer& position_container,
                           ost::seq::SequenceHandle& seq,
                           geom::Vec3List& positions) {

  std::pair<uint64_t, uint64_t> data_range = indexer.GetDataRange(entry_name,
                                                                  chain_name);
  String atomseq;
  atomseq_container.GetCharacters(data_range, atomseq);
  geom::Vec3List extracted_positions;
  position_container.GetPositions(data_range, extracted_positions);

  positions.clear();
  std::vector<char> valid_seq;
  for(uint i = 0; i < atomseq.size(); ++i) {
    if(atomseq[i] != '-') {
      positions.push_back(extracted_positions[i]);
      valid_seq.push_back(atomseq[i]);
    }
  }

  seq = ost::seq::CreateSequence("extracted_seq", String(valid_seq.begin(), 
                                                         valid_seq.end()));
}

void ExtractTemplateData(const String& entry_name, const String& chain_name,
                         const ost::seq::AlignmentHandle& aln,
                         LinearIndexer& indexer,
                         LinearCharacterContainer& seqres_container,
                         LinearCharacterContainer& atomseq_container,
                         LinearPositionContainer& position_container,
                         std::vector<int>& residue_numbers,
                         geom::Vec3List& ca_positions) {

  std::pair<uint64_t, uint64_t> data_range = indexer.GetDataRange(entry_name,
                                                                  chain_name);

  String template_seqres = aln.GetSequence(1).GetGaplessString();
  data_range.first += aln.GetSequence(1).GetOffset();
  data_range.second = data_range.first + template_seqres.size();

  // check, whether the the template seqres is consistent with what
  // we find in seqres_container
  String expected_template_seqres;
  seqres_container.GetCharacters(data_range, expected_template_seqres);
  if(expected_template_seqres != template_seqres) {
    throw std::runtime_error("Template sequence in input alignment is "
                             "inconsistent with sequence in SEQRES container!");
  }

  String template_atomseq;
  atomseq_container.GetCharacters(data_range, template_atomseq);
  geom::Vec3List extracted_positions;
  position_container.GetPositions(data_range, extracted_positions);

  uint current_rnum = aln.GetSequence(0).GetOffset() + 1;
  uint current_template_pos = 0;
  String seqres_seq = aln.GetSequence(0).GetString();
  String template_seq = aln.GetSequence(1).GetString();

  // prepare output
  uint template_atomseq_size = template_atomseq.size();
  ca_positions.clear();
  residue_numbers.clear();
  ca_positions.reserve(template_atomseq_size);
  residue_numbers.reserve(template_atomseq_size);  

  for(int i = 0; i < aln.GetLength(); ++i) {

    if(seqres_seq[i] != '-' && template_seq[i] != '-') {
      if(template_atomseq[current_template_pos] != '-') {
        // it is aligned and we have a valid position!
        residue_numbers.push_back(current_rnum);
        ca_positions.push_back(extracted_positions[current_template_pos]);
      }
    }

    if(seqres_seq[i] != '-') {
      ++current_rnum;
    }

    if(template_seq[i] != '-') {
      ++current_template_pos;
    }
  }
}

}} //ns
