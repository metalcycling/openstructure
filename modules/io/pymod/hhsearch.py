"""
Author: Marco Biasini
"""

import re
from ost import seq

class HHSearchHit:
  def __init__(self, summary, alignment):
    self.summary=summary
    self.alignment=alignment

class HitSummary:
  def __init__(self, pdb_id, chain, prob, e_value, query_start, query_end, 
               template_start, template_end):
   self.prob=prob
   self.pdb_id=pdb_id
   self.chain=chain
   self.e_value=e_value   
   self.query_start=query_start
   self.query_end=query_end
   self.template_start=template_start
   self.template_end=template_end

class HHSearchResult:
  """
  Read HHSearch result file. The result is stored in a list of hh search hits.
  
  Usage:
  
  result=HHSearchResult('output.hhr')
  for hit in result.hits:
    print hit.pdb_id, hit.chain
    print hit.alignment.ToString(80)
  """
  def __init__(self, filename, pipe_separated=False):
    self.pipe_separated=pipe_separated    
    self._Read(filename)
  def _Read(self, filename):
    ifile=open(filename)
    self._ReadHeader(ifile)
    summaries=self._ReadHitSummaries(ifile)
    self.hits=self._ReadHitDetails(ifile, summaries)
    
  def _ReadHeader(self, ifile):
    header_lines=[]
    for line in ifile:
      stripped_line=line.strip()
      if stripped_line=='':
        break
      header_lines.append(stripped_line)
    for header_line in header_lines:
      var, value=re.split('\s+', header_line, 1)
      setattr(self, var.lower(), value)
      
  def _ReadHitSummaries(self, ifile):
    summary_lines=[]
    skip_header=True
    for line in ifile:
      if skip_header==True:
        skip_header=False
        continue

      stripped_line=line.strip()
      if stripped_line=='':
        break
      summary_lines.append(line)
    summaries=[]
    for summary_line in summary_lines:
      pdb_id, chain=(None, None)
      if self.pipe_separated:
        parts=summary_line[4:37].split('|')
        pdb_id=parts[1][:4]
        chain=parts[1][4]
      else:
        pdb_id=summary_line[4:8]
        chain=summary_line[9]
      prob=float(summary_line[36:40])
      e_value=0.0
      query_range=summary_line[76:84].split('-')
      query_start=int(query_range[0].strip())
      query_end=int(query_range[1].strip())
      template_range=summary_line[86:94].split('-')
      template_start=int(template_range[0].strip())
      template_end=int(template_range[1].strip())
      summaries.append(HitSummary(pdb_id, chain, prob, e_value, query_start,
                                  query_end, template_start, template_end))
    return summaries

  def _ReadHitDetails(self, ifile, summaries):
    hits=[]
    for summary in summaries:
      alignment=self._ReadHitDetail(ifile)
      hits.append(HHSearchHit(summary, alignment))
    return hits
  def _ReadHitDetail(self, ifile):
    skip_header=True
    q_seq, t_seq=('', '')
    for line in ifile:
      if skip_header==True:
        if line.startswith('>'):
          skip_header=False
        continue
      if line.startswith('No'):
        break
      if line.strip()=='':
        continue        
      if line.startswith('Q ss_pred') or line.startswith('Q Consensus'):
        continue
      if line.startswith('T ss_pred') or line.startswith('T Consensus'):
        continue
      if line.startswith(' '):
        continue
      if line.startswith('Q'):
        q_seq+=re.split('\s+', line)[3]
      if line.startswith('T'):
        t_seq+=re.split('\s+', line)[3]        
    ali=seq.AlignmentHandle()
    ali.AddSequence(seq.Sequence.FromString('query', q_seq))
    ali.AddSequence(seq.Sequence.FromString('target', t_seq))    
    return ali
