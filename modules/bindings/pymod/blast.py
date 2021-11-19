from ost import settings
import subprocess
from xml.dom import minidom
from ost import io, seq
import ost
import os, re, sys

class AlignedPatch:
  """
  An aligned patch, aka. HSP

  .. attribute:: aln

    The local alignment. Sequence offset of both sequences in the alignment are
    set to the starting position in the query and target sequence, respectively.

    :type: :class:`~ost.seq.AlignmentHandle`

  .. attribute:: bit_score

    The bit score of the HSP

  .. attribute:: score

    The score of the HSP

  .. attribute:: evalue

    The E-value of the HSP

  .. attribute:: seqid

    The sequence identity of the HSP
  """
  def __init__(self, aln, bit_score, score, evalue, seqid):
    self.aln=aln
    self.bit_score=bit_score
    self.score=score
    self.evalue=evalue
    self.seqid=seqid

class BlastHit:
  """
  A positive match found by BLAST.
  
  Each blast hit consists of one or more HSPs, encoded by the
  :class:`AlignedPatch` class.
  
  .. attribute:: identifier

    The identifier of the matching sequence

  .. attribute:: aligned_patches

    list of :class:`AlignedPatch` instances holding the actual HSPs.
  """
  def __init__(self, identifier, aligned_patches):
    self.identifier=identifier
    self.aligned_patches=aligned_patches


def ParseBlastOutput(string, seqid_thres=0, evalue_thres=float("infinity")):
  """
  Parses the blast output and returns a list of BlastHits
  setting no seqid_thres or evalue_thres, restores default behaviour without filtering
  """
  def _GetText(node):
    rc=''
    for child in node.childNodes:
        if child.nodeType==child.TEXT_NODE:
            rc+=child.data
    return rc

  def _GetValue(node, tag_name):
    tag=node.getElementsByTagName(tag_name)
    assert len(tag)==1
    return _GetText(tag[0])

  def _GetInt(node, tag_name):
    return int(_GetValue(node, tag_name))

  def _ParseHsp(query_id, hit_id, hsp, tot_query_len, seqid_thres=0, evalue_thres=float("infinity")):
    bit_score=float(_GetValue(hsp, 'Hsp_bit-score'))
    score=float(_GetValue(hsp, 'Hsp_score'))
    evalue=float(_GetValue(hsp, 'Hsp_evalue'))
    try:
      identity=float(_GetValue(hsp, 'Hsp_identity'))
    except AssertionError:
      # The Hsp_identity tag is not a 'must' in the BLAST XML format. It
      # describes the number of matching characters. Hence we assume, if it is
      # missing, there are 0 matches.
      identity=0
    hsp_align_len=float(_GetValue(hsp, 'Hsp_align-len'))
    seqid=identity/hsp_align_len
    query_offset=_GetInt(hsp, 'Hsp_query-from')-1
    hit_offset=_GetInt(hsp, 'Hsp_hit-from')-1
    query_seq=seq.CreateSequence(str(query_id), str(_GetValue(hsp, 'Hsp_qseq')))
    query_seq.offset=query_offset
    hit_seq=seq.CreateSequence(str(hit_id), str(_GetValue(hsp, 'Hsp_hseq')))
    hit_seq.offset=hit_offset
    try:
      if seqid > float(seqid_thres) and evalue < evalue_thres:
        aln=seq.CreateAlignment(query_seq, hit_seq)
        return AlignedPatch(aln, bit_score, score, evalue, seqid)

    except Exception as e:
      print(str(e), query_seq, hit_seq)

  try:
    doc=minidom.parseString(string)
  except Exception as e:
    ost.LogError('Error while parsing BLAST output: %s' % str(e))
    return None
  hits=[]
  query_id=_GetValue(doc, 'BlastOutput_query-def')
  tot_query_len=_GetValue(doc, 'BlastOutput_query-len')
  for hit in doc.getElementsByTagName('Hit'):
    hit_id=_GetValue(hit, 'Hit_def')
    hsps=hit.getElementsByTagName('Hsp')
    aligned_patches=[_ParseHsp(query_id, hit_id, hsp, tot_query_len, seqid_thres, evalue_thres) for hsp in hsps]
    hits.append(BlastHit(hit_id, aligned_patches))
  return hits



class BlastError(RuntimeError):
  def __init__(self, brief, details):
    self.brief=brief
    self.details=details

  def __str__(self):
    if self.details:
      return '%s\n%s' % (self.brief, self.details)
    else:
      return self.brief

def CreateDB(infasta, dbout, mkdb_cmd=None):
  """
  Create a blast DB from a fasta file

  :param infasta: the pdb fasta from which the database will be created
  :type infasta: :class:`string`

  :param dbout: output location for blastDB file
  :type dbout: :class:`string`


  """
  if mkdb_cmd==None:
    try:
      exe=settings.Locate('formatdb')
      args=[exe, '-i', infasta, '-n', dbout]
    except:
      try:
        exe=settings.Locate('makeblastdb')
        args=[exe, '-in', infasta, '-out', dbout, '-dbtype', 'prot']
      except:
        raise RuntimeError('could not find makeblastdb/formatdb executable')
  else:
    if os.path.basename(mkdb_cmd)=='makeblastdb':
      exe=settings.Locate('makeblastdb',explicit_file_name=mkdb_cmd)
      args=[exe, '-in', infasta, '-out', dbout, '-dbtype', 'prot']
    elif os.path.basename(mkdb_cmd)=='formatdb':
        exe=settings.Locate('formatdb',explicit_filename=mkdb_cmd)
        args=[exe, '-i', infasta, '-n', dbout]
    else:
      raise IOError('mkdb command must either be the path to formatdb or makeblastdb!')

  ost.LogInfo('creating blast DB (%s)' % ' '.join(args))
  blast_pipe=subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
  blast_pipe.communicate()

def BlastVersion(blast_location=None):
  """
  Returns the version of the BLAST executable, e.g. 2.2.24 as a string
  """

  try:
    blast_exe=settings.Locate('blastp',explicit_file_name=blast_location)
  except:
    try:
      blast_exe=settings.Locate('blastall', explicit_file_name=blast_location)
    except:
      raise RuntimeError('could not find blast executable')

  if os.path.basename(blast_exe)=='blastall':
    args=[blast_exe]
    pattern=re.compile(r'\s*blastall (\d+\.\d+\.\d+)\s+arguments:\s*')

  else:
    args=[blast_exe, '-version']
    pattern=re.compile(r'\s*Package: blast (\d+\.\d+\.\d+),\s+')

  blast_pipe=subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
  stdout, _ = blast_pipe.communicate()
  lines=stdout.decode().splitlines()

  for line in lines:
    m=pattern.match(line)
    if m:
      return m.group(1)
  raise IOError("could not determine blast version for '%s'" % blast_exe)



def Blast(query, database, gap_open=11, gap_ext=1, matrix='BLOSUM62',
         blast_location=None, outfmt=0, filter_low_complexity=True):
  """
  Runs a protein vs. protein blast search. The results are returned
  according to the value of the ``outfmt`` parameter.

  :param query: the query sequence
  :type query: :class:`seq.ConstSequenceHandle`

  :param database: The filename of the sequence database. Make sure that
      formatdb has been run on the database and the <database>.pin file exists.
  :param matrix: The substitution matrix to be used. Must be one of 'BLOSUM45',
     'BLOSUM62', 'BLOSUM80', 'PAM30', 'PAM70'.
  :param gap_open: Gap opening penalty. Note that only a subset of gap opening
     penalties is supported for each substitutition matrix. Consult the blast
     docs for more information.
  :param gap_ext: Gap extension penalty. Only a subset of gap extension
     penalties are supported for each of the substitution matrices. Consult the
     blast docs for more information.
  :param outfmt: output format, where '0' corresponds to default output (parsed 
     blast output and 1 to raw string output).
  :param filter_low_complexity: Mask off segments of the query sequence that 
     have low compositional complexity, as determined by the SEG program of 
     Wootton & Federhen (Computers and Chemistry, 1993)
  :rtype: :class:`BlastHit` (with ``outfmt=0``) or :class:`str` 
     (with ``outfmt=1``)
  """
  subst_mats=('BLOSUM45', 'BLOSUM62', 'BLOSUM80', 'PAM30', 'PAM70',)
  if matrix not in subst_mats:
    raise ValueError('matrix must be one of %s' % ', '.join(subst_mats))
  if not os.path.exists('%s.pin' % database) and not os.path.exists('%s.pal' % database):
    raise IOError("Database %s does not exist" % database)
  if blast_location!=None and not os.path.exists(blast_location):
    ost.LogScript('Could not find %s' %blast_location)

  if blast_location==None:
    try:
      blast_exe=settings.Locate('blastp')
    except:
      try:
        blast_exe=settings.Locate('blastall')
      except:
        raise RuntimeError('could not find blast executable')
  else:
    blast_exe=settings.Locate(os.path.basename(blast_location),explicit_file_name=blast_location)

  if os.path.basename(blast_exe)=='blastall':
    args=[blast_exe, '-d', database, '-p', 'blastp',
          '-m', '7', '-M', matrix, '-G', str(gap_open), '-E', str(gap_ext)]
    if filter_low_complexity==False:
      args.append('-F')
      args.append('F')
          
  else:
    complexity_opt='-seg'
    if filter_low_complexity==True:
      complexity_arg='yes'
    else:
      complexity_arg='no'
    args=[blast_exe, '-db', database, '-matrix', matrix,
          '-gapopen', str(gap_open), '-gapextend', str(gap_ext), '-outfmt', '5', complexity_opt, complexity_arg ]

  ost.LogInfo('running BLAST (%s)' % ' '.join(args))
  blast_pipe=subprocess.Popen(args, stderr=subprocess.PIPE,
                              stdout=subprocess.PIPE, stdin=subprocess.PIPE)
  if isinstance(query, str):
    stdout, stderr=blast_pipe.communicate(query.encode())
  else:
    stdout, stderr=blast_pipe.communicate(io.SequenceToString(query, 'fasta').encode())

  if len(stderr)>0:
     pattern=re.compile(r'^\[.*\]\s+ERROR:\s+(.*)')
     lines=stderr.decode().split('\n')
     error_message=pattern.match(lines[0])
     if error_message:
       raise BlastError(error_message.group(1), '\n'.join(lines[1:]))
  if outfmt==0:
    return ParseBlastOutput(stdout.decode())
  else:
    return stdout.decode()
