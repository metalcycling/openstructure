'''HHblits wrapper classes and functions.
'''

import subprocess
import datetime
import os
import shutil
import tempfile

import ost
from ost import settings, seq
from ost.bindings import utils

class HHblitsHit:
    """
    A hit found by HHblits

    .. attribute:: hit_id

      String identifying the hit

      :type: :class:`str`

    .. attribute:: aln

      Pairwise alignment containing the aligned part between the query and the
      target. First sequence is the query, the second sequence the target.

      :type: :class:`~ost.seq.AlignmentHandle`

    .. attribute:: score

      The alignment score

      :type: :class:`float`

    .. attribute:: ss_score

      The secondary structure score

      :type: :class:`float`

    .. attribute:: evalue

      The E-value of the alignment

      :type: :class:`float`

    .. attribute:: pvalue

      The P-value of the alignment

      :type: :class:`float`

    .. attribute:: prob

      The probability of the alignment (between 0 and 100)

      :type: :class:`float`
    """
    def __init__(self, hit_id, aln, score, ss_score, evalue, pvalue, prob):
        self.hit_id = hit_id
        self.aln = aln
        self.score = score
        self.ss_score = ss_score
        self.evalue = evalue
        self.prob = prob
        self.pvalue = pvalue

class HHblitsHeader:
    """Stats from the beginning of search output.

    .. attribute:: query

      The name of the query sequence

      :type: :class:`str`

    .. attribute:: match_columns

      Total of aligned Match columns

      :type: :class:`int`

    .. attribute:: n_eff

      Value of the ``-neff`` option

      :type: :class:`float`

    .. attribute:: searched_hmms

      Number of profiles searched

      :type: :class:`int`

    .. attribute:: date

      Execution date

      :type: :class:`datetime.datetime`

    .. attribute:: command

      Command used to run

      :type: :class:`str`
    """
    def __init__(self):
        self.query = ''
        self.match_columns = 0
        self.n_eff = 0
        self.searched_hmms = 0
        self.date = None
        self.command = ''

def ParseHeaderLine(line):
    '''Fetch header content.

    First, we seek the start of the identifier, that is, the first whitespace
    after the hit number + 1. Since the identifier may contain whitespaces
    itself, we cannot split the whole line

    :param line: Line from the output header.
    :type line: :class:`str`

    :return: Hit information and query/template offsets
    :rtype: (:class:`HHblitsHit`, (:class:`int`, :class:`int`))
    '''
    for i in range(0, len(line)):
        if line[i].isdigit():
            break
    for i in range(i, len(line)):
        if line[i] == ' ':
            break
    assert len(line)-i >= 31 and line[i+1] != ' '
    hit_id = line[i+1:i+31].strip()
    fields = line[i+32:].split()
    prob = float(fields[0])
    evalue = float(fields[1])
    pvalue = float(fields[2])
    score = float(fields[3])
    ss_score = float(fields[4])
    offsets = (int(fields[6].split('-')[0]), int(fields[7].split('-')[0]))
    return (HHblitsHit(hit_id, None, score, ss_score, evalue, pvalue, prob),
            offsets)

def ParseHHblitsOutput(output):
    """
    Parses the HHblits output as produced by :meth:`HHblits.Search` and returns
    the header of the search results and a list of hits.

    :param output: Iterable containing the lines of the HHblits output file
    :type output: iterable (e.g. an open file handle)

    :return: a tuple of the header of the search results and the hits
    :rtype: (:class:`HHblitsHeader`, :class:`list` of :class:`HHblitsHit`)
    """
    lines = iter(output)
    def _ParseHeaderSection(lines):
        value_start_column = 14
        date_pattern = '%a %b %d %H:%M:%S %Y'
        header = HHblitsHeader()
        line = next(lines)
        assert line.startswith('Query')
        header.query = line[value_start_column:].strip()
        line = next(lines)
        assert line.startswith('Match_columns')
        header.match_columns = int(line[value_start_column:].strip())

        line = next(lines)
        assert line.startswith('No_of_seqs')

        line = next(lines)
        assert line.startswith('Neff')
        header.n_eff = float(line[value_start_column:].strip())

        line = next(lines)
        assert line.startswith('Searched_HMMs')
        header.searched_hmms = int(line[value_start_column:].strip())

        line = next(lines)
        assert line.startswith('Date')
        value = line[value_start_column:].strip()
        header.date = datetime.datetime.strptime(value, date_pattern)

        line = next(lines)
        assert line.startswith('Command')
        header.command = line[value_start_column:].strip()

        line = next(lines)
        assert len(line.strip()) == 0
        return header

    def _ParseTableOfContents(lines):
        line = next(lines)
        assert line.startswith(' No Hit')
        hits = []
        while True:
            line = next(lines)
            if len(line.strip()) == 0:
                return hits
            hits.append(ParseHeaderLine(line))
        return hits

    def _ParseResultBody(query_id, hits, lines):
        entry_index = None
        query_str = ''
        templ_str = ''
        def _MakeAln(query_id, hit_id, query_string, templ_string,
                     q_offset, t_offset):
            s1 = seq.CreateSequence(query_id, query_string)
            s1.offset = q_offset-1
            s2 = seq.CreateSequence(hit_id, templ_string)
            s2.offset = t_offset-1
            return seq.CreateAlignment(s1, s2)
        try:
            while True:
                # Lines which we are interested in:
                # - "Done!" -> end of list
                # - "No ..." -> next item in list
                # - "T <hit_id> <start> <data> <end>"
                # - "Q <query_id> <start> <data> <end>"
                # -> rest is to be skipped
                line = next(lines)
                if len(line.strip()) == 0:
                    continue
                if line.startswith('Done!'):
                    if len(query_str) > 0:
                        hits[entry_index][0].aln = _MakeAln(\
                            query_id, hits[entry_index][0].hit_id,
                            query_str, templ_str, *hits[entry_index][1])
                    return [h for h, o in hits]
                if line.startswith('No '):
                    if len(query_str) > 0:
                        hits[entry_index][0].aln = _MakeAln(\
                            query_id, hits[entry_index][0].hit_id,
                            query_str, templ_str, *hits[entry_index][1])
                    entry_index = int(line[3:].strip())-1
                    line = next(lines)
                    hits[entry_index][0].hit_id = line[1:].strip()
                    query_str = ''
                    templ_str = ''
                    # skip the next line. It doesn't contain information we
                    # don't already know
                    next(lines)
                    continue
                assert entry_index != None
                # Skip all "T ..." and "Q ..." lines besides the one we want
                if line[1:].startswith(' Consensus'):
                    continue
                if line[1:].startswith(' ss_pred'):
                    continue
                if line[1:].startswith(' ss_conf'):
                    continue
                if line[1:].startswith(' ss_dssp'):
                    continue
                if line.startswith('T '):
                    for start_pos in range(22, len(line)):
                        if line[start_pos].isalpha() or line[start_pos] == '-':
                            break
                    end_pos = line.find(' ', start_pos)
                    # this can fail if we didn't skip all other "T ..." lines
                    if end_pos == -1:
                        error_str = "Unparsable line '%s' for entry No %d" \
                                    % (line.strip(), entry_index + 1)
                        raise AssertionError(error_str)
                    templ_str += line[start_pos:end_pos]
                if line.startswith('Q '):
                    for start_pos in range(22, len(line)):
                        if line[start_pos].isalpha() or line[start_pos] == '-':
                            break
                    end_pos = line.find(' ', start_pos)
                    # this can fail if we didn't skip all other "Q ..." lines
                    if end_pos == -1:
                        error_str = "Unparsable line '%s' for entry No %d" \
                                    % (line.strip(), entry_index + 1)
                        raise AssertionError(error_str)
                    query_str += line[start_pos:end_pos]
        except StopIteration:
            if len(query_str) > 0:
                hits[entry_index][0].aln = _MakeAln(query_id,
                                                    hits[entry_index][0].hit_id,
                                                    query_str, templ_str,
                                                    *hits[entry_index][1])
            return [h for h, o in hits]
    header = _ParseHeaderSection(lines)
    # parse the table of contents. This is neccessary as some of the properties
    # (i.e. start of alignment) we need are only given there. From the TOC we
    # create a list of hits that is then further filled with data when we parse
    # the actual result body
    hits = _ParseTableOfContents(lines)
    return header, _ParseResultBody(header.query, hits, lines)

def ParseA3M(a3m_file):
    '''
    Parse secondary structure information and the multiple sequence alignment 
    out of an A3M file as produced by :meth:`HHblits.BuildQueryMSA`.
    
    :param a3m_file: Iterable containing the lines of the A3M file
    :type a3m_file: iterable (e.g. an open file handle)
    
    :return: Dictionary containing "ss_pred" (:class:`list`), "ss_conf"
             (:class:`list`) and "msa" (:class:`~ost.seq.AlignmentHandle`).
             If not available, "ss_pred" and "ss_conf" entries are set to None.
    '''
    profile_dict = dict()
    state = 'NONE'
    pred_seq_txt = ''
    conf_seq_txt = ''
    msa_seq = list()
    msa_head = list()
    for line in a3m_file:
        if len(line.rstrip()) == 0:
            continue
        elif line.startswith('>ss_pred'):
            state = 'sspred'
            continue
        elif line.startswith('>ss_conf'):
            state = 'ssconf'
            continue
        elif line[0] == '>':
            msa_seq.append('')
            msa_head.append(line[1:].rstrip())
            state = 'msa'
            continue

        if state == 'sspred':
            pred_seq_txt += line.rstrip()
        elif state == 'ssconf':
            conf_seq_txt += line.rstrip()
        elif state == 'msa':
            msa_seq[len(msa_seq)-1] += line.rstrip()

    if len(pred_seq_txt) > 0:
        profile_dict['ss_pred'] = list()
        profile_dict['ss_conf'] = list()
        for i in range(0, len(pred_seq_txt)):
            profile_dict['ss_pred'].append(pred_seq_txt[i])
            profile_dict['ss_conf'].append(int(conf_seq_txt[i]))
    else:
        profile_dict['ss_pred'] = None
        profile_dict['ss_conf'] = None    

    # post processing
    # MSA
    profile_dict['msa'] = None
    if len(msa_seq) > 1:
        t = msa_seq[0]
        al = seq.AlignmentList()
        for i in range(1, len(msa_seq)):
            qs = ''
            ts = ''
            k = 0
            for c in msa_seq[i]:
                if c.islower():
                    qs += '-'
                    ts += c.upper()
                else:
                    qs += t[k]
                    ts += c
                    k += 1
            nl = seq.CreateAlignment(seq.CreateSequence(msa_head[0], qs), 
                                     seq.CreateSequence(msa_head[i], ts))
            al.append(nl)
        profile_dict['msa'] = seq.alg.MergePairwiseAlignments(\
            al, seq.CreateSequence(msa_head[0], t))
    return profile_dict


def ParseHHM(profile):
    '''
    Parse secondary structure information and the MSA out of an HHM profile as
    produced by :meth:`HHblits.A3MToProfile`.

    :param profile: Opened file handle holding the profile.
    :type profile: :class:`file`

    :return: Dictionary containing "ss_pred" (:class:`list`), "ss_conf"
             (:class:`list`), "msa" (:class:`~ost.seq.AlignmentHandle`) and
             "consensus" (:class:`~ost.seq.SequenceHandle`).
             If not available, "ss_pred" and "ss_conf" entries are set to None.
    '''
    profile_dict = dict()
    state = 'NONE'
    pred_seq_txt = ''
    conf_seq_txt = ''
    consensus_txt = ''
    msa_seq = list()
    msa_head = list()
    for line in profile:
        if len(line.rstrip()) == 0:
            continue
        if line.rstrip() == '>ss_pred PSIPRED predicted secondary structure':
            state = 'sspred'
            continue
        elif line.rstrip() == '>ss_conf PSIPRED confidence values':
            state = 'ssconf'
            continue
        elif line.rstrip() == '>Consensus':
            state = 'consensus'
            continue
        elif line[0] == '>':
            if state == 'consensus' or state == 'msa':
                msa_seq.append('')
                msa_head.append(line[1:].rstrip())
            else:
                raise IOError('Profile file "%s" is missing ' % profile.name+
                              'the "Consensus" section')
            state = 'msa'
            continue
        elif line[0] == '#':
            state = 'NONE'
            continue

        if state == 'sspred':
            pred_seq_txt += line.rstrip()
        elif state == 'ssconf':
            conf_seq_txt += line.rstrip()
        elif state == 'msa':
            msa_seq[len(msa_seq)-1] += line.rstrip()
        elif state == 'consensus':
            consensus_txt += line.rstrip()

    if len(pred_seq_txt) > 0:
        profile_dict['ss_pred'] = list()
        profile_dict['ss_conf'] = list()
        for i in range(0, len(pred_seq_txt)):
            profile_dict['ss_pred'].append(pred_seq_txt[i])
            profile_dict['ss_conf'].append(int(conf_seq_txt[i]))
    else:
        profile_dict['ss_pred'] = None
        profile_dict['ss_conf'] = None

    # post processing
    # MSA
    profile_dict['msa'] = None
    if len(msa_seq):
        t = msa_seq[0]
        al = seq.AlignmentList()
        for i in range(1, len(msa_seq)):
            qs = ''
            ts = ''
            k = 0
            for c in msa_seq[i]:
                if c.islower():
                    qs += '-'
                    ts += c.upper()
                else:
                    qs += t[k]
                    ts += c
                    k += 1
            nl = seq.CreateAlignment(seq.CreateSequence(msa_head[0], qs),
                                     seq.CreateSequence(msa_head[i], ts))
            al.append(nl)
        profile_dict['msa'] = seq.alg.MergePairwiseAlignments(\
            al, seq.CreateSequence(msa_head[0], t))
    # Consensus
    profile_dict['consensus'] = seq.CreateSequence('Consensus', consensus_txt)

    return profile_dict


class HHblits:
    """
    Initialise a new HHblits "search" for the given query. Query may either
    be a :class:`~ost.seq.SequenceHandle` or a string. In the former case, the
    query is the actual query sequence, in the latter case, the query is the
    filename to the file containing the query.

    :param query: Query sequence as file or sequence.
    :type query: :class:`~ost.seq.SequenceHandle` or :class:`str`
    :param hhsuite_root: Path to the top-level directory of your hhsuite
                         installation.
    :type hhsuite_root: :class:`str`
    :param hhblits_bin: Name of the hhblits binary. Will only be used if
                        :attr:`hhsuite_root`:file:`/bin/hhblits` does not exist.
    :type hhblits_bin: :class:`str`
    :param working_dir: Directory for temporary files. Will be created if not
                        present but **not** automatically deleted.
    :type working_dir: :class:`str`
    """
    OUTPUT_PREFIX = 'query_hhblits'
    def __init__(self, query, hhsuite_root, hhblits_bin=None, working_dir=None):
        self.query = query
        self.hhsuite_root = hhsuite_root
        if os.path.exists(os.path.join(self.hhsuite_root, 'bin/hhblits')):
            self.bin_dir = os.path.join(self.hhsuite_root, 'bin')
            self.hhblits_bin = os.path.join(self.hhsuite_root, 'bin/hhblits')
        else:
            self.hhblits_bin = settings.Locate('hhblits',
                                               explicit_file_name=hhblits_bin)
            self.bin_dir = os.path.dirname(self.hhblits_bin)
            # guess root folder (note: this may fail in future)
            self.hhsuite_root = os.path.dirname(self.bin_dir)
                                             
        if working_dir:
            self.needs_cleanup = False
            self.working_dir = working_dir
            if not os.path.exists(working_dir):
                os.mkdir(working_dir)
            if isinstance(query, str):
                self.filename = os.path.abspath(os.path.join(
                    self.working_dir, os.path.basename(query)))
                if self.filename != os.path.abspath(query):
                    shutil.copy(query, self.filename)
            else:
                self.filename = os.path.join(self.working_dir,
                                             '%s.fasta' % HHblits.OUTPUT_PREFIX)
                ost.io.SaveSequence(query, self.filename)
        else:
            self.needs_cleanup = True
            if isinstance(query, str):
                self.working_dir = tempfile.mkdtemp()
                self.filename = os.path.abspath(os.path.join(
                    self.working_dir, os.path.basename(query)))
                shutil.copy(query, self.filename)
            else:
                tmp_dir = utils.TempDirWithFiles((query,))
                self.working_dir = tmp_dir.dirname
                self.filename = tmp_dir.files[0]

    def BuildQueryMSA(self, nrdb, options={}, a3m_file=None, assign_ss=True):
        """Builds the MSA for the query sequence.

        The produced A3M file can be parsed by :func:`ParseA3M`. If the file was
        already produced, hhblits is not called again and the existing file path
        is returned (neglecting the *assign_ss* flag!!!).

        :param nrdb: Database to be align against; has to be an hhblits database
        :type nrdb: :class:`str`

        :param options: Dictionary of options to *hhblits*, one "-" is added in
                        front of every key. Boolean True values add flag without
                        value. Merged with default options 
                        {'cpu': 1, 'n': 1, 'e': 0.001}, where 'n' defines the 
                        number of iterations and 'e' the E-value cutoff for 
                        inclusion of sequences in result alignment.
        :type options: :class:`dict`

        :param a3m_file: a path of a3m_file to be used, optional
        :type a3m_file: :class:`str`

        :param assign_ss: HHblits does not assign predicted secondary structure 
                          by default. You can optionally assign it with the
                          addss.pl script provided by the HH-suite. However, 
                          your HH-suite installation requires you to specify
                          paths to PSIRED etc. We refer to the HH-suite user
                          guide for further instructions. Assignment is done
                          by calling :func:`HHblits.AssignSSToA3M`

        :type assign_ss:  :class:`bool`

        :return: The path to the A3M file containing the MSA
        :rtype: :class:`str`
        """
        if a3m_file is None:
            a3m_file = '%s.a3m' % os.path.splitext(self.filename)[0]
        else:
            a3m_file = os.path.abspath(a3m_file)
        if os.path.exists(a3m_file):
            ost.LogInfo('Reusing already existing query alignment (%s)' % a3m_file)
            return a3m_file
        ost.LogInfo('Using hhblits from "%s"' % self.hhsuite_root)
        full_nrdb = os.path.join(os.path.abspath(os.path.split(nrdb)[0]),
                                 os.path.split(nrdb)[1])
        # create MSA
        opts = {'cpu' : 1,     # no. of cpus used
                'n'   : 1,     # no. of iterations
                'e'   : 0.001} # evalue threshold
        opts.update(options)
        opt_cmd, _ = _ParseOptions(opts)
        hhblits_cmd = '%s -i %s -oa3m %s -d %s %s' % \
                      (self.hhblits_bin, self.filename, a3m_file, full_nrdb,
                       opt_cmd)

        p = subprocess.run(hhblits_cmd, shell=True, cwd=self.working_dir,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        lines = p.stdout.decode().splitlines()
        for line in lines:
            ost.LogVerbose(line.strip())

        lines = p.stderr.decode().splitlines()
        for line in lines:
            ost.LogError(line.strip())

        if not os.path.exists(a3m_file):
            raise RuntimeError('Building query profile failed, no output')

        if assign_ss:
            return self.AssignSSToA3M(a3m_file)
        else:
            return a3m_file

    def AssignSSToA3M(self, a3m_file):
        """
        HHblits does not assign predicted secondary structure by default. 
        You can optionally assign it with the addss.pl script provided by the 
        HH-suite. However, your HH-suite installation requires you to specify
        paths to PSIRED etc. We refer to the HH-suite user guide for further 
        instructions.

        :param a3m_file:  Path to file you want to assign secondary structure to
        :type a3m_file:  :class:`str`
        """

        a3m_file = os.path.abspath(a3m_file)
        addss_cmd = "perl %s %s" % (os.path.join(self.hhsuite_root, 
                                                 'scripts/addss.pl'), 
                                    a3m_file)
        env = dict(os.environ)
        env.update({'PERL5LIB' : os.path.join(self.hhsuite_root, 'scripts'),
                    'HHLIB' : os.path.join(self.hhsuite_root),
                    'PATH' : '%s:%s' % (os.path.join(self.hhsuite_root, 'bin'),
                                        os.environ['PATH'])})

        p = subprocess.run(addss_cmd, shell=True, cwd=self.working_dir,
                           env=env, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)

        lines = p.stdout.decode().splitlines()
        for line in lines:
            ost.LogVerbose(line.strip())
            if 'error' in line.lower() or 'bad interpreter' in line.lower():
                raise RuntimeError('Predicting secondary structure for MSA '+
                                   '(%s) failed, on command: %s' % (a3m_file, line))

        lines = p.stderr.decode().splitlines()
        for line in lines:
            ost.LogError(line.strip())
            if 'error' in line.lower() or 'bad interpreter' in line.lower():
                raise RuntimeError('Predicting secondary structure for MSA '+
                                   '(%s) failed, on command: %s' % (a3m_file, line))

        return a3m_file

    def A3MToProfile(self, a3m_file, hhm_file=None):
        """
        Converts the A3M alignment file to a hhm profile. If hhm_file is not
        given, the output file will be set to <:attr:`a3m_file`-basename>.hhm.

        The produced HHM file can be parsed by :func:`ParseHHM`.

        If the file was already produced, the existing file path is returned
        without recomputing it.

        :param a3m_file: Path to input MSA as produced by :meth:`BuildQueryMSA`
        :type a3m_file: :class:`str`

        :param hhm_file: Desired output file name 
        :type hhm_file: :class:`str`

        :return: Path to the profile file
        :rtype: :class:`str`
        """
        hhmake = os.path.join(self.bin_dir, 'hhmake')
        if not hhm_file:
            hhm_file = '%s.hhm' % os.path.splitext(a3m_file)[0]
        if os.path.exists(hhm_file):
            return hhm_file
        ost.LogVerbose('converting %s to %s' % (a3m_file, hhm_file))
        p = subprocess.run('%s -i %s -o %s' % (hhmake, a3m_file, hhm_file),
                           shell=True, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
        lines = p.stdout.decode().splitlines()
        for line in lines:
            ost.LogVerbose(line.strip())
        lines = p.stderr.decode().splitlines()
        for line in lines:
            ost.LogError(line.strip())

        if p.returncode != 0:
            raise IOError('could not convert a3m to hhm file')

        if not os.path.exists(hhm_file):
            raise RuntimeError('could not convert a3m to hhm file, no output')

        return hhm_file

    def A3MToCS(self, a3m_file, cs_file=None, options={}):
        """
        Converts the A3M alignment file to a column state sequence file. If
        cs_file is not given, the output file will be set to
        <:attr:`a3m_file`-basename>.seq219.

        If the file was already produced, the existing file path is returned
        without recomputing it.

        :param a3m_file: Path to input MSA as produced by :meth:`BuildQueryMSA`
        :type a3m_file: :class:`str`

        :param cs_file: Output file name (may be omitted)
        :type cs_file: :class:`str`

        :param options: Dictionary of options to *cstranslate*, one "-" is added
                        in front of every key. Boolean True values add flag
                        without value.
        :type options: :class:`dict`

        :return: Path to the column state sequence file
        :rtype: :class:`str`
        """
        cstranslate = os.path.join(self.bin_dir, 'cstranslate')
        if not cs_file:
            cs_file = '%s.seq219' % os.path.splitext(a3m_file)[0]
        if os.path.exists(cs_file):
            return cs_file
        opt_cmd, _ = _ParseOptions(options)
        cs_cmd = '%s -i %s -o %s %s' % (
            cstranslate,
            os.path.abspath(a3m_file),
            os.path.abspath(cs_file),
            opt_cmd)
        ost.LogVerbose('converting %s to %s' % (a3m_file, cs_file))
        p = subprocess.run(cs_cmd, shell=True, cwd=self.working_dir,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        if not os.path.exists(cs_file):
            raise RuntimeError('Creating column state sequence file failed, ' +
                               'no output')

        if b'Wrote abstract state sequence to' in p.stdout:
            return cs_file
        else:
            raise RuntimeError('Creating column state sequence file failed')

    def Cleanup(self):
        """Delete temporary data.

        Delete temporary data if no working dir was given. Controlled by
        :attr:`needs_cleanup`.
        """
        if self.needs_cleanup and os.path.exists(self.working_dir):
            shutil.rmtree(self.working_dir)

    def CleanupFailed(self):
        '''In case something went wrong, call to make sure everything is clean.

        This will delete the working dir independently of :attr:`needs_cleanup`.
        '''
        store_needs_cleanup = self.needs_cleanup
        self.needs_cleanup = True
        self.Cleanup()
        self.needs_cleanup = store_needs_cleanup

    def Search(self, a3m_file, database, options={}, prefix=''):
        """
        Searches for templates in the given database. Before running the search,
        the hhm file is copied. This makes it possible to launch several hhblits
        instances at once. Upon success, the filename of the result file is
        returned. This file may be parsed with :func:`ParseHHblitsOutput`.

        :param a3m_file: Path to input MSA as produced by :meth:`BuildQueryMSA`
        :type a3m_file: :class:`str`

        :param database: Search database, needs to be the common prefix of the
                         database files
        :type database: :class:`str`

        :param options: Dictionary of options to *hhblits*, one "-" is added in
                        front of every key. Boolean True values add flag without
                        value. Merged with default options {'cpu': 1, 'n': 1},
                        where 'n' defines the number of iterations.
        :type options: :class:`dict`

        :param prefix: Prefix to the result file
        :type prefix: :class:`str`

        :return: The path to the result file
        :rtype: :class:`str`
        """
        opts = {'cpu' : 1, # no. of cpus used
                'n'   : 1}   # no. of iterations
        opts.update(options)
        opt_cmd, opt_str = _ParseOptions(opts)
        base = os.path.basename(os.path.splitext(a3m_file)[0])
        hhr_file = '%s%s_%s.hhr' % (prefix, base, opt_str)
        hhr_file = os.path.join(self.working_dir, hhr_file)
        search_cmd = '%s %s -e 0.001 -Z 10000 -B 10000 -i %s -o %s -d %s' % (
            self.hhblits_bin,
            opt_cmd,
            os.path.abspath(a3m_file),
            os.path.abspath(hhr_file),
            os.path.join(os.path.abspath(os.path.split(database)[0]),
                         os.path.split(database)[1]))
        ost.LogInfo('searching %s' % database)
        ost.LogVerbose(search_cmd)
        p = subprocess.run(search_cmd, shell=True, cwd=self.working_dir,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        lines = p.stdout.decode().splitlines()
        for line in lines:
            ost.LogVerbose(line.strip())
        lines = p.stderr.decode().splitlines()
        for line in lines:
            ost.LogError(line.strip())

        if p.returncode != 0:
            raise RuntimeError('Sequence search failed')

        if not os.path.exists(hhr_file):
            raise RuntimeError('Sequence search failed, no output')

        return hhr_file


def _ParseOptions(opts):
    """
    :return: Tuple of strings (opt_cmd, opt_str), where opt_cmd can be
             passed to command ("-" added in front of keys, options
             separated by space) and opt_str (options separated by "_")
             can be used for filenames.
    :param opts: Dictionary of options, one "-" is added in front of every
                 key. Boolean True values add flag without value.
    """
    opt_cmd = list()
    opt_str = list()
    for k, val in opts.items():
        if type(val) == type(True):
            if val == True:
                opt_cmd.append('-%s' % str(k))
                opt_str.append(str(k))
        else:
            opt_cmd.append('-%s %s' % (str(k), str(val)))
            opt_str.append('%s%s' % (str(k), str(val)))
    opt_cmd = ' '.join(opt_cmd)
    opt_str = '_'.join(opt_str)
    return opt_cmd, opt_str


__all__ = ['HHblits', 'HHblitsHit', 'HHblitsHeader',
           'ParseHHblitsOutput', 'ParseA3M', 'ParseHHM',
           'ParseHeaderLine']

#  LocalWords:  HHblits MSA hhblits hhtools PSIPRED addss param nrdb str
#  LocalWords:  cpu hhm func ParseHHblitsOutput ss pred conf msa hhsuite dir
#  LocalWords:  attr basename rtype cstranslate tuple HHblitsHeader meth aln
#  LocalWords:  HHblitsHit iterable evalue pvalue neff hmms datetime
#  LocalWords:  whitespace whitespaces
