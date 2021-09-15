import os
from enum import Enum

from ost import GetSharedDataPath

def _StrToFloat(str_item):
    """Returns None if *str_item* looks like invalid number, casts to 
    float otherwise"""
    x = str_item.strip()
    if x.lower() in ["na", "-", "none"]:
        return None
    return float(x)


class AnnoType(str, Enum):
    """Possible types of aaindex entries"""
    SINGLE = 'SINGLE' #: Single amino acid annotation
    PAIR = 'PAIR' #: Pairwise amino acid annotation, substitution/pairwise scores


class AAIndexData:
    """Data object representing an annotation in aaindex, preferably 
    constructed from it's static :func:`Parse` method. The following 
    attributes are available:

    * key: aaindex accession number (e.g. ANDN920101)
    * desc: descriptive title
    * ref: Reference to article if available
    * authors: Authors of article if available
    * title: Title of article if available
    * journal: Journal of article if available
    * anno_type: Enum (:class:`AnnoType`) specifying whether we're dealing 
                 with a single or pairwise amino acid annotation/score.
    * anno: :class:`dict` with annotation. If *anno_type* is SINGLE,
            keys are amino acid one letter codes (single character strings).
            If *anno_type* is PAIR, keys are two one letter codes added 
            together (two character strings). Even when the thing is 
            symmetric, both keys exist. I.e. 'AB' AND 'BA'.
            Values are of type :class:`float` (None if not available).
    """
    def __init__(self):
        self.key = None
        self.desc = None
        self.ref = None
        self.authors = None
        self.title = None
        self.journal = None
        self.anno_type = AnnoType.SINGLE
        self.anno = dict()
    
    @staticmethod 
    def Parse(data):
        """Creates :class:`AAIndexData` from data.

        :param data: Iterable with strings in data format described for aaindex.
        :returns: :class:`AAIndexData`, if iterable contains several entries, 
                   parsing stops at separation sequence ('//'). None is returned
                   if nothing could be parsed.
        :raises: descriptive error in case of corrupt data
        """
        key = ""
        desc = ""
        ref = ""
        authors = ""
        title = ""
        journal = ""
        anno = dict()
        anno_type = None

        # temp variables, used for parsing 
        values = list()
        pair_rows = None
        pair_cols = None

        current_data_type = None
        stuff_read = False
        for line in data:
            if line.startswith("//"):
                break # we're done
            elif line.strip() == "":
                continue # nothing to read
            elif line[0] in ["H", "D", "R", "A", "T", "J", "I", "M"]:
                stuff_read = True
                current_data_type = line[0] # something we're interested in
            elif line.startswith(" "):
                pass # continuation of previous stuff
            else:
                current_data_type = None # unkown data, skip...
            
            if current_data_type == "H":
                key = line[2:].strip()
            elif current_data_type == "D":
                desc += line[2:]
            elif current_data_type == "R":
                ref += line[2:]
            elif current_data_type == "A":
                authors += line[2:]
            elif current_data_type == "T":
                title += line[2:]
            elif current_data_type == "J":
                journal += line[2:]
            elif current_data_type == "I":
                if anno_type == AnnoType.PAIR:
                    raise RuntimeError("Observed single AA and pairwise "
                                       "features in the same aaindex entry")
                anno_type = AnnoType.SINGLE
                if line.startswith("I"):
                    # header, must match expected aa ordering
                    aakeys = [item.strip() for item in line[1:].split()]
                    exp_aa_keys = ["A/L", "R/K", "N/M", "D/F", "C/P", "Q/S", 
                                   "E/T", "G/W", "H/Y", "I/V"]
                    if aakeys != exp_aa_keys:
                        raise RuntimeError(f"Keys in single AA AAIndex entry " 
                                           "are expected to be "
                                           "I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V "
                                           "got {line}")
                else:
                    # its numbers... will be added to anno dict below
                    values += [_StrToFloat(x) for x in line.split()]
            elif current_data_type == "M":
                if anno_type == AnnoType.SINGLE:
                    raise RuntimeError("Observed single AA and pairwise "
                                       "features in the same aaindex entry")
                anno_type = AnnoType.PAIR
                if line.startswith("M"):
                    # header, don't expect strict one letter code ordering here
                    # also because some pair entries include gaps ('-').
                    # expect something like: "M rows = <x>, cols = <x>"
                    split_line = line[1:].split(',')
                    split_line = sorted([item.strip() for item in split_line])
                    # after sorting we should have exactly two elements, the
                    # first starting with cols and the second with rows
                    if len(split_line) != 2 or \
                        not split_line[0].startswith("cols") or \
                        not split_line[1].startswith("rows"):
                        raise RuntimeError(f"Expect value header in pair "
                                           "AAIndex entry to be of form: "
                                           "\"M rows = <x>, cols = <x>\" got: "
                                           "{line}")
                    pair_cols = split_line[0].split("=")[1].strip()
                    pair_rows = split_line[1].split("=")[1].strip()
                    if len(pair_cols) != len(pair_cols):
                        raise RuntimeError(f"Expect rows and cols to have same "
                                           "number of elements when parsing "
                                           "pair AAIndex entry got {line}")
                else:
                    # its numbers... will be added to anno dict below
                    values += [_StrToFloat(x) for x in line.split()]

        if not stuff_read:
            return None

        if key == "":
            raise RuntimeError("Cannot parse AAIndex entry without key...")

        if anno_type == AnnoType.SINGLE:
            olcs = "ARNDCQEGHILKMFPSTWYV"
            if len(olcs) != len(values):
                raise RuntimeError(f"Expected {len(olcs)} values in single AA "
                                   "AAIndex entry, got {len(values)}")
            for olc, value in zip(olcs, values):
                anno[olc] = value
        elif anno_type == AnnoType.PAIR:
            # number of values can differ as the provided matrix can either be
            # diagonal (symmetric A -> B == B -> A) or rectangular (non-symmetric 
            # A -> B != B -> A)
            # For the former, number of columns and rows must be equal, no such
            # requirement for non-symmetric case
            n_values_match = False
            n_cols = len(pair_cols)
            n_rows = len(pair_rows)
            n_nonsym = n_cols * n_rows
            if len(values) == n_nonsym:
                n_values_match = True
                value_idx = 0
                for a in pair_rows:
                    for b in pair_cols:
                        anno[a+b] = values[value_idx]
                        value_idx += 1

            if n_cols == n_rows:
                n_values_match = True
                N = n_cols
                n_sym = (N*N - N) / 2 # number of elements below diagonal
                n_sym += N # add diagonal elements again
                if len(values) == n_sym:
                    value_idx = 0
                    for row_idx, row in enumerate(pair_rows):
                        for col in pair_cols[: row_idx+1]:
                            anno[row+col] = values[value_idx]
                            anno[col+row] = values[value_idx]
                            value_idx += 1
            if not n_values_match:
                raise RuntimeError(f"Number of parsed values doesn't match "
                                   "parsed rows and cols descriptors")
        else:
            raise RuntimeError("Cannot parse AAIndex entry without values...")

        data = AAIndexData()
        data.key = key
        data.title = title
        data.ref = ref
        data.authors = authors
        data.title = title
        data.journal = journal
        data.anno_type = anno_type
        data.anno = anno
        return data

    def GetScore(self, olc):
        """Score/Annotation getter

        :param olc: One letter code of amino acid
        :type olc: :class:`string`
        :returns: Annotation/score for *olc*
        :raises: :class:`ValueError` if *olc* is not known or 
                 :class:`RuntimeError` if anno_type of this 
                 :class:`AAIndexData` object is not AnnoType.SINGLE.
        """
        if self.anno_type == AnnoType.SINGLE:
            if olc in self.anno:
                return self.anno[olc]
            else:
                raise ValueError(f"OLC not in AAIndex: {olc}")
        raise RuntimeError("Cannot return score for single amino acid with "
                           "AAIndex of type PAIR")

    def GetPairScore(self, olc_one, olc_two):
        """Score/Annotation getter

        :param olc_one: One letter code of first amino acid
        :type olc_one: :class:`string`
        :param olc_two: One letter code of second amino acid
        :type olc_two: :class:`string`
        :returns: Pairwise annotation/score for *olc_one*/*olc_two*
        :raises: :class:`ValueError`  if key constructed from *olc_one* and 
                 *olc_two* is not known or 
                 :class:`RuntimeError` if anno_type of this 
                 :class:`AAIndexData` object is not AnnoType.PAIR.
        """
        if self.anno_type == AnnoType.PAIR:
            anno_key = olc_one + olc_two
            if anno_key in self.anno:
                return self.anno[anno_key]
            else:
                raise ValueError(f"Cannot find annotation for following pairs "
                                 "of olcs: {olc_one}, {olc_two}")
        raise RuntimeError("Cannot return score for pair of amino acid "
                           "with AAIndex of type SINGLE")


class AAIndex:
    """Provides access to data from the amino acid index database (aaindex):

    Kawashima, S. and Kanehisa, M.; AAindex: amino acid index database. 
    Nucleic Acids Res. 28, 374 (2000). 

    Files are available `here <https://www.genome.jp/aaindex/>`_

    :param aaindex_files: Paths to aaindex files. If not given, the files
                          aaindex1, aaindex2 and aaindex3 from the specified
                          source are used (Release 9.2, Feb 2017).
    :type aaindex_files: :class:`list` of :class:`str`
    """    
    def __init__(self, aaindex_files = None):
        if aaindex_files is None:
            aaindex_dir = os.path.join(GetSharedDataPath(), 'aaindex')
            self.files_to_load = [os.path.join(aaindex_dir, "aaindex1"),
                                  os.path.join(aaindex_dir, "aaindex2"),
                                  os.path.join(aaindex_dir, "aaindex3")]
        else:
            self.files_to_load = list(aaindex_files)
        self.aaindex_entries = dict()

    def keys(self):
        """Emulate dict like behvaiour and returns all available keys, accession
        numbers respectively. 

        :returns: keys (or accession numbers) of all available aaindex entries
        """
        self._LoadAll()
        return self.aaindex_entries.keys()

    def values(self):
        """Emulate dict like behvaiour and returns all available entries.

        :returns: iterable of entries (type :class:`AAIndexData`)
        """
        self._LoadAll()
        return self.aaindex_entries.values()

    def __getitem__(self, key):
        """Getter by aaindex accession number (e.g. ANDN920101)

        :param key: aaindex accession number
        :type key: :class:`str`
        :returns: :class:`AAIndexData` object
        """
        while True:
            if key in self.aaindex_entries:
                return self.aaindex_entries[key]
            # key is not available let's load the next aaindex file
            if not self._Load():
                # all files loaded, there is entry with *key* in any of them
                break
        raise KeyError(f"{key} does not exist in provided aa_index files")

    def _LoadAll(self):
        """Loads all remaining files specified in self.files_to_load
        """
        while self._Load():
            pass

    def _Load(self):
        """Loads and removes first element in self.files_to_load. Returns False
        if there is no file to load anymore, True if one is successfully loaded.
        """
        if len(self.files_to_load) == 0:
            return False
        path = self.files_to_load.pop(0)
        if not os.path.exists(path):
            raise RuntimeError(f"Tried to load {path} but file doesnt exist.")
        with open(path, 'r') as fh:
            data = fh.readlines()
            data_it = iter(data)
            while True: # read as long as it spits out AAIndexData entries
                entry = AAIndexData.Parse(data_it)
                if entry is None:
                    break # nothing to read anymore...
                else:
                    self.aaindex_entries[entry.key] = entry
        return True
