import sys
import os
import subprocess
import tempfile
import shutil

from ost import io
from ost import mol

def _Setup(mdl, ref, mdl_ch1, mdl_ch2, ref_ch1, ref_ch2):
    """ Performs parameter checks and dumps files for DockQ

    In case of dimeric interfaces the respective chains are selected from
    mdl/trg, renamed to A/B and dumped to disk.

    In case of interfaces with more chains involved, we simply select the
    specified chains and do no renaming before dumping to disk.
    """
    if isinstance(mdl_ch1, str):
        mdl_ch1 = [mdl_ch1]
    if isinstance(mdl_ch2, str):
        mdl_ch2 = [mdl_ch2]
    if isinstance(ref_ch1, str):
        ref_ch1 = [ref_ch1]
    if isinstance(ref_ch2, str):
        ref_ch2 = [ref_ch2]

    if len(mdl_ch1) == 0:
        raise RuntimeError("mdl_ch1 is empty")
    if len(mdl_ch2) == 0:
        raise RuntimeError("mdl_ch2 is empty")

    if len(mdl_ch1) != len(ref_ch1):
        raise RuntimeError("mdl_ch1/ref_ch1 inconsistent in size")
    if len(mdl_ch2) != len(ref_ch2):
        raise RuntimeError("mdl_ch2/ref_ch2 inconsistent in size")

    for cname in mdl_ch1:
        ch = mdl.FindChain(cname)
        if not ch.IsValid():
            raise RuntimeError(f"Chain {cname} specified in mdl_ch1 not "
                               f"present in mdl")

    for cname in mdl_ch2:
        ch = mdl.FindChain(cname)
        if not ch.IsValid():
            raise RuntimeError(f"Chain {cname} specified in mdl_ch2 not "
                               f"present in mdl")

    for cname in ref_ch1:
        ch = ref.FindChain(cname)
        if not ch.IsValid():
            raise RuntimeError(f"Chain {cname} specified in ref_ch1 not "
                               f"present in ref")

    for cname in ref_ch2:
        ch = ref.FindChain(cname)
        if not ch.IsValid():
            raise RuntimeError(f"Chain {cname} specified in ref_ch2 not "
                               f"present in ref")

    mdl_to_dump = mdl.CreateFullView()
    ref_to_dump = ref.CreateFullView()

    if len(mdl_ch1) == 1 and len(mdl_ch2) == 1:
        # Dimer processing of mdl => Create new entity only containing 
        # the two specified chains and rename them to A, B
        mdl_to_dump = mol.CreateEntityFromView(mdl_to_dump, True)
        tmp = mol.CreateEntity()
        ed = tmp.EditXCS()
        ch1 = mdl_to_dump.FindChain(mdl_ch1[0])
        ed.InsertChain("A", ch1, deep=True)
        ch2 = mdl_to_dump.FindChain(mdl_ch2[0])
        ed.InsertChain("B", ch2, deep=True)
        mdl_ch1 = ["A"]
        mdl_ch2 = ["B"]
        mdl_to_dump = tmp

        # Same for ref
        ref_to_dump = mol.CreateEntityFromView(ref_to_dump, True)
        tmp = mol.CreateEntity()
        ed = tmp.EditXCS()
        ch1 = ref_to_dump.FindChain(ref_ch1[0])
        ed.InsertChain("A", ch1, deep=True)
        ch2 = ref_to_dump.FindChain(ref_ch2[0])
        ed.InsertChain("B", ch2, deep=True)
        ref_ch1 = ["A"]
        ref_ch2 = ["B"]
        ref_to_dump = tmp
    else:
        # Interface with more chains...
        raise NotImplementedError("DockQ computations beyond two interacting "
                                  "chains has not been properly tested...")
        #mdl_chain_names = mdl_ch1 + mdl_ch2
        #ref_chain_names = ref_ch1 + ref_ch2
        #mdl_to_dump = mdl_to_dump.Select(f"cname={','.join(mdl_chain_names)}")
        #ref_to_dump = ref_to_dump.Select(f"cname={','.join(ref_chain_names)}")

    # first write structures to string, only create a tmpdir and the actual
    # files if this succeeds
    mdl_str = io.EntityToPDBStr(mdl_to_dump)
    ref_str = io.EntityToPDBStr(ref_to_dump)

    tmp_dir = tempfile.mkdtemp()
    with open(os.path.join(tmp_dir, "mdl.pdb"), 'w') as fh:
        fh.write(mdl_str)
    with open(os.path.join(tmp_dir, "ref.pdb"), 'w') as fh:
        fh.write(ref_str)

    return (tmp_dir, mdl_ch1, mdl_ch2, ref_ch1, ref_ch2)

class DockQResult:
    """ DockQ result object
    """
    def __init__(self, Fnat, Fnonnat, native_contacts, model_contacts, iRMS,
                 LRMS, DockQ):
        self._Fnat = Fnat
        self._Fnonnat = Fnonnat
        self._native_contacts = native_contacts
        self._model_contacts = model_contacts
        self._iRMS = iRMS
        self._LRMS = LRMS
        self._DockQ = DockQ

    @property
    def Fnat(self):
        """ DockQ - Fnat output

        :type: :class:`float`
        """
        return self._Fnat
    
    @property
    def Fnonnat(self):
        """ DockQ - Fnonnat output

        :type: :class:`float`
        """
        return self._Fnonnat

    @property
    def native_contacts(self):
        """ DockQ - number native contacts

        :type: :class:`int`
        """
        return self._native_contacts

    @property
    def model_contacts(self):
        """ DockQ - number model contacts

        :type: :class:`int`
        """
        return self._model_contacts

    @property
    def iRMS(self):
        """ DockQ - iRMS output

        :type: :class:`float`
        """
        return self._iRMS

    @property
    def LRMS(self):
        """ DockQ - LMRS output 

        :type: :class:`float`
        """
        return self._LRMS

    @property
    def DockQ(self):
        """ DockQ - DockQ output

        :type: :class:`float`
        """
        return self._DockQ

    def JSONSummary(self):
        """ Returns JSON serializable summary
        """
        return {"Fnat": self.Fnat,
                "Fnonnat": self.Fnonnat,
                "native_contacts": self.native_contacts,
                "model_contacts": self.model_contacts,
                "iRMS": self.iRMS,
                "LRMS": self.LRMS,
                "DockQ": self.DockQ}

    @staticmethod
    def FromDockQOutput(output):
        """ Static constructor from raw DockQ output

        :param output: Raw output from DockQ executable
        :type output: :class:`str`
        :returns: Object of type :class:`DockQResult`
        """
        Fnat = None
        Fnonnat = None
        native_contacts = None
        model_contacts = None
        iRMS = None
        LRMS = None
        DockQ = None

        for line in output.splitlines():
            if line.startswith('*'):
                continue
            if line.startswith("Fnat"):
                Fnat = float(line.split()[1])
                native_contacts = int(line.split()[5])
            elif line.startswith("Fnonnat"):
                Fnonnat = float(line.split()[1])
                model_contacts = int(line.split()[5])
            elif line.startswith("iRMS"):
                iRMS = float(line.split()[1])
            elif line.startswith("LRMS"):
                LRMS = float(line.split()[1])
            elif line.startswith("DockQ"):
                DockQ = float(line.split()[1])

        return DockQResult(Fnat, Fnonnat, native_contacts, model_contacts,
                           iRMS, LRMS, DockQ)


def DockQ(dockq_exec, mdl, ref, mdl_ch1, mdl_ch2, ref_ch1,
          ref_ch2):
    """ Computes DockQ for specified interface 

    DockQ is available from https://github.com/bjornwallner/DockQ - 
    For this binding to work, DockQ must be properly installed and its
    dependencies must be available (numpy, Biopython).

    :param dockq_exec: Path to DockQ.py script from DockQ repository
    :type dockq_exec: :class:`str`
    :param mdl: Model structure
    :type mdl: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
    :param ref: Reference structure, i.e. native structure
    :type ref: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
    :param mdl_ch1: Specifies chain(s) in model constituting first part of
                    interface
    :type mdl_ch1: :class:`str`/:class:`list` of :class:`str`
    :param mdl_ch2: Specifies chain(s) in model constituting second part of
                    interface
    :type mdl_ch2: :class:`str`/:class:`list` of :class:`str`
    :param ref_ch1: ref equivalent of mdl_ch1
    :type ref_ch1: :class:`str`/:class:`list` of :class:`str`
    :param ref_ch2: ref equivalent of mdl_ch2
    :type ref_ch2: :class:`str`/:class:`list` of :class:`str`
    :returns: Result object of type :class:`DockQResult`
    """
    if not os.path.exists(dockq_exec):
        raise RuntimeError(f"DockQ executable ({dockq_exec}) does not exist")

    tmp_dir, mdl_ch1, mdl_ch2, ref_ch1, ref_ch2 = \
    _Setup(mdl, ref, mdl_ch1, mdl_ch2, ref_ch1, ref_ch2)

    cmd = [sys.executable, dockq_exec, os.path.join(tmp_dir, "mdl.pdb"),
           os.path.join(tmp_dir, "ref.pdb")]

    # add mdl/ref chains
    cmd.append("-model_chain1")
    cmd += mdl_ch1
    cmd.append("-model_chain2")
    cmd += mdl_ch2
    cmd.append("-native_chain1")
    cmd += ref_ch1
    cmd.append("-native_chain2")
    cmd += ref_ch2

    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    shutil.rmtree(tmp_dir) # cleanup, no matter if DockQ executed successfully 

    if proc.returncode != 0:
        raise RuntimeError("DockQ run failed - returncode: " + \
                           str(proc.return_code))

    if proc.stderr.decode() != "":
        raise RuntimeError("DockQ run failed - stderr: " + proc.stderr.decode())

    return DockQResult.FromDockQOutput(proc.stdout.decode())
