''' Unit tests for the HHblits wrapper
'''

import unittest
import os
import sys
import tempfile
import shutil
import filecmp
import datetime

import ost
from ost import seq
from ost import settings
from ost.bindings import hhblits3

class _UnitTestHHblitsLog(ost.LogSink):
  """Dedicated logger to hide some expected warning/ error messages.
  """
  def __init__(self):
    ost.LogSink.__init__(self)
    self.lcwd = os.getcwd()

  def LogMessage(self, message, severity):
    message = message.strip()
    dnem = "could not open file '%s'" % os.path.join(self.lcwd,
                                                     'doesnotexist.a3m')
    if message.endswith(dnem):
        return
    print(message)

def setUpModule():
  ost.PushLogSink(_UnitTestHHblitsLog())

class TestHHblitsBindings(unittest.TestCase):
    def setUp(self):
        # expects default hhblits installation...
        hhblits_bin_dir = os.path.dirname(settings.Locate('hhblits'))
        self.hhroot = os.path.join(hhblits_bin_dir, os.pardir)

    def tearDown(self):
        if hasattr(self, 'hh'):
            self.hh.Cleanup()
        if hasattr(self, 'tmpfile'):
            os.remove(self.tmpfile)

    def testOSTSequenceNoWDir(self):
        # OST sequence as input, no working dir should work
        query_seq = seq.CreateSequence('1AKE.B',
                                       'MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRA'+
                                       'AVKSGSELGKQAKDIMDAGKLVTDELVIALVKERIAQ'+
                                       'EDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDV'+
                                       'PDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDV'+
                                       'TGEELTTRKDDQEETVRKRLVEYHQMTAPLIGYYYYS'+
                                       'KEAEAGNTKYAKVDGTKPVAEVRADLEKILG')
        self.hh = hhblits3.HHblits(query_seq, self.hhroot)
        # this works only as long as utils.TempDirWithFiles() names the first
        # fasta file seq01.fasta
        self.assertEqual(self.hh.filename[-11:], 'seq01.fasta')

    def testOSTSequenceWDir(self):
        # OST sequence as input, with working dir should work
        query_seq = seq.CreateSequence('1AKE.B',
                                       'MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRA'+
                                       'AVKSGSELGKQAKDIMDAGKLVTDELVIALVKERIAQ'+
                                       'EDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDV'+
                                       'PDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDV'+
                                       'TGEELTTRKDDQEETVRKRLVEYHQMTAPLIGYYYYS'+
                                       'KEAEAGNTKYAKVDGTKPVAEVRADLEKILG')
        tmpdir = tempfile.mkdtemp()
        self.hh = hhblits3.HHblits(query_seq, self.hhroot, working_dir=tmpdir)
        self.assertEqual(self.hh.filename, os.path.join(tmpdir,
                                                        'query_hhblits.fasta'))
        self.assertEqual(self.hh.working_dir, tmpdir)
        self.hh.needs_cleanup = True

    def testFileNoWDir(self):
        # we have a sequence but no working dir
        self.hh = hhblits3.HHblits('testfiles/seq2.fasta', self.hhroot)
        self.assertEqual(self.hh.filename[-10:], 'seq2.fasta')

    def testFileWDir(self):
        # we have a sequence and a working dir
        tmpdir = tempfile.mkdtemp()
        self.hh = hhblits3.HHblits('testfiles/seq2.fasta', self.hhroot,
                                   working_dir=tmpdir)
        self.assertEqual(self.hh.filename, os.path.join(tmpdir, 'seq2.fasta'))
        self.assertEqual(self.hh.working_dir, tmpdir)
        self.hh.needs_cleanup = True

    def testDoCleanup(self):
        # test that cleanup cleans up
        tmpdir = tempfile.mkdtemp()
        self.hh = hhblits3.HHblits('testfiles/seq2.fasta', self.hhroot,
                                   working_dir=tmpdir)
        self.assertTrue(os.path.isdir(tmpdir))
        self.hh.needs_cleanup = True
        self.hh.Cleanup()
        self.assertFalse(os.path.isdir(tmpdir))

    def testDoNotCleanup(self):
        # test that cleanup cleans up
        tmpdir = tempfile.mkdtemp()
        self.hh = hhblits3.HHblits('testfiles/seq2.fasta', self.hhroot,
                                   working_dir=tmpdir)
        self.assertTrue(os.path.isdir(tmpdir))
        self.hh.needs_cleanup = False
        self.hh.Cleanup()
        self.assertTrue(os.path.isdir(tmpdir))
        shutil.rmtree(tmpdir)
        self.assertFalse(os.path.isdir(tmpdir))

    def testCleanupFailed(self):
        # test that cleanup cleans up
        tmpdir = tempfile.mkdtemp()
        self.hh = hhblits3.HHblits('testfiles/seq2.fasta', self.hhroot,
                                   working_dir=tmpdir)
        self.assertTrue(os.path.isdir(tmpdir))
        self.hh.needs_cleanup = False
        self.hh.CleanupFailed()
        self.assertFalse(os.path.isdir(tmpdir))

    def testFileNoExists(self):
        # test on calling with a non-existing fasta file
        with self.assertRaises(IOError) as ioe:
            hhblits3.HHblits('notthere.fasta', self.hhroot)
        self.assertEqual(ioe.exception.errno, 2)
        self.assertEqual(ioe.exception.strerror, 'No such file or directory')

    def testProfile(self):
        # try to create a search profile
        query_seq = seq.CreateSequence('Test', 'VLSPADKTNVKAAWGKVGAHAGEYGAEA'+
                                       'LERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVAD'+
                                       'ALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKL'+
                                       'LSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVL'+
                                       'TSKYR')
        self.hh = hhblits3.HHblits(query_seq, self.hhroot)
        a3m = self.hh.BuildQueryMSA('testfiles/hhblitsdb/unittestdb', 
                                    assign_ss = False)
        self.assertTrue(filecmp.cmp(a3m, "testfiles/testali_two_no_ss.a3m"))

    def testA3mToProfileFileName(self):
        # test A3mToProfile to work with a given hhmake_file name
        query_seq = seq.CreateSequence('Test', 'VLSPADKTNVKAAWGKVGAHAGEYGAEA'+
                                       'LERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVAD'+
                                       'ALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKL'+
                                       'LSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVL'+
                                       'TSKYR')
        self.hh = hhblits3.HHblits(query_seq, self.hhroot)
        _, self.tmpfile = tempfile.mkstemp(suffix='.hmm')
        os.remove(self.tmpfile)
        hhfile = self.hh.A3MToProfile("testfiles/testali.a3m",
                                      hhm_file=self.tmpfile)
        with open(hhfile) as tfh:
          tlst = tfh.readlines()
        with open("testfiles/test.hmm") as efh:
          elst = efh.readlines()
        self.assertEqual(len(elst), len(tlst))
        for i in range(0, len(elst)):
            if not elst[i].startswith(('FILE', 'COM', 'DATE')):
                self.assertEqual(elst[i], tlst[i])

    def testA3mToProfileWithoutFileName(self):
        # test A3mToProfile to work without a given hhmake_file name
        query_seq = seq.CreateSequence('Test', 'VLSPADKTNVKAAWGKVGAHAGEYGAEA'+
                                       'LERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVAD'+
                                       'ALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKL'+
                                       'LSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVL'+
                                       'TSKYR')
        self.hh = hhblits3.HHblits(query_seq, self.hhroot)
        hhfile = self.hh.A3MToProfile("testfiles/testali.a3m")
        with open(hhfile) as tfh:
          tlst = tfh.readlines()
        with open("testfiles/test.hmm") as efh:
          elst = efh.readlines()
        self.assertEqual(len(elst), len(tlst))
        for i in range(0, len(elst)):
            if not elst[i].startswith(('FILE', 'COM', 'DATE')):
                self.assertEqual(elst[i], tlst[i])
        os.remove(hhfile)

    def testA3mToProfileWithExistingFile(self):
        # test A3mToProfile to work with an existing hmm file
        query_seq = seq.CreateSequence('Test', 'VLSPADKTNVKAAWGKVGAHAGEYGAEA'+
                                       'LERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVAD'+
                                       'ALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKL'+
                                       'LSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVL'+
                                       'TSKYR')
        self.hh = hhblits3.HHblits(query_seq, self.hhroot)
        hhfile = self.hh.A3MToProfile("testfiles/testali.a3m",
                                      hhm_file="testfiles/test.hmm")
        # when the hmm file already exists, its not touched, so files should be
        # completely equal.
        self.assertTrue(filecmp.cmp(hhfile, "testfiles/test.hmm"))

    def testA3mToProfileWithoutA3m(self):
        # test A3mToProfile to fail if we do not have an a3m file.
        query_seq = seq.CreateSequence('Test', 'VLSPADKTNVKAAWGKVGAHAGEYGAEA'+
                                       'LERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVAD'+
                                       'ALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKL'+
                                       'LSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVL'+
                                       'TSKYR')
        self.hh = hhblits3.HHblits(query_seq, self.hhroot)
        with self.assertRaises(IOError) as ioe:
            self.hh.A3MToProfile("doesnotexist.a3m")
        self.assertEqual(ioe.exception.errno, None)
        self.assertEqual(ioe.exception.args[0],
                         "could not convert a3m to hhm file")

    def testA3mToCSFileName(self):
        # test A3mToCS to work with a given hhmake_file name
        query_seq = seq.CreateSequence('Test', 'VLSPADKTNVKAAWGKVGAHAGEYGAEA'+
                                       'LERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVAD'+
                                       'ALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKL'+
                                       'LSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVL'+
                                       'TSKYR')
        self.hh = hhblits3.HHblits(query_seq, self.hhroot)
        _, self.tmpfile = tempfile.mkstemp(suffix='.seq219')
        os.remove(self.tmpfile)
        csfile = self.hh.A3MToCS("testfiles/testali.a3m",
                                 cs_file=self.tmpfile, options={'-alphabet' :
                                                os.path.join(self.hhroot,
                                                             'data',
                                                             'cs219.lib')})
        self.assertTrue(filecmp.cmp(csfile, 'testfiles/test.seq219'))

    def testA3mToCSNoFileName(self):
        # test A3mToCS to work with a given hhmake_file name
        query_seq = seq.CreateSequence('Test', 'VLSPADKTNVKAAWGKVGAHAGEYGAEA'+
                                       'LERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVAD'+
                                       'ALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKL'+
                                       'LSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVL'+
                                       'TSKYR')
        self.hh = hhblits3.HHblits(query_seq, self.hhroot)

        csfile = self.hh.A3MToCS("testfiles/testali.a3m",
                                 options={'-alphabet' :
                                          os.path.join(self.hhroot,
                                                       'data',
                                                       'cs219.lib')})
        self.assertTrue(filecmp.cmp(csfile, 'testfiles/test.seq219'))
        self.tmpfile = "testfiles/testali.seq219"

    def testA3mToCSFileExists(self):
        # test A3mToCS to work with a given hhmake_file name
        query_seq = seq.CreateSequence('Test', 'VLSPADKTNVKAAWGKVGAHAGEYGAEA'+
                                       'LERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVAD'+
                                       'ALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKL'+
                                       'LSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVL'+
                                       'TSKYR')
        self.hh = hhblits3.HHblits(query_seq, self.hhroot)
        csfile = self.hh.A3MToCS("testfiles/testali.a3m",
                                 cs_file='testfiles/test.seq219',
                                 options={'-alphabet' :
                                          os.path.join(self.hhroot,
                                                       'data',
                                                       'cs219.lib')})
        self.assertTrue(filecmp.cmp(csfile, 'testfiles/test.seq219'))

    def testSearchWorking(self):
        # successful search
        query_seq = seq.CreateSequence('Test', 'VLSPADKTNVKAAWGKVGAHAGEYGAEA'+
                                       'LERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVAD'+
                                       'ALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKL'+
                                       'LSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVL'+
                                       'TSKYR')
        self.hh = hhblits3.HHblits(query_seq, self.hhroot)
        search_file = self.hh.Search("testfiles/testali.a3m",
                                     'testfiles/hhblitsdb/unittestdb')

        with open(search_file) as f:
              header, hits = hhblits3.ParseHHblitsOutput(f)
        self.assertEqual(len(hits), 4)
        exp_scores = [199.3, 12.3, 11.5, 8.1]
        for hit_idx in range(4):
            self.assertAlmostEqual(hits[hit_idx].score, exp_scores[hit_idx])

    def testSearchNotWorking(self):
        query_seq = seq.CreateSequence('Test', 'VLSPADKTNVKAAWGKVGAHAGEYGAEA'+
                                       'LERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVAD'+
                                       'ALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKL'+
                                       'LSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVL'+
                                       'TSKYR')
        self.hh = hhblits3.HHblits(query_seq, self.hhroot)
        with self.assertRaises(RuntimeError) as rte:
            self.hh.Search("doesnotexist.a3m", 'testfiles/hhblitsdb/unittestdb')

    def testParseHHMWorking(self):
        # get info from an HHM file
        with open("testfiles/test.hmm") as hhm_fh:
            prof = hhblits3.ParseHHM(hhm_fh)
        self.assertEqual(''.join([str(x) for x in prof['ss_conf']]),
                         '999999999999998873391557999999998639441123987788888'+
                         '856788999999999998735477789999999887650299989899889'+
                         '999999997536679989999999999999999984329')
        self.assertEqual(''.join(prof['ss_pred']), 'CCCHHHHHHHHHHHHHHHCCCHHHH'+
                         'HHHHHHHHHHCCCCCCCCCCCCCCCCCHHHHHHHHHHHHHHHHHHHCCCCH'+
                         'HHHHHHHHHHHHHHCCCCHHHHHHHHHHHHHHHHHHCCCCCCHHHHHHHHH'+
                         'HHHHHHHHHHHHCC')
        self.assertEqual(str(prof['consensus']), 'xltxxxxxxixxsWxxvxxxxxxxgxx'+
                         'xxxxlfxxxPxxxxxFxxxxxxxxxxxxxxhxxxvxxxlxxxixxldxxxx'+
                         'xlxxlxxxHxxxxgvxxxxxxxxxxxlxxxlxxxxgxxxxxxxxxAWxxxx'+
                         'xxixxxmxxxyx')
        self.assertEqual(prof['msa'].GetCount(), 7)

    def testParseHHMNotWorking(self):
        # get info from an HHM file
        with self.assertRaises(IOError) as ioe:
            with open('testfiles/testali.a3m') as f:
              hhblits3.ParseHHM(f)
        self.assertEqual(ioe.exception.args[0],
                         'Profile file "testfiles/testali.a3m" is missing '+
                         'the "Consensus" section')

    def testParseA3MWorking(self):
        # get info from an HHM file
        with open("testfiles/testali.a3m") as a3m_fh:
            prof = hhblits3.ParseA3M(a3m_fh)
        self.assertEqual(''.join([str(x) for x in prof['ss_conf']]),
                         '999999999999998873391557999999998639441123987788888'+
                         '856788999999999998735477789999999887650299989899889'+
                         '999999997536679989999999999999999984329')
        self.assertEqual(''.join(prof['ss_pred']), 'CCCHHHHHHHHHHHHHHHCCCHHHH'+
                         'HHHHHHHHHHCCCCCCCCCCCCCCCCCHHHHHHHHHHHHHHHHHHHCCCCH'+
                         'HHHHHHHHHHHHHHCCCCHHHHHHHHHHHHHHHHHHCCCCCCHHHHHHHHH'+
                         'HHHHHHHHHHHHCC')
        self.assertEqual(prof['msa'].GetCount(), 253)

    def fastParseHeader(self):
        header_line = '  1 814cbc1899f35c872169524af30fc2 100.0  5E-100' + \
                      '  5E-104  710.5  34.1  277    3-293     2-280 (281)'
        hit, offset = hhblits3.ParseHeaderLine(header_line)
        self.assertEqual(hit.hit_id, '814cbc1899f35c872169524af30fc2')
        self.assertAlmostEqual(hit.evalue, 0)
        self.assertAlmostEqual(hit.prob, 100.0)
        self.assertAlmostEqual(hit.pvalue, 0)
        self.assertAlmostEqual(hit.score, 710.5)
        self.assertAlmostEqual(hit.ss_score, 34.1)

    def testParseHHblitsOutput(self):
        with open("testfiles/test.hhr") as f:
          header, hits = hhblits3.ParseHHblitsOutput(f)
        self.assertEqual(header.query, 'Test')
        self.assertEqual(header.match_columns, 141)
        self.assertEqual(header.n_eff, 9.4)
        self.assertEqual(header.searched_hmms, 5)
        self.assertEqual(header.date,
                         datetime.datetime.strptime('Mon Aug 31 16:45:45 2015',
                                                    '%a %b %d %H:%M:%S %Y'))
        self.assertEqual(header.command, '/import/bc2/apps/HH-suite/2.0.16-go'+
                         'olf-1.4.10/bin/hhblits -cpu 1 -n 1 -e 0.001 -Z 1000'+
                         '0 -B 10000 -i /import/bc2/home/schwede/bienert/git/'+
                         'ost_newenv.git/modules/bindings/tests/testfiles/tes'+
                         'tali.a3m -o /scratch/tmp9E5dor/testali_cpu1_n1.hhr '+
                         '-d /import/bc2/home/schwede/bienert/git/ost_newenv.'+
                         'git/modules/bindings/tests/testfiles/hhblitsdb/hhbl'+
                         'itsdb')
        self.assertEqual(len(hits), 4)
        # hit 1
        self.assertEqual(hits[0].hit_id, '3e7b90809bd446a538f9eb1a1ca0e551')
        self.assertEqual(hits[0].score, float('222.3'))
        self.assertEqual(hits[0].ss_score, float('16.4'))
        self.assertEqual(hits[0].evalue, float('2.7e-42'))
        self.assertEqual(hits[0].pvalue, float('6.9E-45'))
        self.assertEqual(hits[0].prob, float('100.0'))
        self.assertEqual(str(hits[0].aln),
                         'Test                VLSPADKTNVKAAWGKVGAHAGEYGAEALER'+
                         'MFLSFPTTKTYFPHF-DL-S----HGSAQ\n3e7b90809bd446a5... '+
                         'VLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKT'+
                         'EAEMKASED\n\nTest                VKGHGKKVADALTNAVAH'+
                         'VDDMPNALSALSDLHAHK-LRVDPVNFKLLSHCLLVTLAAHL\n3e7b908'+
                         '09bd446a5... LKKHGVTVLTALGAILKKKGHHEAELKPLAQSHA-TKH'+
                         'KIPIKYLEFISEAIIHVLHSRH\n\nTest                PAEFT'+
                         'PAVHASLDKFLASVSTVLTSKYR\n3e7b90809bd446a5... PGDFGA'+
                         'DAQGAMNKALELFRKDIAAKYK\n')
        # hit 2
        self.assertEqual(hits[1].hit_id, 'af828e69a5f2d0fd42a2a213e09eb64f')
        self.assertEqual(hits[1].score, float('215.7'))
        self.assertEqual(hits[1].ss_score, float('17.2'))
        self.assertEqual(hits[1].evalue, float('1.5e-41'))
        self.assertEqual(hits[1].pvalue, float('3.9E-44'))
        self.assertEqual(hits[1].prob, float('100.0'))
        self.assertEqual(str(hits[1].aln),
                         'Test                VLSPADKTNVKAAWGKVGAHAGEYGAEALER'+
                         'MFLSFPTTKTYFPHFDLSHGSAQVKGHGK\naf828e69a5f2d0fd... '+
                         'VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHG'+
                         'SAQVKGHGK\n\nTest                KVADALTNAVAHVDDMPN'+
                         'ALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPA\naf828e6'+
                         '9a5f2d0fd... KVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNF'+
                         'KLLSHCLLVTLAAHLPAEFTPA\n\nTest                VHASL'+
                         'DKFLASVSTVLTSKYR\naf828e69a5f2d0fd... VHASLDKFLASVS'+
                         'TVLTSKYR\n')
        # hit 3
        self.assertEqual(hits[2].hit_id, '9287755aa6aa27583da6be3b2408bfcc')
        self.assertEqual(hits[2].score, float('215.1'))
        self.assertEqual(hits[2].ss_score, float('16.5'))
        self.assertEqual(hits[2].evalue, float('3.8e-41'))
        self.assertEqual(hits[2].pvalue, float('9.7E-44'))
        self.assertEqual(hits[2].prob, float('100.0'))
        self.assertEqual(str(hits[2].aln),
                         'Test                VLSPADKTNVKAAWGKVGAHAGEYGAEALER'+
                         'MFLSFPTTKTYFPHF-DLS-----HGSAQ\n9287755aa6aa2758... '+
                         'HLTPEEKSAVTALWGKVN--VDEVGGEALGRLLVVYPWTQRFFESFGDLST'+
                         'PDAVMGNPK\n\nTest                VKGHGKKVADALTNAVAH'+
                         'VDDMPNALSALSDLHAHK-LRVDPVNFKLLSHCLLVTLAAHL\n9287755'+
                         'aa6aa2758... VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHC-DKL'+
                         'HVDPENFRLLGNVLVCVLAHHF\n\nTest                PAEFT'+
                         'PAVHASLDKFLASVSTVLTSKYR\n9287755aa6aa2758... GKEFTP'+
                         'PVQAAYQKVVAGVANALAHKYH\n')
        # hit 4
        self.assertEqual(hits[3].hit_id, 'e69e1ac0a4b2554df25f2c2183b0fba0')
        self.assertEqual(hits[3].score, float('12.6'))
        self.assertEqual(hits[3].ss_score, float('2.6'))
        self.assertEqual(hits[3].evalue, float('64.0'))
        self.assertEqual(hits[3].pvalue, float('0.16'))
        self.assertEqual(hits[3].prob, float('0.6'))
        self.assertEqual(str(hits[3].aln),
                         'Test                VDPVNFKLLSHCLLVTLAAHL\ne69e1ac0'+
                         'a4b2554d... ATPEQAQLVHKEIRKIVKDTC\n')

if __name__ == "__main__":
    try:
        from ost.bindings import hhblits
        hhblits_version_string = hhblits.GetHHblitsVersionString()
    except:
        print("No HHblits installation found: skip unit tests")
        sys.exit(0)
    if hhblits_version_string is None:
        print("Could not identify HHblits version: skip unit tests")
        sys.exit(0)
    hhblits_version_major = int(hhblits_version_string.split('.')[0])
    if hhblits_version_major != 3:
        print("Found hhblits with major version %i. Unit tests only run with "\
              "version 3: skip unit tests"%(hhblits_version_major))
        sys.exit(0)

    from ost import testutils
    testutils.RunTests()
