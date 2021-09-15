import subprocess
import ost
from ost import settings

def createdb(infasta, resultDB, exe_path=None):
    """
    Convert fasta files containing query and/or target sequences into a mmseqs2 database.
    
    :param infasta: The fasta file from which the mmseqs2 database will be created.
    :type infasta: :class:`string`

    :param resultDB: The output location for mmseqs2 database.
    :type resultDB: :class:`string`  

    :param exe_path: The path where mmseqs2 executable is located
    :type exe_path: :class:`string`

    """

    mmseqs2_exe = settings.Locate('mmseqs2', explicit_file_name=exe_path)
    args=[mmseqs2_exe, 'createdb', infasta, resultDB]

    ost.LogInfo(f"running MMseqs2 {' '.join(args)}")
    mmseqs2_pipe=subprocess.run(args)



def create_index(trg_db, exe_path=None, directory=None):
    """

    An index file of the targetDB is computed for a fast read-in.
    It is recommended to compute the index if the targetDB is reused for several searches.   
    A directory for temporary files is generated. 
    It is recommended to create this temporary folder on a local drive. 
    
    :param trg_db: The target database mmseqs2 file.
     (You need to initially convert a fasta file into a mmseqs2 database using createdb).
    :type trg_db: :class:`string`

    :param exe_path: The path where mmseqs2 executable is located.
    :type exe_path: :class:`string'

    :param directory: The directory for temperary files.
    :type directory: :class:`string`

    """

    mmseqs2_exe = settings.Locate('mmseqs2', explicit_file_name=exe_path)
    args=[mmseqs2_exe, 'createindex', trg_db, directory]

    ost.LogInfo(f"running MMseqs2 {' '.join(args)}")
    mmseqs2_pipe=subprocess.run(args)



def alignment(query_db, trg_db, resultDB, directory, resultDB_m8, sen=None, exe_path=None,
              start_sens=None, sens_steps=None, fmt=None):
    """

    The alignment consists of two steps the prefilter and alignment.

    :param query_db: The query database mmseqs2 file.
     (You need to initially convert a fasta file into a mmseqs2 database using createdb).
    :type query_db: :class:`string`
    
    :param trg_db: The target database mmseqs2 file.
     (You need to initially convert a fasta file into a mmseqs2 database using createdb).
    :type trg_db: :class:`string`

    :param resultDB: The output location.
     (Output of createdb)
    :type resultDB: :class:`string`  

    :param exe_path: The path where mmseqs2 executable is located.
    :type exe_path: :class:`string`

    :param directory: The directory for temperary files.
    :type directory: :class:`string`

    :param sen: It controls the speed and sensitivity of the search.
                A very fast search would use a sensitivity of 1.0,
                while a very sensitive search would use a sensitivity of up to 7.0.  
    :type sen: :class:`float`

    :param start_sens: Best hit fast. The lowest sensitivity is defined with --start-sens.                      
    :type start_sens: :class:`int`

    :param sens_steps: Best hit fast. 
           The number of steps to reach the highest sensitivity can be defined with --sens-steps.
    :type sens_steps: :class:`int`

    Convert the result database into a BLAST tab formatted file.
    The file is formatted as a tab-separated list with 12 columns: 
    (1,2) identifiers for query and target sequences/profiles, 
    (3) sequence identity, 
    (4) alignment length, 
    (5) number of mismatches, 
    (6) number of gap openings,
    (7-8, 9-10) domain start and end-position in query and in target, 
    (11) E-value, 
    and (12) bit score.

    The option --format-output defines a custom output format. 
    The fields that are supported can be found in the following link:
    https://github.com/soedinglab/mmseqs2/wiki#custom-alignment-format-with-convertalis

    :param resultDB_m8: The output location
    :type resultDB_m8: :class:`string`

    :param fmt: Format output type, if the default is not used.
    :type fmt: :class:`string`

    """
    mmseqs2_exe = settings.Locate('mmseqs2', explicit_file_name=exe_path)
    command=[mmseqs2_exe, 'search', query_db, trg_db, resultDB, directory, '-a']

    if sen:
        sen=str(sen)
        command.append('-s')
        command.append(sen)

    if start_sens and sens_steps:
        start_sens=str(start_sens)
        command.append('--start-sens')
        command.append(start_sens)
        sens_steps=str(sens_steps)
        command.append('--sens-steps')
        command.append(sens_steps)


    ost.LogInfo(f"running MMseqs2 {' '.join(command)}")

    mmseqs2_pipe=subprocess.run(command)


    args=[mmseqs2_exe, 'convertalis', query_db, trg_db, resultDB, resultDB_m8]
   
    if fmt:
        args.append('--format-output')
        args.append(fmt)


    ost.LogInfo(f"running MMseqs2 (' '.join(args))")

    mmseqs2_pipe=subprocess.run(args)
