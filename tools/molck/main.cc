#include <unistd.h>
#include <boost/program_options.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <ost/base.hh>
#include <ost/boost_filesystem_helper.hh>
#include <ost/platform.hh>
#include <ost/conop/conop.hh>
#include <ost/io/mol/pdb_writer.hh>
#include <ost/io/mol/pdb_reader.hh>
#include <ost/io/mol/mmcif_reader.hh>
#include <ost/io/io_exception.hh>
#include <ost/mol/alg/molck.hh>
#include <ost/log.hh>
#if defined(__APPLE__)
#include <mach-o/dyld.h>
#endif
using namespace ost;
using namespace ost::mol;
using namespace ost::io;

namespace po=boost::program_options;
namespace fs=boost::filesystem;

const char* USAGE=
"This is molck - the molecule checker\n"
"Usage: molck [options] file1.pdb [file2.pdb [...]]\n"
"Options\n"
"  --complib=path       Location of the compound library file. If not provided,\n"
"                       the following locations are searched in this order:\n" 
"                       1. Working directory,\n"
"                       2. OpenStructure standard library location (if the\n"
"                          executable is part of a standard OpenStructure\n"
"                          installation)\n"
"  --rm=<a>,<b>         Remove atoms and residues matching some criteria:\n"
"                       - zeroocc - Remove atoms with zero occupancy\n"
"                       - hyd - Remove hydrogen atoms\n"
"                       - oxt - Remove terminal oxygens\n"
"                       - nonstd - Remove all residues not one of the\n"
"                                  20 standard amino acids\n"
"                       - unk - Remove unknown atoms not following\n"
"                               the nomenclature\n"
"                       Default: hyd\n"
"  --fix-ele            Clean up element column\n"
"  --stdout             Write cleaned file(s) to stdout\n"
"  --out=filename       Write cleaned file(s) to disk. % characters in the\n"
"                       filename are replaced with the basename of the input\n"
"                       file without extension. Default: %-molcked.pdb\n"
"  --color=auto|on|off  Whether output should be colored. Delault: auto\n"
"  --map-nonstd         Maps modified residues back to the parent amino acid,\n"
"                       for example: MSE -> MET, SEP -> SER.\n";

void usage()
{
  LOG_INFO(USAGE);
  exit(0);
}	

EntityHandle load_x(const String& file, const IOProfile& profile)
{
  try {
    EntityHandle ent = CreateEntity();
    const size_t fl = file.length();
    if (   (fl > 6 && file.compare(fl - 6, 6, ".mmcif") == 0)
        || (fl > 4 && file.compare(fl - 4, 4, ".cif") == 0)) {
      MMCifReader reader(file, ent, profile);
      reader.Parse();
      return ent;
    } else {
      // if its not mmcif, we assume the file to be in pdb format without even
      // looking at the filename
      PDBReader reader(file, profile);
      if (reader.HasNext()) {
        reader.Import(ent);
        return ent;
      }
      std::stringstream ss;
      ss << "ERROR: '" << file << "' does not contain any ATOM records. "
      << "Are you sure this is a PDB file?";
      LOG_INFO(ss.str());
      return EntityHandle();
    }
  } catch (std::exception& e) {
    LOG_ERROR("ERROR: " << e.what());
    return EntityHandle();
  }
}

// load compound library, exiting if it could not be found...
ost::conop::CompoundLibPtr load_compound_lib(const String& custom_path)
{
  if (custom_path!="") {
    if (fs::exists(custom_path)) {  
      return ost::conop::CompoundLib::Load(custom_path);
    } else {
      LOG_INFO("Could not find compounds.chemlib at the provided location, trying other options");
    }
  } 
  if (fs::exists("compounds.chemlib")) {
    return ost::conop::CompoundLib::Load("compounds.chemlib");
  }
  char result[ 1024 ]; 
  ost::conop::CompoundLibPtr lib;
  String exe_path; 
  #if defined(__APPLE__)
  uint32_t size=1023;
  if (!_NSGetExecutablePath(result, &size)) {
    exe_path=String(result); 
  }
  #else 
  ssize_t count = readlink( "/proc/self/exe", result, 1024 );
  exe_path = std::string( result, (count > 0) ? count : 0 );
  #endif
  if (exe_path.empty()) { 
    LOG_INFO("Could not determine the path of the molck executable. Will only "
             "look for compounds.chemlib in the current working directory");
  } else {
    fs::path path_and_exe(exe_path);
    fs::path path_only=path_and_exe.branch_path();
    fs::path share_path = path_only.branch_path();
    share_path = share_path / "share" / "openstructure" / "compounds.chemlib";

    String share_path_string=BFPathToString(share_path);
      
    if (fs::exists(share_path_string)) {
      return ost::conop::CompoundLib::Load(share_path_string);
    }  
  }
  if (!lib) {
    LOG_ERROR("Could not load compounds.chemlib");
    exit(-1);
  }
  return ost::conop::CompoundLibPtr();
}


int main(int argc, char *argv[])
{
  // setup logging
  ost::Logger::Instance().PushVerbosityLevel(ost::Logger::INFO);

  if (argc<2) {
    usage();
  }
  IOProfile prof;
  prof.fault_tolerant = true;
  // set all settings to false by default (set to true if args given)
  ost::mol::alg::MolckSettings settings(false, false, false, false, false,
                                        false, false, false);
  String rm;
  String color;

  bool write_to_stdout = false;
  bool write_to_file = false;
  String output_blueprint_string;
  String custom_path="";

  po::options_description desc("Options");
  desc.add_options()
    ("rm", po::value<String>(&rm)->default_value("hyd"), "atoms to be removed")
    ("color", po::value<String>(&color)->default_value("auto"),
     "whether the output should be colored.")
    ("files", po::value< std::vector<String> >(), "input file(s)")
    ("stdout", "write cleaned file(s) to stdout")
    ("out,o", po::value<String>(&output_blueprint_string)->default_value("%-molcked.pdb"),
     "write cleaned file to output using blueprint to determine path")
    ("map-nonstd", "map non standard residues back to standard ones (e.g.: MSE->MET,SEP->SER,etc.)")
    ("fix-ele", "insert element")
    ("complib", po::value<String>(&custom_path)->default_value(""),"location of the compound library file")
  ;
  po::positional_options_description p;
  p.add("files", -1);
  std::vector<String> files;
  po::variables_map vm;
  
  try {
    po::store(po::command_line_parser(argc, argv).
                options(desc).positional(p).run(),
              vm);
  } catch (std::exception& e) {
    LOG_ERROR(e.what());
    usage();
    exit(-1);
  }
  po::notify(vm); 
  if (vm.count("files")) {
    files = vm["files"].as<std::vector<String> >();
  } else {
    usage();
    exit(-1);
  }
  if (vm.count("complib")) {
    custom_path = vm["complib"].as<String>();
  }    
  if (vm.count("stdout")) {
    write_to_stdout = true;
  } else {
    write_to_file = true;
    output_blueprint_string = vm["out"].as<String>();
  }
  if (vm.count("map-nonstd")) {
    settings.map_nonstd_res = true;
  }
  if (vm.count("fix-ele")) {
    settings.assign_elem = true;
  }
  
  std::vector<StringRef> rms=StringRef(rm.c_str(), rm.size()).split(',');
  for (size_t i=0; i<rms.size(); ++i) {
    if (rms[i] == StringRef("unk", 3)) {
      settings.rm_unk_atoms = true;
    } else if (rms[i] == StringRef("nonstd", 6)) {
      settings.rm_non_std = true;
    } else if (rms[i] == StringRef("hyd", 3)) {
      settings.rm_hyd_atoms = true;
    } else if (rms[i] == StringRef("oxt", 3)) {
      settings.rm_oxt_atoms = true;
    } else if (rms[i] == StringRef("zeroocc", 7)) {
      settings.rm_zero_occ_atoms = true;
    } else {
      std::stringstream ss;
      ss << "unknown value to remove '" << rms[i] << "'";
      LOG_ERROR(ss.str());
      usage();
      exit(-1);
    }
  }
  if (color=="auto") {
    settings.colored = isatty(STDERR_FILENO);
  } else if (color == "on" || color == "1" || color == "yes") {
    settings.colored = true;
  } else if (color == "off" || color == "0" || color == "no") {
    settings.colored = false;	
  } else {
    usage();
    exit(-1);
  }
  ost::conop::CompoundLibPtr lib=load_compound_lib(custom_path);  
  for (unsigned int i = 0; i < files.size(); ++i) {
    EntityHandle ent=load_x(files[i], prof);
    if (!ent.IsValid()) {
      continue;
    }
    
    ost::mol::alg::Molck(ent, lib, settings);
 
    if (write_to_stdout) {
      PDBWriter writer(std::cout, prof);
      writer.Write(ent);
    }
    if (write_to_file) {
      fs::path input_file_path(files[i]);
      fs::path input_filename = BFPathStem(input_file_path);
      String input_filename_string=BFPathToString(input_filename);
      size_t replstart =output_blueprint_string.find('%');	
      String output_blueprint_string_copy = output_blueprint_string;
      if (replstart != String::npos) {
        output_blueprint_string_copy.replace(replstart,1,input_filename_string); 
      } 
      try {
        fs::path out_path(output_blueprint_string_copy);
        if (out_path.has_parent_path() && !exists(out_path.parent_path())) {
          LOG_ERROR("Output path does not exist: " + output_blueprint_string_copy);
          exit(-1);
        }
      } catch (std::exception& e) {

        size_t perden = String(e.what()).find("Permission denied");	

        if (perden != String::npos) {
          LOG_ERROR("Cannot write into output directory: " + output_blueprint_string_copy);
          exit(-1);
        } else {
          LOG_ERROR(e.what());
          exit(-1);
        }
      }
      LOG_INFO("Writing out file: " + output_blueprint_string_copy);
      PDBWriter writer(output_blueprint_string_copy, prof);
      writer.Write(ent);
    }
  }
 
  return 0;
}
