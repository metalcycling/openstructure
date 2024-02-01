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
// along with this library; if not, write to the  Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//------------------------------------------------------------------------------
#include <iomanip>
#if defined (_MSC_VER)
  #define BOOST_ALL_DYN_LINK 1
  #include <BaseTsd.h>
  typedef SSIZE_T ssize_t;
#endif
#include <boost/program_options.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/convenience.hpp>
#include <ost/base.hh>
#include <ost/boost_filesystem_helper.hh>
#include <ost/mol/chain_view.hh>
#include <ost/mol/alg/local_dist_diff_test.hh>
#include <ost/mol/alg/filter_clashes.hh>
#include <ost/io/mol/pdb_reader.hh>
#include <ost/io/io_exception.hh>
#include <ost/io/stereochemical_params_reader.hh>
#include <ost/conop/conop.hh>
#include <ost/conop/compound_lib.hh>
#include <ost/string_ref.hh>
#include <ost/conop/amino_acids.hh>
#include <ost/conop/rule_based.hh>
#include <ost/platform.hh>
#include <ost/log.hh>
#include <ost/mol/alg/consistency_checks.hh>
#include <ost/version.hh>

#if defined(__APPLE__)
#include <mach-o/dyld.h>
#endif
#include <ost/dyn_cast.hh>

using namespace ost;
using namespace ost::io;
using namespace ost::conop;
using namespace ost::mol;
using namespace ost::mol::alg;

namespace po=boost::program_options;
namespace fs=boost::filesystem;

// loads a file
EntityHandle load(const String& file, const IOProfile& profile)
{
  try {
    PDBReader reader(file, profile);
    if (reader.HasNext()) {
      EntityHandle ent=CreateEntity();
      reader.Import(ent);
      if (profile.processor) {
        profile.processor->Process(ent);
      }

      if (ent.GetChainList().size()!=1) {
        std::cout << "WARNING: File " << file << " has more than one chain" << std::endl;
      }    
      return ent;
    }
    std::cerr << "ERROR: '" << file << "' does not contain any ATOM records. "
              << "Are you sure this is a PDB file?" << std::endl;
    return EntityHandle();
  } catch (io::IOException& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return EntityHandle();
  }
}

// prints usage output
void usage()
{
  std::cerr << "usage: lddt [options] <mod1> [mod1 [mod2]] <re1>[,ref2,ref3]" << std::endl;
  std::cerr << "   -s <sel>   selection performed on ref" << std::endl;
  std::cerr << "   -c         use Calphas only" << std::endl;
  std::cerr << "   -f         perform structural checks and filter input data" << std::endl;
  std::cerr << "   -t         fault tolerant parsing" << std::endl;
  std::cerr << "   -p <file>  use specified parmeter file. Mandatory" << std::endl;
  std::cerr << "   -v <level> verbosity level (0=results only,1=problems reported, 2=full report)" << std::endl;
  std::cerr << "   -b <value> tolerance in stddevs for bonds" << std::endl;
  std::cerr << "   -a <value> tolerance in stddevs for angles" << std::endl;
  std::cerr << "   -r <value> distance inclusion radius" << std::endl;
  std::cerr << "   -i <value> sequence separation" << std::endl;
  std::cerr << "   -e         print version" << std::endl;
  std::cerr << "   -x         ignore residue name consistency checks" << std::endl;
  std::cerr << "   -l <path>  location of the compound library file" << std::endl;
}

CompoundLibPtr load_compound_lib(const String& custom_path)
{
  if (custom_path!="") {
    if (fs::exists(custom_path)) {  
      return CompoundLib::Load(custom_path);
    } else {
      std::cerr << "Could not find compounds.chemlib at the provided location, trying other options" << std::endl;
    }
  } 
  if (fs::exists("compounds.chemlib")) {
    return CompoundLib::Load("compounds.chemlib");
  }
  char result[ 1024 ]; 
  CompoundLibPtr lib;
  String exe_path; 
  #if defined(__APPLE__)
  uint32_t size=1023;
  if (!_NSGetExecutablePath(result, &size)) {
    exe_path=String(result); 
  }
  #else
    #if defined (_MSC_VER)
      // todo find exe path on Windows
    #else
      ssize_t count = readlink( "/proc/self/exe", result, 1024 );
      exe_path = std::string( result, (count > 0) ? count : 0 );
    #endif
  #endif
  if (exe_path.empty()) { 
    std::cerr << "Could not determine the path of the lddt executable. Will only "
       "look for compounds.chemlib in the current working directory" << std::endl;
  } else {
    fs::path path_and_exe(exe_path);
    fs::path path_only=path_and_exe.branch_path();
    fs::path share_path = path_only.branch_path();
    share_path = share_path / "share" / "openstructure" / "compounds.chemlib";

    String share_path_string=BFPathToString(share_path);
      
    if (fs::exists(share_path_string)) {
      return CompoundLib::Load(share_path_string);
    }  
  }
  if (!lib) {
    std::cerr << "Could not load compounds.chemlib" << std::endl;
    exit(-1);
  }
  return CompoundLibPtr();
}

int main (int argc, char **argv)
{
  // sets some default values for parameters
  String version = OST_VERSION_STRING;
  lDDTSettings settings;
  String parameter_file_path;
  bool structural_checks = false;
  bool ignore_consistency_checks = false;
  Real bond_tolerance = 12.0;
  Real angle_tolerance = 12.0;
  String sel;
  // creates the required loading profile
  IOProfile profile;
  // parses options
  String custom_path;
  po::options_description desc("Options");
  desc.add_options()
    ("calpha,c", "consider only calpha atoms")
    ("sel,s", po::value<String>(&sel)->default_value(""), "selection performed on reference structure")
    ("tolerant,t", "fault tolerant mode")
    ("structural-checks,f", "perform stereo-chemical and clash checks")
    ("ignore-consistency-checks,x", "ignore residue name consistency checks")
    ("version,e", "version")
    ("parameter-file,p", po::value<String>(), "stereo-chemical parameter file")
    ("verbosity,v", po::value<int>(), "verbosity level")
    ("bond_tolerance,b", po::value<Real>(), "tolerance in stddev for bonds")
    ("complib,l", po::value<String>(&custom_path)->default_value(""), "location of the compound library file")
    ("angle_tolerance,a", po::value<Real>(), "tolerance in stddev for angles")
    ("inclusion_radius,r", po::value<Real>(), "distance inclusion radius")
    ("sequence_separation,i", po::value<int>(), "sequence separation")
    ("files", po::value< std::vector<String> >(), "input file(s)")
    ("reference",po::value<String>(),"reference(s)")
  ;
  po::positional_options_description p;
  p.add("files", -1);
  po::variables_map vm;
  try {
    po::store(po::command_line_parser(argc, argv).
                options(desc).positional(p).run(),
              vm);
  } catch (std::exception& e) {
    std::cout << e.what() << std::endl;
    usage();
    exit(-1);
  }
  po::notify(vm);
  conop::CompoundLibPtr lib = load_compound_lib(custom_path);
  if (!lib) {
    exit(0);
  }
  profile.processor = conop::RuleBasedProcessorPtr(new conop::RuleBasedProcessor(lib));
  profile.processor->SetCheckBondFeasibility(false);
  if (vm.count("version")) {
    std::cout << "Version: " << version << std::endl;
    exit(0);
  }
  std::vector<String> files;
  if (vm.count("files")) {
    files=vm["files"].as<std::vector<String> >();
  } else {
    usage();
    exit(-1);
  }
  if (vm.count("calpha")) {
    profile.calpha_only=true;
  }
  if (vm.count("structural-checks")) {
    structural_checks=true;
  }
  if (vm.count("ignore-consistency-checks")) {
    ignore_consistency_checks=true;
  }
  if (vm.count("tolerant")) {
    profile.fault_tolerant=true;
  }

  if (vm.count("parameter-file")) {
      parameter_file_path=vm["parameter-file"].as<String>();
    } else if (structural_checks==true) {
    std::cout << "Please specify a stereo-chemical parameter file" << std::endl;
    exit(-1);
  }
  int verbosity_level=0;
  int ost_verbosity_level=2;
  if (vm.count("verbosity")) {
    verbosity_level=vm["verbosity"].as<int>();
    if (verbosity_level==0) {
      ost_verbosity_level=2;    
    } else if (verbosity_level==1) {
      ost_verbosity_level=3;
    } else if (verbosity_level==2) {
      ost_verbosity_level=4;
    } else {
      std::cout << "Verbosity level " << verbosity_level << " is not available" << std::endl;
      exit(-1);
    }
  }

  Logger::Instance().PushVerbosityLevel(ost_verbosity_level);

  if (vm.count("bond_tolerance")) {
    bond_tolerance=vm["bond_tolerance"].as<Real>();
  }
  if (vm.count("angle_tolerance")) {
    angle_tolerance=vm["angle_tolerance"].as<Real>();
  }
  if (vm.count("inclusion_radius")) {
    settings.radius=vm["inclusion_radius"].as<Real>();
  }
  if (vm.count("sequence_separation")) {
    settings.sequence_separation=vm["sequence_separation"].as<int>();
  }

  std::vector<EntityView> ref_list;  
    
  // loads the reference file and creates the list of distances to check in lddt    
  // if the reference file is a comma-separated list of files, switches to multi-
  // reference mode
  GlobalRDMap glob_dist_list; 
  String ref_file=files.back();
  ost::StringRef ref_file_sr(ref_file.c_str(),ref_file.length());
  std::vector<StringRef> ref_file_split_sr=ref_file_sr.split(',');
  for (std::vector<StringRef>::const_iterator ref_file_split_sr_it = ref_file_split_sr.begin();
       ref_file_split_sr_it != ref_file_split_sr.end();++ref_file_split_sr_it) {
    String ref_filename = ref_file_split_sr_it->str();
    EntityHandle ref=load(ref_filename, profile);
    if (!ref) {
      exit(-1);
    }
    if (sel != ""){
      std::cout << "Performing \"" << sel << "\" selection on reference " << ref_filename << std::endl;
      try {
        ref_list.push_back(ref.Select(sel));
      } catch (const ost::mol::QueryError& e) {
        std::cerr << "Provided selection argument failed." << std::endl << e.GetFormattedMessage() << std::endl;
        exit(-1);
      }
    }
    else {
      ref_list.push_back(ref.CreateFullView());
    }
  }
  CleanlDDTReferences(ref_list);
  if (ref_list.size()==1) {
    std::cout << "Multi-reference mode: Off" << std::endl;
  } else {
    std::cout << "Multi-reference mode: On" << std::endl;
  }
  glob_dist_list = PreparelDDTGlobalRDMap(ref_list,
                                          settings.cutoffs,
                                          settings.sequence_separation,
                                          settings.radius);
  files.pop_back();

  // prints out parameters used in the lddt calculation
  std::cout << "Verbosity level: " << verbosity_level << std::endl;
  settings.PrintParameters();
  if (structural_checks) {
    LOG_INFO("Log entries format:");
    // verbosity level/format must match filter_clashes::CheckStereoChemistry
    LOG_INFO("BOND INFO FORMAT:  Chain  Residue  ResNum  Bond  Min  Max  Observed  Z-score  Status");
    LOG_INFO("ANGLE INFO FORMAT:  Chain  Residue  ResNum  Angle  Min  Max  Observed  Z-score  Status");
    // verbosity level/format must match filter_clashes::FilterClashes
    LOG_INFO("CLASH INFO FORMAT:  Chain1  Residue1  ResNum1  Atom1  Chain2  Residue2  ResNum2  Atom2  Threshold  Observed  Difference  Status");
  }
  // verbosity level/format must match local_dist_diff_test::calc_overlap1
  LOG_VERBOSE("LDDT INFO FORMAT:  Chain1  Residue1  ResNum1  Atom1  Chain2  Residue2  ResNum2  Atom2  ModelDist  TargetDistMin  TargetDistMax  Tolerance  Status");

  // error if the reference structure is empty
  if (glob_dist_list.size()==0) {
    std::cout << "ERROR: No valid distance to check in the reference structure(s). Please check that the first chain of the reference structure(s) contains a peptide." << std::endl;
    exit(-1);
  }

  // cycles through the models to evaluate 
  for (size_t i=0; i<files.size(); ++i) {
    EntityHandle model=load(files[i], profile);
    if (!model) {
      if (!profile.fault_tolerant) {
        exit(-1);
      }
      continue;
    }
    EntityView model_view = model.GetChainList()[0].Select("peptide=true");

    boost::filesystem::path pathstring(files[i]);
    String filestring=BFPathToString(pathstring);
    std::cout << "File: " << files[i] << std::endl;

    if (structural_checks) {
      StereoChemicalProps stereochemical_params;
      try {
        stereochemical_params = ost::io::ReadStereoChemicalPropsFile(parameter_file_path, true);
      } catch (std::exception& e) {
        std::cout << e.what() << std::endl;
        exit(-1);
      }
      
      try {
        CheckStructure(model_view,
                       stereochemical_params.bond_table,
                       stereochemical_params.angle_table,
                       stereochemical_params.nonbonded_table,
                       bond_tolerance,
                       angle_tolerance);
      } catch (std::exception& e) {
        std::cout << e.what() << std::endl;
        exit(-1);
      }
    }

    // Check consistency
    for (std::vector<EntityView>::const_iterator ref_list_it = ref_list.begin();
         ref_list_it != ref_list.end(); ++ref_list_it) {
      bool cons_check = ResidueNamesMatch(model_view,*ref_list_it, ignore_consistency_checks);
      if (cons_check==false) {
        if (ignore_consistency_checks==false) {
          throw std::runtime_error("Residue names in model and in reference structure(s) are inconsistent.");
        } else {
          LOG_WARNING("Residue names in model and in reference structure(s) are inconsistent.");
        }
      }
    }

    // computes the lddt score   
    LocalDistDiffTest(model_view, ref_list, glob_dist_list, settings);

    // prints the residue-by-residue statistics
    std::vector<lDDTLocalScore> local_scores;
    EntityView outview = model.GetChainList()[0].Select("peptide=true");
    local_scores = GetlDDTPerResidueStats(outview, glob_dist_list, structural_checks, settings.label);
    PrintlDDTPerResidueStats(local_scores, structural_checks, settings.cutoffs.size());
  }
  return 0;
}
