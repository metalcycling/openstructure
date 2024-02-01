#include <ost/seq/alg/hmm_pseudo_counts.hh>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include <limits>
#include <cmath>

#include <ost/string_ref.hh>


namespace {

  // to mimic the HHblits behaviour, R is based on the Gonnet substitution matrix
  // every entry R[a][b] corresponds to P(a|b) 
  const Real R[20][20] = {
    {0.133435, 0.0861429, 0.0716694, 0.0767892, 0.0451931, 0.0861612, 0.0638644, 0.0638589, 0.0700223, 0.0582424, 0.065352, 0.0716654, 0.0822722, 0.0733284, 0.0668668, 0.0989164, 0.0881558, 0.0785747, 0.0335151, 0.0462615},
    {0.0212019, 0.266857, 0.00903972, 0.009469, 0.0157131, 0.01193, 0.0140111, 0.0146727, 0.00991927, 0.013372, 0.0153475, 0.0124885, 0.00924679, 0.0108782, 0.0113849, 0.0193327, 0.0168459, 0.018897, 0.0149936, 0.0168408},
    {0.0504801, 0.0258694, 0.15964, 0.100729, 0.0191763, 0.0553552, 0.0593027, 0.0225536, 0.0606905, 0.0215294, 0.0270972, 0.0897624, 0.0460354, 0.0665564, 0.0504718, 0.0606802, 0.0540841, 0.027733, 0.0163566, 0.0283808},
    {0.0598873, 0.0300042, 0.111533, 0.1372, 0.0243854, 0.0498144, 0.0656566, 0.0321571, 0.0789513, 0.0314283, 0.0377994, 0.0736845, 0.0533622, 0.0885886, 0.0656581, 0.0627186, 0.0585313, 0.0386611, 0.0222899, 0.0321599},
    {0.0206017, 0.0291031, 0.0124111, 0.0142537, 0.175334, 0.0105647, 0.0341859, 0.0440577, 0.0163643, 0.0554499, 0.0505966, 0.0171498, 0.0145874, 0.0192437, 0.0167459, 0.0183671, 0.021082, 0.0358121, 0.0801796, 0.113219},
    {0.0848077, 0.0477099, 0.0773563, 0.0628701, 0.0228112, 0.34552, 0.054741, 0.0268108, 0.0586764, 0.0274347, 0.0337462, 0.0828826, 0.0522808, 0.0600294, 0.0600437, 0.0828881, 0.058672, 0.0353579, 0.0300674, 0.0300801},
    {0.0306873, 0.0273537, 0.0404564, 0.0404524, 0.0360342, 0.0267232, 0.146872, 0.0222367, 0.0423625, 0.0238214, 0.0273249, 0.0486327, 0.0286452, 0.0486303, 0.042362, 0.0352323, 0.0344236, 0.0232737, 0.0307088, 0.0222177},
    {0.0420777, 0.0392812, 0.0210988, 0.027169, 0.0636824, 0.017948, 0.0304931, 0.127064, 0.0311844, 0.0963867, 0.0899444, 0.0265473, 0.0278066, 0.0326653, 0.0291056, 0.0334299, 0.0440638, 0.10328, 0.0334349, 0.0430405},
    {0.0544335, 0.0313295, 0.0669828, 0.0786964, 0.0279058, 0.0463415, 0.0685348, 0.0367906, 0.124704, 0.0368031, 0.0432189, 0.0717651, 0.0519718, 0.0842986, 0.111139, 0.0610879, 0.0610786, 0.040354, 0.0266998, 0.0368013},
    {0.0759224, 0.0708227, 0.0398451, 0.0525312, 0.158562, 0.0363335, 0.0646247, 0.190686, 0.0617143, 0.251356, 0.190682, 0.0501533, 0.0589235, 0.0692529, 0.0602971, 0.0617101, 0.074195, 0.151465, 0.0852309, 0.100056},
    {0.0187229, 0.0178647, 0.0110218, 0.0138856, 0.0317983, 0.00982238, 0.0162919, 0.0391074, 0.0159279, 0.0419077, 0.0592039, 0.0132612, 0.0126454, 0.0174664, 0.0148745, 0.0159425, 0.0191539, 0.0317932, 0.0174792, 0.0210003},
    {0.0375111, 0.0265585, 0.066705, 0.0494529, 0.0196915, 0.0440748, 0.052976, 0.0210883, 0.0483208, 0.0201381, 0.0242281, 0.096418, 0.0326617, 0.0472207, 0.0430638, 0.0494368, 0.0450912, 0.0242096, 0.0175593, 0.0291164},
    {0.0486405, 0.0222116, 0.0386411, 0.0404524, 0.0189187, 0.0314024, 0.0352449, 0.0249495, 0.039526, 0.026724, 0.0260953, 0.0368921, 0.261161, 0.0433597, 0.0369035, 0.0497801, 0.0464422, 0.0299902, 0.0143521, 0.0222177},
    {0.0312223, 0.0188189, 0.0402341, 0.0483655, 0.0179742, 0.0259677, 0.0430922, 0.0211081, 0.0461724, 0.0226204, 0.0259586, 0.0384126, 0.0312272, 0.0608874, 0.0461829, 0.0342238, 0.0326925, 0.0231499, 0.0175593, 0.0220909},
    {0.0447524, 0.0309584, 0.0479587, 0.0563456, 0.0245857, 0.0408271, 0.059004, 0.0295632, 0.0956849, 0.0309579, 0.0347482, 0.0550639, 0.0417761, 0.072593, 0.151669, 0.049072, 0.049074, 0.0324263, 0.0355997, 0.0339353},
    {0.0601482, 0.0477629, 0.0523859, 0.0489009, 0.0244998, 0.0512063, 0.0445856, 0.0308503, 0.0477837, 0.0287859, 0.0338373, 0.057432, 0.0511994, 0.0488754, 0.0445844, 0.077481, 0.065948, 0.0370921, 0.0218089, 0.0301309},
    {0.0817285, 0.0634542, 0.0711878, 0.0695787, 0.0428747, 0.0552625, 0.0664169, 0.0619975, 0.072842, 0.0527675, 0.061982, 0.0798664, 0.0728268, 0.0711834, 0.067978, 0.100547, 0.12659, 0.0711838, 0.0318313, 0.0459572},
    {0.0744872, 0.0727841, 0.0373259, 0.0469937, 0.0744727, 0.0340536, 0.0459161, 0.148588, 0.0492103, 0.110149, 0.105201, 0.0438468, 0.0480877, 0.0515413, 0.0459295, 0.0578264, 0.0727877, 0.159228, 0.0400096, 0.0564827},
    {0.00545379, 0.00991306, 0.0037789, 0.00465085, 0.0286213, 0.00497084, 0.0103997, 0.0082571, 0.00558903, 0.0106396, 0.00992804, 0.00545903, 0.0039503, 0.00671079, 0.00865564, 0.00583628, 0.00558714, 0.00686789, 0.328817, 0.0321092},
    {0.0237983, 0.0351993, 0.0207284, 0.0212132, 0.127766, 0.0157211, 0.0237862, 0.0336026, 0.0243534, 0.0394855, 0.0377084, 0.0286163, 0.0193322, 0.02669, 0.0260839, 0.0254908, 0.025501, 0.0306509, 0.101507, 0.237902}
  };

}

namespace ost{ namespace seq{ namespace alg{


void ContextProfileDB::Save(const String& filename) const {

  if(profiles_.empty()) {
    throw Error("Cannot save empty ContextProfileDB");
  }

  std::ofstream out_stream(filename.c_str(), std::ios::binary);

  if (!out_stream){
    std::stringstream ss;
    ss << "the file '" << filename << "' could not be opened.";
    throw Error(ss.str());
  }

  std::vector<Real> data;
  for(auto it = profiles_.begin(); it != profiles_.end(); ++it) {
    const std::vector<Real>& p_data = it->GetData();
    data.insert(data.end(), p_data.begin(), p_data.end());
  }

  uint32_t magic_number = 424242;
  out_stream.write(reinterpret_cast<char*>(&magic_number), sizeof(uint32_t));

  uint8_t version = 1;
  out_stream.write(reinterpret_cast<char*>(&version), sizeof(uint8_t));

  uint32_t length = profiles_[0].GetLength();
  out_stream.write(reinterpret_cast<char*>(&length),sizeof(uint32_t));

  uint32_t n_profiles = profiles_.size();
  out_stream.write(reinterpret_cast<char*>(&n_profiles),sizeof(uint32_t));

  uint32_t data_size = data.size();
  out_stream.write(reinterpret_cast<char*>(&data_size),sizeof(uint32_t));

  out_stream.write(reinterpret_cast<char*>(&data[0]), data_size*sizeof(Real));
}


ContextProfileDBPtr ContextProfileDB::Load(const String& filename) {

  std::ifstream in_stream(filename.c_str(), std::ios::binary);
  if (!in_stream){
    std::stringstream ss;
    ss << "the file '" << filename << "' could not be opened.";
    throw Error(ss.str());
  }

  uint32_t magic_number;
  in_stream.read(reinterpret_cast<char*>(&magic_number), sizeof(uint32_t));

  if(magic_number != 424242) {
    std::stringstream ss;
    ss << "Could not read magic number in " << filename<<". Either the file ";
    ss << "is corrupt or does not contain a ContextProfileDB.";
    throw Error(ss.str());
  }

  uint8_t version;
  in_stream.read(reinterpret_cast<char*>(&version), sizeof(uint8_t));
  if(version != 1) {
    std::stringstream ss;
    ss << "ContextProfileDB in " << filename << " is of version " << version;
    ss << " but only version 1 can be read.";
    throw Error(ss.str());
  }

  uint32_t length;
  in_stream.read(reinterpret_cast<char*>(&length), sizeof(uint32_t));

  uint32_t n_profiles;
  in_stream.read(reinterpret_cast<char*>(&n_profiles),sizeof(uint32_t));

  uint32_t data_size;
  in_stream.read(reinterpret_cast<char*>(&data_size), sizeof(uint32_t));

  std::vector<Real> data(data_size, 0.0);
  in_stream.read(reinterpret_cast<char*>(&data[0]), data_size*sizeof(Real));

  ContextProfileDBPtr db(new ContextProfileDB);

  int data_loc = 0;
  int profile_data_size = ContextProfile::DataSize(length);
  for(uint i = 0; i < n_profiles; ++i, data_loc+=profile_data_size) {
    db->profiles_.push_back(ContextProfile(length, &data[data_loc]));
  }

  return db;
}

ContextProfileDBPtr ContextProfileDB::FromCRF(const String& filename) {

  // open it up
  boost::iostreams::filtering_stream<boost::iostreams::input> in;
  boost::filesystem::ifstream stream(filename);
  if(!stream) {
    throw Error("Could not open " + filename);
  }
  // add unzip if necessary
  if(boost::iequals(".gz", boost::filesystem::extension(filename))) {
    in.push(boost::iostreams::gzip_decompressor());
  }
  in.push(stream);

  // tmp. storage
  std::string line;
  ost::StringRef sline;
  std::vector<ost::StringRef> chunks;

  int size = -1;
  int length = -1;

  while(std::getline(in, line)) {

    sline = ost::StringRef(line.c_str(), line.length());

    if(sline.length()>=8 &&
       sline.substr(0, 8) == ost::StringRef("CrfState", 8)) {
      throw Error("Require to read SIZE and LENG before first ContextProfile");
    }

    if(sline.length()>4 && 
        sline.substr(0, 4) == ost::StringRef("SIZE", 4)) {
      chunks = sline.split();
      if(chunks.size() != 2) {
        throw Error("Badly formatted line: " + line);
      }
      std::pair<bool, int> s = chunks[1].to_int();
      if (!s.first) {
        throw Error("Badly formatted line: " + line);
      }
      size = s.second;
    }

    if(sline.length()>4 && 
        sline.substr(0, 4) == ost::StringRef("LENG", 4)) {
      chunks = sline.split();
      if(chunks.size() != 2) {
        throw Error("Badly formatted line: " + line);
      }
      std::pair<bool, int> l = chunks[1].to_int();
      if (!l.first) {
        throw Error("Badly formatted line: " + line);
      }
      length = l.second;
    }

    if(size != -1 && length != -1) {
      break;
    }
  }

  ContextProfileDBPtr db(new ContextProfileDB);
  bool in_crf_state=false;
  ContextProfile current_context_profile(length);

  // stuff that we read for every context profile
  Real cp_bias = std::numeric_limits<Real>::quiet_NaN();
  int cp_leng = -1;
  int cp_alph = -1;
  std::vector<char> cp_olcs;
  std::vector<std::vector<Real> > cp_weights;
  std::vector<Real> cp_pc;

  while(std::getline(in, line)) {

    sline = ost::StringRef(line.c_str(), line.length());

    if(!in_crf_state) {
      if(sline.length()>=8 &&
         sline.substr(0, 8) == ost::StringRef("CrfState", 8)) {
        in_crf_state = true;
      }
      continue;
    }

    if(sline.length()>4 && 
        sline.substr(0, 4) == ost::StringRef("BIAS", 4)) {
      chunks = sline.split();
      if(chunks.size() != 2) {
        throw Error("Badly formatted line: " + line);
      }
      std::pair<bool, Real> b = chunks[1].to_float();
      if (!b.first) {
        throw Error("Badly formatted line: " + line);
      }
      cp_bias = b.second;
      continue;
    }

    if(sline.length()>4 && 
        sline.substr(0, 4) == ost::StringRef("LENG", 4)) {
      chunks = sline.split();
      if(chunks.size() != 2) {
        throw Error("Badly formatted line: " + line);
      }
      std::pair<bool, int> l = chunks[1].to_int();
      if (!l.first) {
        throw Error("Badly formatted line: " + line);
      }
      cp_leng = l.second;
      continue;
    }

    if(sline.length()>4 && 
        sline.substr(0, 4) == ost::StringRef("ALPH", 4)) {

      chunks = sline.split();
      if(chunks.size() != 2) {
        throw Error("Badly formatted line: " + line);
      }
      std::pair<bool, int> a = chunks[1].to_int();
      if (!a.first) {
        throw Error("Badly formatted line: " + line);
      }
      cp_alph = a.second;
      continue;
    }

    if(sline.length()>7 && 
        sline.substr(0, 7) == ost::StringRef("WEIGHTS", 7)) {

      if(cp_alph == -1 || cp_leng == -1) {
        throw Error("Require LENG and ALPH before reading WEIGHTS in entry");
      }

      chunks = sline.split();
      if(chunks.size() != static_cast<uint>(cp_alph+1)) {
        throw Error("Badly formatted line: " + line);
      }
      for(int i = 0; i < cp_alph; ++i) {
        cp_olcs.push_back(chunks[i+1][0]);
      }

      cp_weights = std::vector<std::vector<Real> >(cp_leng, std::vector<Real>());

      // the next cp_leng lines should be the weights
      for(int i = 0; i < cp_leng; ++i) {
        if(!std::getline(in, line)) {
          throw Error("Failed to load all WEIGHTS in entry");
        }
        sline = ost::StringRef(line.c_str(), line.length());

        if(sline.length()>=2 &&
           sline.substr(0, 2) == ost::StringRef("//", 2)) {
          throw Error("Arrived at end of entry before reading all weights");
        }        

        chunks = sline.split();

        if(chunks.size() != static_cast<uint>(cp_alph + 1)) {
          throw Error("Failed to load all WEIGHTS in entry");
        }

        // read the position of the weights. if it can't be parsed to int, 
        // something is fishy (e.g. already another key word)
        std::pair<bool, int> p = chunks[0].to_int();
        if(!p.first) {
          throw Error("Badly formatted line: " + line);
        }

        if(p.second < 1 || p.second > cp_leng) {
          throw Error("Badly formatted line: " + line);
        }

        for(int j = 0; j < cp_alph; ++j) {
          std::pair<bool, Real> w = chunks[j+1].to_float();
          if(!w.first) {
            throw Error("Badly formatted line: " + line);
          }
          cp_weights[p.second-1].push_back(0.001 * w.second);
        }
      }
      continue;
    }

    if(sline.length()>2 && 
        sline.substr(0, 2) == ost::StringRef("PC", 2)) {

      chunks = sline.split();
      if(chunks.size() != static_cast<uint>(cp_alph+1)) {
        throw Error("Badly formatted line: " + line);
      }
      for(int i = 0; i < cp_alph; ++i) {
        std::pair<bool, Real> w = chunks[i+1].to_float();
        if(!w.first) {
          throw Error("Badly formatted line: " + line);
        }
        cp_pc.push_back(0.001*w.second);
      }
    }

    if(sline.length()>=2 &&
       sline.substr(0, 2) == ost::StringRef("//", 2)) {

      //check if data is OK
      if(std::isnan(cp_bias)) {
        throw Error("Oberved entry without BIAS value");
      }
      if(cp_leng == -1) {
        throw Error("Observed entry without LENG value");
      }
      if(cp_leng != length) {
        throw Error("Require all entries to be of same length as specified in file header");
      }
      if(cp_alph == -1) {
        throw Error("Observed Entry without ALPH value");
      }
      if(cp_alph != 20) {
        throw Error("Expect ALPH to be 20 for all entries");
      }
      if(cp_olcs.empty()) {
        throw Error("Observed Entry without WEIGHTS value");          
      }
      if(cp_olcs.size() != 20) {
        throw Error("Expect exactly 20 items after WEIGHTS key word");
      }
      for(int i = 0; i < length; ++i) {
        if(cp_weights[i].size() != 20) {
          throw Error("Observed Entry with != 20 weight values for each pos");
        }
      }
      if(cp_pc.size() != 20) {
        throw Error("Observed Entry with != 20 PC values");
      }

      // fill data in context profile
      current_context_profile = ContextProfile(length);
      current_context_profile.SetBias(cp_bias);
      for(int i = 0; i < length; ++i) {
        for(int j = 0; j < 20; ++j) {
          current_context_profile.SetWeight(i, cp_olcs[j], cp_weights[i][j]);
        }
      }

      // magic rescaling of context profiles as implemented in hhblits
      Real max = -std::numeric_limits<Real>::max();
      for(int i = 0; i < 20; ++i) {
        max = std::max(max, cp_pc[i]);
      }

      Real sum = 0.0;
      for(int i = 0; i < 20; ++i) {
        sum += std::exp(cp_pc[i]-max);
      }

      Real tmp = max + std::log(sum);
      for(int i = 0; i < 20; ++i) {
        current_context_profile.SetPseudoCount(cp_olcs[i], std::exp(cp_pc[i]-tmp));
      }

      db->AddProfile(current_context_profile);

      // invalidate all variables for a next profile
      cp_bias = std::numeric_limits<Real>::quiet_NaN();
      cp_leng = -1;
      cp_alph = -1;
      cp_olcs = std::vector<char>();
      cp_weights = std::vector<std::vector<Real> >();
      cp_pc = std::vector<Real>();
      in_crf_state = false;
    }      
  }

  if(db->size() != static_cast<size_t>(size)) {
    throw Error("Number of read entries does not correspond to what was promised in the header");
  }

  return db;
}


void AddTransitionPseudoCounts(ost::seq::ProfileHandle& profile,
                               Real gapb, Real gapd, Real gape) {

  // a priori probabilities estimated with default values of HHblits
  Real pM2D = gapd * 0.0286;
  Real pM2I = pM2D; 
  Real pM2M = 1 - pM2D - pM2I;
  Real pI2I = 1.0 * gape / (gape - 1 + 1.0 / 0.75);;
  Real pI2M = 1 - pI2I;
  Real pD2D = pI2I;
  Real pD2M = 1 - pD2D;

  for (size_t col_idx = 0; col_idx < profile.size(); ++col_idx) {

    ost::seq::HMMDataPtr data = profile[col_idx].GetHMMData();

    // Transitions from M state
    Real neff = data->GetNeff();
    Real p0 = (neff - 1) * data->GetProb(ost::seq::HMM_M2M) + gapb*pM2M;
    Real p1 = (neff - 1) * data->GetProb(ost::seq::HMM_M2D) + gapb*pM2D;
    Real p2 = (neff - 1) * data->GetProb(ost::seq::HMM_M2I) + gapb*pM2I;
    Real sum = p0 + p1 + p2;
    data->SetProb(ost::seq::HMM_M2M, p0/sum);
    data->SetProb(ost::seq::HMM_M2D, p1/sum);
    data->SetProb(ost::seq::HMM_M2I, p2/sum);

    // Transitions from I state
    Real neff_i = data->GetNeff_I();
    p0 = neff_i * data->GetProb(ost::seq::HMM_I2M) + gapb*pI2M;
    p1 = neff_i * data->GetProb(ost::seq::HMM_I2I) + gapb*pI2I;
    sum = p0 + p1;
    data->SetProb(ost::seq::HMM_I2M, p0/sum);
    data->SetProb(ost::seq::HMM_I2I, p1/sum);

    // Transitions from D state
    Real neff_d = data->GetNeff_D();
    p0 = neff_d * data->GetProb(ost::seq::HMM_D2M) + gapb*pD2M;
    p1 = neff_d * data->GetProb(ost::seq::HMM_D2D) + gapb*pD2D;
    sum = p0 + p1;
    data->SetProb(ost::seq::HMM_D2M, p0/sum);
    data->SetProb(ost::seq::HMM_D2D, p1/sum);
  }
}


void AddAAPseudoCounts(ost::seq::ProfileHandle& profile,
                       Real a, Real b, Real c) {

  Real full_admixture [20];
  for(size_t col_idx = 0; col_idx < profile.size(); ++col_idx) {
    Real* col_freq = profile[col_idx].freqs_begin();
    for(int i = 0; i < 20; ++i) {
      full_admixture[i] = 0.0;
      for(int j = 0; j < 20; ++j) {
        full_admixture[i] += R[i][j] * col_freq[j];
      }
    }
    // tau estimated as in hhblits in diversity dependent mode:
    // tau = a/(1+((Neff[i]-1)/b)^c) with default values a=1.0, b=1.5, c=1.0
    // this is the equation they write in HHblits when you display the help
    // (Neff[i]-1) got rid of the -1. Well, that's how HHblits implements it
    Real neff = profile[col_idx].GetHMMData()->GetNeff();
    Real tau = std::min(1.0, a / (1.0 + std::pow((neff) / b, c)));
    for (int i = 0; i < 20; ++i) {
      col_freq[i] = (1. - tau) * col_freq[i] + tau * full_admixture[i];
    }
  }
}


void AddAAPseudoCounts(ost::seq::ProfileHandle& profile,
                       const ContextProfileDB& db,
                       Real a, Real b, Real c) {

  std::vector<Real> cp_scores(db.size(), 0.0);
  int cp_length = db.profile_length();
  if(cp_length % 2 != 1) {
    throw Error("Length of profiles in db must be an odd number");
  } 
  // extension from center to both directions
  int cp_ext = (cp_length - 1) / 2;

  std::vector<std::vector<Real> > 
  count_profile(profile.size(), std::vector<Real>(20, 0.0));
  std::vector<std::vector<Real> > 
  context_profile(profile.size(), std::vector<Real>(20, 0.0));

  // fill counts profile
  for(size_t col_idx = 0; col_idx < profile.size(); ++col_idx) {
    HMMDataPtr hmm_data = profile[col_idx].GetHMMData(); 
    Real neff = hmm_data->GetNeff();
    Real* col_freq = profile[col_idx].freqs_begin();
    std::vector<Real>& counts = count_profile[col_idx];
    for(int i = 0; i < 20; ++i) {
      counts[i] = col_freq[i] * neff;
    }
  }

  // process columns
  for(size_t col_idx = 0; col_idx < count_profile.size(); ++col_idx) {
    int min = std::max(0, static_cast<int>(col_idx) - cp_ext);
    int max = std::min(col_idx + cp_ext, count_profile.size() - 1); 
    int cp_min = static_cast<int>(col_idx) < cp_ext ? 
                 std::abs(static_cast<int>(col_idx) - cp_ext) : 0; 
    Real max_score = -std::numeric_limits<Real>::max();
    // estimate score for each context profile in db
    for(size_t cp_idx = 0; cp_idx < db.size(); ++cp_idx) {
      const ContextProfile& cp = db[cp_idx];
      Real score = cp.GetBias();
      for(int i = min, j = cp_min; i<=max; ++i, ++j) {
        const Real* weights = cp.GetWeights(j);
        const std::vector<Real>& counts = count_profile[i];
        for(int k = 0; k < 20; ++k) {
          score += weights[k] * counts[k];
        }
      }
      max_score = std::max(max_score, score);
      cp_scores[cp_idx] = score;
    }

    // same in hhblits code: log-sum-exp trick to avoid overflows
    Real summed_exp_score = 0.0;
    for(size_t i = 0; i < cp_scores.size(); ++i) {
      summed_exp_score += std::exp(cp_scores[i]-max_score);
    }
    Real tmp = max_score - std::log(summed_exp_score);

    std::vector<Real>& col_freq = context_profile[col_idx];
    for(size_t cp_idx = 0; cp_idx < db.size(); ++cp_idx) {
      Real w = std::exp(cp_scores[cp_idx]-tmp);
      const Real* cp_pc = db[cp_idx].GetPseudoCounts();
      for(int i = 0; i < 20; ++i) {
        col_freq[i] += w*cp_pc[i];
      }
    }
    // normalize
    Real sum = 0.0;
    for(int i = 0; i < 20; ++i) {
      sum += col_freq[i];
    }
    Real norm_factor = 1.0 / sum;
    for(int i = 0; i < 20; ++i) {
      col_freq[i] *= norm_factor;
    }
  }

  // mix together count and context profile to get final frequencies
  for(size_t col_idx = 0; col_idx < profile.size(); ++col_idx) {
    // tau estimated as in hhblits in diversity dependent mode:
    // tau = a/(1+((Neff[i]-1)/b)^c) with default values a=0.9, b=4.0, c=1.0
    // this is the equation they write in HHblits when you display the help
    // (Neff[i]-1) got rid of the -1. Well, that's how HHblits implements it
    Real neff = profile[col_idx].GetHMMData()->GetNeff();
    Real tau = std::min(1.0, a / (1.0 + std::pow((neff) / b, c)));
    Real* col_freq = profile[col_idx].freqs_begin();
    const std::vector<Real>& counts = count_profile[col_idx];
    const std::vector<Real>& context = context_profile[col_idx];
    for(int i = 0; i < 20; ++i) {
      col_freq[i] = tau*context[i] + (1.-tau)*counts[i]/neff;
    }
  }
}


void AddNullPseudoCounts(ost::seq::ProfileHandle& profile) {
  Real mixing_factor = 100.0 / profile.GetNeff();
  const Real* current_null_freq = profile.GetNullModel().freqs_begin();
  Real null_freq[20];
  for(int i = 0; i < 20; ++i) {
    null_freq[i] = current_null_freq[i] * mixing_factor;
  } 
  for(size_t i = 0; i < profile.size(); ++i) {
    const Real* freq = profile[i].freqs_begin();
    for(int j = 0; j < 20; ++j) {
      null_freq[j] += freq[j]; 
    }
  }
  // normalize
  Real summed_p = 0.0;
  for(int i = 0; i < 20; ++i) {
    summed_p += null_freq[i];
  }
  Real factor = 1.0/summed_p;
  for(int i = 0; i < 20; ++i) {
    null_freq[i] *= factor;
  }
  // create new nullmodel and set it
  ost::seq::ProfileColumn new_null_model;
  Real* new_null_freq = new_null_model.freqs_begin();
  for(int i = 0; i < 20; ++i) {
    new_null_freq[i] = null_freq[i];
  }
  profile.SetNullModel(new_null_model);
}


}}} // ns
