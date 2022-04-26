#include "worm.hpp"

void worm::print_params(std::ostream& os) const {
#ifdef UNISYS
  os << "# config : UNIFORM                                  : yes\n";
#else
  os << "# config : UNIFORM                                  : no\n";
#endif
#ifdef CAN_WINDOW
  os << "# config : CWINDOW                                  : yes\n";
#else
  os << "# config : CWINDOW                                  : no\n";
#endif
  Diagram_type::print_name(os);
  os << "# run time limit                                    : " << runtimelimit << "\n";
  os << "# reset statistics when restored (requires hack)    : " << reset_statistics << "\n";
  os << "# maximum number of sweeps                          : " << total_sweeps << "\n";
  os << "# maximum number of sweeps for thermalization       : " << thermalization_sweeps << "\n";
  os << "# inverse temperature                               : " << beta << "\n";
  os << "# canonical measurement value                       : " << canonical << "\n";
  os << "# random number seed                                : " << parameters["seed"].as<size_t>() << "\n";
  os << "# Ntest                                             : " << Ntest<< "\n";
  os << "# Nsave                                             : " << Nsave<< "\n";
  os << "# Nmeasure O(1)                                     : " << Nmeasure << "\n";
  os << "# Nmeasure O(N)                                     : " << Nmeasure2 << "\n";
  os << "# energy shift                                      : " << E_off << "\n";
  os << "# weight worm configurations                        : " << C_worm << "\n";
  os << "# update probability insert worm                    : " << update_prob[update_tag::insertworm]    << "\t" << update_prob_cuml[update_tag::insertworm]    << "\n";
  os << "# update probability move worm                      : " << update_prob[update_tag::moveworm]      << "\t" << update_prob_cuml[update_tag::moveworm]    << "\n";
  os << "# update probability insert kink                    : " << update_prob[update_tag::insertkink]    << "\t" << update_prob_cuml[update_tag::insertkink]    << "\n";
  os << "# update probability delete kink                    : " << update_prob[update_tag::deletekink]    << "\t" << update_prob_cuml[update_tag::deletekink]    << "\n";
  os << "# update probability glue worm                      : " << update_prob[update_tag::glueworm]      << "\t" << update_prob_cuml[update_tag::glueworm]      << "\n";
  
  MyLatt->print_name(os);
  MyLatt->print_params(os);
  MyModel->print_params(os);
  
}


boost::multi_array<double,2> worm::serialize_op_string() const {
  size_t maxsize = 0;
  for (auto const& vertexlist : operator_string){
    if (vertexlist.size() > maxsize)
      maxsize = vertexlist.size();
  }
  boost::multi_array<double, 2> op_string_arr(
    boost::extents[operator_string.size()][maxsize*(zcmax + 5) + 1]);

  auto row_it = op_string_arr.begin();

  for (auto const& vertexlist : operator_string){
    auto col_it = (row_it++)->begin();
    *(col_it++) = vertexlist.size();
    for (auto const& vertex : vertexlist){
      *(col_it++) = vertex.before();
      *(col_it++) = vertex.after();
      *(col_it++) = vertex.link();
      *(col_it++) = vertex.time();
      *(col_it++) = vertex.color();
      if (vertex.link() != -1) {
        for (size_t d = 0; d < zcmax; ++d, ++col_it)
          *col_it = vertex.get_assoc_time(d);
      }
    }
  }
  return op_string_arr;
}


void worm::save(alps::hdf5::archive & ar) const {
  // Most of the save logic is already implemented in the base class
  alps::mcbase::save(ar);

  // random number engine
  std::ostringstream engine_ss;
  engine_ss << MyGenerator;

  std::vector<double> wind_vec(mWinding.data(), mWinding.data() + mWinding.size());

  ar["checkpoint/random"] << engine_ss.str();
  ar["checkpoint/sweeps"] << sweeps;
  ar["checkpoint/counter"] << counter;
  ar["checkpoint/update_statistics"] << update_statistics;
  ar["checkpoint/configuration/operator_string"] << serialize_op_string();
  ar["checkpoint/configuration/mWinding"] << wind_vec;
  ar["checkpoint/configuration/nrvertex"] << nrvertex;
  ar["checkpoint/configuration/nr_of_particles"] << number_of_particles;
  ar["checkpoint/configuration/state"] << state;
#ifdef UNISYS
  ar["checkpoint/configuration/hist_densmat"] << hist_densmat;
#endif
#ifdef CAN_WINDOW
  ar["checkpoint/configuration/hist_gt"] << hist_gt;
#endif

}


void worm::load(alps::hdf5::archive & ar) {
  alps::mcbase::load(ar);

  worm_diag = true;
  worm_at_stop = 0;
  worm_passes_nb_kink = 0;
  worm_dtime = 0.;
  new_measurement = true;
  worm_meas_densmat = false;


  auto find_elem = [&](Diagram_type::iterator begin, Diagram_type::iterator end, 
                       double findTime) -> Diagram_type::iterator
  {
    for(auto it = begin; it != end; ++it) {
        if (it->time() == findTime)
            return it;
    }
    return end;
  };

  std::string engine_str;
  std::vector<double> wind_vec;
  boost::multi_array<double, 2> op_string_arr;

  ar["checkpoint/random"] >> engine_str;
  ar["checkpoint/sweeps"] >> sweeps;
  ar["checkpoint/counter"] >> counter;
  ar["checkpoint/update_statistics"] >> update_statistics;
  ar["checkpoint/configuration/operator_string"] >> op_string_arr;
  ar["checkpoint/configuration/mWinding"] >> wind_vec;
  ar["checkpoint/configuration/nrvertex"] >> nrvertex;
  ar["checkpoint/configuration/nr_of_particles"] >> number_of_particles;
  ar["checkpoint/configuration/state"] >> state;
#ifdef UNISYS
  ar["checkpoint/configuration/hist_densmat"] >> hist_densmat;
#endif
#ifdef CAN_WINDOW
  ar["checkpoint/configuration/hist_gt"] >> hist_gt;
#endif

  std::istringstream engine_ss(engine_str);
  engine_ss >> MyGenerator;

  operator_string.clear();
  for (auto const& row : op_string_arr){
    Diagram_type vertexlist;
    auto col_it = row.begin() + 1;
    auto end = col_it + *row.begin() * (zcmax + 5);
    while (col_it != end) {
      auto before = StateType(*(col_it++));
      auto after = StateType(*(col_it++));
      auto link = SiteIndex(*(col_it++));
      auto time = *(col_it++);
      auto color = int(*(col_it++));
      vertexlist.insert(vertexlist.end(), Element(before, after, link, time, color));
      col_it += zcmax;
    }
    operator_string.push_back(vertexlist);
  }

  for (SiteIndex i = 0; i < Nsites+1; i++) {
    dummy_it[i] = find_elem(operator_string[i].begin(), operator_string[i].end(), beta);
    Nprtcls += dummy_it[i]->before() * beta;
  }

  std::cout << "# Restoring existing associations ... ";
  for (SiteIndex i = 0; i < Nsites; i++) {
    for (size_t j = 0; j < zcmax; j++) {
      if (MyLatt->has_neighbor(i,j))
        dummy_it[i]->set_assoc(j, dummy_it[nb[i][j]]);
      else
        dummy_it[i]->set_assoc(j, *(dummy_it.end() - 1));
    }
  }

  for (SiteIndex i = 0; i < Nsites; i++) {
    size_t elem_idx = 0;
    for (auto elem_it = operator_string[i].begin(); elem_it != operator_string[i].end(); ++elem_it, ++elem_idx) {

      std::vector<double> assoc_vec(zcmax);
      auto row_it = op_string_arr.begin() + i;
      auto col_it = row_it->begin() + 1 + elem_idx*(zcmax + 5) + 5;
      std::copy(col_it, col_it + zcmax, assoc_vec.begin());

      for (size_t j = 0; j < zcmax; j++) {
        if (MyLatt->has_neighbor(i,j)) {
            auto assoc_elem_it = find_elem(operator_string[nb[i][j]].begin(), operator_string[nb[i][j]].end(), 
                                           assoc_vec.at(j));
            if (assoc_elem_it != operator_string[nb[i][j]].end())
                elem_it->set_assoc(j, assoc_elem_it);
            else
                throw std::runtime_error("Element not found during reconstruction of associations.");
        }
        else {
            elem_it->set_assoc(j, *(dummy_it.end() - 1));
        }
      }
    }
  }

  mWinding = Eigen::Map<LATTICE::Direction>(wind_vec.data());

  worm_head_it = dummy_it[0];
  worm_tail_it = dummy_it[0];

  std::cout << "...done. Computing potential energies...";
  Epot_tot = calc_potential_energy_nb() + calc_potential_energy_loc();
  std::cout << "...done.\n";
  cout << "# Potential Energy tot : " << Epot_tot << endl;
  std::cout << "\n# Finished loading.\n";
  reset_statistics = parameters["reset_statistics"];
  if (reset_statistics == 1) {
    std::cout << "# Resetting statistics...\n";
    force_reset_statistics();
  }

}

void worm::print_conf(std::ostream& os) const {
  os << "\n\n Printing operator string";
  for (SiteType i = 0; i < Nsites; i++) {
    os << "\nSite : " << i;
    int n = operator_string[i].size();
    int ii = 0;
    for (Diagram_type::const_iterator it = operator_string[i].begin(); it != operator_string[i].end(); ++it, ++ii) {
      it->print();
      if (ii > n) {
        os << "\n Error with list!\n";
        char ch; cin >> ch;
      }
    }
    os << "\n--------------------\n\n\n";
  }
  os << "\n Worm head iterator time : "<< worm_head_it->time()  << "\n";
  os << "\n Worm tail iterator time : "<< worm_tail_it->time()  << "\n";
  os << "\n worm_at_stop = " << worm_at_stop << "\t worm_passes_nb_kink = " << worm_passes_nb_kink << "\t worm_measdensmat = " << worm_meas_densmat << "\n";
}


void worm::print_update_statistics(std::ostream& os) const {
  std::vector<string> name;
  os << setprecision(10);
  name.push_back("INSERT WORM       ");
  name.push_back("MOVE WORM         ");
  name.push_back("INSERT KINK       ");
  name.push_back("DELETE KINK       ");
  name.push_back("GLUE WORM         ");
  os << "\n\n# UPDATE STATISTICS";
  os << "\n" << "# col 1 : all updates"
     << "\n" << "# col 2 : impossible updates"
     << "\n" << "# col 3 : rejected updates"
     << "\n" << "# col 4 : accepted updates"
     << "\n" << "# col 5 : acceptance factor with respect to Metropolis ratio only"
     << "\n" << "# col 6 : acceptance factor with respect to all attempts.\n";
  for (size_t i = 0; i < update_count; i++) {
    os << "\n" << name[i]
       << "\t" << update_statistics[total_attempted][i]
       << "\t" << update_statistics[impossible][i]
       << "\t" << update_statistics[rejected][i]
       << "\t" << update_statistics[accepted][i]
       << "\t" << update_statistics[accepted][i] / (update_statistics[total_attempted][i] - update_statistics[impossible][i])
       << "\t" << update_statistics[accepted][i] / update_statistics[total_attempted][i];
  }
  os << "\n\n";
  double s = 0;
  for (size_t i = 0; i < update_count; i++) {
    s += update_statistics[total_attempted][i];
  }
  os << "# Total number of steps : " << s << "\t or 10^" << log10(s) << "\n";
  
}

