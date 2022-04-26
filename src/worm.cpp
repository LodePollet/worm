#include "worm.hpp"

#include <alps/params/convenience_params.hpp>

std::string worm::code_name() {
    return "Simulation of the Bose-Hubbard or the XXZ model with the worm algorithm";
}

// Defines the parameters for the worm simulation
void worm::define_parameters(parameters_type & parameters) {
    // If the parameters are restored, they are already defined
    if (parameters.is_restored()) {
        return;
    }
    
    // Adds the parameters of the base class
    alps::mcbase::define_parameters(parameters);
    // Adds the convenience parameters (for save/load)
    // followed by the worm specific parameters
    alps::define_convenience_parameters(parameters)
        .description(worm::code_name())
        .define<size_t>("runtimelimit",         60,      "run time limit in seconds")
        .define<int>("reset_statistics",        0,       "force reset statistics when restoring from old configuration")
        .define<double>("beta",                 1.0,     "inverse temperature")
        .define<int>("sweeps",                  1000,    "maximum number of sweeps")
        .define<int>("thermalization",          200,     "number of sweeps for thermalization")
        .define<size_t>("seed",                 91,      "seed for random number generation")
        .define<double>("E_off",                1.0,     "energy offset (technical parameter in move update")
        .define<double>("C_worm",               2.0,     "weight factor of worm configurations vs diagonal configurations (technical parameter in insertworm/glueworm updates")
        .define<int>("canonical",               -1,       "-1 for grand-canonical measurement; a positive number indicates the canonical measurement at this value")
        .define<unsigned long>("Ntest",         10000,    "test configuration after this number of updates "  )
        .define<unsigned long>("Nsave",         100000,   "save configuration after this number of updates"   )
        .define<size_t>("Nmeasure",             1,        "measure O(1) observables after this number of updates"   )
        .define<size_t>("Nmeasure2",            1,        "measure O(N) observables after this number of updates"   )
        .define<double>("p_insertworm",         1.0,      "update probability to insert worm in diagonal configuration")
        .define<double>("p_moveworm",           0.3,      "update probability to move worm head around")
        .define<double>("p_insertkink",         0.2,      "update probability to insert kink at position of worm head")
        .define<double>("p_deletekink",         0.2,      "update probability tp remove kink at position of worm head")
        .define<double>("p_glueworm",           0.3,      "update_probability to glue worm head and tail together and go back to diagonal configuration");

  
  
  LATTICE::Base::define_parameters(parameters);
  model::define_parameters(parameters);

  std::string ModelClassifierName = parameters["model"].as<std::string>();
  if ( ModelClassifierName == "BoseHubbard" ) {
    BoseHubbard::define_custom_model_parameters(parameters);
  }
  else if ( ModelClassifierName == "XXZ" ) {
    XXZ::define_custom_model_parameters(parameters);
  }
  else {
    throw std::runtime_error("Invalid Model \"" + ModelClassifierName + "\"");
  }

}


// Creates a new simulation.
// We always need the parameters and the seed as we need to pass it to
// the alps::mcbase constructor. We also initialize our internal state,
// mainly using values from the parameters.
worm::worm(parameters_type const & parameters, std::size_t seed_offset) : alps::mcbase(parameters, seed_offset)
    , sweeps(0)
    , thermalization_sweeps(int(parameters["thermalization"]))
    , total_sweeps(parameters["sweeps"])
    , runtimelimit(parameters["runtimelimit"])
    , reset_statistics(parameters["reset_statistics"])
    , beta(parameters["beta"])
    , E_off(parameters["E_off"])
    , C_worm(parameters["C_worm"])
    , canonical(parameters["canonical"])
    , Ntest(parameters["Ntest"])
    , Nsave(parameters["Nsave"])
    , Nmeasure(parameters["Nmeasure"])
    , Nmeasure2(parameters["Nmeasure2"])
    , update_statistics(boost::extents[statistics_tag::statistics_count][update_tag::update_count])
    , nbtime(zcmax+1)
    , dtime(zcmax+1)
    , nbval(zcmax)
    , MyGenerator(std::size_t(parameters["seed"]) + seed_offset)
{

  //set dtol
  int beta_int = int(beta);
  int bin_digits = 0;
  for (; beta_int > 0; beta_int >>= 1) bin_digits++;
  dtol = pow(2,-DBL_MANT_DIG+bin_digits);

  initialize_update_prob();

  {
  size_t Lx = parameters["Lx"].as<size_t>();
  size_t Ly = parameters["Ly"].as<size_t>();
  size_t Lz = parameters["Lz"].as<size_t>();
  bool pbcx = parameters["pbcx"].as<bool>();
  bool pbcy = parameters["pbcy"].as<bool>();
  bool pbcz = parameters["pbcz"].as<bool>();
  MyLatt = unique_ptr<LATTICE>(new LATTICE(Lx, Ly, Lz, pbcx, pbcy, pbcz));
#ifdef DEBUGMODE
  MyLatt->print();
#endif
  }

  winding_element = MyLatt->directions();
  Nsites = MyLatt->get_Nsites();
  zc.resize(Nsites);
  nb.resize(boost::extents[Nsites][zcmax]);
  bond_index.resize(boost::extents[Nsites][zcmax]);
  opposite_direction.resize(boost::extents[Nsites][zcmax]);

  for (SiteIndex s = 0; s < Nsites; s++) {
    for (size_t d = 0; d < zcmax; ++d) {
      if (MyLatt->has_neighbor(s,d)) {
        zc[s].push_back(d);
        nb[s][d] = MyLatt->nb(s,d);
        opposite_direction[s][d] = MyLatt->get_opposite_direction(s, nb[s][d]);
        bond_index[s][d] = MyLatt->get_bond_index(s,d);
      }
      else { //-1 means non existent
        nb[s][d] = -1;
        opposite_direction[s][d] = -1;
        bond_index[s][d] = -1;
      }
    }
  }

  ModelClassifierName = parameters["model"].as<std::string>();
  if (ModelClassifierName == "BoseHubbard") {
    MyModel = unique_ptr<BoseHubbard>(new BoseHubbard(parameters));
    hist_dm_fac = 1.;
  }
  else if (ModelClassifierName == "XXZ") {
    MyModel = unique_ptr<XXZ>(new XXZ(parameters));
    hist_dm_fac = 2.;
  }
  else {
    throw std::runtime_error("Invalid Model \"" + ModelClassifierName + "\"");
  }


#ifdef CAN_WINDOW
  if (canonical < 0)
    throw std::runtime_error("canonical value cannot be negative when CAN_WINDOW flag is set.");
#endif
  if (Nsites * MyModel->get_nmax() < canonical)
    throw std::runtime_error("canonical value should not be larger than Nsites * nmax.");


  measurements
    << alps::accumulators::LogBinningAccumulator<double >("Total_Energy")
    << alps::accumulators::LogBinningAccumulator<double >("Kinetic_Energy")
    << alps::accumulators::LogBinningAccumulator<double >("Potential_Energy")
    << alps::accumulators::LogBinningAccumulator<double >("Number_of_particles")
    << alps::accumulators::LogBinningAccumulator<double >("Number_of_particles_squared")
    << alps::accumulators::LogBinningAccumulator<vector<double> >("Density_Distribution")
#ifdef UNISYS
    << alps::accumulators::LogBinningAccumulator<vector<double> >("Density_Matrix")
    //<< alps::accumulators::LogBinningAccumulator<vector<double> >("Density_Matrix2")
    << alps::accumulators::LogBinningAccumulator<vector<double> >("DensDens_CorrFun")
    << alps::accumulators::LogBinningAccumulator<vector<double> >("Winding_number_squared")
#endif
#ifdef CAN_WINDOW
    << alps::accumulators::LogBinningAccumulator<vector<double> >("Greenfun_p0_tau")
#endif
  ;
#ifdef UNISYS
  hist_densmat.resize(Nsites);
  //hist_dd.resize(Nsites);
  for (size_t i=0; i < hist_densmat.size(); i++) {
    hist_densmat[i] = 0;
  }
#endif

#ifdef CAN_WINDOW
  hist_gt.resize(Ntimes_gt);	
  for (size_t i=0; i < hist_gt.size(); i++) hist_gt[i] = 0;
#endif
  
  
  nmax = MyModel->get_nmax();
  nmin = MyModel->get_nmin();
#ifdef UNISYS
  site_weight_diag.resize(nmax - nmin + 1);
  bond_weight_diag.resize(boost::extents[nmax - nmin + 1][nmax - nmin + 1]);
  site_weight_offdiag.resize(boost::extents[nmax - nmin + 1][nmax - nmin + 1]);
  bond_weight_offdiag.resize(boost::extents[nmax - nmin + 1][nmax - nmin + 1][nmax - nmin + 1][nmax - nmin + 1]);
  for (StateType n = 0; n < nmax - nmin + 1; n++ ) {
    site_weight_diag[n] = MyModel->site_weight_diag(0,n + nmin);
    for (StateType m = 0; m < nmax - nmin + 1; m++ ) {
      bond_weight_diag[n][m] = MyModel->bond_weight_diag( 0,n + nmin, m + nmin);
      site_weight_offdiag[n][m] = MyModel->site_weight_offdiag(0, n+nmin, m+nmin);
      for (StateType k = 0; k < nmax - nmin + 1; k++ ) {
        for (StateType l = 0; l < nmax - nmin + 1; l++ ) {
          bond_weight_offdiag[n][m][k][l] = MyModel->bond_weight_offdiag(0, n+nmin, m+nmin,k+nmin,l+nmin);
        }
      }
    }
  }
#endif
  
  state_eval.resize(MyModel->get_nmax() - MyModel->get_nmin() + 1);
  for (size_t i=0; i < state_eval.size(); i++) {
    state_eval[i] = MyModel->stateval(i);
  }
  
  C_NBW = C_worm / (Nsites * beta);
  operator_string.resize(Nsites + 1);
  dummy_it.resize(Nsites + 1);

  state.resize(Nsites);
  av_dns.resize(Nsites);
  av_state.resize(Nsites);
  av_state_sq.resize(Nsites);
    
  MCdiag = 0;
 
  counter.resize(counter_tag::SAMPLE_NUM + 1);
  for (size_t j=0; j < counter.size(); j++) counter[j] = 0;
  
}



void worm::initialize_update_prob() {
  initialize_update_prob(insertworm, moveworm);
  initialize_update_prob(moveworm, update_count);
}


void worm::initialize_update_prob(update_tag begin, update_tag end) {
  double norm = 0.;
  for (size_t upd = begin; upd < end; ++upd) {
    if (parameters.defined("p_" + update_names[upd])) {
      
      update_prob[upd] = double(parameters["p_" + update_names[upd]]);
      if (update_prob[upd] < 0) {
        throw std::runtime_error("Negative update probability: "
                                         + update_names[upd]);
      }
    }
    else {
      update_prob[upd] = 0.;
    }
    norm += update_prob[upd];
  }
  if (abs(norm - 1.) > 1e-10) {
    std::cerr << "Update probabilities do not add up to one.\n"
              << "Renormalizing..." << std::endl;
    for (size_t upd = begin; upd < end; ++upd) {
      update_prob[upd] /= norm;
      parameters["p_" + update_names[upd]] = update_prob[upd];
    }
  }
  norm = 0.;
  for (size_t upd = begin; upd < end; ++upd) {
    norm += update_prob[upd];
    update_prob_cuml[upd] = norm;
  }
}



void worm::initialize() {
  worm_diag = true;
  worm_at_stop = 0;
  worm_passes_nb_kink = 0;
  worm_dtime = 0.;
  new_measurement = true;
  worm_meas_densmat = false;
  mWinding = LATTICE::Direction::Zero();
  nrvertex = 0;
  Epot_measure = 0;

#ifndef UNISYS
  // The user must set the parameters of the Hamiltonian
  // The user can change the code below (eg, read from a user specified file)
  std::size_t Nbonds = MyLatt->get_Nbonds();
  std::size_t Ns = static_cast<std::size_t>(Nsites);
  if (ModelClassifierName == "BoseHubbard") {
    double t_ = 1;
    double U_ = 4;
    double chempot_ = 0.7;
    double V_ = 0;
    std::vector<std::size_t> dims{Nbonds, Ns, Ns, Nbonds};
    std::vector<double> val;
    for (size_t i=0; i < dims[0]; i++) val.push_back(t_);
    for (size_t i=0; i < dims[1]; i++) val.push_back(U_);
    for (size_t i=0; i < dims[2]; i++) val.push_back(chempot_);
    for (size_t i=0; i < dims[3]; i++) val.push_back(V_);
    MyModel->user_set_params(val, dims);
  }
  else {
    double hmag = 0;
    double J_ = 1.;
    double Jz = 1.; 
    std::vector<std::size_t> dims{Nbonds, Nbonds, Ns};
    std::vector<double> val;
    for (size_t i=0; i < dims[0]; i++) val.push_back(J_);
    for (size_t i=0; i < dims[1]; i++) val.push_back(Jz);
    for (size_t i=0; i < dims[2]; i++) val.push_back(hmag);
    MyModel->user_set_params(val, dims);
  }
#endif


#ifdef CAN_WINDOW
  for (SiteType s = 0; s < Nsites; s++) {
      state[s] = canonical / Nsites;
      if (s < canonical % Nsites) ++state[s];
  }
  number_of_particles = canonical;
#else
  for (SiteType s = 0; s < Nsites; s++)
    state[s] = 1;
  number_of_particles = Nsites;
#endif
  
  // insert dummy elements on all sites, keeping the code simpler and making a line of where to measure diag properties
  Nprtcls = 0.;
  for (SiteIndex i = 0; i < Nsites; i++) {
    Element_t new_elem(state[i],state[i], i,  beta, 0);
    dummy_it[i] = operator_string[i].insert(operator_string[i].end(), new_elem);
    Nprtcls += dummy_it[i]->before() * beta;
  }
  //Fictitious site at the end (index -1) for empty associations
  //only used for non-uniform coordination number
  dummy_it[Nsites] = (operator_string[Nsites].insert(operator_string[Nsites].end(), 
                                Element(0, 0, -1, beta, 0)));
  //Assoc is set to dummy_it.end()-1 (i.e. the fictitious site) if neighbor along direction index j doesn't exist
  for (SiteIndex i = 0; i < Nsites; i++) {
    for (size_t j = 0; j < zcmax; j++) {
      if (MyLatt->has_neighbor(i,j))
        dummy_it[i]->set_assoc(j, dummy_it[nb[i][j]]);
      else
        dummy_it[i]->set_assoc(j, *(dummy_it.end() - 1));
    }
  }
  worm_head_it = dummy_it[0];
  worm_tail_it = dummy_it[0];

  //initialize Epot_measure
  for (SiteType j = 0; j < Nsites; j++) {
    int n = dummy_it[j]->before();
    Epot_measure +=  MyModel->site_weight_diag(j, n);
    for (size_t const &iz : zc[j])
      Epot_measure += 0.5 * MyModel->bond_weight_diag(bond_index[j][iz], n, dummy_it[nb[j][iz]]->before());
  }
  
  Epot_tot = calc_potential_energy_nb() + calc_potential_energy_loc();
  std::cout << "# Potential Energy tot : " << Epot_tot << endl;
  std::cout << "# Finished intializing.\n";

}


void worm::measure() {
  ++sweeps;
  if (sweeps < thermalization_sweeps) return;
  if (sweeps == thermalization_sweeps) {
#ifdef UNISYS
    for (size_t i=0; i < hist_densmat.size(); i++) hist_densmat[i] = 0;
#endif
#ifdef CAN_WINDOW
    for (size_t i=0; i < hist_gt.size(); i++) hist_gt[i] = 0;
#endif
    return;
  }

  counter[counter_tag::WRITE]++;
  counter[counter_tag::MEASURE]++;
  counter[counter_tag::MEASURE2]++;
  counter[counter_tag::TEST]++;

  if(counter[counter_tag::MEASURE] >= Nmeasure) {
    counter[counter_tag::MEASURE] = 0;
    double Ek = nrvertex/beta * (-1.);
    vector<double> wsq(latt_dim, 0);
    for (size_t i=0; i < latt_dim; i++) {
      wsq[i] = mWinding[i] * mWinding[i] / (MyLatt->get_Ls(i) * MyLatt->get_Ls(i) * 1.);
    }
    measurements["Kinetic_Energy"] << Ek;
    measurements["Potential_Energy"] << Epot_measure;
    measurements["Total_Energy"] << Epot_measure + Ek;
    measurements["Number_of_particles"] << number_of_particles * 1.;
    measurements["Number_of_particles_squared"] << number_of_particles * number_of_particles * 1.;
#ifdef UNISYS
    measurements["Winding_number_squared"] << wsq;
#endif
  }
    
  if(counter[counter_tag::MEASURE2] >= Nmeasure2) {
    counter[counter_tag::MEASURE2] = 0;
    measure_corrfun();
    Epot_measure = calc_potential_energy_measure();
  }

  if (counter[counter_tag::TEST] >= Ntest) {
    std::cout << "\n# Testing...";
    try {
      test_conf();
      std::cout << "\n# Testing OK...\n";
      counter[counter_tag::TEST] = 0;
    }
    catch (const char* e) {
      std::cerr << e << std::endl;
      exit(1);
    }
    counter[counter_tag::TEST] = 0;
  }

  if (counter[counter_tag::WRITE] >= Nsave) {
    alps::hdf5::archive ar(parameters["outputfile"].as<std::string>(), "w");
    ar["samples/operator_string_" + to_string(counter[counter_tag::SAMPLE_NUM]++)] << serialize_op_string();
    counter[counter_tag::WRITE] = 0;
  }
}

void worm::force_reset_statistics() {
  if (sweeps > thermalization_sweeps) {
    for (size_t i=0; i < counter.size(); i++) counter[i] = 0;
    reset(measurements["Total_Energy"]);
    reset(measurements["Kinetic_Energy"]);
    reset(measurements["Potential_Energy"]);
    reset(measurements["Number_of_particles"]);
    reset(measurements["Number_of_particles_squared"]);
    reset(measurements["Density_Distribution"]);
#ifdef UNISYS
    reset(measurements["Density_Matrix"]);
    //reset(measurements["Density_Matrix2"]);
    reset(measurements["DensDens_CorrFun"]);
    reset(measurements["Winding_number_squared"]);
    for (size_t i=0; i < hist_densmat.size(); i++) {
      hist_densmat[i] = 0;
    }
#endif
#ifdef CAN_WINDOW
    reset(measurements["Greenfun_p0_tau"]);
    for (size_t i=0; i < hist_gt.size(); i++) hist_gt[i] = 0;
#endif 
  }
  for (size_t i=0; i < statistics_tag::statistics_count; i++) {
    for (size_t j=0; j < update_tag::update_count; j++) update_statistics[i][j] = 0;
  }
  sweeps = thermalization_sweeps;
}

void worm::measure_corrfun() {
#ifdef UNISYS
  //Corrfun is only measured for periodic boundary conditions
  if (LATTICE::n_basis > 1) return;
  if (!parameters["pbcx"]) return;
  if (LATTICE::dim > 1 && !parameters["pbcy"]) return;
  if (LATTICE::dim > 2 && !parameters["pbcz"]) return;
  vector<double> hist_dd(hist_densmat.size(), 0.);
  vector<double> dens_distr(Nsites);
  for (auto i=0; i < Nsites; i++) {
    size_t dist = MyLatt->relNumbering(0,i);
    if (hist_dd[dist] == 0) {
        hist_dd[dist] = state_eval[dummy_it[0]->before()] * state_eval[dummy_it[i]->before()];
    }
   dens_distr[i] = state_eval[dummy_it[i]->before()];
  }

  vector<double> hist_dm(hist_densmat.size());
  //vector<double> hist_dm2(hist_densmat.size());
  //double norm = av_dens_Nmeasure2 * 1. / Nmeasure2 * Nmeasure * 1. / Nsites;
  for (size_t i=0; i < hist_dm.size(); i++) {
    hist_dm[i] = hist_densmat[i] * hist_dm_fac / Nmeasure2;
  }
 
  measurements["Density_Distribution"] << dens_distr;
  measurements["Density_Matrix"] << hist_dm;
  //measurements["Density_Matrix2"] << hist_dm2;
  measurements["DensDens_CorrFun"] << hist_dd;
  for (size_t i=0; i < hist_densmat.size(); i++) {
    hist_densmat[i] = 0;
  }
#endif
#ifdef CAN_WINDOW
  vector<double> hist_gft(hist_gt.size());
  double dtau = (2 * can_window * beta) / hist_gt.size();
  for (size_t i=0; i < hist_gft.size(); i++) {
    hist_gft[i] = hist_gt[i] * hist_dm_fac / Nmeasure2 / dtau;
    hist_gt[i] = 0;
  }
  measurements["Greenfun_p0_tau"] << hist_gft;
#endif
}


void worm::measure_density_matrix() {
#ifdef UNISYS
  if (LATTICE::n_basis > 1) return;
  if (!parameters["pbcx"]) return;
  if (LATTICE::dim > 1 && !parameters["pbcy"]) return;
  if (LATTICE::dim > 2 && !parameters["pbcz"]) return;
  size_t dd = MyLatt->relNumbering(worm_head_it->link(), worm_tail_it->link());
  if (dd) {
    hist_densmat[dd] += 1. / (C_worm);
  }
  else {
    if (ModelClassifierName == "XXZ") {
      hist_densmat[dd] += 0.5 / (C_worm);
    }
    else {
      if (worm_at_stop == +1) {
        int outer_dens = worm_tail_it->before();
        int inner_dens = worm_tail_it->after();
        hist_densmat[0] += (outer_dens > inner_dens ? 0.5/C_worm : outer_dens/( inner_dens * 2 * C_worm  ) );
      }
      else {
        int outer_dens = worm_tail_it->after();
        int inner_dens = worm_tail_it->before();
        hist_densmat[0] += (outer_dens > inner_dens ? 0.5/C_worm : outer_dens/( inner_dens * 2*C_worm  ) );
      }
    }
  }
#endif
}

#ifdef CAN_WINDOW
void worm::measure_Gpt() {
  double dt = Nprtcls/beta - canonical;
  if (( dt < -can_window) || (dt > can_window)) {
    cerr << "# dt is out of bounds in measure_gpt " << dt << "\t" << Nprtcls << "\n";
    exit(1);
  }
  size_t index = static_cast<size_t>((can_window + dt) / (can_window * 2) * hist_gt.size());
  if (index >= hist_gt.size()) std::runtime_error("index in measure_gpt is out of bounds ");
  hist_gt[index] += 1./(C_worm);
}
#endif


// Returns a number between 0.0 and 1.0 with the completion percentage
double worm::fraction_completed() const {
  if ( sweeps > thermalization_sweeps && (sweeps - thermalization_sweeps) % (total_sweeps/10) == 0) {
    std::cout << "# Fraction completed : " << (sweeps - thermalization_sweeps) << "/" << total_sweeps << std::endl;
  }
  return (sweeps < thermalization_sweeps ? 0. : ( sweeps - thermalization_sweeps ) / double(total_sweeps));
}



void worm::test_conf() {
  for (SiteIndex i = 0; i < Nsites; i++) {
    for (Diagram_type::iterator it = operator_string[i].begin(); it != operator_string[i].end(); ++it) {
      if (MyModel->range_fail(it->before() )) {
        cerr << "\n# TEST_CONF :  Error with density  (before) << " << i << "\t" << it->time() << "\t" << it->before() << "\t" << it->after() << "\n";
        throw exception();
      }
      if (MyModel->range_fail(it->after() )) {
        cerr << "\n# TEST_CONF : Error with density  (after) << " << i << "\t" << it->time() << "\t" << it->before() << "\t" << it->after() << "\n";
        throw exception();
      }
      Diagram_type::iterator itn = it;
      ++itn; if (itn == operator_string[i].end()) itn = operator_string[i].begin();
      if ( (it->after() == it->before()) && (it->color() != 0)) {
        cerr << "\n# TEST_CONF : error in configuration (before==after) on site  " << i << endl;
        throw exception();
      }
      if ( (it->color() == 0) && ((it->after() != it->before()))) {
        cerr << "\n# TEST_CONF : error in configuration (dummy) on site  " << i << endl;
        throw exception();
      }
      if (it->after() != itn->before() ) {
        cerr << "\n# TEST_CONF : error in configuration (diag) on site  " << i << endl;
        throw exception();
      }
      if (it->time() > itn->time() && !(it == dummy_it[i] ) ) {
        cerr << "\n# TEST_CONF : chronology broken on site  " << i << endl;
        throw exception();
      }
      if (it->color() == 1) {
        bool b = false;
        bool b2 = false;
        for (size_t const& j : zc[i]) {
          if ( b && abs(it->time() - it->get_assoc(j)->time() ) < 1e-16 && it->get_assoc(j)->color() == 1 ) b2 = true;
          if (abs(it->time() - it->get_assoc(j)->time() ) < 1e-16 && it->get_assoc(j)->color() == 1) b = true;
          if ( (it->link() == nb[i][j] ) && (it->get_assoc(j)->link() == i) && (it->get_assoc(j)->time() != it->time())) {
            cerr << "\n# TEST_CONF : site is linked and has different association time; site " << i <<  " link " << it->link() << " at time " << it->time() << " and " <<  it->get_assoc(j)->time() << endl;
            throw exception();
          }
        }
        if (b2) {
          cerr << "\n# TEST_CONF : error in configuration, found two links at same time on site " << i << " for time " << it->time() << "\t links : " << it->get_assoc(0)->time() << "\t" << it->get_assoc(1)->time();
          throw exception();
        }
        if (!b) {
          cerr << "\n# TEST_CONF : error in configuration, found no link at same time on site " << i << " for time " << it->time() << "\t links : " << it->get_assoc(0)->time() << "\t" << it->get_assoc(0)->color()  << "\t" <<  it->get_assoc(1)->time() << "\t" << it->get_assoc(1)->color() << endl;
          throw exception();
        }
      }
      if (it->color() != 0) {         // extensive check on associations
        for (size_t const& j : zc[i]) {
          
          if (it->get_assoc(j)->time() < it->time() ) {
            cerr << "\n# TEST_CONF : assoc has earlier time than current link for current link on site " << i << " at time " << it->time() << " and nb " <<j << endl;
            throw exception();
          }
          Diagram_type::iterator itp = it->get_assoc(j);
          SiteIndex adj_site = nb[i][j];
          if (itp == operator_string[adj_site].begin()) itp = operator_string[adj_site].end();
          --itp;
          if (!(itp == dummy_it[adj_site]) && (itp->time() > it->time())  ) {
            cerr << "\n# TEST_CONF : there is an earlier element which should have been associated to cursite " << i << " at time " << it->time() << " and nb dir  " <<j << " site " <<  adj_site << "\t times it->get_assoc(j) and its predecessor : " << it->get_assoc(j)->time() << "\t" << itp->time() << endl;
            throw exception();
          }
        }
      } // ... if it->color() != 0
    } // ... iterate over configuration list
  } // iterate over all sites
  if (!worm_diag) {
    if (worm_head_it->color() != -1 ) {
      cerr << "\n# TEST_CONF : worm_head_it points to wrong element : " << worm_head_it->color()   << endl;
      throw exception();
    }
    if (worm_tail_it->color() != -1 ) {
      cerr << "\n# TEST_CONF : worm_tail_it points to wrong element : " << worm_tail_it->color()   << endl;
      throw exception();
    }
  }
  if (worm_meas_densmat) {
    if (worm_tail_it->time() != worm_head_it->time() ) {
      cerr << "\n# TEST_CONF : worm_meas_densmat is set but worms are at unequal times : " << worm_tail_it->time()  << "\t" << worm_head_it->time() << endl;
      throw exception();
    }
  }
  if ( (worm_at_stop == +1)  && (!worm_passes_nb_kink) && (!worm_meas_densmat) ) {
    Diagram_type::iterator itp = worm_head_it;
    if (itp == operator_string[worm_head_it->link()].begin()) itp = operator_string[worm_head_it->link()].end();
    --itp;
    if ( ! ( (itp->time() == worm_head_it->time()) || ( abs(itp->time()-beta + worm_head_it->time()) < 1e-10 ) ) )  {
      cerr << "\n# TEST_CONF : worm_at_stop == +1 but times are unequal times : " << worm_head_it->time()  << "\t" << itp->time() << "\t site : " << worm_head_it->link() << endl;
      throw exception();
    }
  }
  if ( (worm_at_stop == -1) && (!worm_passes_nb_kink) && (!worm_meas_densmat) ) {
    Diagram_type::iterator itn = worm_head_it;
    ++itn;
    if (itn == operator_string[worm_head_it->link()].end()) itn = operator_string[worm_head_it->link()].begin();
    if (itn->time() != worm_head_it->time()) {
      cerr << "\n# TEST_CONF : worm_at_stop == -1 but times are unequal times : " << worm_head_it->time()  << "\t" << itn->time() << "\t site : " << worm_head_it->link() << endl;
      throw exception();
    }
  }
  if (worm_passes_nb_kink) {
    bool b = false;
    for (size_t const& j : zc[worm_head_it->link()]) {
      SiteIndex adj_site = nb[worm_head_it->link()][j];
      for (Diagram_type::iterator it = operator_string[adj_site].begin(); it != operator_string[adj_site].end(); ++it) if (it->time() == worm_head_it->time()) b = true;
    }
    if (!b) {
      cerr << "\n# TEST_CONF : worm_passes_nb_kink is set but none is found; worm head time and site : " << worm_head_it->time() << "\t" << worm_head_it->link() << endl;
      throw exception();
    }
  }
  double Ep = calc_potential_energy_loc() + calc_potential_energy_nb();
  if (is_not_close(Ep, Epot_tot, 1e-8)) {
    cerr << "# Potential total energies do not match " << Ep << "\t" << Epot_tot << "\n";
    throw exception();
  }
  else {
#ifdef DEBUGMODE
    cout << "# Potential energies OK : " << Ep << "\t" << Epot_tot << "\t nrvertex : " << nrvertex << "\t Nprtcls : " << Nprtcls / Nsites / beta << "\n";
#endif
  }
}
  


double worm::calc_potential_energy_loc() {
  double En = 0.;
  double dt = 0.;
  for (int isite = 0; isite < Nsites; isite++) {
    double t_old = 0.;
    Diagram_type::iterator it1 = operator_string[isite].begin();
    int n1 = it1->before();
    for (;;) {
      dt = it1->time() - t_old;
#ifdef UNISYS
      En += site_weight_diag[n1 - nmin] * dt;
#else
      En += MyModel->site_weight_diag(isite, n1) * dt;
#endif
      
      n1 = it1->after();
      t_old = it1->time();
      ++it1;
      if (it1 == operator_string[isite].end()) break;
    }
  }
  return (En / beta );
}

double worm::calc_potential_energy_nb() {
  double En = 0.;
  double dt = 0.;
  for (int isite = 0; isite < Nsites; isite++) {
    for (auto const& iz : zc[isite]) {
      Diagram_type::iterator it1 = operator_string[isite].begin();
      int n1 = it1->before();
      SiteType jsite = nb[isite][iz];
      Diagram_type::iterator it2 = operator_string[jsite].begin();
      int n2 = it2->before();
      double t_old = 0.;
      for (;;) {
        for (;;) {
          if (it2 == operator_string[jsite].end()) break;
          if (it2->time() > it1->time()) break;
          dt = it2->time() - t_old;
#ifdef UNISYS
          En += bond_weight_diag[ n1][ n2] * dt;
#else
          En += MyModel->bond_weight_diag(bond_index[isite][ iz], n1, n2) * dt;
#endif
          t_old = it2->time();
          n2 = it2->after();
          ++it2;
        }
        dt = it1->time() - t_old;
#ifdef UNISYS
        En += bond_weight_diag[n1][ n2] * dt;
#else
        En += MyModel->bond_weight_diag(bond_index[isite][ iz], n1, n2) * dt;
#endif
        n1 = it1->after();
        t_old = it1->time();
        ++it1;
        if (it1 == operator_string[isite].end()) break;
      }
    }
  }
  return (En / 2. / beta);
}


double worm::calc_potential_energy_measure(bool with_chem_pot)
{
  double En = 0.;
  for (SiteType j = 0; j < Nsites; j++) {
    int n = dummy_it[j]->before();
    
#ifdef UNISYS
    En += site_weight_diag[n];
    for (size_t const& iz : zc[j]) {
      En += 0.5 * bond_weight_diag[n][dummy_it[nb[j][iz]]->before()];
    }
#else
    En +=  MyModel->site_weight_diag(j, n, with_chem_pot) ;
    for (size_t const &iz : zc[j]) {
      En += 0.5 * MyModel->bond_weight_diag(bond_index[j][iz], n, dummy_it[nb[j][iz]]->before());
    }
#endif
  }
  return (En);
}



void worm::find_assoc_insert(const SiteIndex cursite, Diagram_type::iterator  it, const int shift) {
  // it is the iterator to the newly insterted element
  // memory for its associations has been allocated already but the associations must be set correctly here
  // and we need to check on the neighbors if their associations need to be changed
  
  // part 1 : iterators for the new element located on site cursite
  //std::cout<< "# Welcome to find_assoc_insert " << cursite << "\t" << *it << "\n";
  // we go one element up where we by assumption have a properly associated element
  Diagram_type::iterator ito, itl, itp, itw;
  ito= it;
  ++ito;
  if (ito == operator_string[cursite].end()) ito = operator_string[cursite].begin();
  
  double t0 = it->time();
  it->time( t0 + shift * dtol*2);
 
  for (size_t j = 0; j < zcmax; ++j) {
    SiteIndex s = nb[cursite][j];
    if (s == -1) {
        it->set_assoc(j, *(dummy_it.end() - 1));
        continue;
    }
    itl = ito->get_assoc(j);
    itp = itl;
    if (itp == operator_string[s].begin()) itp = operator_string[s].end();
    --itp;
    it->set_assoc(j, itl);
    while (!t_between(it->time(), itp->time(), itl->time())) {
      itl = itp;
      if (itp == operator_string[s].begin()) itp = operator_string[s].end();
      --itp;
      it->set_assoc(j, itl);
    }
  }
  
  
  // part 2 : iterators on the neighbors might have to be set newly in the "opposite" direction
  for (size_t const& j : zc[cursite]) {
    SiteIndex s = nb[cursite][j];
    size_t oppdir = opposite_direction[cursite][j];
    itl = it->get_assoc(j);
    itp = itl;
    if (itp == operator_string[s].begin()) itp = operator_string[s].end();
    --itp;
    itw = itp->get_assoc(oppdir);
    if (itl->time() == it->time()) {
      itl->set_assoc(oppdir, it);
    }
    while ((itp->time() != itw->time())  && ( t_between(it->time(), itp->time(), itw->time()) )) {
      itp->set_assoc(oppdir, it);
      if (itp == operator_string[s].begin()) itp = operator_string[s].end();
      --itp;
      itw = itp->get_assoc(oppdir);
    }
  }
  
  it->time(t0);
}




void worm::find_assoc_delete(const SiteIndex cursite, Diagram_type::iterator  it) {
  if (operator_string[cursite].size() == 1) {
    std::cerr <<"# How to erase a list of length one when there must be a dummy?\n";
    return;
  }
 
  Diagram_type::iterator itn=it;
  ++itn;
  if (itn == operator_string[cursite].end()) itn = operator_string[cursite].begin();
  Diagram_type::iterator itp = it;
  if (itp == operator_string[cursite].begin()) itp=operator_string[cursite].end();
  --itp;
  // on nb sites, go down until you find an element that points at the current element
  Diagram_type::iterator it_end, it_begin, itt;
  for (size_t const& j : zc[cursite]) {
    //SiteIndex s = nbsite[j];
    SiteIndex s = nb[cursite][j];
    size_t oppdir = opposite_direction[cursite][j];
    it_end=it->get_assoc(j);
    it_begin=itp->get_assoc(j);
    itt=it_begin;
    if ((itt == it_end) && (itt == dummy_it[s])) {
      for (Diagram_type::iterator itt=operator_string[s].begin(); itt != operator_string[s].end(); ++itt) {
        if (itt->get_assoc(oppdir) == it) itt->set_assoc(oppdir, itn);
      }
    }
    else {
      while (itt != it_end) {
        if (itt->get_assoc(oppdir) == it) itt->set_assoc(oppdir, itn);
        ++itt;
        if (itt==operator_string[s].end()) itt=operator_string[s].begin();
      }
      if (it_end->get_assoc(oppdir) == it) it_end->set_assoc(oppdir, itn);
    } // else
  } //for
}
  

