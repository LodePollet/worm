/* $Id: worm.hpp,v 1.1 2006/09/09 8:59:59 pollet Exp $ */

#ifndef worm_HPP
#define worm_HPP

#pragma once


#include "lattice.hpp"
#include "worm.Element.hpp"
#include "model.hpp"
#include <alps/mc/mcbase.hpp>
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/hdf5/multi_array.hpp>
#include <alps/accumulators.hpp>
#include <alps/mc/api.hpp>
#include <alps/mc/stop_callback.hpp>

#include <random>
#include <vector>
#include <list>
#include <numeric>
#include <time.h>
#include <iomanip>
#include <exception>
#include <cstdlib>
#include <fstream>


typedef boost::multi_array<double, 2> DMatrix;
typedef boost::multi_array<double, 2>::index DMatrix_index;

typedef boost::multi_array<double, 4> DTensor;
typedef boost::multi_array<double, 4>::index DTensor_index;


typedef boost::multi_array<int, 2> IMatrix;
typedef boost::multi_array<int, 2>::index IMatrix_index;



typedef Element Element_t;

enum statistics_tag {
    impossible,
    rejected,
    accepted,
    total_attempted,
    statistics_count
};
namespace {
    const std::string statistics_names[] = {
        "impossible",
        "rejected",
        "accepted",
        "total_attempted"
    };
    const std::string update_names[] = {
        "insertworm",
        "moveworm",
        "insertkink",
        "deletekink",
        "glueworm",
    };
}


using namespace std;


class worm : public alps::mcbase {
public :
    static constexpr double tol = 1e-14;
    using StateType = model::StateType;
    using SiteType = LATTICE::Base::SiteType;
    using SiteIndex = LATTICE::Base::SiteIndex;
    using BondType = LATTICE::Base::BondType;
    using BondIndex = LATTICE::Base::BondIndex;
    static const size_t zcmax = LATTICE::zcmax;
    static const size_t latt_dim = LATTICE::Base::dim;
  
  worm(parameters_type const & parms, std::size_t seed_offset = 0);
    
  static std::string code_name();
    
  // ALPS overloads:
  static void define_parameters(parameters_type & parameters);

  virtual void update();
  virtual void measure();
  void measure_corrfun();
  virtual double fraction_completed() const override;
  
  using alps::mcbase::save;
  virtual void save(alps::hdf5::archive & ar) const;

  using alps::mcbase::load;
  virtual void load(alps::hdf5::archive & ar);
  
  bool is_not_close(const double val1, const double val2, const double tol=1e-10) const {
    bool answ = false;
    if (fabs(val1) > 1.) {
      if (fabs(val2) < 0.001) return true;
      if (fabs(1. - val1/val2) > tol) return true;
    }
    else {
      if (fabs(val1 - val2) > tol) return true;
    }
    return answ;
  }

  
  bool Metropolis(const double& x) {
 //#ifdef DEBUGMODE
     if (x < 0 || std::isnan(x) ) {
       cerr << "x < 0 in Metropolis? " << x << "\n";
       char ch; cin >> ch;
     }
 //#endif
     if (x >= 1) return (true);
     if (rnd(MyGenerator) < x) return (true);
     //if (random() < x) return (true);
     return (false);
   }

  double Diag_energy(const SiteIndex s, const StateType n, const size_t dir, const size_t m) {
#ifdef UNISYS
    double d = site_weight_diag[n-nmin];
    if (nb[s][dir] != -1)
        d  += bond_weight_diag[n-nmin][m-nmin];
    return d;
#else
    double d = MyModel->site_weight_diag(s,n);
    if (nb[s][dir] != -1)
        d += MyModel->bond_weight_diag( bond_index[s][dir],n, m);
    return d;
#endif
  }

  double Diag_energy(const SiteIndex s, const StateType n, const vector<StateType>& snb) {
#ifdef UNISYS
    double d = site_weight_diag[n - nmin];
    for (size_t m=0 ; m < snb.size(); m++) {
        if (nb[s][m] != -1)
            d += bond_weight_diag[n-nmin][snb[m]-nmin];
    }
    return d;
#else
    double d = MyModel->site_weight_diag(s,n);
    for (size_t i=0 ; i < snb.size(); i++) {
        if (nb[s][i] != -1)
            d+= MyModel->bond_weight_diag( bond_index[s][i],n, snb[i]);
    }
    return d;
#endif
  }
  
 
  //void reset_measurements();
  void force_reset_statistics();
  
  double calc_potential_energy_measure(bool with_chem_pot=true);
  double calc_potential_energy_nb();
  double calc_potential_energy_loc();
  
  int INSERTWORM();
  int MOVEWORM();
  int INSERTKINK();
  int DELETEKINK();
  int GLUEWORM();
  int PASSINTERACTION(const int, Diagram_type::iterator);
  void PASS_DUMMY(const int, const SiteType);
  int UPDATE_DENSITYMATRIX();
  
  void output_final(const alps::results_type<worm>::type& results, alps::hdf5::archive& ar);
  void initialize();
  void test_conf();
  void print_update_statistics(std::ostream&) const;
  void print_params(std::ostream& os) const;
  void print_conf(ostream&) const;
  boost::multi_array<double,2> serialize_op_string() const;

  void measure_density_matrix();
  void update_hist();
#ifdef CAN_WINDOW
  void measure_Gpt();
#endif
  
protected:
  unsigned long sweeps;
  unsigned long thermalization_sweeps;
  unsigned long total_sweeps;
  size_t Nloop;
  size_t runtimelimit;
  int reset_statistics;
  
  int MCdiag;
  
  enum update_tag {insertworm, moveworm, insertkink, deletekink, glueworm, update_count};
  enum counter_tag {WRITE, MEASURE, MEASURE2, TEST, SAMPLE_NUM};
  double update_prob_cuml[update_count];
  double update_prob[update_count];
private:
  void initialize_update_prob(update_tag begin, update_tag end);
  void initialize_update_prob();

    
private :
  double dtol;

  std::unique_ptr<LATTICE> MyLatt;
  std::unique_ptr<model> MyModel;
  std::string ModelClassifierName;
  
  double beta;            // inverse temperature
  double E_off;           // energy offset, internal parameter of the code
  double C_worm, C_NBW;
  int canonical;          // canonical number of particles
  
  unsigned long Ntest;
  unsigned long Nsave;
  size_t Nmeasure, Nmeasure2;
  
  // internal variables
  unsigned long long nrvertex;
  long long number_of_particles;
  double Epot_measure, Epot_tot, Nprtcls;
  bool passed_dummy;
  std::vector<int> state; // occupation per site
  std::vector<Diagram_type> operator_string;
 
  // worm variables
  Diagram_type::iterator worm_tail_it;  // iterator pointing to the tail
  Diagram_type::iterator worm_head_it;  // iterator pointing to the tail
  vector<Diagram_type::iterator > dummy_it; // iterator pointing at dummy elements useful for measuring diag prop

  bool worm_diag;
  int worm_at_stop;
  int worm_passes_nb_kink;
  bool worm_meas_densmat;
  double worm_dtime;
  bool new_measurement;
  
  void find_assoc_insert(const SiteIndex, Diagram_type::iterator, const int );
  void find_assoc_delete(const SiteIndex, Diagram_type::iterator );
  bool shift_right_kink(double&, double, double&);
  bool shift_left_kink(double&, double, double&);
  bool t_between(const double, const double, const double) const;
  void worm_left();
  void worm_right();
  void worm_pass_next();
  void worm_pass_prev();
  void worm_changedir_toleft();
  void worm_changedir_toright();
  void get_diaginfo_nb_pos(const SiteIndex s);
  void get_diaginfo_nb_neg(const SiteIndex s);
  
  // measurement variables
  double mZ_dns, mZ_state;
  double mZ, mG;
  double mZ_green;
  vector<double> av_dns;
  vector<double> av_state;
  vector<double> av_state_sq;
  vector<double> av_green;
  vector<double> time_green;

  LATTICE::Direction mWinding;
  std::array<LATTICE::Direction, LATTICE::zcmax> winding_element;
  
  DMatrix update_statistics;
  vector<double> nbtime;
  vector<double> dtime;
  vector<StateType> nbval;
  
  // lattice variables
  vector<vector<size_t>> zc;
  SiteIndex Nsites;
  IMatrix nb;
  IMatrix bond_index;
  IMatrix opposite_direction;

  StateType nmin,nmax;
  std::vector<double> state_eval;
#ifdef UNISYS
  std::vector<double> site_weight_diag;
  DMatrix bond_weight_diag;
  DMatrix site_weight_offdiag;
  DTensor bond_weight_offdiag;
  bool range_check(const StateType n) {
    return ( (n >= nmin) && (n <= nmax));
  }
#endif

  double hist_dm_fac;
#ifdef UNISYS
  std::vector<double> hist_densmat;
#endif
  std::vector<size_t> counter;
#ifdef CAN_WINDOW
  static constexpr size_t Ntimes_gt = 200;
  static constexpr double can_window = 0.1;  // fraction of the imaginary time above and below the worm tail in which the worm head can move
  std::vector<double> hist_gt;
#endif
  
  mt19937 MyGenerator;
  uniform_real_distribution<double> rnd;
  
};


inline void worm::get_diaginfo_nb_pos(const SiteIndex s) {
  Diagram_type::iterator its, itp;
  for (size_t const& i : zc[s]) {
    SiteType snb = nb[s][i];
    its = worm_head_it->get_assoc(i);
    itp = its;
    if (itp == operator_string[snb].begin()) itp=operator_string[snb].end();
    --itp;
    if ( worm_head_it->time() == its->time() ) {
      ++its;
      if (its == operator_string[snb].end() ) its = operator_string[snb].begin();
      dtime[i] = its->time() - worm_head_it->time();
      if (dtime[i] <= 0.) dtime[i] += beta;
      nbtime[i] = its->time();
      nbval[i] = its->before();
    }
    else {
      nbtime[i] = its->time();
      dtime[i] = its->time() - worm_head_it->time();
      if (dtime[i] <= 0.) dtime[i] += beta;
      nbval[i] = its->before();
    }
  }
}

inline void worm::get_diaginfo_nb_neg(const SiteIndex s) {
  Diagram_type::iterator its, itp;
  for (size_t const& i : zc[s]) {
    SiteType  snb =  nb[s][i];
    its = worm_head_it->get_assoc(i);
    itp = its;
    if (itp == operator_string[snb].begin()) itp=operator_string[snb].end();
    --itp;
#ifdef DEBUGMODE
    cout << "\n" << i << "\t" << snb << "\t" << worm_head_it->time() << "\t" << itp->time() << "\t" << its->time() << "\n";
#endif
    nbtime[i] = itp->time();
    dtime[i] = worm_head_it->time() - itp->time();
    if (dtime[i] <= 0.) dtime[i] += beta;
    nbval[i] = itp->after();
#ifdef DEBUGMODE
    cout << i << "\n-\t" << dtime[i] << "\t" << its->time() << "\t" << itp->time()   << "\toccup : " << itp->after() << "\n";
#endif
  }
}


inline bool worm::t_between(const double t0, const double t1, const double t2) const
// check if t0 is between t1 and t2, given that t1 happens before t2 (mod beta)
{
  if (t2 >t1)
  {
    if ((t2>=t0) && (t0 > t1)) return (true);
  }
  else
  {
    if ((t0 > t1) || (t0 <= t2)) return (true);
  }
  return (false);
}

#endif
