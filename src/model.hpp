#pragma once

#include "lattice.hpp"
#include <math.h>
#include <vector>
#include <string>
#include <algorithm>
#include <iterator>

class model {
public:
  using StateType = int;
  using SiteIndex = typename LATTICE::Base::SiteIndex;
  using BondIndex = typename LATTICE::Base::BondIndex;

  model(alps::params const& params) {
    model_name = params["model"].as<std::string>();
  };

  //first argument only needed in case of non-uniform model, e.g. random on-site potential
  virtual double site_weight_diag(const SiteIndex, const StateType, const bool b = true) const = 0;
  virtual double site_weight_offdiag(const SiteIndex, const StateType, const StateType) const = 0;

  //first argument only needed in case of anisotropic model, e.g. kitaev model
  virtual double bond_weight_diag(const BondIndex, const StateType, const StateType) const = 0;
  virtual double bond_weight_offdiag(const BondIndex, const StateType, const StateType, const StateType, const StateType) const = 0;

#ifndef UNISYS
  virtual void user_set_params(const std::vector<double>&, const std::vector<std::size_t>& ) = 0; 
#endif
  
  virtual double stateval(const StateType n) const = 0;

  virtual void print_params(std::ostream&) const = 0;

  StateType get_nmax() const { return nmax;}
  StateType get_nmin() const { return nmin;}
  
  static void define_parameters(alps::params& params) {
    params.define<std::string>("model", "BoseHubbard", "name of the Hamiltonian. Choose either BoseHubbard or XXZ");
  }

  static bool parameters_supplied(alps::params const& params) {
    return params.supplied("model");
  }

  bool range_fail(const StateType n) const {
    return n < nmin || n > nmax;
  }
  std::string get_name() const { return model_name;}
  
protected:
  std::string model_name;
  StateType nmax;
  StateType nmin;
};


class BoseHubbard : public model {
public:
  BoseHubbard(alps::params const& params) 
    : model(params)
#ifdef UNISYS
    ,  U_on(params["U_on"])
    ,  V_nn(params["V_nn"])
    ,  mu(params["mu"])
    ,  t_hop(params["t_hop"])
#endif
  {
    nmin = 0;
    nmax = params["nmax"].as<StateType>();
  }

  static void define_custom_model_parameters(alps::params & params) {
    params
#ifdef UNISYS
      .define<double>("t_hop",                1.0,    "hopping amplitude in the Bose-Hubbard model t b^{\\dagger}_i b_j")
      .define<double>("U_on",                 4.0,    "on-site density-density repulsion in Bose-Hubbard model U n_i(n_i-1)")
      .define<double>("V_nn",                 0.0,    "nearest neighbor density density repulsion in Bose-Hubbard model V n_i n_j")
      .define<double>("mu",                   0.7,    "chemical potential - mu n_i")
#endif
      .define<int>("nmax",                    10,     "maximum occupation number");
  }
  
  double stateval(const StateType n) const override {
    return n*1.;
  }

  double site_weight_offdiag(const SiteIndex s, const StateType n1, const StateType n2) const override {
    if (range_fail(n1) || range_fail(n2)) return 0;
    if (n1 == n2 + 1) {
      return sqrt(n1 * 1.);
    }
    else if (n2 == n1 + 1) {
      return sqrt(n2 * 1.);
    }
    return 0;
  }

  double site_weight_diag(const SiteIndex s, const StateType n, const bool b = true) const override {
#ifdef UNISYS
    double w = U_on * n * (n-1.)/2. ;
    if (b) w += -mu * n;
    return w;
#else
    double w = U_on[s] * n * (n-1.)/2. ;
    if (b) w += -mu[s] * n;
    return w;
#endif
  }

  double bond_weight_diag(const BondIndex b, const StateType n, const StateType m) const override {
#ifdef UNISYS
    return V_nn * n * m;
#else
    return V_nn[b] * n * m;
#endif
  }

  double bond_weight_offdiag(const BondIndex b, const StateType n1, const StateType n2, 
                             const StateType m1, const StateType m2) const override {
    if (range_fail(n1) || range_fail(n2) || range_fail(m1) || range_fail(m2))  return 0;
    double w = 0.;
    if ((n1 == n2 + 1) && (m2 == m1 + 1)) {
      w = sqrt(n1 * m2 * 1.);
    }
    else if  ((n2 == n1 + 1) && (m1 == m2 + 1)) {
      w = sqrt(n2 * m1 * 1.);
    }
#ifdef UNISYS
    return (t_hop * w);
#else
    return (t_hop[b] * w);
#endif
  }


#ifndef UNISYS
  void user_set_params(const std::vector<double>& p, const std::vector<std::size_t>& dims) override {
    std::copy(p.begin(), p.begin() + dims[0], std::back_inserter(t_hop));
    std::copy(p.begin() + dims[0], p.begin() + dims[0] + dims[1], std::back_inserter(U_on));
    std::copy(p.begin() + dims[0] + dims[1], p.begin() + dims[0] + dims[1] + dims[2], std::back_inserter(mu));
    std::copy(p.begin() + dims[0] + dims[1] + dims[2], p.end(), std::back_inserter(V_nn));
  }
#endif
  
  
  void print_params(std::ostream& os) const override {
#ifdef UNISYS
    os << "# model : name                                      : " << model_name << "\n";
    os << "# model : nmin, nmax                                : " << nmin << "\t" << nmax << "\n";
    os << "# model : U_on                                      : " << U_on    << "\n";
    os << "# model : V_nn                                      : " << V_nn    << "\n";
    os << "# model : mu                                        : " << mu      << "\n";
    os << "# model : hopping amplitude                         : " << t_hop   << "\n";
#else
    os << "# model : name                                      : " << model_name << "\n";
    os << "# model : nmin, nmax                                : " << nmin << "\t" << nmax << "\n";
    os << "# User specified inhomogeneous parameters:\n";
    os << "# model : U_on                                      : "; for (size_t i=0; i < U_on.size(); ++i) os  << U_on[i] << " "; os << "\n";
    os << "# model : V_nn                                      : "; for (size_t i=0; i < V_nn.size(); ++i) os  << V_nn[i] << " "; os << "\n";
    os << "# model : mu                                        : "; for (size_t i=0; i < mu.size()  ; ++i) os  << mu[i]   << " "; os << "\n";
    os << "# model : t_hop                                     : "; for (size_t i=0; i < t_hop.size(); ++i) os << t_hop[i] << " "; os << "\n";
#endif
  }

private:
#ifdef UNISYS
  double U_on;
  double V_nn;
  double mu;
  double t_hop;
#else
  std::vector<double> U_on;
  std::vector<double> V_nn;
  std::vector<double> mu;
  std::vector<double> t_hop;
#endif
};


//the Hamiltonian is H_XXZ = 
//         + \sum_{\langle i, j\rangle}( Jzz*(S_i^z*S_j^z) - 1/2 * Jpm*(S_i^{+}S_j^{-} + S_i^{-}S_j^{+}) ) + 
//         - \sum_i h * S_i^z
class XXZ : public model {
public:
  XXZ(alps::params const& params) 
    : model(params)
#ifdef UNISYS
    ,  Jpm(params["Jpm"])
    ,  Jzz(params["Jzz"])
    ,  h(params["h"])
#endif
  {
    nmin = 0;
    nmax = params["nmax"].as<StateType>();
    Sspin = nmax/2.;
  }

  static void define_custom_model_parameters(alps::params & params) {
    params
#ifdef UNISYS
      .define<double>("Jpm",            1.0,    "XX coupling")
      .define<double>("Jzz",            1.0,    "Z coupling")
      .define<double>("h",              0.2,    "ext. field coupling")
#endif
      .define<int>("nmax",                1,    "2S (two times spin)");
  }
  
  double site_weight_offdiag(const SiteIndex s, const StateType n1, const StateType n2) const override {
    if (range_fail(n1) || range_fail(n2)) return 0;
    if (n2 == n1 + 1) return sqrt( Sspin * (Sspin + 1.) - stateval(n1) * stateval(n2));
    if (n1 == n2 + 1) return sqrt( Sspin * (Sspin + 1.) - stateval(n1) * stateval(n2));
    return 0;
  }
  
  double stateval(const StateType n) const override{
    return (n - Sspin);
  }

  double site_weight_diag(const SiteIndex s, const StateType n, const bool b = true) const override {
#ifdef UNISYS
     return b ? - h * stateval(n) : 0;
#else
     return b ? - h[s] * stateval(n) : 0;
#endif
  }

  double bond_weight_diag(const BondIndex b, const StateType n, const StateType m) const override {
#ifdef UNISYS
    return Jzz * stateval(n)  * stateval(m);
#else
    return Jzz[b] * stateval(n)  * stateval(m);
#endif
  }

  double bond_weight_offdiag(const BondIndex b, const StateType n1, const StateType n2, 
                             const StateType m1, const StateType m2) const override {
    if (range_fail(n1) || range_fail(n2) || range_fail(m1) || range_fail(m2))  return 0;
#ifdef UNISYS
    if (n2 == n1 + 1 && m1 == m2 + 1) {
        return 0.5 * Jpm * sqrt(Sspin * (Sspin+1.) - stateval(n1)*stateval(n2)) * sqrt(Sspin * (Sspin+1.) - stateval(m1)*stateval(m2));
    }
    else if (n1 == n2 + 1 && m2 == m1 + 1) {
        return 0.5 * Jpm * sqrt(Sspin * (Sspin+1.) - stateval(n1)*stateval(n2)) * sqrt(Sspin * (Sspin+1.) - stateval(m1)*stateval(m2));
    }
#else
    if (n2 == n1 + 1 && m1 == m2 + 1) {
        return 0.5 * Jpm[b] * sqrt(Sspin * (Sspin+1.) - stateval(n1)*stateval(n2)) * sqrt(Sspin * (Sspin+1.) - stateval(m1)*stateval(m2));
    }
    else if (n1 == n2 + 1 && m2 == m1 + 1) {
        return 0.5 * Jpm[b] * sqrt(Sspin * (Sspin+1.) - stateval(n1)*stateval(n2)) * sqrt(Sspin * (Sspin+1.) - stateval(m1)*stateval(m2));
    }
#endif
    return 0;
  }

#ifndef UNISYS  
  void user_set_params(const std::vector<double>& p ,const std::vector<std::size_t>& dims) override {
    std::copy(p.begin(), p.begin() + dims[0], std::back_inserter(Jpm));
    std::copy(p.begin()+dims[0], p.begin() + dims[0] + dims[1], std::back_inserter(Jzz));
    std::copy(p.begin() + dims[0] + dims[1], p.end(), std::back_inserter(h));
  }
#endif
  
  void print_params(std::ostream& os) const override {
    os << "# model : name                                 : " << model_name << "\n";
    os << "# model : nmin, nmax                           : " << nmin << "\t" << nmax << "\n";
#ifdef UNISYS
    os << "# model : Jpm                                  : " << Jpm    << "\n";
    os << "# model : Jzz                                  : " << Jzz    << "\n";
    os << "# model : h                                    : " << h    << "\n";
#else
    os << "# User specified inhomogeneous parameters:\n";
    os << "# model : Jpm                                      : "; for (size_t i=0; i < Jpm.size(); ++i) os  << Jpm[i] << " "; os << "\n";
    os << "# model : Jzz                                      : "; for (size_t i=0; i < Jzz.size(); ++i) os  << Jzz[i] << " "; os << "\n";
    os << "# model : h                                        : "; for (size_t i=0; i < h.size()  ; ++i) os  << h[i]   << " "; os << "\n";
#endif
  }

private:
  double Sspin;
#ifdef UNISYS
  double Jpm;
  double Jzz;
  double h;
#else
  std::vector<double> Jpm;
  std::vector<double> Jzz;
  std::vector<double> h;
#endif
};

