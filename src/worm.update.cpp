#include "worm.hpp"


void worm::update() {
  
  size_t Nupd = 0;
  do {
    double q = rnd(MyGenerator);
    int a;
    if (worm_diag) {
      if (q < update_prob_cuml[insertworm]) {
        a = INSERTWORM();
        update_statistics[statistics_tag::total_attempted][update_tag::insertworm] += 1.;
        update_statistics[a][update_tag::insertworm] += 1.;
      }
    }
    else {
      if (q < update_prob_cuml[moveworm]) {
        a = MOVEWORM();
        update_statistics[statistics_tag::total_attempted][update_tag::moveworm] += 1.;
        update_statistics[a][update_tag::moveworm] += 1.;
      }
      else if (q < update_prob_cuml[insertkink]) {
        a = INSERTKINK();
        update_statistics[statistics_tag::total_attempted][update_tag::insertkink] += 1.;
        update_statistics[a][update_tag::insertkink] += 1.;
      }
      else if (q < update_prob_cuml[deletekink]) {
        a = DELETEKINK();
        update_statistics[statistics_tag::total_attempted][update_tag::deletekink] += 1.;
        update_statistics[a][update_tag::deletekink] += 1.;
      }
      else if (q < update_prob_cuml[glueworm]) {
        a = GLUEWORM();
        update_statistics[statistics_tag::total_attempted][update_tag::glueworm] += 1.;
        update_statistics[a][update_tag::glueworm] += 1.;
      }
    }
#ifdef DEBUGMODE
    print_conf(std::cout);
    try {
      test_conf();
    }
    catch (const exception& e) {
      cerr << e.what() << endl;
      print_conf(std::cout);
      exit(1);
    }
#endif
    
    if (worm_diag) {
      mZ += 1.;
      MCdiag += 1;
    }
    else {
      mG += 1.;
#ifdef UNISYS
      if (worm_meas_densmat) measure_density_matrix();
#endif
#ifdef CAN_WINDOW
   if (worm_at_stop == 0) measure_Gpt();
#endif
    }
    Nupd++;
    if (Nupd == 10000) {
      Nupd = 0;
    }
    
  } while (!worm_diag);

}


int worm::INSERTWORM() {
#ifdef DEBUGMODE
  if (!worm_diag) {
    cerr << "\nNon diagonal configuration in INSERTWORM ? \n" ; 
  }
  std::cout << "# Welcome to INSERTWORM\n";
  print_conf(std::cout);
#endif

  worm_dtime = 0;
  mZ_dns += 1.;
  SiteIndex isite=static_cast<SiteIndex>(rnd(MyGenerator)*Nsites);
  double start_time=rnd(MyGenerator)*beta;
#ifdef DEBUGMODE
  cout << "INSERTWORM : isite " << isite << "\t start_time "<< start_time << "\n";
#endif
  // find suitable place to make a worm pair
  Diagram_type::iterator cycle_it = dummy_it[isite];
  Diagram_type::iterator prev_it = cycle_it;
  if (operator_string[isite].size() != 1) // dummy element always there
  {
    
    if (prev_it ==  operator_string[isite].begin()) prev_it = operator_string[isite].end();
    --prev_it;

    // you cannot  create a worm at the same place where the previous one was removed...
    // you really need to go around in space-time. maybe for insulating phases 
    // the next couple of lines are time consuming
    while (!t_between(start_time, prev_it->time(), cycle_it->time()))
    {
      cycle_it = prev_it;
      if (prev_it == operator_string[isite].begin()) prev_it = operator_string[isite].end();
      --prev_it;
    }
  }

  
  int worm_dir = ( (rnd(MyGenerator) < 0.5) ? 1 : -1 );
  StateType n_out = cycle_it->before();
  StateType n_mid = n_out + ((rnd(MyGenerator) < 0.5) ? 1 : -1);
#ifdef UNISYS
  if (!range_check(n_mid)) return impossible;
#else
  if ( MyModel->range_fail(n_mid) ) return impossible;
#endif
  
#ifdef UNISYS
  double wgt = site_weight_offdiag[n_out][ n_mid];
#else
  double wgt = MyModel->site_weight_offdiag(isite, n_out, n_mid);
#endif
  
  double ratio =  2* C_NBW * Nsites * beta * wgt * wgt * update_prob[glueworm] / update_prob[insertworm];

  
#ifdef DEBUGMODE
  cout << "INSERTWORM ratio " <<  ratio << "\t" <<  2* C_NBW * Nsites * beta << "\t" <<  wgt << "\t" <<  update_prob[glueworm] / update_prob[insertworm] << "\n";
  //char ch; cin >> ch;
#endif
  if (Metropolis(ratio)) {
#ifdef DEBUGMODE
      cout << "INSERTWORM accepted " << cycle_it->time() << "\t" << passed_dummy << "\t ratio : " << ratio << "\n";
#endif
    Element_t lower_kink(n_out, n_mid, isite, start_time, -1);

    Element_t upper_kink(n_mid, n_out, isite, start_time, -1);
    worm_head_it = operator_string[isite].insert(cycle_it, upper_kink);
    find_assoc_insert(isite, worm_head_it,0);
    worm_tail_it = operator_string[isite].insert(worm_head_it, lower_kink);
    find_assoc_insert(isite, worm_tail_it, -1);
    
    if (worm_dir == 1) {
      worm_at_stop = 1;
      worm_passes_nb_kink = 0;
    }
    else {
      std::swap(worm_head_it, worm_tail_it);
      worm_at_stop = -1;
      worm_passes_nb_kink = 0;
    }
#ifdef DEBUGMODE
    cout << "\nInitial head (it): " << *worm_head_it;
    cout << "\nInitial tail (it): " << *worm_tail_it;
#endif
    worm_diag = 0;
    new_measurement = true;
    worm_meas_densmat = true;
    return accepted;
  }
  else {
    return rejected;
  }
}




int worm::MOVEWORM() {
#ifdef DEBUGMODE
  if (worm_diag) cerr << "\n# Configuration is diagonal in MOVEWORM ? " << "\n";
  std::cout << "# Welcome to MOVEWORM.\n";
#endif
  
  int dir = (rnd(MyGenerator) < 0.5 ? 1 : -1);
  
  SiteIndex isite = worm_head_it->link();
  Diagram_type::iterator itn = worm_head_it;
  ++itn;
  if (itn == operator_string[isite].end() ) itn = operator_string[isite].begin();
  Diagram_type::iterator itp = worm_head_it;
  if (itp == operator_string[isite].begin() ) itp = operator_string[isite].end();
  --itp;
  
      
   if (dir == +1) {
     if (worm_at_stop == -1) {
       // several cases : dummy, worm_tail, density matrix, interaction
       if (worm_meas_densmat) {
         if ( itn->color() < 0 ) return impossible; // cannot pass worm tail
         
         for (size_t const& j : zc[isite]) {
            if (nb[isite][j] == worm_tail_it->link() ) {
              size_t oppdir = opposite_direction[isite][j];
              Diagram_type::iterator it_nb = worm_head_it->get_assoc(j);
              Diagram_type::iterator it_up = it_nb;
              ++it_up; if (it_up == operator_string[worm_tail_it->link()].end() ) it_up = operator_string[worm_tail_it->link()].begin();
              worm_head_it->set_assoc(j, it_up);
              it_nb->set_assoc(oppdir, worm_head_it);
              break;
            }
         }
         
         worm_at_stop = +1;
       }
       else if (worm_passes_nb_kink) {
         Diagram_type::iterator it_nb, it_up;
         for (auto const& zci : zc[worm_head_it->link()]) {
           SiteIndex adj_site = nb[worm_head_it->link()][zci];
           it_nb = worm_head_it->get_assoc(zci) ;
           it_up = it_nb;
           ++it_up; if (it_up == operator_string[adj_site].end()) it_up == operator_string[adj_site].begin();
           size_t oppdir = opposite_direction[worm_head_it->link()][zci];
           if (!is_not_close(it_nb->time(), worm_head_it->time(), 2*tol)) {
                worm_head_it->set_assoc(zci, it_up);
                it_nb->set_assoc(oppdir, worm_head_it);
           }
         }

         worm_at_stop = +1;
       }
       // when measuring density matrix on other site : no problem
       else if (itn->color() == 0) { // dummy
         PASS_DUMMY(dir, isite ); // no further problem
         return accepted;
       }
       else if (itn->color() == 1) {
         return impossible; // cannot pass interaction
       }
       else {
         cerr << "\n# Error in MOVEWORM : " << dir << "\t" << itn->color() << "\t" << itn->time() << "\t" << worm_head_it->time() << "\t" << worm_tail_it->time() << "\n";
         exit(1);
       }
     }
     
    
     get_diaginfo_nb_pos(isite);
     
     double t_loc = ( t_between(worm_tail_it->time(), worm_head_it->time(), itn->time() ) ) ? worm_tail_it->time() : itn->time() ; // for density matrix
     if (worm_meas_densmat) t_loc = itn->time(); // exceptional case
     nbtime[zcmax] = t_loc;
     double dt_loc = t_loc - worm_head_it->time();
     if (dt_loc < 0) dt_loc += beta;
     dtime[zcmax] = dt_loc;
     
     size_t nb_next = zcmax; double dt_min = dtime[zcmax];
     double tmax = nbtime[zcmax];
     if (t_between(worm_tail_it->time() , worm_head_it->time(), tmax) ) {
       tmax = worm_tail_it->time();
       dt_min = worm_tail_it->time() - worm_head_it->time();
       if (dt_min < 0) dt_min += beta;
     }
     if (worm_meas_densmat) {
       tmax = nbtime[zcmax];
       dt_min = dtime[zcmax];
     }
     for (size_t const& i : zc[isite]) if (dtime[i] < dt_min) {
       dt_min = dtime[i];
       nb_next = i;
       tmax = nbtime[i];
     }
    

#ifdef DEBUGMODE 
     cout << "# MOVEWORM dir == +1 tmax " << tmax << "\tworm_meas_densmat " << ( worm_meas_densmat ? "true" : "false") << "\n";
     for (auto const& i : zc[isite]) cout << "# MOVEWORM dir == +1 nb, dtime, nbtime " << i << "\t" << dtime[i] << "\t" << nbtime[i] << "\n";
     cout << "# MOVEWORM dir == +1 nb, dtime, nbtime " << zcmax << "\t" << dtime[zcmax] << "\t" << nbtime[zcmax] << "\n";
#endif

     double tau = ( (tmax > worm_head_it->time() ) ? tmax - worm_head_it->time() : beta + tmax - worm_head_it->time() );
     StateType n_R = itn->before();
     double E_R = Diag_energy(isite, n_R, nbval);
     StateType n_L = itp->after();
     double E_L = Diag_energy(isite, n_L, nbval);
     double ratio, pexp, en_nom, en_denom;
     if (E_R > E_L) {
       pexp = -log(rnd(MyGenerator))/E_off;
       en_denom = (tau > pexp || abs(pexp-tau) < dtol) ?  E_off : 1.;
       en_nom = (worm_at_stop == 0) ? E_R - E_L + E_off : 1.;
     }
     else {
       pexp = -log(rnd(MyGenerator)) / (E_L - E_R + E_off);
       en_denom = (tau > pexp || abs(pexp-tau) < dtol) ? (E_L - E_R + E_off) : 1.;
       en_nom = (worm_at_stop == 0) ? E_off : 1.;
     }
     if (pexp < dtol) return impossible;
     ratio = en_nom / en_denom;
#ifdef DEBUGMODE
     cout << "MOVEWORM : worm head site-time : " << worm_head_it->link() << "\t" << worm_head_it->time()   << "\n";
     cout << "MOVEWORM : dir " << dir << "\t tau " << tau << "\t pexp : " << pexp << "\t n_R " << n_R << "\t n_L " << n_L << "\t E_R : " << E_R << "\t E_L : " << E_L << "\t ratio : " << ratio << "\n";
#endif
     if (Metropolis(ratio)) {
       if (pexp > tau || abs(pexp-tau) < dtol) {
#ifdef CAN_WINDOW
         if (abs( (Nprtcls + (n_L - n_R) * tau)/beta - canonical) >= can_window) return impossible;
#endif
         double t_w = tmax; //worm_head.time() + tau;
        
         worm_head_it->time(t_w);
         worm_dtime += tau;
         Epot_tot += tau * (E_L - E_R) / beta;
         Nprtcls += tau * (n_L - n_R);
         worm_at_stop = -1;
         worm_passes_nb_kink = (nb_next == zcmax ? 0 : 1);
         worm_meas_densmat =  (worm_head_it->time() == worm_tail_it->time() ? true : false);
         return accepted;
       }
       else {
         double t_w = worm_head_it->time() + pexp;
         if (t_w > beta) t_w -= beta;
#ifdef CAN_WINDOW
         if (abs( (Nprtcls + (n_L - n_R) * pexp)/beta - canonical) >= can_window) return impossible;
#endif
         worm_head_it->time(t_w);
         worm_dtime += pexp;
         Epot_tot += pexp * (E_L - E_R) / beta;
         Nprtcls += pexp * (n_L - n_R);
         worm_at_stop = 0;
         worm_passes_nb_kink = 0;
         worm_meas_densmat = false;
         return accepted;
       }
     }
     else {
       return rejected;
     }
   } // ... dir == +1
   else { // dir == -1
     if (worm_at_stop == +1) {
       // several cases : dummy, worm_tail, density matrix, interaction
       if ( worm_meas_densmat ) {
         if (itp->color() < 0) return impossible; // cannot pass worm tail
         
         Diagram_type::iterator it_up = worm_head_it;
         ++it_up; if (it_up == operator_string[isite].end()) it_up = operator_string[isite].begin();
         for (size_t const& j : zc[isite]) {
           if (nb[isite][j] == worm_tail_it->link()) {
             size_t oppdir = opposite_direction[isite][j];
             worm_head_it->set_assoc(j, worm_tail_it);
             worm_tail_it->set_assoc(oppdir, it_up);
             break;
          }
        }
         
         worm_at_stop = -1;
         // else do nothing
       }
       else if (worm_passes_nb_kink) {
         // do nothing ?
         Diagram_type::iterator it_up = worm_head_it;
         ++it_up; if (it_up == operator_string[isite].end()) it_up = operator_string[isite].begin();
         Diagram_type::iterator it_nb;
         for (auto const& zci : zc[isite]) {
           SiteIndex adj_site =nb[isite][zci];
           it_nb = worm_head_it->get_assoc(zci);
           if (it_nb == operator_string[adj_site].begin()) it_nb = operator_string[adj_site].end();
           --it_nb;
           size_t oppdir = opposite_direction[worm_head_it->link()][zci];
           if (!is_not_close(it_nb->time(), worm_head_it->time(), 2*tol)) {
                worm_head_it->set_assoc(zci, it_nb);
                it_nb->set_assoc(oppdir, it_up );
           }
         }
         worm_at_stop = -1;
       }
       // when measuring density matrix on other side : no problem
       else if (itp->color() == 0) { // dummy
         PASS_DUMMY(dir, isite ); // no further problem
         return accepted;
       }
       else if (itp->color() == 1) {
         return impossible; // cannot pass interaction          // XXXX WHY NOT PASSINTERACTION???
       }
       else  {
         cerr << "\n# Error in MOVEWORM : " << dir << "\t" << itp->color() << "\t" << itp->time() << "\t" << worm_head_it->time() << "\t" << worm_tail_it->time() << "\n";
         exit(1);
       }
     }
    
     get_diaginfo_nb_neg(isite);
     
     double t_loc = ( t_between(worm_tail_it->time(), itp->time(), worm_head_it->time() ) ) ? worm_tail_it->time() : itp->time() ; // for density matrix
     if (worm_meas_densmat) t_loc = itp->time(); // exceptional case
     nbtime[zcmax] = t_loc;
     double dt_loc =  worm_head_it->time() - t_loc;
     if (dt_loc < 0) dt_loc += beta;
     dtime[zcmax] = dt_loc;
    
     double tmax = itp->time();
    
     size_t nb_next = zcmax; double dt_min = dtime[zcmax];
     if (t_between(worm_tail_it->time(), tmax, worm_head_it->time()) ) { // for density matrix
       tmax = worm_tail_it->time();
       dt_min = worm_head_it->time() - worm_tail_it->time();
       if (dt_min < 0) dt_min += beta;
     }
     if (worm_meas_densmat) {
       tmax = nbtime[zcmax];
       dt_min = dtime[zcmax];
     }
     for (size_t const& i : zc[isite]) if (dtime[i] < dt_min) {
       dt_min = dtime[i];
       nb_next = i;
       tmax = nbtime[i];
     }
#ifdef DEBUGMODE
     cout << "# MOVEWORM dir == -1, tmax " << tmax << "\tworm_meas_densmat " << ( worm_meas_densmat ? "true" : "false") << "\n";
     for (size_t const& i : zc[isite]) cout << "# MOVEWORM dir == -1, nb, dtime, nbtime " << i << "\t" << dtime[i] << "\t" << nbtime[i] << "\n";
     cout << "# MOVEWORM dir == -1, nb, dtime, nbtime " << zcmax << "\t" << dtime[zcmax] << "\t" << nbtime[zcmax] << "\n";
     //char ch; cin >> ch;
#endif
     double tau = ( (tmax < worm_head_it->time() ) ? worm_head_it->time() - tmax : beta + worm_head_it->time() - tmax );
    
     StateType n_R = worm_head_it->after();
     StateType n_L = worm_head_it->before();
     double E_L = Diag_energy(isite, n_L, nbval );
     double E_R = Diag_energy(isite, n_R, nbval );
     double pexp, ratio, en_nom, en_denom;
     if (E_R > E_L) {
       pexp = -log(rnd(MyGenerator)) / (E_R - E_L + E_off);
       //ratio = E_off/ (E_R - E_L + E_off);
       en_denom = (tau > pexp || abs(pexp-tau) < dtol) ? E_R - E_L + E_off : 1.;
       en_nom = (worm_at_stop == 0) ? E_off : 1.;
     }
     else {
       pexp = -log(rnd(MyGenerator)) / (E_off);
       en_denom = (tau > pexp || abs(pexp-tau) < dtol) ? E_off : 1.;
       en_nom = (worm_at_stop == 0) ? E_L - E_R + E_off : 1.;
     }
     if (pexp < dtol) return impossible;
     ratio = en_nom / en_denom;
#ifdef DEBUGMODE
     cout << "MOVEWORM : worm head site-time : " << worm_head_it->link() << "\t" << worm_head_it->time()  << "\n";
     cout << "MOVEWORM : dir " << dir << "\t tau " << tau << "\t pexp : " << pexp << "\tn_R " << n_R << "\tn_L " << n_L << "\tE_R : " << E_R << "\t E_L : " << E_L << "\t ratio : " << ratio << "\n";
#endif
     if (Metropolis(ratio)) {
       if (pexp > tau || abs(pexp-tau) < dtol) {
#ifdef CAN_WINDOW
         if (abs( (Nprtcls + (n_R - n_L) * tau)/beta - canonical) >= can_window) return impossible;
#endif
         double t_w = tmax;
         if (tmax == beta) t_w = 0;
         
         worm_head_it->time(t_w);
         worm_dtime -= tau;
         worm_at_stop = +1;
         Epot_tot += (E_R - E_L) * tau / beta;
         Nprtcls += tau * (n_R - n_L);
         worm_passes_nb_kink = (nb_next == zcmax ? 0 : 1);
#ifdef DEBUGMODE
         cout << "MOVEWORM : worm head site-time : " << worm_head_it->link() << "\t" << worm_head_it->time()  << "\n";
         cout << "MOVEWORM : dir " << dir << "\t tau " << tau << "\t pexp : " << pexp << "\tn_R " << n_R << "\tn_L " << n_L << "\tE_R : " << E_R << "\t E_L : " << E_L << "\t ratio : " << ratio << "\n";
         cout << "MOVEWORM : " << worm_head_it->time() << "\t" << worm_tail_it->time() << "\n";
#endif
         worm_meas_densmat =  (worm_head_it->time() == worm_tail_it->time() ? true : false);
         return accepted;
       }
       else {
         double t_w = worm_head_it->time() - pexp;
         if (t_w < 0) t_w += beta;
#ifdef CAN_WINDOW
         if (abs( (Nprtcls + (n_R - n_L) * pexp)/beta - canonical) >= can_window) return impossible;
#endif
         worm_head_it->time(t_w);
         worm_dtime -= pexp;
         worm_at_stop = 0;
         Epot_tot += pexp * (E_R - E_L) / beta;
         Nprtcls  += pexp * (n_R - n_L);
         worm_passes_nb_kink = 0;
         worm_meas_densmat = false;
         return accepted;
       }
     }
     else {
       return rejected;
    }
  } // ... dir == -1
}



int worm::INSERTKINK() {
 
#ifdef DEBUGMODE
  if (worm_diag) cerr << "Diagonal configuration in INSERTKINK ? \n\n";
#endif
  if (worm_at_stop != 0) return impossible;
  if (worm_meas_densmat) return impossible;
  if (worm_passes_nb_kink == 1) return impossible;
  
  SiteIndex isite = worm_head_it->link();

  // choose a direction (temporal)
  int dir = ( (rnd(MyGenerator) < 0.5) ? -1 : 1 );
  
  // choose a neighbor (spatial)
  size_t dirIndex = static_cast<size_t> (rnd(MyGenerator) * zc[isite].size());
  size_t nbs = zc[isite][dirIndex];
  SiteType adj_site = nb[isite][nbs];

#ifdef DEBUGMODE
  cout << "INSERTKINK dir : "<< dir << "\tnbs " << nbs <<"\tadj_site " << adj_site << "\n";
#endif
  
  StateType n_up = worm_head_it->after();
  StateType n_down = worm_head_it->before();
  StateType n_diff = n_up - n_down;
  
  // let's set the iterator on the adjacent site right
  Diagram_type::iterator it_adj = worm_head_it->get_assoc(nbs);
  Diagram_type::iterator itp_adj = it_adj;

  if (itp_adj == operator_string[adj_site].begin()) itp_adj = operator_string[adj_site].end();
  --itp_adj;
  while (!t_between(worm_head_it->time(), itp_adj->time(), it_adj->time() ) ) {
    it_adj = itp_adj;
    if (itp_adj == operator_string[adj_site].begin()) itp_adj = operator_string[adj_site].end();
    --itp_adj;
  }  
  StateType nbs_out = it_adj->before();
  StateType nbs_mid = nbs_out - dir*n_diff;

#ifdef DEBUGMODE
  cout << "INSERTKINK : n_up" << n_up << "\tn_down" << n_down << "\tnbs_out " << nbs_out << "\tnbs_mid " << nbs_mid << "\n";
#endif
  
#ifdef UNISYS
  if (!range_check(nbs_mid)) return impossible;
#else
  if (MyModel->range_fail(nbs_mid)) return impossible;
#endif
    
    
#ifdef UNISYS
  double weight_old = site_weight_offdiag[n_down][ n_up];
  double weight_new = (dir == 1 ? bond_weight_offdiag[n_down][ n_up][ nbs_out][ nbs_mid] : bond_weight_offdiag[ n_down][ n_up][ nbs_mid][ nbs_out] );
  weight_new *= site_weight_offdiag[nbs_mid][ nbs_out];
#else
  BondIndex b = bond_index[isite][ nbs];
  double weight_old = MyModel->site_weight_offdiag(isite, n_down, n_up);
  double weight_new = (dir == 1 ? MyModel->bond_weight_offdiag(b, n_down, n_up, nbs_out, nbs_mid) : MyModel->bond_weight_offdiag(b, n_down, n_up, nbs_mid, nbs_out) );
  weight_new *= MyModel->site_weight_offdiag(adj_site, nbs_mid, nbs_out);
#endif
    
  double ratio = 2 * zc[isite].size() * weight_new /weight_old * update_prob[deletekink] / (update_prob[insertkink]);
  
  
#ifdef DEBUGMODE
  cout << "INSERTKINK : ratio " << ratio << "\t" << weight_new << "\t" << weight_old << "\n";
#endif
  
  
  if (Metropolis(ratio)) {
#ifdef DEBUGMODE
    cout << "INSERTKINK accepted. worm_head_it " << worm_head_it->time() << "\tratio : " << ratio << "\n";
    cout << "INSERTKINK adj_site : " << adj_site << "\tit_adj " << it_adj->time() << "\n";
#endif
   
    // insert the interaction : iw on the worm site, change the properties of the worm_head_it and insert a new kink on the adj site + a new worm there

    Element_t new_elem1(nbs_mid, nbs_out, isite,     worm_head_it->time(),  1);
    Element_t new_elem2(nbs_out, nbs_mid, adj_site,  worm_head_it->time(), -1);
    worm_head_it->color(1);
    worm_head_it->link(adj_site);
    
    if (dir ==1 ) {
      new_elem1.before(nbs_out);
      new_elem1.after(nbs_mid);
      new_elem2.before(nbs_mid);
      new_elem2.after(nbs_out);
      worm_head_it = operator_string[adj_site].insert(it_adj, new_elem2);
      find_assoc_insert(adj_site, worm_head_it, +1);
      it_adj = operator_string[adj_site].insert(worm_head_it, new_elem1);
      find_assoc_insert(adj_site, it_adj, 0);
      worm_at_stop = +1;
    }
    else {
      it_adj = operator_string[adj_site].insert(it_adj, new_elem1);
      find_assoc_insert(adj_site, it_adj, 0);
      Diagram_type::iterator new_it = operator_string[adj_site].insert(it_adj, new_elem2);
      find_assoc_insert(adj_site, new_it, -1);
      worm_head_it->set_assoc(nbs, it_adj);
      worm_head_it = new_it;
      worm_at_stop = -1;
    }
        
    // worm jumps
    nrvertex++;
    worm_passes_nb_kink = 0;
    worm_meas_densmat = false;
    if (MyLatt->relNumbering(isite, adj_site)) {
        mWinding += winding_element[nbs]  *  (n_down - n_up);
    }
    return accepted;
  }
  else {
    return rejected;
  }
}


int worm::DELETEKINK() {
  if (worm_meas_densmat ) return impossible; // no density matrix measurement
  
#ifdef DEBUGMODE
  cout << "# Welcome to DELETEKINK : \n";
#endif

  Diagram_type::iterator it = worm_head_it;                                    // it should become the iterator pointing to the kink to be deleted; needs to be changed depending on worm_at_stop == +/- 1
  SiteIndex isite = it->link();                                                   // worm is linked to itself
  StateType n_out;
#ifdef DEBUGMODE
  StateType n_mid;                                                         // n_mid is the occupation between worm and interaction, n_out is the "outside"
#endif
  
  if (worm_at_stop == +1 && worm_passes_nb_kink == 0) {                           // ie, the worm is just above the kink
    if (it == operator_string[isite].begin()) it = operator_string[isite].end();
    --it;
    n_out = it->before();
    if (it->color() != +1) return impossible;
#ifdef DEBUGMODE
    n_mid = it->after();
#endif
    if (n_out != worm_head_it->after() ) return ( PASSINTERACTION(-1, it) );      // this is a situation where n_out = 1, n_mid = 2, and the occup after the worm is 3.
  } 
  else if (worm_at_stop == -1 && worm_passes_nb_kink == 0) { // dir==+1           // ie, the worm is just below the kink
    ++it; if (it == operator_string[isite].end()) it = operator_string[isite].begin();
    if (it->color() != +1 ) return impossible;
#ifdef DEBUGMODE
    n_mid = it->before();
#endif
    n_out = it->after();
    if (n_out != worm_head_it->before() ) return (PASSINTERACTION(+1, it));       // this is a situation where n_out = 1, n_mid = 2, and the occup after the worm is 3.
  }
  else {
    return impossible;
  }
  
  SiteIndex adj_site = it->link();                                    // site to wich the current site is linked
  size_t linkdir = 0;
    for (auto const& zci : zc[isite]) {
      linkdir = zci;
      if (nb[isite][linkdir] == adj_site) break;
    }
  Diagram_type::iterator itlink = it->get_assoc(linkdir);
  StateType n_A_link = itlink->after();                               // occupancy on kink after  the hopping event on the linked site
  StateType n_B_link = itlink->before();                              // occupancy on kink before the hopping event on the linked site
  
  // weight of tunneling element
#ifdef UNISYS
  double weight_old = bond_weight_offdiag[it->before()][ it->after()][ n_B_link][ n_A_link] * site_weight_offdiag[worm_head_it->before()][ worm_head_it->after() ];
  double weight_new = site_weight_offdiag[n_B_link][ n_A_link];
#else
  BondIndex b = bond_index[isite][linkdir];               // bond index of the kink
  double weight_old = MyModel->bond_weight_offdiag(b, it->before(), it->after(), n_B_link, n_A_link) * MyModel->site_weight_offdiag(isite, worm_head_it->before(), worm_head_it->after() );
  double weight_new = MyModel->site_weight_offdiag(adj_site, n_B_link, n_A_link);
#endif
    
  double ratio = weight_new / weight_old * update_prob[insertkink] / (update_prob[deletekink] * 2* zc[isite].size());
  
  
#ifdef DEBUGMODE
   cout << "DELETEKINK : n_out " << n_out << "\tn_mid " << n_mid  << "\n";
   cout << "DELETEKINK it->link() " << it->link() << "\tit->time()" << it->time() << "\n";
#endif

  
#ifdef DEBUGMODE
  
   cout << "DELETEKINK : ratio " << ratio << "\t" << 1./ratio << "\n";
  cout << "DELETEKINK old weight, new weight : " << "\t" << weight_old << "\t" << weight_new << "\n";
#endif

  if (Metropolis(ratio)) {
    if (worm_at_stop == +1) {
      find_assoc_delete(isite, worm_head_it);
      operator_string[isite].erase(worm_head_it);
      find_assoc_delete(isite, it);
      operator_string[isite].erase(it);
    }
    else {
      find_assoc_delete(isite, it);
      operator_string[isite].erase(it);
      find_assoc_delete(isite, worm_head_it);
      operator_string[isite].erase(worm_head_it);
    }
    worm_head_it = itlink;
    worm_head_it->color(-1);
    worm_head_it->link(adj_site);
    nrvertex--;
    worm_at_stop = 0;
    worm_passes_nb_kink = 0;
    worm_meas_densmat = false;
    if (MyLatt->relNumbering(isite, adj_site)) {
        mWinding += winding_element[linkdir] * (n_B_link - n_A_link);
    }
    return accepted;
  }
  else { // Metropolis rejection
    return rejected;
  }
  
}


void worm::PASS_DUMMY(const int dir, const SiteType isite) {
#ifdef DEBUGMODE
  cout << "\n# PASSING DUMMY... dir = " << dir << "\n";
#endif

  auto update_internal_variables = [&] (StateType old_dens, StateType new_dens) {
#ifdef UNISYS
    double Ep_substract = site_weight_diag[old_dens];
    double Ep_add = site_weight_diag[new_dens];
    for (size_t const& iz : zc[isite]) {
      Ep_substract += bond_weight_diag[old_dens][dummy_it[nb[isite][iz]]->before()];
      Ep_add += bond_weight_diag[new_dens][dummy_it[nb[isite][iz]]->before()];
    }
#else
    double Ep_substract = MyModel->site_weight_diag(isite, old_dens);
    double Ep_add = MyModel->site_weight_diag(isite, new_dens);
    for (size_t const &iz : zc[isite]) {
      Ep_substract += MyModel->bond_weight_diag(bond_index[isite][iz], old_dens, dummy_it[nb[isite][iz]]->before());
      Ep_add += MyModel->bond_weight_diag(bond_index[isite][iz], new_dens, dummy_it[nb[isite][iz]]->before());
    }
#endif

    Epot_measure += Ep_add;
    Epot_measure -= Ep_substract;

    number_of_particles += new_dens - old_dens;
    state[isite] = new_dens;
  };

  if (dir == +1) {
    StateType new_dens = worm_head_it->before();
    StateType old_dens = dummy_it[isite]->before();
    update_internal_variables(old_dens, new_dens);

    dummy_it[isite]->before(new_dens);
    dummy_it[isite]->after(new_dens);
    Element_t new_elem = *worm_head_it;
    find_assoc_delete(isite, worm_head_it);
    operator_string[isite].erase(worm_head_it);
    new_elem.time(0);
    worm_head_it = operator_string[isite].insert(operator_string[isite].begin(), new_elem);
    find_assoc_insert(isite, worm_head_it, 0);
    worm_at_stop = +1;
    worm_passes_nb_kink = 0;
    worm_meas_densmat = false;
  }
  else {
    StateType new_dens = worm_head_it->after();
    StateType old_dens = dummy_it[isite]->before();
    update_internal_variables(old_dens, new_dens);

    dummy_it[isite]->before(new_dens);
    dummy_it[isite]->after(new_dens);
    Element_t new_elem = *worm_head_it;
    find_assoc_delete(isite, worm_head_it);
    operator_string[isite].erase(worm_head_it);
    new_elem.time(beta);
    worm_head_it = operator_string[isite].insert(dummy_it[isite], new_elem);
    find_assoc_insert(isite, worm_head_it, -1);
    worm_at_stop = -1;
    worm_passes_nb_kink = 0;
    worm_meas_densmat = false;
  }
#ifdef DEBUGMODE
  print_conf(std::cout);
  std::cout << "# ... done with PASSDUMMY\n";
#endif
  
  return;
}

int worm::PASSINTERACTION(const int dir, const Diagram_type::iterator it) {
  // here : it is the iterator to the kink that will be passed, not the worm_head_it
  SiteIndex isite = worm_head_it->link();                             // current site
  StateType n_A = it->after();                                        // occupancy on kink after the hopping event on the current site  (before and after as always in the view of positive time [0, beta[)
  StateType n_B = it->before();                                       // occupancy on kink before the hopping event on the current site
  SiteIndex adj_site = it->link();                                    // site to wich the current site is linked
  size_t linkdir = 0;
    for (auto const& zci : zc[isite]) {
      linkdir = zci;
      if (nb[isite][linkdir] == adj_site) break;
    }
  // TODO : assert linkdir is not out of bounds
  StateType n_A_link = it->get_assoc(linkdir)->after();               // occupancy on kink after  the hopping event on the linked site
  StateType n_B_link = it->get_assoc(linkdir)->before();              // occupancy on kink before the hopping event on the linked site
  
  if (dir == 1) {
    
    //n_O = 2*n_B - n_A;
    StateType n_O = worm_head_it->before();
    
#ifdef DEBUGMODE
    cout << "PASSINTERACTION in positive direction. it->time " << it->time() << "\n";
    cout << "PASSINTERACTION accepted.\n";
#endif
    
#ifdef UNISYS
    double wgt_old = site_weight_offdiag[n_O][ n_B] * bond_weight_offdiag[n_B][ n_A][ n_B_link][ n_A_link];
    double wgt_new = site_weight_offdiag[n_B][ n_A] * bond_weight_offdiag[n_O][ n_B][ n_B_link][ n_A_link];
#else
    BondIndex b = bond_index[isite][linkdir];               // bond index of the kink
    double wgt_old = MyModel->site_weight_offdiag(isite, n_O, n_B) * MyModel->bond_weight_offdiag(b, n_B, n_A, n_B_link, n_A_link);
    double wgt_new = MyModel->site_weight_offdiag(isite, n_B, n_A) * MyModel->bond_weight_offdiag(b, n_O, n_B, n_B_link, n_A_link);
#endif
    
    double ratio = wgt_new / wgt_old;
    if (Metropolis(ratio)) {
      it->before(n_O);
      it->after(n_B);
      Element_t new_elem = *worm_head_it;
      find_assoc_delete(isite, worm_head_it);
      operator_string[isite].erase(worm_head_it);
      new_elem.time(it->time());
      new_elem.before(n_B);
      new_elem.after(n_A);
      Diagram_type::iterator itn = it;
      ++itn; if (itn == operator_string[isite].end() ) itn = operator_string[isite].begin();
      worm_head_it = operator_string[isite].insert(itn, new_elem);
      find_assoc_insert(isite, worm_head_it, +1);
      worm_at_stop = +1;
      worm_passes_nb_kink = 0;
      worm_meas_densmat = false;
      return accepted;
    }
    else {
      std::cerr << "# How can Metropolis be rejected in PASSINTERACTION, dir == +1 ???\n";
      return rejected;
    }
  }
  else { // dir == -1
    StateType n_O = worm_head_it->after();
#ifdef UNISYS
    double wgt_old = site_weight_offdiag[n_A][ n_O] * bond_weight_offdiag[n_B][ n_A][ n_B_link][ n_A_link];
    double wgt_new = site_weight_offdiag[n_B][ n_A] * bond_weight_offdiag[n_A][ n_O][ n_B_link][ n_A_link];
    
#else
    BondIndex b = bond_index[isite][linkdir];               // bond index of the kink
    double wgt_old = MyModel->site_weight_offdiag(isite, n_A, n_O) * MyModel->bond_weight_offdiag(b, n_B, n_A, n_B_link, n_A_link);
    double wgt_new = MyModel->site_weight_offdiag(isite, n_B, n_A) * MyModel->bond_weight_offdiag(b, n_A, n_O, n_B_link, n_A_link);
#endif
    double ratio = wgt_new / wgt_old;
  
    if (Metropolis(ratio)) {
      
#ifdef DEBUGMODE
      cout << "PASSINTERACTION in negative direction. it->time " << it->time() << "\n";
      cout << "PASSINTERACTION accepted. " << n_B << "\t" << n_O << "\t" << n_A << "\n";
#endif
      it->before(n_A);
      it->after(n_O);
      Element_t new_elem = *worm_head_it;
      find_assoc_delete(isite, worm_head_it);
      operator_string[isite].erase(worm_head_it);
      new_elem.time(it->time());
      new_elem.before(n_B);
      new_elem.after(n_A);
      worm_head_it = operator_string[isite].insert(it, new_elem);
      find_assoc_insert(isite, worm_head_it, -1);
      worm_at_stop = -1;
      worm_passes_nb_kink = 0;
      worm_meas_densmat = false;
      return accepted;
    }
    else {
      std::cerr << "# How can Metropolis be rejected in PASSINTERACTION, dir == -1 ???\n";
      return rejected;
    }
  } // ...dir == -1
}



int worm::GLUEWORM() {
#ifdef DEBUGMODE
  cout << "# Welcome to GLUEWORM  : worm_at_stop : " << worm_at_stop << "\n";
  if (worm_diag) cerr << "Diagonal configuration in GLUEWORM ? \n\n";
#endif
  size_t isite = worm_head_it->link();
  if (worm_at_stop == 0) return impossible;
  if (worm_passes_nb_kink) return impossible;
  if (worm_tail_it->time() != worm_head_it->time()) return impossible;
  Diagram_type::iterator it_lower = worm_head_it;
  Diagram_type::iterator it_upper = worm_head_it;
  if (worm_at_stop == -1) {
    ++it_upper; if (it_upper == operator_string[isite].end()) it_upper = operator_string[isite].begin();
    if (!(it_upper == worm_tail_it)) return impossible;
  }
  else {
    if (it_lower == operator_string[isite].begin()) it_lower = operator_string[isite].end();
    --it_lower;
    if (!(it_lower == worm_tail_it)) return impossible;
  }
  StateType n_out = it_upper->after();
  if (n_out != it_lower->before()) {
    return impossible;
  }
  StateType n_mid = it_upper->before();
  if (n_mid != it_lower->after()) {
    std::cerr << "Really strange n_mid in GLUEWORM : " << it_upper->before() << "\t" << it_lower->after() << "\n";
    return impossible;
  }
#ifdef UNISYS
  double wgt = site_weight_offdiag[n_out][ n_mid];
#else
  double wgt = MyModel->site_weight_offdiag(isite, n_out, n_mid);
#endif
  
  double ratio = 2 * C_NBW * Nsites * beta * wgt * wgt*  update_prob[glueworm] / update_prob[insertworm];
  
  
#ifdef DEBUGMODE
  cout << "GLUEWORM  ratio " << ratio << "\t" << 1./ratio << "\n";
#endif
  if (Metropolis(1./ratio)) {
#ifdef DEBUGMODE
    cout << "GLUEWORM accepted. ratio " << ratio << "\n";
#endif
    find_assoc_delete(isite, it_upper);
    operator_string[isite].erase(it_upper);
    find_assoc_delete(isite, it_lower);
    operator_string[isite].erase(it_lower);
    
    worm_head_it = dummy_it[isite];
    worm_tail_it = dummy_it[isite];
    worm_diag = 1;
    worm_at_stop = 0;
    worm_passes_nb_kink = 0;
    worm_meas_densmat = false;
    return accepted;
  } // Metropolis
  else { // Metropolis
#ifdef DEBUGMODE
    cout << "GLUEWORM rejected. ratio " << ratio << "\n";
#endif
    return rejected;
  } // else...Metropolis
}

