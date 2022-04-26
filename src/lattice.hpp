#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <numeric>
#include <map>
#include <utility>
#include <vector>

#include <alps/params.hpp>

#include <Eigen/Dense>

//NOTE: crystal_dir 0,1,2 corresponds to directions x,y,z respectively
//this is important for boundary conditions

template <size_t DIM, size_t N_BASIS>
struct lattice {
    static const size_t dim = DIM;
    static const size_t n_basis = N_BASIS;

    //for readability 
    using SiteType = int;       //must allow value -1 for pseudo site which is used for non-uniform coordination number
    using SiteIndex = int;      //same here
    using BondIndex = size_t;
    using BondType = std::pair<size_t, size_t>; 
    using DirectionIndex = size_t;
    using Direction = Eigen::Matrix<double, DIM, 1>;
    using NeighborType = std::map<DirectionIndex, SiteIndex>;

    static void define_parameters(alps::params& params) {
        params.define<std::string>("lattice", "chain", "name of the lattice")
              .define<size_t>("Lx", 4, "number of unitcells along x axis of the lattice")
              .define<size_t>("Ly", 1, "number of unitcells along y axis of the lattice")
              .define<size_t>("Lz", 1, "number of unitcells along z axis of the lattice")
              .define<bool>("pbcx", true, "periodic boundary condition along x axis (PBC:1, OBC:0)")
              .define<bool>("pbcy", true, "periodic boundary condition along y axis (PBC:1, OBC:0)")
              .define<bool>("pbcz", true, "periodic boundary condition along z axis (PBC:1, OBC:0)");
    }

    static bool parameters_supplied(alps::params const& params) {
        return params.supplied("Lx") &&
               params.supplied("Ly") &&
               params.supplied("Lz") &&
               params.supplied("pbcx") &&
               params.supplied("pbcy") &&
               params.supplied("pbcz");
    }

    //interface functions
    size_t get_dim() const {return DIM;}
    size_t get_nbasis() const {return N_BASIS;}
    size_t get_zcoord(SiteIndex s) const {return sites.at(s).size();}
    size_t get_Nsites() const {return sites.size();}
    size_t get_Nbonds() const {return BondIndices.size();}
    size_t get_N_nb(SiteIndex s) const {return sites.at(s).size();}
    size_t get_Ls(int d) const {
      if ((d < 0) || (d > 2)) throw std::runtime_error("call to get_Ls() with invalid argument");
      if (d == 0) { return Lx;}
      else if (d == 1) { return Ly; }
      else { return Lz; }
    };

    SiteType operator[](SiteIndex s) const {return s;}
    SiteType get_site(SiteIndex s) const {return s;}

    bool has_neighbor(SiteIndex s, DirectionIndex dir) const {
        auto search = sites.at(s).find(dir);
        return search != sites.at(s).end();
    }

    SiteType nb(SiteIndex s, DirectionIndex dir) const {return sites.at(s).at(dir);}

    DirectionIndex get_direction(SiteIndex s1, SiteIndex s2) const {
        return BondDirections.at(std::make_pair(s1,s2));
    }
    DirectionIndex get_opposite_direction(SiteIndex s1, SiteIndex s2) const {
        return BondDirections.at(std::make_pair(s2,s1));
    }
    BondIndex get_bond_index(SiteIndex s, DirectionIndex dir) const {
        SiteIndex other = sites.at(s).at(dir);
        BondType bond = s < other ? BondType(s, other) : BondType(other, s);
        return BondIndices.at(bond);
    }

    // Relative numbering of site s2 w.r.t site s1. 
    // This is NOT symmetric in s1, s2, i.e. relNumbering(s1,s2) != relNumbering(s2,s1)
    size_t relNumbering(SiteType s1, SiteType s2) const {
        size_t s1x = (s1/N_BASIS) % Lx;
        size_t s2x = (s2/N_BASIS) % Lx;
        size_t nx = (s2x < s1x ? s2x + Lx - s1x  : s2x - s1x);
        if (DIM > 1) {
            size_t s1y = (s1/N_BASIS) / Lx % Ly;
            size_t s2y = (s2/N_BASIS) / Lx % Ly;
            size_t ny = (s2y < s1y ? s2y + Ly - s1y  : s2y - s1y);
            if (DIM > 2) {
                size_t s1z = (s1/N_BASIS) / Lx / Ly % Lz;
                size_t s2z = (s2/N_BASIS) / Lx / Ly % Lz;
                size_t nz = (s2z < s1z ? s2z + Lz - s1z  : s2z - s1z);
                return nx + Lx*ny + Lx*Ly*nz;
            }
            return nx + Lx*ny;
        }
        return nx;
    }

    virtual void print_name(std::ostream& os) const = 0;

    void print() const {
        std::cout << "# Nsites, dim, zcoord " << sites.size() << " " << DIM << " " 
                  << sites.at(0).size() << "\n";
        for (auto const& s : sites){
            std::cout << "site: " << s.first << "\t\t" << "has neighbors: ";
            for (auto const& s2 : s.second)
                std::cout << s2.second << "(dir " << s2.first << ")\t";
            std::cout << std::endl;
        }
        for (auto const& b : BondIndices){
            std::cout << "bond: " << b.first.first << "," << b.first.second << "\t\t" << "has index: "
                        << b.second << std::endl;
        }
    }

    void print_params(std::ostream& os) const {
        os << "# lattice: dim                                      : " << DIM << "\n";
        os << "# lattice: pbcx                                     : " << pbcx << "\n";
        os << "# lattice: pbcy                                     : " << pbcy << "\n";
        os << "# lattice: pbcz                                     : " << pbcz << "\n";
        os << "# lattice: length                                   : " 
                << "Lx, Ly, Lz  " << Lx << "\t" << Ly << "\t" << Lz << "\t"
                << "num Sites  " << sites.size() << "\t"
                << "num Bonds  "  << BondIndices.size() << "\n";
    }

    
protected:
    constexpr SiteType next(size_t crystal_dir, SiteType bidx) const {
        if(bidx == -1) return bidx;
        size_t lower = plengths[crystal_dir];
        size_t upper = plengths[crystal_dir+1];
        bidx += lower * N_BASIS;
        if (bidx / N_BASIS % upper / lower == 0) {
            if ((crystal_dir == 0 && pbcx) ||
                (crystal_dir == 1 && pbcy) ||
                (crystal_dir == 2 && pbcz))
                bidx -= upper * N_BASIS;
            else
                return -1;
        }
        return bidx;
    }

    constexpr SiteType previous(size_t crystal_dir, SiteType bidx) const {
        if(bidx == -1) return bidx;
        size_t lower = plengths[crystal_dir];
        size_t upper = plengths[crystal_dir+1];
        if (bidx / N_BASIS % upper / lower == 0) {
            if ((crystal_dir == 0 && pbcx) ||
                (crystal_dir == 1 && pbcy) ||
                (crystal_dir == 2 && pbcz))
                bidx += upper * N_BASIS;
            else
                return -1;
        }
        return bidx -= lower * N_BASIS;
    }

    template <typename action>
    void generate_lattice(action&& define_unitcell) {
        for (size_t bs = 0; bs < N_BASIS * Lx * Ly * Lz; bs += N_BASIS)
            define_unitcell(bs);
        size_t bondIndex = 0;
        for (auto& s : sites){
            for (NeighborType::iterator nb_it = (s.second).begin(); nb_it != (s.second).end();) {
                //if any pseudosites present because of open boundary conditions, remove them
                if (nb_it->second == -1) {
                    nb_it = (s.second).erase(nb_it);
                }
                else {
                    BondType bond = s.first <= nb_it->second ?
                                    BondType(s.first, nb_it->second) : BondType(nb_it->second, s.first);
                    if (BondIndices.find(bond) == BondIndices.end())
                        BondIndices.emplace(std::make_pair(bond, bondIndex++));
                    BondDirections.emplace(std::make_pair(BondType(s.first, nb_it->second), nb_it->first));
                    ++nb_it;
                }
            }
        }
    }
   
    lattice(size_t Lx, size_t Ly = 1, size_t Lz = 1,
            bool pbcx = true, bool pbcy = true, bool pbcz = true)
        : Lx(Lx),
          Ly(Ly),
          Lz(Lz),
          pbcx(pbcx),
          pbcy(pbcy),
          pbcz(pbcz)
    {
        plengths = { 1, Lx, Ly, Lz};
        std::partial_sum(plengths.begin(), plengths.end(), plengths.begin(),
                         std::multiplies<size_t>{});
    }

    std::map<SiteType, NeighborType> sites;

private:
    std::map<BondType, BondIndex> BondIndices; //without double counting
    std::map<BondType, DirectionIndex> BondDirections; //with double counting in opposite direction

    size_t const Lx;
    size_t const Ly;
    size_t const Lz;
    std::vector<size_t> plengths;
    bool const pbcx;
    bool const pbcy;
    bool const pbcz;
};


struct cubic : lattice<3,1> {
    using Base = lattice;
    static const size_t zcmax = 6;

    cubic(size_t Lx, size_t Ly, size_t Lz, bool pbcx=true, bool pbcy=true, bool pbcz=true)
        : Base(Lx, Ly, Lz, pbcx, pbcy, pbcz)
    {
        if (((Lx < 3) && pbcx) || ((Ly < 3) && pbcy) || ((Lz < 3) && pbcz))
            throw std::runtime_error("PBC with linear length < 3 leads to double counting of bonds");
        generate_lattice([&](size_t bs) {return define_unitcell(bs);});
    }
    
    void define_unitcell (size_t base) {
        auto mp = [](Base::DirectionIndex d, size_t site) {return std::make_pair(d, site);};
        NeighborType nb0 = {mp(0,next(0,base)),
                            mp(3,previous(0,base)), 
                            mp(1,next(1,base)), 
                            mp(4,previous(1,base)), 
                            mp(2,next(2,base)),
                            mp(5,previous(2,base)),
                           };

        sites.emplace(std::make_pair(base, nb0));
    }

    std::array<Direction,zcmax> directions() {
        return {
            (Direction() <<  1,  0,  0).finished(),
            (Direction() <<  0,  1,  0).finished(),
            (Direction() <<  0,  0,  1).finished(),
            (Direction() << -1,  0,  0).finished(),
            (Direction() <<  0, -1,  0).finished(),
            (Direction() <<  0,  0, -1).finished(),
        };
    }

    void print_name(std::ostream& os) const {
        os << "# lattice: Name                                     : cubic\n";
    };
};


struct square : lattice<2,1> {
    using Base = lattice;
    static const size_t zcmax = 4;

    square(size_t Lx, size_t Ly, size_t Lz, bool pbcx=true, bool pbcy=true, bool pbcz=false)
        : Base(Lx, Ly, 1, pbcx, pbcy, false)
    {
        if (((Lx < 3) && pbcx) || ((Ly < 3) && pbcy))
            throw std::runtime_error("PBC with linear length < 3 leads to double counting of bonds");
        generate_lattice([&](size_t bs) {return define_unitcell(bs);});
    }
    
    void define_unitcell (size_t base) {
        auto mp = [](Base::DirectionIndex d, size_t site) {return std::make_pair(d, site);};
        NeighborType nb0 = {mp(0,next(0,base)),
                            mp(2,previous(0,base)),
                            mp(1,next(1,base)),
                            mp(3,previous(1,base)),
                           };
        sites.emplace(std::make_pair(base, nb0));
    }

    std::array<Direction,zcmax> directions() {
        return {
            (Direction() <<  1, 0).finished(),
            (Direction() <<  0, 1).finished(),
            (Direction() << -1, 0).finished(),
            (Direction() <<  0,-1).finished()
        };
    }

    void print_name(std::ostream& os) const {
        os << "# lattice: Name                                     : square\n";
    };
};


struct chain : lattice<1,1> {
    using Base = lattice;
    static const size_t zcmax = 2;

    chain(size_t Lx, size_t Ly, size_t Lz, bool pbcx=true, bool pbcy=false, bool pbcz=false)
        : Base(Lx, 1, 1, pbcx, false, false)
    {
        if ((Lx < 3) && pbcx)
            throw std::runtime_error("PBC with linear length < 3 leads to double counting of bonds");
        generate_lattice([&](size_t bs) {return define_unitcell(bs);});
    }

    void define_unitcell (size_t base) {
        auto mp = [](Base::DirectionIndex d, size_t site) {return std::make_pair(d, site);};
        NeighborType nb0 = {mp(0,next(0,base)), mp(1,previous(0,base))};

        sites.emplace(std::make_pair(base, nb0));
    }

    std::array<Direction,zcmax> directions() {
        return {
            (Direction() <<  1).finished(),
            (Direction() << -1).finished(),
        };
    }

    void print_name(std::ostream& os) const {
        os << "# lattice: Name                                     : chain\n";
    };
};


//treated as special case of square with Ly=2 and pbcy=false
struct ladder : square {

    ladder(size_t Lx, size_t Ly, size_t Lz, bool pbcx=true, bool pbcy=false, bool pbcz=false)
        : square(Lx, 2, 1, pbcx, false, false)
    {
        if ((Lx < 3) && pbcx)
            throw std::runtime_error("PBC with linear length < 3 leads to double counting of bonds");
        generate_lattice([&](size_t bs) {return define_unitcell(bs);});
    }

    void print_name(std::ostream& os) const {
        os << "# lattice: Name                                     : ladder\n";
    };
};


struct honeycomb : lattice<2,2> {
    using Base = lattice;
    static size_t const zcmax = 6;

    honeycomb(size_t Lx, size_t Ly, size_t Lz, bool pbcx=true, bool pbcy=true, bool pbcz=false)
        : Base(Lx, Ly, 1, pbcx, pbcy, false)
    {
        if (((Lx < 2) && pbcx) || ((Ly < 2) && pbcy))
            throw std::runtime_error("PBC with linear length < 2 leads to double counting of bonds");
        generate_lattice([&](size_t bs) {return define_unitcell(bs);});
    }

    void define_unitcell (size_t base) {
        auto mp = [](Base::DirectionIndex d, size_t site) {return std::make_pair(d, site);};

        NeighborType
            nb0 = {mp(0,base+1), mp(1,previous(0,base+1)), mp(2,previous(1,base+1))},
            nb1 = {mp(3,base+0), mp(5,next(0,base+0)), mp(4,next(1,base+0))};

        sites.emplace(std::make_pair(base+0, nb0));
        sites.emplace(std::make_pair(base+1, nb1));
    }

    //Directions are set to produce the correct Windingnumbers in crystal dir 0 and 1
    std::array<Direction,zcmax> directions() {
        return {
            (Direction() <<   0,  0).finished(),
            (Direction() <<  -1,  0).finished(),
            (Direction() <<   0, -1).finished(),
            (Direction() <<   0,  0).finished(),
            (Direction() <<   0,  1).finished(),
            (Direction() <<   1,  0).finished(),
        };
    }

    void print_name(std::ostream& os) const {
        os << "# lattice: Name                                     : honeycomb\n";
    };
};


struct triangular : lattice<2,1> {
    using Base = lattice;
    static size_t const zcmax = 6;

    triangular(size_t Lx, size_t Ly, size_t Lz, bool pbcx=true, bool pbcy=true, bool pbcz=false)
        : Base(Lx, Ly, 1, pbcx, pbcy, false)
    {
        if (((Lx < 3) && pbcx) || ((Ly < 3) && pbcy))
            throw std::runtime_error("PBC with linear length < 3 leads to double counting of bonds");
        generate_lattice([&](size_t bs) {return define_unitcell(bs);});
    }

    void define_unitcell (size_t base) {
        auto mp = [](Base::DirectionIndex d, size_t site) {return std::make_pair(d, site);};

        NeighborType
            nb0 = {mp(0,next(0,base)), mp(1,next(1,base)), mp(2,next(1,previous(0,base))),
                   mp(3,previous(0,base)), mp(4,previous(1,base)), mp(5,next(0,previous(1,base)))};

        sites.emplace(std::make_pair(base, nb0));
    }

    //Directions are set to produce the correct Windingnumbers in crystal dir 0 and 1
    std::array<Direction,zcmax> directions() {
        return {
            (Direction() <<   1,  0).finished(),
            (Direction() <<   0,  1).finished(),
            (Direction() <<  -1,  1).finished(),
            (Direction() <<  -1,  0).finished(),
            (Direction() <<   0, -1).finished(),
            (Direction() <<   1, -1).finished(),
        };
    }

        void print_name(std::ostream& os) const {
        os << "# lattice: Name                                     : triangular\n";
    };
};

