/* $Id: worm.hpp,v 1.1 2006/09/09 9:21:44 pollet Exp $ */

#ifndef worm_Element_HPP
#define worm_Element_HPP

#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <list>
#include <stdint.h>

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>

#include "lattice.hpp"
#include "model.hpp"

class Element;  // forward declaration

#if defined DSTRUC_AVL
#include "avl_tree.hpp"

class AVL_diagram : public AVL::tree<Element> {
public:

    
    static void print_name(std::ostream& os) {
        os << "# config : DSTRUC                                   : AVL\n";

    }

    using AVL::tree<Element>::insert;
    inline iterator insert(const_iterator pos, const Element& v) {
        return insert(v);
    }
    inline iterator insert(const_iterator pos, Element&& v) {
        return insert(std::move(v));
    }
};

typedef AVL_diagram Diagram_type;


#elif defined DSTRUC_LIST
#include <list>

class list_diagram : public std::list<Element > {
public:

    
    static void print_name(std::ostream& os) {
        os << "# config : DSTRUC                                   : LIST\n";
    }

    iterator at(size_t i) {
        if (i > size())
            return end();
        iterator it = begin();
        for (; i > 0; --i)
            ++it;
        return it;
    }
    const_iterator at(size_t i) const {
        if (i > size())
            return end();
        const_iterator it = begin();
        for (; i > 0; --i)
            ++it;
        return it;
    }
};

typedef list_diagram Diagram_type;

#elif defined DSTRUC_LIST_STACK
#include "list_stack.hpp"

class list_stack_diagram : public list_stack<Element> {
public:

    static void print_name(std::ostream& os) {
        os << "# config : DSTRUC                                   : LIST_STACK\n";
    }
};

typedef list_stack_diagram Diagram_type;
#endif



class Element
{
 public:
    using StateType = model::StateType;
    using SiteIndex = LATTICE::Base::SiteIndex;
    static const size_t ZC = LATTICE::zcmax;

  Element() {};
  Element(const StateType n0, const StateType n1, const SiteIndex s0, 
          const double t, const int color, const Diagram_type::iterator v[]) {
    mTime = t;
    mBefore = n0;
    mAfter = n1;
    mLink = s0;
    mColor = color;
    for (std::size_t i = 0; i < ZC; i++) mAssoc[i] = v[i];
  }
  Element(const StateType n0, const StateType n1, const SiteIndex s0, const double t, const int color) {
    mTime = t;
    mBefore = n0;
    mAfter = n1;
    mLink = s0;
    mColor = color;
  }

  void init(const StateType n0, const StateType n1, const SiteIndex s0, const double t, const int col) {
    mTime = t;
    mBefore = n0;
    mAfter = n1;
    mLink = s0;
    mColor = col;
  }
  void init(const StateType n0, const StateType n1, const SiteIndex s0, 
            const double t, const int col, Diagram_type::iterator v[] ) {
    mTime = t;
    mBefore = n0;
    mAfter = n1;
    mLink = s0;
    mColor = col;
    for (std::size_t i = 0; i < ZC; i++) mAssoc[i] = v[i];
  }
  ~Element() {}
  Element(const Element& src) 
    : mTime(src.mTime), mBefore(src.mBefore), mAfter(src.mAfter), mLink(src.mLink),  mColor(src.mColor) {
    for (std::size_t i = 0; i < ZC; i++) mAssoc[i] = src.mAssoc[i];
  }
  bool operator<(const Element& rhs) const {
      // this is an extra function needed for the AVL tree data structure to work, which only compares elements based on <
      // we can use the time function for that,
      // however we have equal times occurruing which are differentiated by the worm_at_stop variable (and can be infinite greater or smaller), which is not part of the element data structure.
      // my fix is to look at the occupancy in that case, which should be continuous. So if the "after" occupation of the current element coincides with the "before" of the rhs element,
      // then we know that the current element is before the rhs element and vice verso. 
      // the fix only works for worm variables (and not for dummies) but worm_at_stop = +1 or -1 and on-site, equal time variables only makes sense for worms anyway.
      if (time() == rhs.time()) {
        return after() == rhs.before();
      }
      else {
        return time() < rhs.time();
      }
    } 
  Element& operator=(const Element rhs) {
    if (this == &rhs) return (*this);
    mTime = rhs.mTime;
    mBefore = rhs.mBefore;
    mAfter = rhs.mAfter;
    mLink = rhs.mLink;
    mColor = rhs.mColor;
    for (std::size_t i = 0; i < ZC; i++) mAssoc[i] = rhs.mAssoc[i];
    return *this;
  }
  friend std::ostream &operator<<( std::ostream &os, const Element& rhs) {
    os  << rhs.color() << "\t"
        << rhs.time() << "\t"
        << rhs.link() << "\t"
        << rhs.before() << "\t"
        << rhs.after();
    return (os);
  }
  friend std::istream& operator>>(std::istream& is, Element& e) {
    is >> e.mColor>> e.mTime >> e.mLink >> e.mBefore >> e.mAfter;
    return (is);
  }
  friend bool operator== (const Element& lhs, const Element& rhs) {
    return  ( (lhs.mColor == rhs.mColor) && (lhs.mLink == rhs.mLink) && (lhs.mTime == rhs.mTime) && (lhs.mBefore == rhs.mBefore) && (lhs.mAfter == rhs.mAfter) );
  }


  double time() const {return (mTime);}
  StateType before() const {return (mBefore);}
  StateType after() const {return (mAfter);}
  int color() const {return (mColor);}
  SiteIndex link() const { return mLink;}
  
  void time(const double t) {mTime=t;};
  void before(const StateType n) {mBefore = n;};
  void after(const StateType n) {mAfter = n;};
  void link(const SiteIndex l)  { mLink = l;}
  void color(const int col) {mColor = col;}
  void assoc(const Diagram_type::iterator   v[ZC]) {
    for (size_t i = 0; i < ZC; i++) mAssoc[i] = v[i];
  }

  void print() const {
    std::cout << "\n" << this << "\tcolor : " << mColor << "\ttime : " 
              << mTime << "\tlink : " << mLink << "\tbefore : " << mBefore << "\tafter : " << mAfter;
    for (size_t j=0; j < ZC; j++) {
        if(mAssoc[j]->mLink == -1)
            std::cout << "\t" << "none";
        else
            std::cout << "\t" << mAssoc[j]->time() * mAssoc[j]->color();
    }
  }


  Diagram_type::iterator& get_assoc(SiteIndex s) { return mAssoc[s];}
  double get_assoc_time(SiteIndex s) const {return mAssoc[s]->time();}
  void set_assoc(const SiteIndex j, const Diagram_type::iterator v) {
    mAssoc[j] = v;
  }

 private:
  double mTime;  // time of the interaction
  StateType mBefore;     // occupation on the site before the interaction
  StateType mAfter;      // occupation on the site after the interaction
  SiteIndex mLink;   // site to which element is linked
  Diagram_type::iterator mAssoc[ZC];   // associations; i.e. iterator to the kink on nb sites that are equal or just greater in time
  int mColor;      // (1) red boson interaction (2) blue boson interaction (-1) red boson worm (-2) blue boson worm (0) dummy etc
};






#endif
