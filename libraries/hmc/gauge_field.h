#ifndef GAUGE_FIELD_H
#define GAUGE_FIELD_H


#include "datatypes/sun.h"
#include "datatypes/representations.h"
#include "integrator.h"


/// Define a standard base gauge class. Gauge field types (represented, 
/// smeared, etc) inherit from this
template<typename sun>
class gauge_field_base{
  public:
    /// The base type of the matrix (double, int,...)
    using basetype = typename sun::base_type;
    /// The matrix type
    using gauge_type = sun;
    /// The size of the matrix
    static constexpr int N = sun::size;

    /// A matrix field for each direction
    field<sun> gauge[NDIM];
    /// Also create a momentum field. This is only
    /// allocated if necessary
    field<sun> momentum[NDIM];

    /// Recalculate represented or smeared field
    virtual void refresh(){}
    /// Set the field to unity
    virtual void set_unity(){}
    /// Draw a random gauge field
    virtual void random(){}

    /// Update the momentum by given force
    virtual void add_momentum(field<SquareMatrix<N,cmplx<basetype>>> *force){}
    /// Draw gaussian random momentum
    virtual void draw_momentum(){}
    /// Set the momentum to zero
    virtual void zero_momentum(){}
    /// Make a copy of the gauge field for HMC
    virtual void backup(){}
    /// Restore the gauge field from the backup.
    /// Used when an HMC trajectory is rejected.
    virtual void restore_backup(){}

    /// If the base type is double, this will be the 
    /// corresponding floating point type.
    using gauge_type_flt = sun;
};

/// Specialize double precision SU(N) matrix types
template<template <int,typename> class M, int N>
class gauge_field_base<M<N,double>>{
  public:
    /// The basetype is double
    using basetype = double;
    /// The matrix type
    using gauge_type =  M<N,double>;

    /// A matrix field for each direction
    field<gauge_type> gauge[NDIM];
    /// Also create a momentum field. This is only
    /// allocated if necessary
    field<gauge_type> momentum[NDIM];

    /// Recalculate represented or smeared field
    virtual void refresh(){}
    /// Set the field to unity
    virtual void set_unity(){}
    /// Draw a random gauge field
    virtual void random(){}

    /// Update the momentum by given force
    virtual void add_momentum(field<SquareMatrix<N,cmplx<basetype>>> *force){}
    /// Draw gaussian random momentum
    virtual void draw_momentum(){}
    /// Set the momentum to zero
    virtual void zero_momentum(){}
    /// Make a copy of the gauge field for HMC
    virtual void backup(){}
    /// Restore the gauge field from the backup.
    /// Used when an HMC trajectory is rejected.
    virtual void restore_backup(){}

    /// This is the single precision type
    using gauge_type_flt = M<N,float>;
    /// Return a single precision copy of the gauge field
    gauge_field_base<M<N,float>> get_single_precision(){
      gauge_field_base<M<N,float>> gauge_flt;
      foralldir(dir){
        gauge_flt.gauge[dir] = gauge[dir];
      }
      return gauge_flt;
    }
};




/// Calculate the Polyakov loop for a given gauge field.
template<int N>
double polyakov_loop(direction dir, field<SU<N>> (&gauge)[NDIM]){
  // This implementation uses the onsites() to cycle through the
  // NDIM-1 dimensional planes. This is probably not the most
  // efficient implementation.
  coordinate_vector vol = lattice->size();
  field<SU<N>> polyakov; polyakov[ALL] = 1;
  for(int t=0; t<vol[dir]; t++){
    onsites(ALL){
      if(X.coordinates()[dir]==(t+1)%vol[dir]){
        polyakov[X] = polyakov[X] * gauge[dir][X-dir];
      }
    }
  }

  double poly=0;
  onsites(ALL){
    if(X.coordinates()[dir]==0){
      poly += polyakov[X].trace().re;
    }
  }

  double v3 = lattice->volume()/vol[dir];
  return poly/(N*v3);
}





/// Calculate the sum of staples in direction dir2 
/// connected to links in direction dir1
/// This version takes two different fields for the
/// different directions and is necessary for HEX
/// smearing
template<typename SUN>
field<SUN> calc_staples(field<SUN> *U1, field<SUN> *U2, direction dir1, direction dir2)
{
  field<SUN> staple_sum;
  static field<SUN> down_staple;
  staple_sum[ALL] = 0;
  //Calculate the down side staple.
  //This will be communicated up.
  down_staple[ALL] = U2[dir2][X+dir1].conj()
                   * U1[dir1][X].conj()
                   * U2[dir2][X];
  // Forward staple
  staple_sum[ALL]  = staple_sum[X]
                   + U2[dir2][X+dir1]
                   * U1[dir1][X+dir2].conj()
                   * U2[dir2][X].conj();
  // Add the down staple
  staple_sum[ALL] = staple_sum[X] + down_staple[X-dir2];
  return staple_sum;
}


/// Calculate the sum of staples connected to links in direction dir
template<typename SUN>
field<SUN> calc_staples(field<SUN> *U, direction dir)
{
  field<SUN> staple_sum;
  static field<SUN> down_staple;
  staple_sum[ALL] = 0;
  foralldir(dir2) if(dir2!=dir) {
    //Calculate the down side staple.
    //This will be communicated up.
    down_staple[ALL] = U[dir2][X+dir].conj()
                     * U[dir][X].conj()
                     * U[dir2][X];
    // Forward staple
    staple_sum[ALL]  = staple_sum[X]
                     + U[dir2][X+dir]
                     * U[dir][X+dir2].conj()
                     * U[dir2][X].conj();
    // Add the down staple
    staple_sum[ALL] = staple_sum[X] + down_staple[X-dir2];
  }
  return staple_sum;
}




/// Measure the plaquette
template<int N, typename radix>
double plaquette_sum(field<SU<N,radix>> *U){
  double Plaq=0;
  foralldir(dir1) foralldir(dir2) if(dir2 < dir1){
    onsites(ALL){
      element<SU<N,radix>> temp;
      temp = U[dir1][X] * U[dir2][X+dir1]
           * U[dir1][X+dir2].conj()
           * U[dir2][X].conj();
      Plaq += 1-temp.trace().re/N;
    }
  }
  return Plaq;
}


/// The plaquette measurement for square matrices
template<int N, typename radix>
double plaquette_sum(field<Matrix<N,N,radix>> *U){
  double Plaq=0;
  foralldir(dir1) foralldir(dir2) if(dir2 < dir1){
    onsites(ALL){
      element<SU<N,radix>> temp;
      temp = U[dir1][X] * U[dir2][X+dir1]
           * U[dir1][X+dir2].conj()
           * U[dir2][X].conj();
      Plaq += 1-temp.trace()/N;
    }
  }
  return Plaq;
}













/// A gauge field contains a SU(N) matrix in each
/// direction for the gauge field and for the momentum.
/// Defines methods for HMC to update the field and the 
/// momentum.
template<typename matrix>
class gauge_field : public gauge_field_base<matrix> {
  public:
  /// The matrix type
  using gauge_type = matrix;
  /// The fundamental gauge type. In the standard case
  /// it is the same as the matrix type
  using fund_type = matrix;
  /// The base type (double, float, int...)
  using basetype = typename matrix::base_type;
  /// The size of the matrix
  static constexpr int N = matrix::size;
  /// Storage for a backup of the gauge field
  field<matrix> gauge_backup[NDIM];

  /// Set the gauge field to unity
  void set_unity(){
    foralldir(dir){
      onsites(ALL){
        this->gauge[dir][X] = 1;
      }
    }
  }

  /// Draw a random gauge field
  void random(){
    foralldir(dir){
      onsites(ALL){
        if(disable_avx[X]==0){};
        this->gauge[dir][X].random();
      }
    }
  }

  /// Gaussian random momentum for each element
  void draw_momentum(){
    foralldir(dir) {
      onsites(ALL){
        if(disable_avx[X]==0){};
        this->momentum[dir][X].gaussian_algebra();
      }
    }
  }

  /// Set the momentum to zero
  void zero_momentum(){
    foralldir(dir) {
      this->momentum[dir][ALL]=0;
    }
  }

  /// Update the gauge field with time step eps
  void gauge_update(double eps){
    foralldir(dir){
      onsites(ALL){
        element<matrix> momexp = eps*this->momentum[dir][X];
        momexp.exp();
        this->gauge[dir][X] = momexp*this->gauge[dir][X];
      }
    }
  }

  /// Project a force term to the algebra and add to the
  /// momentum
  void add_momentum(field<SquareMatrix<N,cmplx<basetype>>> *force){
    foralldir(dir){
      onsites(ALL){
        force[dir][X] = this->gauge[dir][X]*force[dir][X];
        project_antihermitean(force[dir][X]);
        this->momentum[dir][X] = this->momentum[dir][X] + force[dir][X];
      }
    }
  }

  /// Make a copy of fields updated in a trajectory
  void backup(){
    foralldir(dir) gauge_backup[dir] = this->gauge[dir];
  }

  /// Restore the previous backup
  void restore_backup(){
    foralldir(dir) this->gauge[dir] = gauge_backup[dir];
  }

  /// Read the gauge field from a file
  void read_file(std::string filename){
    std::ifstream inputfile;
    inputfile.open(filename, std::ios::in | std::ios::binary);
    foralldir(dir){
      read_fields(inputfile, this->gauge[dir]);
    }
    inputfile.close();
  }

  /// Write the gauge field to a file
  void write_file(std::string filename){
    std::ofstream outputfile;
    outputfile.open(filename, std::ios::out | std::ios::trunc | std::ios::binary);
    foralldir(dir){
      write_fields(outputfile, this->gauge[dir]);
    }
    outputfile.close();
  }
  
  
  /// Calculate the plaquette
  double plaquette(){
    return plaquette_sum(this->gauge)/(lattice->volume()*NDIM*(NDIM-1)/2);
  }

  /// Calculate the polyakov loop
  double polyakov(int dir){
    return polyakov_loop(dir, this->gauge);
  }

  /// Return a reference to the momentum field
  field<gauge_type> & get_momentum(int dir){
    return this->momentum[dir];
  }
  /// Return a reference to the gauge field
  field<gauge_type> & get_gauge(int dir){
    return this->gauge[dir];
  }
};




/// A gauge field, similar to the standard gauge_field class above,
/// but with the gauge field projected into a higher representation.
template<class repr>
class represented_gauge_field : public gauge_field_base<repr> {
  public: 
  /// The matrix type
  using gauge_type = repr;
  /// The type of a fundamental representation matrix
  using fund_type = typename repr::sun;
  /// The base type (double, float, int...)
  using basetype = typename repr::base_type;
  /// The size of the matrix
  static constexpr int Nf = fund_type::size;
  /// The size of the representation
  static constexpr int N = repr::size;
  /// Reference to the fundamental gauge field
  gauge_field<fund_type> &fundamental;

  /// Construct from a fundamental field
  represented_gauge_field(gauge_field<fund_type>  &f) : fundamental(f){
    gauge_field_base<repr>();
  }
  /// Copy constructor
  represented_gauge_field(represented_gauge_field &r)
    : fundamental(r.fundamental){
      gauge_field_base<repr>();
  }


  /// Represent the fields
  void refresh(){
    foralldir(dir){
      this->gauge[dir].check_alloc();
      onsites(ALL){
        if(disable_avx[X]==0){};
        this->gauge[dir][X].represent(fundamental.gauge[dir][X]);
      }
    }
  }

  /// Set the gauge field to unity. This will set the
  /// underlying fundamental field
  void set_unity(){
    fundamental.set_unity();
    refresh();
  }

  /// Draw a random gauge field. This will set the
  /// underlying fundamental field
  void random(){
    fundamental.random();
    refresh();
  }


  /// Project a force term to the algebra and add to the
  /// momentum
  void add_momentum(field<SquareMatrix<N,cmplx<basetype>>> (&force)[NDIM]){
    foralldir(dir){
      onsites(ALL){
        if(disable_avx[X]==0){};
        element<fund_type> fforce;
        fforce = repr::project_force(this->gauge[dir][X]*force[dir][X]);
        fundamental.momentum[dir][X] = fundamental.momentum[dir][X] + fforce;
      }
    }
  }

  /// This gets called if there is a represented gauge action term.
  /// If there is also a fundamental term, it gets called twice... 
  void draw_momentum(){
    fundamental.draw_momentum();
  }

  /// Set the momentum to zero
  void zero_momentum(){
    fundamental.zero_momentum();
  }

  /// Make a backup of the fundamental gauge field
  /// Again, this may get called twice.
  void backup(){
    fundamental.backup();
  }

  /// Restore the previous backup
  void restore_backup(){
    fundamental.restore_backup();
  }


  /// Return a reference to the momentum field
  field<fund_type> & get_momentum(int dir){
    return fundamental.get_momentum(dir);
  }
  /// Return a reference to the gauge field
  field<fund_type> & get_gauge(int dir){
    return fundamental.get_gauge(dir);
  }
};



/// Shortcuts for represented gauge fields
template<int N, typename radix>
using symmetric_gauge_field = represented_gauge_field<symmetric<N,radix>>;
/// Shortcuts for represented gauge fields
template<int N, typename radix>
using antisymmetric_gauge_field = represented_gauge_field<antisymmetric<N,radix>>;
/// Shortcuts for represented gauge fields
template<int N, typename radix>
using adjoint_gauge_field = represented_gauge_field<adjointRep<N,radix>>;

















/*******************
 * Action terms 
 *******************/

/// The action of the canonical momentum of a gauge field.
/// Momentum actions are a special case. It does not contain
/// a force_step()-function. It contains a step()-function,
/// which updates the momentum itself. It can be used as
/// the lowest level of an integrator.
template<typename gauge_field>
class gauge_momentum_action : public action_base, public integrator_base {
  public:
    /// The underlying gauge field type
    using gauge_field_type = gauge_field;
    /// The gauge matrix type
    using gauge_mat = typename gauge_field::gauge_type;
    /// The size of the gauge matrix
    static constexpr int N = gauge_mat::size;
    /// The type of the momentum field
    using momtype = SquareMatrix<N, cmplx<typename gauge_mat::base_type>>;

    /// A reference to the gauge field
    gauge_field &gauge;

    /// construct from a gauge field
    gauge_momentum_action(gauge_field &g) 
    : gauge(g){}
    /// construct a copy
    gauge_momentum_action(gauge_momentum_action &ga)
    : gauge(ga.gauge) {}

    /// The gauge action
    double action(){
      double Sa = 0;
      foralldir(dir) {
        onsites(ALL){
          Sa += gauge.momentum[dir][X].algebra_norm();
        }
      }
      return Sa;
    }

    /// Gaussian random momentum for each element
    void draw_gaussian_fields(){
      gauge.draw_momentum();
    }


    /* The following allow using a gauge action as the lowest level
       of an integrator. */
    /// Make a copy of fields updated in a trajectory
    void backup_fields(){
      gauge.backup();
    }

    /// Restore the previous backup
    void restore_backup(){
      gauge.restore_backup();
    }

    /// A momentum action is also the lowest level of an
    /// integrator hierarchy and needs to define the an step
    /// to update the gauge field using the momentum 
    
    /// Update the gauge field with momentum
    void step(double eps){
      gauge.gauge_update(eps);
    }

};



/// The Wilson plaquette action of a gauge field.
/// Action terms contain a force_step()-function, which
/// updates the momentum of the gauge field. To do this,
/// it needs to have a reference to the momentum field.
template<typename gauge_field>
class gauge_action : public action_base {
  public:
    /// The underlying gauge field type
    using gauge_field_type = gauge_field;
    /// The gauge matrix type
    using gauge_mat = typename gauge_field::gauge_type;
    /// The size of the gauge matrix
    static constexpr int N = gauge_mat::size;
    /// The type of the momentum field
    using momtype = SquareMatrix<N, cmplx<typename gauge_mat::base_type>>;

    /// A reference to the gauge field
    gauge_field &gauge;
    /// The coupling
    double beta;

    /// Construct out of a gauge field
    gauge_action(gauge_field &g, double b) 
    : gauge(g), beta(b){}
    /// Construct a copy
    gauge_action(gauge_action &ga)
    : gauge(ga.gauge), beta(ga.beta) {}

    /// The gauge action
    double action(){
      double Sg = beta*plaquette_sum(gauge.gauge);
      return Sg;
    }

    /// Update the momentum with the gauge field force
    void force_step(double eps){
      field<gauge_mat> staple;
      field<momtype> force[NDIM];
      foralldir(dir){
        staple = calc_staples(gauge.gauge, dir);
        onsites(ALL){
          force[dir][X] = (-beta*eps/N)*staple[X];
        }
      }
      gauge.add_momentum(force);
    }
};





#endif
