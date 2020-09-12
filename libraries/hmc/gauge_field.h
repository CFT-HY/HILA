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
    using basetype = typename sun::base_type;
    static constexpr int N = sun::size;
    field<sun> gauge[NDIM];
    field<sun> momentum[NDIM];

    virtual void refresh(){}
    virtual void set_unity(){}
    virtual void random(){}


    virtual void add_momentum(field<squarematrix<N,cmplx<basetype>>> *force){}
    virtual void draw_momentum(){}
    virtual void zero_momentum(){}
    virtual void backup(){}
    virtual void restore_backup(){}
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
  down_staple[ALL] = U2[dir2][X+dir1].conjugate()
                   * U1[dir1][X].conjugate()
                   * U2[dir2][X];
  // Forward staple
  staple_sum[ALL]  = staple_sum[X]
                   + U2[dir2][X+dir1]
                   * U1[dir1][X+dir2].conjugate()
                   * U2[dir2][X].conjugate();
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
    down_staple[ALL] = U[dir2][X+dir].conjugate()
                     * U[dir][X].conjugate()
                     * U[dir2][X];
    // Forward staple
    staple_sum[ALL]  = staple_sum[X]
                     + U[dir2][X+dir]
                     * U[dir][X+dir2].conjugate()
                     * U[dir2][X].conjugate();
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
           * U[dir1][X+dir2].conjugate()
           * U[dir2][X].conjugate();
      Plaq += 1-temp.trace().re/N;
    }
  }
  return Plaq;
}

template<int N, typename radix>
double plaquette_sum(field<squarematrix<N,radix>> *U){
  double Plaq=0;
  foralldir(dir1) foralldir(dir2) if(dir2 < dir1){
    onsites(ALL){
      element<SU<N,radix>> temp;
      temp = U[dir1][X] * U[dir2][X+dir1]
           * U[dir1][X+dir2].conjugate()
           * U[dir2][X].conjugate();
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
  using gauge_type = matrix;
  using fund_type = matrix;
  using basetype = typename matrix::base_type;
  static constexpr int N = matrix::size;
  field<matrix> gauge_backup[NDIM];

  // Set the gauge field to unity
  void set_unity(){
    foralldir(dir){
      onsites(ALL){
        this->gauge[dir][X] = 1;
      }
    }
  }

  void random(){
    foralldir(dir){
      onsites(ALL){
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
  /// mometum
  void add_momentum(field<squarematrix<N,cmplx<basetype>>> *force){
    foralldir(dir){
      onsites(ALL){
        force[dir][X] = this->gauge[dir][X]*force[dir][X];
        project_antihermitean(force[dir][X]);
        this->momentum[dir][X] = this->momentum[dir][X] + force[dir][X];
      }
    }
  }

  // Make a copy of fields updated in a trajectory
  void backup(){
    foralldir(dir) gauge_backup[dir] = this->gauge[dir];
  }

  // Restore the previous backup
  void restore_backup(){
    foralldir(dir) this->gauge[dir] = gauge_backup[dir];
  }

  // Read the gauge field from a file
  void read_file(std::string filename){
    std::ifstream inputfile;
    inputfile.open(filename, std::ios::in | std::ios::binary);
    foralldir(dir){
      read_fields(inputfile, this->gauge[dir]);
    }
    inputfile.close();
  }

  // Write the gauge field to a file
  void write_file(std::string filename){
    std::ofstream outputfile;
    outputfile.open(filename, std::ios::out | std::ios::trunc | std::ios::binary);
    foralldir(dir){
      write_fields(outputfile, this->gauge[dir]);
    }
    outputfile.close();
  }

  
  // Simple measurables that only depend on the gauge field
  double plaquette(){
    return plaquette_sum(this->gauge)/(lattice->volume()*NDIM*(NDIM-1)/2);
  }

  double polyakov(int dir){
    return polyakov_loop(dir, this->gauge);
  }


  field<gauge_type> & get_momentum(int dir){
    return this->momentum[dir];
  }
  field<gauge_type> & get_gauge(int dir){
    return this->gauge[dir];
  }
};





template<typename repr>
class represented_gauge_field : public gauge_field_base<repr> {
  public: 
  using gauge_type = repr;
  using fund_type = typename repr::sun;
  using basetype = typename repr::base_type;
  static constexpr int Nf = fund_type::size;
  static constexpr int N = repr::size;
  gauge_field<fund_type> &fundamental;

  represented_gauge_field(gauge_field<fund_type>  &f) : fundamental(f){
    gauge_field_base<repr>();
  }
  represented_gauge_field(represented_gauge_field &r)
    : fundamental(r.fundamental){
      gauge_field_base<repr>();
    }


  // Represent the fields
  void refresh(){
    foralldir(dir){
      this->gauge[dir].check_alloc();
      onsites(ALL){
        if(disable_avx[X]==0){};
        this->gauge[dir][X].represent(fundamental.gauge[dir][X]);
      }
    }
  }

  void set_unity(){
    fundamental.set_unity();
    refresh();
  }

  void random(){
    fundamental.random();
    refresh();
  }


  void add_momentum(field<squarematrix<N,cmplx<basetype>>> (&force)[NDIM]){
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
  void zero_momentum(){
    fundamental.zero_momentum();
  }

  // Make a backup of the fundamental gauge field
  // Again, this may get called twice.
  void backup(){
    fundamental.backup();
  }

  // Restore the previous backup
  void restore_backup(){
    fundamental.restore_backup();
  }


  field<fund_type> & get_momentum(int dir){
    return fundamental.get_momentum(dir);
  }
  field<fund_type> & get_gauge(int dir){
    return fundamental.get_gauge(dir);
  }
};


/* Shortcuts for represented gauge fields */
template<int N, typename radix>
using symmetric_gauge_field = represented_gauge_field<symmetric<N,radix>>;
template<int N, typename radix>
using antisymmetric_gauge_field = represented_gauge_field<antisymmetric<N,radix>>;
template<int N, typename radix>
using adjoint_gauge_field = represented_gauge_field<adjoint<N,radix>>;

















/*******************
 * Action terms 
 *******************/

/// The Wilson plaquette action of a gauge field
/// The gauge action is a bit special, other action terms
/// only contain the force step for the MC integrator.
/// The gauge action also contains an update step that
/// updates the gauge field. This is the lowest level of the
/// integrator
template<typename gauge_field>
class gauge_action : public action_base, public integrator_base {
  public:
    using gauge_field_type = gauge_field;
    using gauge_mat = typename gauge_field::gauge_type;
    static constexpr int N = gauge_mat::size;
    using momtype = squarematrix<N, cmplx<typename gauge_mat::base_type>>;

    gauge_field &gauge;
    field<gauge_mat> gauge_copy[NDIM];
    double beta;

    gauge_action(gauge_field &g, double b) 
    : gauge(g), beta(b){}

    gauge_action(gauge_action &ga)
    : gauge(ga.gauge), beta(ga.beta) {}

    //The gauge action
    double action(){
      double Sa = 0;
      double Sg = beta*plaquette_sum(gauge.gauge);
      foralldir(dir) {
        onsites(ALL){
          Sa += gauge.momentum[dir][X].algebra_norm();
        }
      }
      
      return Sg+Sa;
    }

    /// Gaussian random momentum for each element
    void draw_gaussian_fields(){
      gauge.draw_momentum();
    }

    // Update the momentum with the gauge field
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

    // Draw a random gauge field
    void random(){
      foralldir(dir){
        onsites(ALL){
          gauge.gauge[dir][X].random();
        }
      }
    }


    /* The following are functions an integrator must have */

    // Make a copy of fields updated in a trajectory
    void backup_fields(){
      gauge.backup();
    }

    // Restore the previous backup
    void restore_backup(){
      gauge.restore_backup();
    }

    // Momentum step and step are required for an integrator.
    // A gauge action is also the lowest level of an
    // integrator hierarchy.
    
    // Update the gauge field with momentum
    void momentum_step(double eps){
      gauge.gauge_update(eps);
    }

};









#endif