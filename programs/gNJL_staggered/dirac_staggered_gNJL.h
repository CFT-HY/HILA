

template<int N>
class dirac_staggered_gNJL {
  private:
    dirac_staggered<SU_vector<N>, SU<N>> D;
    field<double> &sigma, &pi;

  public:

    using vector_type = SU_vector<N>;

    // Constructor: initialize mass, gauge and eta
    dirac_staggered_gNJL(dirac_staggered_gNJL &d) : D(d.D), sigma(d.sigma), pi(d.pi) {}
  
    // Constructor: initialize mass, gauge and eta
    dirac_staggered_gNJL(double m, field<SU<N>> (&g)[NDIM], field<double> &s, field<double> &p) :
    D(m, g), sigma(s), pi(p) {}


    // Applies the operator to in
    void apply( const field<SU_vector<N>> & in, field<SU_vector<N>> & out){
      D.apply(in, out);

      //out[ALL]  = out[X] + sigma[X]*in[X];
      //out[EVEN] = out[X] + cmplx(0,1)*pi[X]*in[X];
      //out[ODD]  = out[X] - cmplx(0,1)*pi[X]*in[X];
    }

    // Applies the conjugate of the operator
    void dagger( const field<SU_vector<N>> & in, field<SU_vector<N>> & out){
      D.dagger(in, out);

      //out[ALL]  = out[X] + sigma[X]*in[X];
      //out[EVEN] = out[X] - cmplx(0,1)*pi[X]*in[X];
      //out[ODD]  = out[X] + cmplx(0,1)*pi[X]*in[X];
    }

    // Applies the derivative of the Dirac operator with respect
    // to the gauge field
    void force( const field<SU_vector<N>> & psi, const field<SU_vector<N>> & chi, field<SU<N>> (&force)[NDIM]){
      D.force(psi, chi, force);
    }
};

