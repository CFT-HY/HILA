#include "test.h"
#include "plumbing/fft.h"
//#include "plumbing/FFT.h"

constexpr double Pi = 3.14159265358979;

int main(int argc, char **argv) {

    // using T = Matrix<2,2,Complex<double>>;
    using T = Complex<double>;

    test_setup(argc, argv);

    Field<T> f, f2, p, p2;
    double sum = 0;

    for (int iter = 0; iter < 3; iter++) {

        // Start with unit field
        f = 1;

        // After one FFT the field is 0 except at coord 0
        p2 = 0;

        p2[{0,0,0,0}] = lattice.volume();

        hila::out0 << "Start fft\n";

        FFT_field(f, p);

        sum = 0;
        onsites(ALL) { sum += (p[X] - p2[X]).squarenorm(); }
        hila::out0 << "Sum " << sum << '\n';
        if (fabs(sum) < 1e-10) {
            hila::out0 << "First FFT passed" << std::endl;
        } else {
            hila::out0 << "First FFT FAILED; limit 1e-10" << std::endl;
        }

        // After two applications the field should be back to a constant * volume
        f2[ALL] = lattice.volume();

        FFT_field(p, f, fft_direction::back);

        sum = 0;
        double tnorm = 0;
        onsites(ALL) {
            sum += (f[X] - f2[X]).squarenorm();
            tnorm += f[X].squarenorm();
        }
        hila::out0 << "Norm " << sum / tnorm << '\n';
        if (fabs(sum / tnorm) < 1e-10) {
            hila::out0 << "2nd FFT passed\n";
        } else {
            hila::out0 << "Second FFT FAILED, limit 1e-10" << std::endl;
        }


        onsites(ALL) {
            double d = X.coordinate(e_x)*2.0*Pi/lattice.size(e_x);
            f[X] = Complex<double>(cos(d),sin(d));
        }

        FFT_field(f, p);

        p2 = 0;
        p2[{1,0,0,0}] = lattice.volume();


        sum = 0;
        onsites(ALL) { sum += (p[X] - p2[X]).squarenorm(); }
        hila::out0 << "Wave sum " << sum << '\n';
        if (fabs(sum) > 1e-10) {
            hila::out0 << "Wave FFT FAILED" << std::endl;
        } 
    }


    // write_fields("test_config_filename", p, f);
    // read_fields("test_config_filename", p, f2);

    // sum=0;
    // onsites(ALL) {
    //  sum += (f2[X]-f[X]).squarenorm();
    //}

    // assert(sum==0 && "Write and read field");

    hila::finishrun();
}
