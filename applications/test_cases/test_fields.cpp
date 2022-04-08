#include "hila.h"
#include "test.h"

/////////////////////
/// test_case 2
/// Coverage:
/// - directions, onsites and foralldir environments
/// - operations between fields
/// - foralldir env inside onsites
/// - referring to an array of fields in a loop
/// - calling a function with const parameters
///   - requiring communication of a const field
/// - calling a function from a loop
/////////////////////

template <typename A, typename B, typename C>
void sum_test_function(A &a, const B &b, const C &c) {
    onsites(ALL) { a[X] = b[X] + c[X + e_x]; }
}

template <typename T> T test_template_function(T a) { return 2 * a; }

element<Complex<double>> test_nontemplate_function(element<Complex<double>> a) {
    element<Complex<double>> b = a;
    return 2 * a;
}

#include "dirac/staggered.h"

int main(int argc, char **argv) {

    // check that you can increment a Direction correctly
    Direction d = e_x;
    Direction d2 = (Direction)(NDIRS - 2);
#if NDIM > 1
    d = next_direction(d);
    d2 = next_direction(d2);
    assert(d == e_y);
    assert(e_x == 0);
    assert(d2 == -e_x);
#endif
    double sum = 0;

    test_setup(argc, argv);

    Field<Complex<double>> s1, s2, s3;
    Field<Complex<double>> s4[3];

    // Test field assingment
    s1 = 0.0;
    s2 = 1.0;
    s3 = 1.0;

    // Test sum and move constructor
    s1 = s2 + s3;

    sum = 0;
    onsites(ALL) { sum += s1[X].re; }
    assert(sum == 2 * (double)lattice->volume() && "onsites reduction");
    s1 = 0;
    s2 = 0;
    s3 = 0;
    sum = 0;

    // Test field-parity expressions
    s1[ALL] = 0.0;
    s2[EVEN] = 1.0;
    s3[ODD] = 1.0;

    s1[ALL] = s2[X] + s3[X];

    onsites(ALL) { sum += s1[X].re; }
    assert(sum == (double)lattice->volume() && "test setting field with parity");

    // Test communication functions
    s1[ALL] = 0;

    assert(s1.is_allocated());
    assert(s1.is_initialized(EVEN));
    assert(s1.is_initialized(ODD));

    s1.mark_changed(ALL);

    // Check initial state of communication functions
    foralldir(d) {
        assert(!s1.is_gathered(d, EVEN));
        assert(!s1.is_gathered(d, ODD));
        assert(!s1.is_gathered(d, ALL));
        assert(!s1.is_move_started(d, EVEN));
        assert(!s1.is_move_started(d, ODD));
        assert(!s1.is_move_started(d, ALL));
        assert(s1.move_not_done(d, EVEN) && "move not done initially");
        assert(s1.move_not_done(d, ODD) && "move not done initially");
        assert(s1.move_not_done(d, ALL) && "move not done initially");
    }

    // Test marking move started and gathered
    foralldir(d) {
        s1.mark_move_started(d, EVEN);
        assert(s1.is_move_started(d, EVEN));
        assert(!s1.is_gathered(d, EVEN));
        assert(!s1.move_not_done(d, EVEN) && "move not done after starting");
        assert(!s1.is_gathered(d, ODD));
        assert(!s1.is_gathered(d, ALL));
        assert(!s1.is_move_started(d, ODD));
        assert(!s1.is_move_started(d, ALL));
        assert(s1.move_not_done(d, ODD));
        assert(s1.move_not_done(d, ALL));

        s1.mark_gathered(d, EVEN);
        assert(s1.is_gathered(d, EVEN));
        assert(!s1.is_move_started(d, EVEN));
        assert(!s1.move_not_done(d, EVEN));

        s1.mark_changed(ALL);

        s1.mark_move_started(d, ODD);
        assert(s1.is_move_started(d, ODD));
        assert(!s1.is_gathered(d, ODD));
        assert(!s1.move_not_done(d, ODD) && "move not done after starting");

        s1.mark_gathered(d, ODD);
        assert(s1.is_gathered(d, ODD));
        assert(!s1.is_move_started(d, ODD));
        assert(!s1.move_not_done(d, ODD));

        s1.mark_changed(ALL);

        s1.mark_move_started(d, ALL);
        assert(s1.is_move_started(d, ALL));
        assert(!s1.is_gathered(d, ALL));
        assert(!s1.move_not_done(d, ALL) && "move not done after starting");

        s1.mark_gathered(d, ALL);
        assert(s1.is_gathered(d, ALL));
        assert(!s1.is_move_started(d, ALL));
        assert(!s1.move_not_done(d, ALL));

        s1.mark_changed(ALL);
    }

    // Try setting an element on node 0
    CoordinateVector coord;
    foralldir(d) { coord[d] = 0; }
    s1.set_element(Complex<double>(1), coord);
    Complex<double> elem = s1.get_element(coord);
    assert(elem.re == 1 && elem.im == 0 && "set_element");

    // Now try setting on a different node, if the lattice is split
    foralldir(d) { coord[d] = nd[d] - 1; }
    s1.set_element(Complex<double>(1), coord);
    elem = s1.get_element(coord);
    assert(elem.re == 1 && elem.im == 0 && "set_element on other node");

    // Now try actually moving the data
    foralldir(d) { coord[d] = 0; }
    foralldir(d) {
        // Move data up in Direction d
        s2[ALL] = s1[X - d];

        // Should still be on this node
        CoordinateVector coord2 = coord;
        coord2[d] += 1;
        Complex<double> moved = s2.get_element(coord2);
        assert(elem.re == 1 && elem.im == 0);

        // Move data down
        s2[ALL] = s1[X + d];

        // Now it may be on a different node
        coord2 = coord;
        coord2[d] = (coord[d] - 1 + nd[d]) % nd[d];
        moved = s2.get_element(coord2);
        if (elem.re != 1 || elem.im != 0) {
            output0 << "Problem in communicating to Direction " << d << "\n";
            output0 << "Received " << moved << "\n";
            assert(elem.re == 1 && elem.im == 0);
        }
    }

    // Check boundary conditions by moving the data
    {
        Field<double> s1, s2;
        CoordinateVector coord1 = 0;
        CoordinateVector coord2 = 0;
        coord1[0] = (nd[0] - 1) % nd[0];

        // Move data up accross the X boundary
        s1 = 0;
        s2 = 0;
        s1.set_element(1, coord1);
        s2[ALL] = s1[X - Direction(0)];
        double moved_element = s2.get_element(coord2);
        assert(moved_element == 1 && "moved up");

        // The other way
        s1 = 0;
        s2 = 0;
        s1.set_element(1, coord2);
        s2[ALL] = s1[X + Direction(0)];
        moved_element = s2.get_element(coord1);
        assert(moved_element == 1 && "moved down");

        s1.set_boundary_condition(Direction(0), BoundaryCondition::ANTIPERIODIC);
        s2.set_boundary_condition(Direction(0), BoundaryCondition::ANTIPERIODIC);

        // Now try antiperiodic boundaries
        s1 = 0;
        s2 = 0;
        s1.set_element(1, coord1);
        s2[ALL] = s1[X - Direction(0)];
        moved_element = s2.get_element(coord2);
        assert(moved_element == -1 && "moved up antiperiodic");

        // The other way
        s1 = 0;
        s2 = 0;
        s1.set_element(1, coord2);
        s2[ODD] = s1[X + Direction(0)];
        moved_element = s2.get_element(coord1);
        assert(moved_element == -1 && "moved down antiperiodic");

        s1 = 0;
        s1[EVEN] = s2[X - Direction(0)];
        moved_element = s1.get_element(coord2);
        assert(moved_element == 1 && "moved back antiperiodic");
    }

    // Communication and copy test with full field
    foralldir(d) {
        s1 = 1.0;
        s2 = 1.0;
        s3 = 1.0;
        double sum = 0;
        double sum2 = 0;
        onsites(EVEN) {
            double a = s1[X + d].re;
            double b = s2[X].re;
            sum += a - b;
        }
        assert(sum == 0 && "Test communicating a filled field");
    }

    // Test interleaving by calculating a sum of coordinates
    // when a field is communicated
    foralldir(d) {
        s1 = 1;
        s2 = 1;
        sum = 0;
        onsites(ALL) {
            sum += X.coordinates()[0];
            s1[X];
        }

        double result = lattice->volume() * (nd[0] - 1) / 2;
        assert(sum == result && "Reproduce write problem 1");

        s1 = 1;
        s2 = 1;
        sum = 0;
        onsites(ALL) {
            sum += X.coordinates()[0];
            s1[X + d];
        }

        assert(sum == result && "Reproduce write problem 2");
    }

    // Check that communicated and changed fields get marked
    // correctly
    CoordinateVector c(0);
    foralldir(dir) {
        for (int i = 0; i < 2; i++) {
            s1[EVEN] = 0;
            s1[c] = i;

            s2[ALL] = 0;
            s2[ODD] = s2[X] - s1[X + dir];
            s2[EVEN] = s2[X] + s2[X - dir];

            double r = s2[c].re;
            sum = (r + i) * (r + i);
            if (hila::myrank() == 0) {
                assert(sum < 1e-8 && "Even-odd comm test");
            }
        }
    }

    // Test starting communication manually

    s1[EVEN] = 1.0;
    s2[EVEN] = 1.0;
    s2[ODD] = -s1[X + e_x];
    s2.start_gather(e_x, ODD);

    sum = 0;
    onsites(ALL) { sum += s2[X].re; }
    assert(sum == 0);

    // Test referring to an array of fields

    s4[0] = s1;
    s4[1] = s1;

    s4[2][ALL] = s4[0][X] - s4[1][X];

    sum = 0;
    onsites(ALL) { sum += (s4[2][X] * s4[2][X]).re; }
    assert(sum == 0);

    // Test function call outside loop
    s1[ALL] = 0.0;
    s2[ALL] = 1.0;
    sum_test_function(s3, s1, s2); // s3 = s1 + s2
    onsites(ALL) {
        element<double> diff = s3[X].re - 1.0;
        sum += diff * diff;
    }
    assert(sum == 0);

    // Test function calls in loop
    s1[ALL] = 1.0;
    s2[ALL] = 1.0;
    onsites(ALL) {
        s1[X] = test_template_function(s1[X]);
        s2[X] = test_nontemplate_function(s2[X]);
    }
    onsites(ALL) {
        element<double> diff1 = s1[X].re - 2.0;
        element<double> diff2 = s2[X].re - 2.0;
        sum += diff1 * diff1 + diff2 * diff2;
    }
    assert(sum == 0);

    // Test array reduction
    Field<double> dfield;
    dfield[ALL] = 1;

#if NDIM == 4
    VectorReduction<double> arraysum(nd[e_t]);
    arraysum = 0;

    onsites(ALL) {
        CoordinateVector l = X.coordinates();
        int t = l[e_t];

        arraysum[t] += dfield[X];
    }

    for (int t = 0; t < nd[e_t]; t++) {
        assert(arraysum[t] == nd[e_x] * nd[e_y] * nd[e_z]);
    }
#endif

    hila::finishrun();
}
