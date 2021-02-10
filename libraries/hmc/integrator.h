#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <sys/time.h>
#include <ctime>

/// Define the standard action term class.
/// Action terms are used in the HMC algorithm and
/// implement calculating the action itself and
/// updating the underlying fields
class action_base {
  public:
    /// Calculate and return the action
    virtual double action() { return 0; }

    /// Draw any fields with a gaussian distribution,
    /// including the momentum
    virtual void draw_gaussian_fields() {}

    /// Update the momentum with the derivative
    /// of the action term
    virtual void force_step(double eps) {}

    /// Make a copy of fields updated in a trajectory
    virtual void backup_fields() {}

    /// Restore the previous backup
    virtual void restore_backup() {}
};

/// Represents a sum of two action terms. Useful for adding them
/// to the same integrator level.
class action_sum : public action_base {
  public:
    /// left hand side
    action_base &a1;
    /// right hand side
    action_base &a2;

    /// Construct as sum of two actions
    action_sum(action_base &_a1, action_base &_a2) : a1(_a1), a2(_a2) {}

    /// Copy
    action_sum(action_sum &asum) : a1(asum.a1), a2(asum.a2) {}

    /// The action
    double action() { return a1.action() + a2.action(); }

    /// Gaussian random momentum for each element
    void draw_gaussian_fields() {
        a1.draw_gaussian_fields();
        a2.draw_gaussian_fields();
    }

    /// Update the momentum with the gauge field
    void force_step(double eps) {
        a1.force_step(eps);
        a2.force_step(eps);
    }

    /// Make a copy of fields updated in a trajectory
    void backup_fields() {
        a1.backup_fields();
        a2.backup_fields();
    }

    /// Restore the previous backup
    void restore_backup() {
        a1.restore_backup();
        a2.restore_backup();
    }
};

/// Sum operator for creating an action_sum object
action_sum operator+(action_base a1, action_base a2) {
    action_sum sum(a1, a2);
    return sum;
}

/// A base for an integrator. An integrator updates the gauge and
/// momentum fields in the HMC trajectory, approximately conserving
/// the action
class integrator_base {
  public:
    /// Return the sum of the action terms at this integrator and
    /// all levels below
    virtual double action() { return 0; }

    /// Refresh fields that can be drawn from a gaussian distribution
    /// This is needed at the beginning of a trajectory
    virtual void draw_gaussian_fields() {}

    /// Make a copy of fields updated in a trajectory
    virtual void backup_fields() {}

    /// Restore the previous backup
    virtual void restore_backup() {}

    /// Update the momentum with the gauge field
    virtual void force_step(double eps) {}

    /// Run a lower level integrator step
    virtual void step(double eps) {}
};

/// Build integrator hierarchically by adding a force step on
/// top of an existing integrator
class action_term_integrator : public integrator_base {
  public:
    /// The action term used to update the momentum on
    /// this level
    action_base &action_term;
    /// Lower level integrator, updates the momentum
    integrator_base &lower_integrator;

    /// Constructor from action and lower level integrator.
    /// also works with momentum actions as long as it inherits
    /// the integrator_base.
    action_term_integrator(action_base &a, integrator_base &i)
        : action_term(a), lower_integrator(i) {}

    /// The current total action of fields updated by this
    /// integrator. This is kept constant up to order eps^3.
    double action() { return action_term.action() + lower_integrator.action(); }

    /// Refresh fields that can be drawn from a gaussian distribution
    /// This is needed at the beginning of a trajectory
    void draw_gaussian_fields() {
        action_term.draw_gaussian_fields();
        lower_integrator.draw_gaussian_fields();
    }

    /// Make a copy of fields updated in a trajectory
    void backup_fields() {
        action_term.backup_fields();
        lower_integrator.backup_fields();
    }

    /// Restore the previous backup
    void restore_backup() {
        action_term.restore_backup();
        lower_integrator.restore_backup();
    }

    /// Update the momentum with the gauge field
    void force_step(double eps) { action_term.force_step(eps); }

    /// Update the gauge field with momentum
    void momentum_step(double eps) { lower_integrator.step(eps); }
};

/// Define an integration step for a Molecular Dynamics
/// trajectory.
class leapfrog_integrator : public action_term_integrator {
  public:
    int n = 1;
    leapfrog_integrator(action_base &a, integrator_base &i, int steps)
        : action_term_integrator(a, i), n(steps) {}
    leapfrog_integrator(action_base &a, integrator_base &i)
        : action_term_integrator(a, i) {}

    // Run the integrator update
    void step(double eps) {
        for (int i = 0; i < n; i++)
            this->lower_integrator.step(0.5 * eps / n);
        force_step(eps);
        for (int i = 0; i < n; i++)
            this->lower_integrator.step(0.5 * eps / n);
    }
};

/// Define an integration step for a Molecular Dynamics
/// trajectory.
class O2_integrator : public action_term_integrator {
  public:
    int n = 1;

    O2_integrator(action_base &a, integrator_base &i, int steps)
        : action_term_integrator(a, i), n(steps) {}
    O2_integrator(action_base &a, integrator_base &i) : action_term_integrator(a, i) {}

    // Run the integrator update
    void step(double eps) {
        double zeta = eps * 0.1931833275037836;
        double middlestep = eps - 2 * zeta;
        for (int i = 0; i < n; i++) {
            this->lower_integrator.step(zeta / n);
        }
        force_step(0.5 * eps);
        for (int i = 0; i < n; i++) {
            this->lower_integrator.step(middlestep / n);
        }
        force_step(0.5 * eps);
        for (int i = 0; i < n; i++) {
            this->lower_integrator.step(zeta / n);
        }
    }
};

#endif
