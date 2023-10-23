#ifndef ODEPP_ODE_FORWARD_EULER_HPP
#define ODEPP_ODE_FORWARD_EULER_HPP

#include "odepp/integrators/forward_euler.hpp"
#include "odepp/ode_solve.hpp"

namespace odepp {
template <class State, typename Fn, class... FnArgs>
constexpr auto ode_forward_euler(Fn&& f, const RealType t0, const RealType t1,
                                 const RealType dt, const State& y0,
                                 FnArgs&&... args) {
  return ode_solve(integrator::forward_euler, f, t0, t1, dt, y0,
                   std::forward<FnArgs>(args)...);
}
}  // namespace odepp
#endif