#ifndef ODEPP_ODE_SOLVE_HPP
#define ODEPP_ODE_SOLVE_HPP
#include <type_traits>

#include "odepp/integrators/forward_euler.hpp"
#include "odepp/types.hpp"
namespace odepp {

template <class Integrator, class Fn, class State, class... FnArgs>
  requires std::is_invocable_r_v<State, Fn, RealType&, const State&,
                                 FnArgs&&...>
auto ode_solve(Integrator&& integrate, Fn&& f, const RealType t0,
               const RealType t1, const RealType dt, const State y0,
               FnArgs&&... args) {
  using OutputType = std::pair<std::vector<RealType>, std::vector<State>>;

  auto t = t0;
  auto y = y0;
  // Create output with the initial timestep.
  auto output = OutputType{{t}, {y}};
  while (t < t1) {
    // Do we want one_step to update t?
    y = integrate(f, dt, t, y, std::forward<FnArgs>(args)...);
    // We can add a check to determine if we want to log this timestep or not.
    output.first.emplace_back(t);
    output.second.emplace_back(y);
  }
  return output;
}
}  // namespace odepp
#endif