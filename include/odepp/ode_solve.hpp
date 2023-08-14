#ifndef ODEPP_ODE_SOLVE_HPP
#define ODEPP_ODE_SOLVE_HPP
#include <type_traits>

#include "odepp/integrators/forward_euler.hpp"
#include "odepp/types.hpp"
namespace odepp {

template <class Integrator, class Fn, class State, class... Args>
  requires std::is_invocable_r_v<State, Fn, RealType&, const State&, Args&&...>
auto ode_solve(Integrator&& integrate, const RealType t0, const RealType t1,
               const RealType dt, const State y0, Fn&& f, Args&&... args) {
  using OutputType = std::pair<std::vector<RealType>, std::vector<State>>;

  auto t = t0;
  auto y = y0;
  // Create the empty output
  auto output = OutputType{{t}, {y}};
  while (t < t1) {
    // Do we want one_step to update t?
    y = integrate(t, dt, y, f, std::forward<Args>(args)...);
    // We can add a check to determine if we want to log this timestep or not.
    output.first.emplace_back(t);
    output.second.emplace_back(y);
  }
  return output;
}
}  // namespace odepp
#endif