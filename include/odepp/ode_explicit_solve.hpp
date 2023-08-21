#ifndef ODEPP_ODE_SOLVE_HPP
#define ODEPP_ODE_SOLVE_HPP
#include <type_traits>

#include "odepp/types.hpp"
namespace odepp {
/**
 * @brief
 *
 * @tparam State
 */
template <class State>
using ODEOutput = std::pair<std::vector<RealType>, std::vector<State>>;

/**
 * @brief Explicit Ordinary differential equations must be invocable in the form
 * of f(t,y,arguments) and return a type that is equivalent to the type of y.
 * @tparam Fn
 * @tparam State
 * @tparam FnArgs
 */
template <typename Fn, typename State, typename... FnArgs>
concept ExplicitOdeFn =
    std::is_invocable_r_v<State, Fn, const RealType, const State&, FnArgs&&...>;

/**
 * @brief
 *
 * @tparam Integrator
 * @tparam State
 * @tparam FnArgs
 * @tparam Fn
 * @param integrate
 * @param f
 * @param t0
 * @param t1
 * @param dt
 * @param y0
 * @param args
 * @return ODEOutput<State>
 */
template <class Integrator, class State, class... FnArgs,
          ExplicitOdeFn<State, FnArgs...> Fn>
auto ode_explicit_solve(Integrator&& integrate, Fn&& f, const RealType t0,
                        const RealType t1, const RealType dt, const State& y0,
                        FnArgs&&... args) -> ODEOutput<State> {
  auto t = t0;
  auto y = y0;
  // Create output with the initial timestep.
  auto output = ODEOutput<State>{{t}, {y}};
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