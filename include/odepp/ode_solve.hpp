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
using OdeOutput = std::pair<std::vector<RealType>, std::vector<State>>;

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
 * @brief Explicit Ordinary differential equations must be solved using an
 * explicit integrator.
 *
 * @tparam Fn
 * @tparam State
 * @tparam FnArgs
 */
template <typename Integrator, typename Fn, typename State, typename... FnArgs>
concept ExplicitIntegrator =
    std::is_invocable_r_v<State, Integrator, Fn, const RealType, RealType&,
                          const State&, FnArgs&&...>;

/**
 * @brief
 *
 * @tparam State
 * @tparam FnArgs
 * @tparam Fn
 * @tparam Integrator
 * @param integrate
 * @param f
 * @param t0
 * @param t1
 * @param dt
 * @param y0
 * @param args
 * @return OdeOutput<State>
 */
template <class State, class... FnArgs, ExplicitOdeFn<State, FnArgs...> Fn,
          ExplicitIntegrator<Fn, State, FnArgs...> Integrator>
auto ode_solve(Integrator&& integrate, Fn&& f, const RealType t0,
               const RealType t1, const RealType dt, const State& y0,
               FnArgs&&... args) -> OdeOutput<State> {
  auto t = t0;
  auto y = y0;
  // Create output with the initial timestep.
  auto output = OdeOutput<State>{{t}, {y}};
  while (t < t1) {
    // Do we want one_step to update t?
    y = integrate(f, dt, t, y, std::forward<FnArgs>(args)...);
    // We can add a check to determine if we want to log this timestep or not.
    output.first.emplace_back(t);
    output.second.emplace_back(y);
  }
  return output;
}

/**
 * @brief Implicit Ordinary differential equations must be invocable in the form
 * of f(t,y,yp,arguments) and return a type that is equivalent to the type of y.
 *
 * @tparam Fn
 * @tparam State
 * @tparam FnArgs
 */
template <typename Fn, typename State, typename... FnArgs>
concept ImplicitOdeFn =
    std::is_invocable_r_v<State, Fn, const RealType, const State&, const State&,
                          FnArgs&&...>;

/**
 * @brief Implicit Ordinary Differential Equations must be solved using an
 * Implicit integrator.
 *
 * @tparam Fn
 * @tparam State
 * @tparam FnArgs
 */
template <typename Integrator, typename Fn, typename State, typename... FnArgs>
concept ImplicitIntegrator =
    std::is_invocable_r_v<State, Integrator, Fn, const RealType, RealType&,
                          const State&, const State&, FnArgs&&...>;

/**
 * @brief
 *
 * @tparam State
 * @tparam FnArgs
 * @tparam Fn
 * @tparam Integrator
 * @param integrate
 * @param f
 * @param t0
 * @param t1
 * @param dt
 * @param y0
 * @param yp0
 * @param args
 * @return OdeOutput<State>
 */
template <class State, class... FnArgs, ImplicitOdeFn<State, FnArgs...> Fn,
          ImplicitIntegrator<Fn, State, FnArgs...> Integrator>
auto ode_solve(Integrator&& integrate, Fn&& f, const RealType t0,
               const RealType t1, const RealType dt, const State& y0,
               const State& yp0, FnArgs&&... args) -> OdeOutput<State> {
  auto t = t0;
  auto y = y0;
  auto yp = yp0;
  // Create output with the initial timestep.
  auto output = OdeOutput<State>{{t}, {y}};
  while (t < t1) {
    // Do we want one_step to update t?
    y = integrate(f, dt, t, y, yp, std::forward<FnArgs>(args)...);
    // We can add a check to determine if we want to log this timestep or
    // not .
    output.first.emplace_back(t);
    output.second.emplace_back(y);
  }
  return output;
}
}  // namespace odepp
#endif