#ifndef BALANCE_HPP
#define BALANCE_HPP

namespace balance {

template <typename V>
double compute_holding(double quantity, double inventory, V demand, V holding) {
  int n = demand.size();
  double total_demand = 0;
  double cost = 0;
  double passive_hold, active_hold;
  for (int t = 0; t < n; ++t) {
    total_demand += demand[t];
    passive_hold = total_demand > inventory ?
                   total_demand - inventory : 0;
    if (quantity < passive_hold) {
      break;
    }
    active_hold = quantity - passive_hold;
    cost += active_hold*holding[t];
  }
  return cost;
}

template <typename V>
double compute_penalty(double quantity, double inventory, V demand, V penalty) {
  double n_back = demand[0] > (inventory + quantity) ?
                  demand[0] - (inventory + quantity) : 0;
  double cost = penalty[0]*n_back;
  return cost;
}

template <typename V>
double compute_ordering(double quantity, V ordering) {
  return quantity > 0 ? ordering[0] : 0;
}

template <typename V>
double pr_cost(double quantity, double inventory, V demand, V holding, V penalty) {
  double holding_cost = compute_holding(quantity, inventory, demand, holding);
  double penalty_cost = compute_penalty(quantity, inventory, demand, penalty);
  return std::max(holding_cost, penalty_cost);
}

template <typename V>
double pr_root(double quantity, double inventory, V demand, V holding, V penalty) {
  double holding_cost = compute_holding(quantity, inventory, demand, holding);
  double penalty_cost = compute_penalty(quantity, inventory, demand, penalty);
  return holding_cost - penalty_cost;
}

template <typename V>
double ls_cost(double quantity, double inventory, V demand, V holding, V penalty, V ordering) {
  double ordering_cost = compute_ordering(quantity, ordering);
  double dual_cost = pr_cost<V>(quantity, inventory, demand, holding, penalty);
  return std::max(ordering_cost, dual_cost);
}

}

#endif//BALANCE_HPP
