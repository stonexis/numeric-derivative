#include <iostream>
#include <stdlib.h>
#include "numerical_differentiation.hpp"

using namespace std;
int main(){
    std::size_t count_nodes_in_h;
    std::size_t count_nodes_in_h_2;
    std::size_t count_nodes_M_viz;
    //Считаем функции аналитически
    auto grid_func_der_analytics_in_h = gen_grid_func_and_analytic_derivative<long double>(count_nodes_in_h); //Создаем сетку с шагом h, вычисляем на ней функцию и производную
    auto grid_func_der_analytics_in_h_2 = gen_grid_func_and_analytic_derivative(count_nodes_in_h_2, grid_func_der_analytics_in_h, count_nodes_in_h, 2); // Измельчаем сетки
    auto grid_func_der_analytics_in_M_viz = gen_grid_func_and_analytic_derivative(count_nodes_M_viz, grid_func_der_analytics_in_h_2, count_nodes_in_h_2, Task_const::M_viz_ratio);
    
    //Распаковываем кортежи
    auto& [grid_h, func_in_h, derivative_analytics_in_h] = grid_func_der_analytics_in_h;
    auto& [grid_h_2, func_in_h_2, derivative_analytics_in_h_2] = grid_func_der_analytics_in_h_2;
    auto& [grid_M_viz, func_in_M_viz, derivative_analytics_in_M_viz] = grid_func_der_analytics_in_M_viz;

    //Считаем производные
    auto* derivative_in_h = gen_derivative_func(func_in_h, Task_const::M, Task_const::H);
    auto* derivative_in_h_2 = gen_derivative_func(func_in_h_2, Task_const::M_2, Task_const::STEP_H_2);
    //Уточняем производную методом Рунге
    auto runge_leading_error = gen_runge_romberg(derivative_in_h_2, derivative_in_h);

    auto& [updated_runge, leading_error] = runge_leading_error;

    //Считаем ошибки
    auto errors_in_h = calculate_norms(derivative_analytics_in_h, derivative_in_h, Task_const::M);
    auto errors_in_h_2 = calculate_norms(derivative_analytics_in_h_2, derivative_in_h_2, Task_const::M_2);
    auto errors_runge = calculate_norms(derivative_analytics_in_h_2, updated_runge, Task_const::M);
    auto* norms_leading_error = calculate_norms(leading_error, Task_const::M);

    write_data_to_file(grid_M_viz, grid_h, grid_h_2, derivative_analytics_in_M_viz, derivative_in_h, derivative_in_h_2, updated_runge, Task_const::M, Task_const::M_2, count_nodes_M_viz); 
    print_error_table(errors_in_h, errors_in_h_2, errors_runge, norms_leading_error);

    system("python plotter.py");

    delete[] grid_h_2;
    delete[] grid_h;
    delete[] grid_M_viz;

    delete[] func_in_h_2;
    delete[] func_in_h;
    delete[] func_in_M_viz;

    delete[] derivative_analytics_in_h_2;
    delete[] derivative_analytics_in_h;
    delete[] derivative_analytics_in_M_viz;

    delete[] derivative_in_h;
    delete[] derivative_in_h_2;

    delete[] leading_error;
    delete[] updated_runge;
    
    delete[] errors_in_h.first;
    delete[] errors_in_h.second;

    delete[] errors_in_h_2.first;
    delete[] errors_in_h_2.second;

    delete[] errors_runge.first;
    delete[] errors_runge.second;

    delete[] norms_leading_error;

    return 0;
}