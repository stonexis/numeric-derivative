#include <iostream>
#include <stdlib.h>
#include "numerical_differentiation.hpp"

using namespace std;
int main(){
    //Создаем сетки
    auto* grid_h_2 = gen_uniform_grid(Task_const::STEP_H_2, Task_const::M_2);
    auto* grid_h = gen_uniform_grid(Task_const::H, Task_const::M);

    //Считаем функции аналитически
    auto* derivative_analytics_in_h_2 = gen_func_arr(false, grid_h_2, Task_const::M_2);
    auto* derivative_analytics_in_h = gen_func_arr(false, grid_h, Task_const::M);
    auto* func_in_h_2 = gen_func_arr(true, grid_h_2, Task_const::M_2);
    auto* func_in_h = gen_func_arr(true, grid_h, Task_const::M);
    //Считаем производные
    auto* derivative_in_h = gen_derivative_func(func_in_h, Task_const::M, Task_const::H);
    auto* derivative_in_h_2 = gen_derivative_func(func_in_h_2, Task_const::M_2, Task_const::STEP_H_2);
    long double* leading_error = nullptr;
    auto* updated_runge = gen_runge_romberg(derivative_in_h_2, derivative_in_h, &leading_error);
    //Создаем сетку M_viz
    std::size_t count_nodes_M_viz;
    auto* grid_M_viz = gen_uniform_grid_in_exist(grid_h_2, Task_const::M_2, count_nodes_M_viz);
    auto* derivative_analytics_in_M_viz = gen_func_arr(false, grid_M_viz, count_nodes_M_viz);
    auto* errors_in_h = calculate_norms(derivative_analytics_in_h, derivative_in_h, Task_const::M);
    auto* errors_in_h_2 = calculate_norms(derivative_analytics_in_h_2, derivative_in_h_2, Task_const::M_2);
    auto* errors_runge = calculate_norms(derivative_analytics_in_h_2, updated_runge, Task_const::M_2);

    write_data_to_file(grid_M_viz, grid_h, grid_h_2, derivative_analytics_in_M_viz, derivative_in_h, derivative_in_h_2, updated_runge, Task_const::M, Task_const::M_2, count_nodes_M_viz); 
    print_error_table(errors_in_h, errors_in_h_2, errors_runge, leading_error, Task_const::M);

    system("python plotter.py");

    delete[] grid_h_2;
    delete[] grid_h;
    delete[] derivative_analytics_in_h_2;
    delete[] derivative_analytics_in_h;
    delete[] func_in_h_2;
    delete[] func_in_h;
    delete[] derivative_in_h;
    delete[] derivative_in_h_2;
    delete[] leading_error;
    delete[] updated_runge;
    delete[] grid_M_viz;
    delete[] derivative_analytics_in_M_viz;
    delete[] errors_in_h;
    delete[] errors_in_h_2;
    delete[] errors_runge;
    return 0;
}